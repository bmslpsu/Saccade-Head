function [] = MakeData_Static_HeadFree_Sacd(Fc,direction)
%% MakeData_Static_HeadFree_Sacd:
%   INPUTS:
%       Fc          : head data low-pass filter cutoff frequency [Hz]
%       direction   : only get saccades in this direction
%
%   OUTPUTS:
%       -
%

% Fc = 40;
% direction = 0; % get saccades in all directions

switch direction
    case 0
        dirlabel = 'ALL';
    case 1
        dirlabel = 'CW';
    case -1
        dirlabel = 'CCW';
    otherwise
        error('direction must be 0,1, or -1')
end
    

% Data location
rootdir = 'H:\EXPERIMENTS\RIGID\Experiment_Static_Wave';

% Output file name
filename = ['NewStatic_HeadFree_SACCD_' dirlabel '_filt=' num2str(Fc)];

% Setup Directories 
PATH.daq = rootdir; % DAQ data location
PATH.head = fullfile(PATH.daq,'\tracked_head'); % tracked head data location
PATH.benifly = fullfile(PATH.daq,'\tracked_head_wing'); % tracked benifly data location

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.head,'*.mat',false,'fly','trial','wave','vel');

%% Get Data %%
disp('Loading...')
showplot = false;
tintp = (0:(1/200):(10 - 1/200))';
SACCADE = [I , table(num2cell(zeros(N.file,1)))]; % store saccade objects
SACCADE.Properties.VariableNames{5} = 'saccade';
ALL_DATA = cell(N.fly,N.wave);
COUNT = cell(N.fly,N.wave);
SACCADE_STATS = []; % store saccade stats
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    load(fullfile(PATH.head, [basename{kk} '.mat']),'hAngles'); % load head angles % time arrays
    % benifly = ImportBenifly(fullfile(PATH.benifly, [basename{kk} '.csv'])); % load head & wing benifly angles
    % hAngles = rad2deg(benifly.Head);
    
    % Sync video with trigger & pattern
    trig.raw_time   = t_p; % DAQ raw times for trigger
    trig.pos        = round(data(:,1)); % trigger values
    trig.diff       = diff(trig.pos);
    [~,trig.locs]   = findpeaks(trig.diff); % where each frame starts
    trig.time       = trig.raw_time(trig.locs+1); % where each frame starts
    trig.time       = trig.time - trig.time(1);
   	
    % Get head data
    Head = process_signal(trig.time, hAngles, Fc, tintp, [4 8 16 32 64]);
    
    % Get Saccade Stats
    peaks = []; % find peaks automatically
    head_saccade = saccade(Head.X(:,1), Head.Time, 3.5, direction, peaks, showplot);
    % figure (1)
    
    COUNT{I.fly(kk),I.wave(kk)} = head_saccade.count;
    
    if head_saccade.count==0
        rep = 1;
        SACCADE{kk,5} = {[]};
        % ALL_DATA{I.fly(kk),I.wave(kk)}(end+1,1) = saccade();
    else
        SACCADE{kk,5} = {head_saccade};
        ALL_DATA{I.fly(kk),I.wave(kk)}(end+1,1) = head_saccade;
        rep = head_saccade.count;
    end
    
    VTable = table(D.wave(kk),'VariableNames',{'Wave'});
    ITable = [I(kk,:),VTable];
    ITable = repmat(ITable,rep,1);
    SACCADE_STATS = [SACCADE_STATS ; [ITable , head_saccade.SACD]];
    
    if showplot
        pause
        close all
    end
end

% Fill in empty saccade trials
empty_idx = cellfun(@(x) isempty(x), ALL_DATA);
ALL_DATA(empty_idx) = {saccade(nan*Head.X(:,1), nan*Head.Time, 350, 0, [], false)};

%% Extract & group saccades & intervals by speed & by fly
fields = {'normpeak_saccade','norm_interval','normstart_interval','normend_interval'};
nfield = length(fields);
norm_fields = {'time','position','velocity'};

center = 0; % normalization cenetr for saccades & inter-saccade intervals
dim = 1; % dimension to center

FLY = []; % all trials per speed per fly
GRAND = []; % all trials per speed
for kk = 1:nfield % for each field in the saccade structure
    if kk < 2 % even padding around center only for saccades (not intervals)
        even = true;
    else
        even = false;
    end
            
    FLY.(fields{kk}) = cellfun(@(x) struct_center([x.(fields{kk})], center, even, dim, norm_fields), ...
                            ALL_DATA, 'UniformOutput', true);
    for jj = 1:N.wave
        GRAND.(fields{kk})(:,jj) = struct_center(FLY.(fields{kk})(:,jj), center, even, dim, norm_fields);
    end
end

%% SAVE %%
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','D','I','U','N','T','direction','-v7.3')
disp('SAVING DONE')
end