function [] = MakeData_Static_HeadFree_WingSacd(Fc,direction)
%% MakeData_Static_HeadFree_WingSacd:
%
%   INPUTS:
%       Fc          : head data low-pass filter cutoff frequency [Hz]
%       direction   : only get saccades in this direction
%
%   OUTPUTS:
%       -
%

Fc = 40;
direction = 0; % get saccades in all directions

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
filename = ['NewStatic_HeadFree_HeadWingSACCD_' dirlabel '_filt=' num2str(Fc)];

% Setup Directories 
PATH.daq = rootdir; % DAQ data location
PATH.head = fullfile(PATH.daq,'\tracked_head'); % tracked head data location
PATH.benifly = fullfile(PATH.daq,'vid_filt\tracked_head_wing'); % tracked benifly data location

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.benifly,'*.csv',false,'fly','trial','wave','vel');

%% Get Data %%
disp('Loading...')
showplot = true;
tintrp = (0:(1/200):10)';
SACCADE = [I , splitvars(table(num2cell(zeros(N.file,2))))]; % store saccade objects
SACCADE.Properties.VariableNames{5} = 'head_saccade';
SACCADE.Properties.VariableNames{5} = 'wing_saccade';

HEAD_DATA = cell(N.fly,1);
WING_DATA = cell(N.fly,1);

HEAD_COUNT = cell(N.fly,1);
WING_COUNT = cell(N.fly,1);

HEAD_STATS = []; % store saccade stats
WING_STATS = []; % store saccade stats

[wing.b,wing.a] = butter(2,20/(200/2),'low');
[head.b,head.a] = butter(2,Fc/(200/2),'low');

wing_amp_cut = 8;
wing_sacd_length = 0.8;
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    % load(fullfile(PATH.head, [basename{kk} '.mat']),'hAngles'); % load head angles % time arrays
    benifly = ImportBenifly(fullfile(PATH.benifly, [basename{kk} '.csv'])); % load head & wing benifly angles
    
    % Sync video with trigger & pattern
    trig.raw_time   = t_p; % DAQ raw times for trigger
    trig.pos        = round(data(:,1)); % trigger values
    trig.diff       = diff(trig.pos);
    [~,trig.locs]   = findpeaks(trig.diff); % where each frame starts
    trig.time       = trig.raw_time(trig.locs+1); % where each frame starts
    trig.time       = trig.time - trig.time(1) - 0.1;
    
    % Get pattern data
    pat.time = t_p;
    pat.pos = 3.75*round((96/10)*data(:,2));
    pat.pos = interp1(pat.time, pat.pos, tintrp, 'nearest');
   	pat.pos = pat.pos - pat.pos(1);
    
    % Get head data
    head.pos = rad2deg(benifly.Head);
    % head.pos = hAngles;
    head.pos = interp1(trig.time, head.pos, tintrp, 'pchip');
    head.pos = filtfilt(head.b, head.a, head.pos);
    
  	% Get wing data
	% wing.left  = data(:,4);
    % wing.right = data(:,5);
	wing.left  = rad2deg(benifly.LWing);
    wing.right = rad2deg(benifly.RWing);
    
    wing.left  = interp1(trig.time, wing.left, tintrp, 'pchip');
    wing.right = interp1(trig.time, wing.right, tintrp, 'pchip');
    
    wing.left  = hampel(tintrp, wing.left);
    wing.right = hampel(tintrp, wing.right);
    wing.left  = filtfilt(wing.b, wing.a, wing.left);
    wing.right = filtfilt(wing.b, wing.a, wing.right);
    
    wing.wba = wing.left - wing.right;
    
    % Get Saccade Stats
    peaks = []; % find peaks automatically
    wing_saccade = saccade(wing.wba, tintrp, 3, direction, peaks, showplot, wing_amp_cut);
	%head_saccade = saccade(head.pos, tintrp, 350, direction, peaks, false);
    figure (1)
    
%     if any(any( abs(wing_saccade.normpeak_saccade.velocity) > 2000 ))
%         wing_saccade = saccade(wing.wba, tintrp, 3, direction, peaks, true, wing_amp_cut);
%         pause
%     end
    
    HEAD_COUNT{I.fly(kk),1}(end+1,1) = head_saccade.count;
   	WING_COUNT{I.fly(kk),1}(end+1,1) = wing_saccade.count;
    
    if head_saccade.count==0
        head_rep = 1;
        SACCADE{kk,5} = {[]};
    else
        SACCADE{kk,5} = {head_saccade};
        HEAD_DATA{I.fly(kk),1}(end+1,1) = head_saccade;
        head_rep = head_saccade.count;
    end
    
    if wing_saccade.count==0
        wing_rep = 1;
        SACCADE{kk,6} = {[]};
    else
        SACCADE{kk,6} = {head_saccade};
        WING_DATA{I.fly(kk),1}(end+1,1) = wing_saccade;
        wing_rep = wing_saccade.count;
    end
    
    VTable = table(D.wave(kk),'VariableNames',{'Wave'});
    
    ITable = [I(kk,:),VTable];
    ITable = repmat(ITable,head_rep,1);
    HEAD_STATS = [HEAD_STATS ; [ITable , head_saccade.SACD]];
    
    ITable = [I(kk,:),VTable];
    ITable = repmat(ITable,wing_rep,1);
    WING_STATS = [WING_STATS ; [ITable , wing_saccade.SACD]];
    
    if showplot
        pause
        close all
    end
end

% Fill in empty saccade trials
empty_idx = cellfun(@(x) isempty(x), HEAD_DATA);
HEAD_DATA(empty_idx) = {saccade(nan*head.pos, nan*head.pos, 3, 0, [], false)};

empty_idx = cellfun(@(x) isempty(x), WING_DATA);
WING_DATA(empty_idx) = {saccade(nan*wing.wba, nan*wing.wba, 3, 0, [], false)};

%% Extract & group saccades & intervals by speed & by fly
fields = {'normpeak_saccade','norm_interval','normstart_interval','normend_interval'};
nfield = length(fields);
norm_fields = {'time','position','velocity'};

center = 0; % normalization cenetr for saccades & inter-saccade intervals
dim = 1; % dimension to center

HEAD_FLY = []; % all trials per speed per fly
HEAD_GRAND = []; % all trials per speed
WING_FLY = []; % all trials per speed per fly
WING_GRAND = []; % all trials per speed
for kk = 1:nfield % for each field in the saccade structure
    if kk < 2 % even padding around center only for saccades (not intervals)
        even = true;
    else
        even = false;
    end
            
    HEAD_FLY.(fields{kk}) = cellfun(@(x) struct_center([x.(fields{kk})], center, even, dim, norm_fields), ...
                            HEAD_DATA, 'UniformOutput', true);
                        
    WING_FLY.(fields{kk}) = cellfun(@(x) struct_center([x.(fields{kk})], center, even, dim, norm_fields), ...
                            WING_DATA, 'UniformOutput', true);
                        
    for jj = 1:1
        HEAD_GRAND.(fields{kk})(:,jj) = struct_center(HEAD_FLY.(fields{kk})(:,jj), ...
            center, even, dim, norm_fields);
        WING_GRAND.(fields{kk})(:,jj) = struct_center(WING_FLY.(fields{kk})(:,jj), ...
            center, even, dim, norm_fields);
    end
end

%% SAVE %%
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','HEAD_COUNT','WING_COUNT','SACCADE','HEAD_STATS','WING_STATS','HEAD_FLY','HEAD_GRAND',...
      'WING_FLY','WING_GRAND','D','I','U','N','T','direction','-v7.3')
disp('SAVING DONE')
end