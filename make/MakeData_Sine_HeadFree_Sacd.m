function [] = MakeData_Sine_HeadFree_Sacd(amp,Fc,direction)
%% MakeData_Ramp_HeadFree_Sacd:
%   INPUTS:
%       amp         :   amplitude of sineusoid stimulus
%       Fc          :   head data low-pass filter cutoff frequency [Hz]
%       direction   :   only get saccades in this direction
%
%   OUTPUTS:
%       -
%

amp = 15;
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
rootdir = ['H:\EXPERIMENTS\RIGID\Experiment_Sinusoid\' num2str(amp)];

% Output file name
filename = ['SS_HeadFree_SACCD_' dirlabel '_filt=' num2str(Fc) '_Amp=' num2str(amp)];

% Setup Directories 
PATH.daq = rootdir; % DAQ data location
PATH.vid = fullfile(PATH.daq,'\Vid'); % video data location
% PATH.ang = fullfile(PATH.daq,'\Vid\tracked'); % tracked kinematic data location
PATH.ang = fullfile(PATH.daq,'\Vid\tracked_head'); % tracked kinematic data location

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.ang,'*.mat',false);

%% Get Data %%
disp('Loading...')
showplot = false;
thresh = 300;
amp_cut = 7;
tintrp = (0:(1/200):(10 - 1/200))';
SACCADE = [I , table(num2cell(zeros(N.file,1)))]; % store saccade objects
SACCADE.Properties.VariableNames{4} = 'saccade'; 
ALL_DATA = cell(N.fly,N.freq);
COUNT = cell(N.fly,N.freq);
SACCADE_STATS = []; % store saccade stats
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    load(fullfile(PATH.ang, [basename{kk} '.mat']),'hAngles'); % load head angles % time arrays
    % benifly = ImportBenifly(fullfile(PATH.ang, FILES{kk})); % load head angles
    
    % Sync video with trigger & pattern
    trig.raw_time   = t_p; % DAQ raw times for trigger
    trig.pos        = round(data(:,1)); % trigger values
    trig.diff       = diff(trig.pos); % trigger derivative (rising edge triggers frame)
    [~,trig.locs]   = findpeaks(trig.diff); % where each frame starts
    trig.time       = [0 ; trig.raw_time(trig.locs+1)]; % where each frame starts
    
  	% Get pattern data
    pat.time = t_p;
    pat.pos = 3.75*round((96/10)*data(:,2));
    pat.pos = interp1(pat.time, pat.pos, tintrp, 'nearest');
   	pat.pos = pat.pos - pat.pos(1);
    
    % Get head data
    % benifly.Head(1) = benifly.Head(2);
    Head = process_signal(trig.time, hAngles, Fc, tintrp, [4 8 16 32 64]);
    
    % Get Saccade Stats
    peaks = [];
    head_saccade = saccade(Head.X(:,1), Head.Time, thresh, amp_cut, direction, peaks, nan, showplot);
    % figure (1) ; suptitle(num2str(D.freq(kk)))
    head_saccade = stimSaccade(head_saccade, pat.pos, false); % with actual pattern position
    
    COUNT{I.fly(kk),I.freq(kk)}(end+1,1) = head_saccade.count;
    
    if head_saccade.count==0
        rep = 1;
        SACCADE{kk,4} = {head_saccade};
        % ALL_DATA{I.fly(kk),I.vel(kk)}(end+1,1) = saccade();
    else
        SACCADE{kk,4} = {head_saccade};
        ALL_DATA{I.fly(kk),I.freq(kk)}(end+1,1) = head_saccade;
        rep = head_saccade.count;
    end
    VTable = table(D.freq(kk),'VariableNames',{'Vel'});
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
ALL_DATA(empty_idx) = {saccade(nan*Head.X(:,1), nan*Head.Time, 300, 0, 0, [], [], false)};

all = cellfun(@(x) x.velocity, SACCADE.saccade, 'UniformOutput', false);
all = cat(1,all{:});

 %% Extract & group saccades & intervals by speed & by fly
fields = {'normpeak_saccade','norm_interval','normstart_interval','normend_interval',...
    'normstart_stimulus','error','int_error'};
nfield = length(fields);
norm_fields = {'time','position','velocity'};

center = 0; % normalization center for saccades & inter-saccade intervals
dim = 1; % dimension to center

clc
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
    for jj = 1:N.freq
        GRAND.(fields{kk})(:,jj) = struct_center(FLY.(fields{kk})(:,jj),center, even , dim, norm_fields);
    end
end

%% SAVE %%
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','ALL_DATA','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end