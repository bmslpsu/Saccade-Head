function [] = MakeData_Ramp_HeadFree_Sacd_HeadWing(wave,Fc)
%% MakeData_Ramp_HeadFree_Sacd_HeadWing:
%   INPUTS:
%       wave    :   spatial wavelength of data
%       Fc      :   head data low-pass filter cutoff frequency [Hz]
%
%   OUTPUTS:
%       -
%

wave = 30;
Fc = 10;

% Data location
rootdir = ['H:\EXPERIMENTS\RIGID\Experiment_Asymmetry_Control_Verification\HighContrast\' num2str(wave)];

% Output file name
filename = ['NewRamp_HeadFree_SACCD_Anti_filt=' num2str(Fc) '_Wave=' num2str(wave)];

% Setup Directories 
PATH.daq  = rootdir; % DAQ data location
PATH.vid  = fullfile(PATH.daq,'Vid','vid_filt'); % video data location
PATH.head = fullfile(PATH.daq,'Vid','tracked_head'); % tracked kinematic data location
PATH.wing = fullfile(PATH.vid,'tracked_head_wing'); % tracked kinematic data location

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.wing,'*.csv',false,'fly','trial','vel','wave');

%% Get Data %%
disp('Loading...')
showplot = false;
amp_cut = 7;

Fs = 200;
tintrp = (0:(1/Fs):(10 - 1/Fs))';
[wing.b,wing.a] = butter(2,10/(Fs/2),'low');
[head.b,head.a] = butter(2,10/(Fs/2),'low');
Vel = U.vel{1};
Stim = (Vel*tintrp')';
SACCADE = [I , splitvars(table(num2cell(zeros(N.file,2))))]; % store saccade objects
SACCADE.Properties.VariableNames(5:6) = {'head_saccade','dWBA'}; 
ALL_DATA = cell(N.fly,N.vel);
COUNT = cell(N.fly,N.vel);
SACCADE_STATS = []; % store saccade stats
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    load(fullfile(PATH.head, [basename{kk} '.mat']),'hAngles'); % load head angles % time arrays
    benifly = ImportBenifly(fullfile(PATH.wing, [basename{kk} '.csv'])); % load head & wing angles
    
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
    % head.pos = rad2deg(benifly.Head);
    head.pos = hAngles;
    head.pos = interp1(trig.time, head.pos, tintrp, 'pchip');
    head.pos = filtfilt(head.b, head.a, head.pos);
    % Head = process_signal(trig.time, , Fc, tintrp, [4 8 16 32 64]);
    
  	% Get wing data
	% wing.left  = data(:,4);
    % wing.right = data(:,5);
    % wing.left  = interp1(pat.time, data(:,4), tintrp, 'nearest');
    % wing.right = interp1(pat.time, data(:,5), tintrp, 'nearest');
  	wing.left  = hampel(tintrp, benifly.LWing);
    wing.right = hampel(tintrp, benifly.RWing);
    wing.left  = interp1(trig.time, wing.left, tintrp, 'pchip');
    wing.right = interp1(trig.time, wing.right, tintrp, 'pchip');

    wing.left  = filtfilt(wing.b, wing.a, wing.left);
    wing.right = filtfilt(wing.b, wing.a, wing.right);
    wing.wba = wing.left - wing.right;
    
    [b_high,a_high] = butter(2,0.5/(200/2), 'high');
    wing.wba = filtfilt(b_high, a_high, wing.wba);
    
    % Get Saccade Stats
    peaks = [];
    direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    % direction = sign(D.vel(kk));
    head_saccade = saccade(head.pos, tintrp, 300, amp_cut, direction, peaks, nan, showplot);
    head_saccade = stimSaccade(head_saccade, Stim(:,I.vel(kk)), false); % with approximate pattern position
    % head_saccade = stimSaccade(head_saccade, pat.pos, false); % with actual pattern position
    % figure (1)
    
    % [acor,timeLag,maxCC,timeDiff] = CrossCorr(head_saccade.velocity, wing.wba, 200, 'normalized');
    maxlag = Fs*0.1;
    [acor,lag] = xcorr(head_saccade.position, wing.wba, maxlag);
    timeLag = lag/Fs;
    [~,idx] = max(abs(acor));
    maxCC = acor(idx);
    timeDiff = timeLag(idx);
    
    disp(timeDiff)
 	figure (1)
    subplot(4,1,1) ; cla
        plot(tintrp, head_saccade.position, 'c', 'LineWidth', 1)
    subplot(4,1,2) ; cla
        plot(tintrp, head_saccade.velocity, 'b', 'LineWidth', 1)
    subplot(4,1,3) ; cla
        plot(tintrp, wing.wba, 'r', 'LineWidth', 1)
    subplot(4,1,4) ; cla ; hold on
        plot(timeLag, acor, 'k', 'LineWidth', 1)
        plot(timeDiff, maxCC, 'ro', 'MarkerSize', 10)
    pause()
    
    COUNT{I.fly(kk),I.vel(kk)}(end+1,1) = head_saccade.count;
    
    if head_saccade.count==0
        rep = 1;
        SACCADE{kk,5} = {head_saccade};
        % ALL_DATA{I.fly(kk),I.vel(kk)}(end+1,1) = saccade();
    else
        SACCADE{kk,5} = {head_saccade};
        ALL_DATA{I.fly(kk),I.vel(kk)}(end+1,1) = head_saccade;
        rep = head_saccade.count;
    end
    SACCADE{kk,6} = {wing.wba};
    VTable = table(D.vel(kk),'VariableNames',{'Vel'});
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
ALL_DATA(empty_idx) = {saccade(nan*head.pos, nan*head.pos, 300, 0, 0, [], [], false)};

%% Extract & group saccades & intervals by speed & by fly
fields = {'normpeak_saccade','norm_interval','normstart_interval','normend_interval',...
    'normstart_stimulus','error','int_error'};
nfield = length(fields);
norm_fields = {'time','position','velocity'};

center = 0; % normalization center for saccades & inter-saccade intervals
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
    for jj = 1:N.vel
        GRAND.(fields{kk})(:,jj) = struct_center(FLY.(fields{kk})(:,jj), center, even , dim, norm_fields);
    end
end

%% SAVE %%
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','ALL_DATA','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end