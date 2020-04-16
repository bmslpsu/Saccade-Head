function [] = Saccade_HeadWing()
%% Saccade_HeadWing:
wave = 30;
Fc = 40;

% Data location
rootdir = ['H:\EXPERIMENTS\RIGID\Experiment_Asymmetry_Control_Verification\HighContrast\' num2str(wave)];

% Output file name
% filename = ['NewRampStim_HeadFree_SACCD_Anti_filt=' num2str(Fc) '_Wave=' num2str(wave)];

% Setup Directories 
PATH.daq = rootdir; % DAQ data location
PATH.vid = fullfile(PATH.daq,'\Vid'); % video data location
PATH.head = fullfile(PATH.daq,'\Vid\tracked_head'); % tracked head kinematic data location
PATH.wing = fullfile(PATH.daq,'\Vid\tracked_head_wing'); % tracked wing kinematic data location

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.head,'*.mat',false,'fly','trial','vel','wave');

%% Get Data %%
disp('Loading...')
showplot = true;
tintrp = (0:(1/200):(10 - 1/200))';
Vel = U.vel{1};
Stim = (Vel*tintrp')';
SACCADE = [I , table(num2cell(zeros(N.file,1)))]; % store saccade objects
SACCADE.Properties.VariableNames{5} = 'saccade'; 
ALL_DATA = cell(N.fly,N.vel);
COUNT = cell(N.fly,N.vel);
SACCADE_STATS = []; % store saccade stats
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    load(fullfile(PATH.head, [basename{kk} '.mat']),'hAngles'); % load head angles % time arrays
    % benifly = ImportBenifly(fullfile(PATH.wing, [basename{kk} '.csv'])); % load head angles
    
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
    
    % Get wing data
	wing.left  = data(:,4);
    wing.right = data(:,5);
%     wing.left  = interp1(pat.time, data(:,4), tintrp, 'nearest');
%     wing.right = interp1(pat.time, data(:,5), tintrp, 'nearest');
%     wing.left  = interp1(trig.time, benifly.LWing, tintrp, 'nearest');
%     wing.right = interp1(trig.time, benifly.RWing, tintrp, 'nearest');
%     wing.left  = hampel(tintrp, wing.left);
%     wing.right = hampel(tintrp, wing.right);
%     [b,a] = butter(2,0.3,'low');
%     wing.left  = filtfilt(b, a, wing.left);
%     wing.right = filtfilt(b, a, wing.right);

    wing.wba   = wing.left - wing.right;
    Wing       = process_signal(t_p, wing.wba , 100, tintrp, [4 8 16 32 64]);
    
    % Get head data
    % benifly.Head(1) = benifly.Head(2);
    Head = process_signal(trig.time, hAngles, Fc, tintrp, [4 8 16 32 64]);
    
    % Get Saccade Stats
    peaks = [];
    direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    % direction = sign(D.vel(kk));
    head_saccade = saccade(Head.X(:,1), Head.Time, 350, direction, peaks, showplot);
    head_saccade = stimSaccade(head_saccade, Stim(:,I.vel(kk)), false); % with approximate pattern position
    % head_saccade = stimSaccade(head_saccade, pat.pos, false); % with actual pattern position
    
    wing_saccade = saccade(Wing.X(:,1), Wing.Time, 3, 0, peaks, showplot);
    
    COUNT{I.fly(kk),I.vel(kk)}(end+1,1) = head_saccade.count;
    
%     if head_saccade.count==0
%         rep = 1;
%         SACCADE{kk,5} = {head_saccade};
%         % ALL_DATA{I.fly(kk),I.vel(kk)}(end+1,1) = saccade();
%     else
%         SACCADE{kk,5} = {head_saccade};
%         ALL_DATA{I.fly(kk),I.vel(kk)}(end+1,1) = head_saccade;
%         rep = head_saccade.count;
%     end
%     VTable = table(D.vel(kk),'VariableNames',{'Vel'});
%     ITable = [I(kk,:),VTable];
%     ITable = repmat(ITable,rep,1);
%     SACCADE_STATS = [SACCADE_STATS ; [ITable , head_saccade.SACD]];
    
    if showplot
        pause
        close all
    end
end


end