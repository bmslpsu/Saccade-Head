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
% filename = ['NewRamp_HeadFree_SACCD_Anti_filt=' num2str(Fc) '_Wave=' num2str(wave)];

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
[wing.b,wing.a] = butter(2,Fc/(Fs/2),'low');
[head.b,head.a] = butter(2,40/(Fs/2),'low');
Vel = U.vel{1};
Stim = (Vel*tintrp')';
SACCADE = [I , splitvars(table(num2cell(zeros(N.file,3))))]; % store saccade objects
SACCADE.Properties.VariableNames(5:7) = {'head_saccade','wing_saccade','head2wing'}; 
Head_Wing = cell(N.fly,N.vel);
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
    
    % Get head data
    head.pos = hAngles;
    head.pos = interp1(trig.time, head.pos, tintrp, 'pchip');
    head.pos = filtfilt(head.b, head.a, head.pos);
    
  	% Get wing data
  	wing.left  = hampel(tintrp, benifly.LWing);
    wing.right = hampel(tintrp, benifly.RWing);
    wing.left  = interp1(trig.time, wing.left, tintrp, 'pchip');
    wing.right = interp1(trig.time, wing.right, tintrp, 'p chip');
    wing.left  = filtfilt(wing.b, wing.a, wing.left);
    wing.right = filtfilt(wing.b, wing.a, wing.right);
    wing.dwba  = rad2deg(wing.left - wing.right); 
    [b_high,a_high] = butter(2,0.5/(200/2), 'high');
    wing.wba = filtfilt(b_high, a_high, wing.dwba);
    
    % Get head saccades
    peaks = [];
    direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    head_saccade = saccade(head.pos, tintrp, 250, amp_cut, direction, peaks, nan, showplot);
    head_saccade = stimSaccade(head_saccade, Stim(:,I.vel(kk)), false); % with approximate pattern position   
    
    % Get wing saccades
    [dwba_pks,~,~] = wingSaccade(wing.dwba, tintrp, direction, showplot);
    wing_saccade = saccade(wing.dwba, tintrp, 0, 0, direction, dwba_pks, 0.5, showplot);
    
    % head2wing
    win = 0.2; % [s]
    head2wing = head_wing_saccade(head_saccade, wing_saccade, win, true);
    pause
    
    if head_saccade.count==0
        % skip
     	SACCADE{kk,5} = {[]};
        SACCADE{kk,6} = {[]};
        SACCADE{kk,7} = {[]};
    else
        Head_Wing{I.fly(kk),I.vel(kk)}(end+1,1) = head2wing;
        SACCADE{kk,5} = {head_saccade};
        SACCADE{kk,6} = {wing_saccade};
        SACCADE{kk,7} = {head2wing};
    end


    if showplot
        pause
        close all
    end
end

% Fill in empty saccade trials
empty_idx = cellfun(@(x) isempty(x), Head_Wing);
Head_Wing(empty_idx) = {saccade(nan*head.pos, nan*head.pos, 300, 0, 0, [], [], false)};

empty_idx = cellfun(@(x) isempty(x), SACCADE.head2wing);
SACCADE = SACCADE(~empty_idx,:);

%% Extract & group saccades & intervals by speed & by fly
Tlag_all = cellfun(@(x) x.timelags, SACCADE.head2wing, 'UniformOutput', false);
Tlag_all = mean(cat(2,Tlag_all{:}),2);

Acor_all = cellfun(@(x) x.acor, SACCADE.head2wing, 'UniformOutput', false);
Acor_all = cat(2,Acor_all{:});

Tlag_sacd_all = cellfun(@(x) x.int_lags, SACCADE.head2wing, 'UniformOutput', false);
Tlag_sacd_all = mean(cat(2,Tlag_sacd_all{:}),2);

Acor_sacd_all = cellfun(@(x) x.int_acor, SACCADE.head2wing, 'UniformOutput', false);
Acor_sacd_all = cat(2,Acor_sacd_all{:});

fig = figure (1);
set(fig, 'Color', 'w')
    ax(1) = subplot(2,1,1) ; cla ; hold on ; ylabel('Cross Correlation')
            plot(Tlag_all, Acor_all, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1)
            plot(Tlag_all, mean(Acor_all,2), 'Color', 'k', 'LineWidth', 2)
    ax(2) = subplot(2,1,2) ; cla ; hold on ; ylabel('Cross Correlation')
            plot(Tlag_sacd_all, Acor_sacd_all, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1)
            plot(Tlag_sacd_all, nanmean(Acor_sacd_all,2), 'Color', 'green', 'LineWidth', 2)
            xlabel('Time Lag (ms)')

 	set(ax, 'LineWidth', 1.5, 'XLim', 0.5*[-1 1], 'YLim', [-1 1])
    
%% SAVE %%
% disp('Saving...')
% save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
%       'PATH','ALL_DATA','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N','T','-v7.3')
% disp('SAVING DONE')
end