function [] = Make_Ramp_HeadFree_head_roll()
%% Make_Ramp_HeadFree_head_roll:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%
warning('off', 'signal:findpeaks:largeMinPeakHeight')

% Data location
rootdir = 'E:\EXPERIMENTS\RIGID\Experiment_Ramp_forRoll';

% Output file name
filename = 'Ramp_HeadFree_head_roll';

% Setup Directories 
PATH.daq  = rootdir; % DAQ data location
PATH.head = fullfile(PATH.daq,'tracked_head_edge'); % tracked kinematic data location

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.head,'*.mat',false,'fly','trial','vel');

%% Get Data %%
clc
disp('Loading...')
Fs = 200; % sampling frequency [s]
tintrp = (0:(1/Fs):(10))'; % time vector for interpolation

% HEAD saccade detection parameters
head.showplot = false;
head.Fc_detect = [10 nan];
head.Fc_ss = [nan nan];
head.amp_cut = 4;
head.dur_cut = inf;
head.thresh = [70 , 2, 1, 0];
head.true_thresh = 250;
head.sacd_length = nan;
head.pks = [];
head.min_pkdist = 0.1;
head.min_pkwidth = 0.03;
head.min_pkprom = 75;
head.min_pkthresh = 0;
head.boundThresh = [0.2 60];
head_carry.Fc = 40;
[head_carry.b, head_carry.a] = butter(3, head_carry.Fc / (Fs/2) ,'low');

Vel = U.vel{1};
Stim = (Vel*tintrp')';
SACCADE = [I , splitvars(table(num2cell(zeros(N.file,6))))];
SACCADE.Properties.VariableNames(4:9) = {'head_yaw','head_roll','yaw_scd_pos','yaw_scd_vel',...
                                            'roll_scd_pos','roll_scd_vel'};
HEAD_DATA = cell(N.fly,N.vel);
HEAD_SACCADE_STATS = [];
hampel_dx = 5;
scd_win = 0.5;
roll_cal = 70;
for kk = 1:N.file
    disp(kk)
    %disp(basename{kk})
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p');
    load(fullfile(PATH.head, [basename{kk} '.mat']),'yaw','roll_idx');
 
  	% Sync frames and get pattern data
	[trig,~] = sync_pattern_trigger(t_p, data(:,2), 10, data(:,1), true, [], false, false);
    
    % Get head data
    head.yaw = hampel(trig.time_sync, yaw, hampel_dx);
    head.roll = hampel(trig.time_sync, roll_idx, hampel_dx);
    
    head.yaw = interp1(trig.time_sync, head.yaw , tintrp, 'pchip');
    head.roll = interp1(trig.time_sync, roll_cal*head.roll, tintrp, 'pchip');
    
    head.yaw_vel = central_diff(head.yaw, 1/Fs);
    head.roll_vel = central_diff(head.roll, 1/Fs);
    
	% Head with filter
 	head.yaw_filt = filtfilt(head_carry.b, head_carry.a, head.yaw);
    head.roll_filt = filtfilt(head_carry.b, head_carry.a, head.roll);
    head.yaw_vel_filt = central_diff(head.yaw_filt, 1/Fs);
    head.roll_vel_filt = central_diff(head.roll_filt, 1/Fs);
    
 	%head.yaw_vel_filt = central_diff(head.yaw, 1/Fs);
    %head.yaw_vel_filt = filtfilt(head_carry.b, head_carry.a, head.yaw_vel_filt);
 	%head.roll_vel_filt = central_diff(head.roll, 1/Fs);
    %head.roll_vel_filt = filtfilt(head_carry.b, head_carry.a, head.roll_vel_filt);
       
    % Extract head yaw saccades
    direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    head_yaw = saccade_all(head.yaw_filt, tintrp, head.thresh, head.true_thresh, head.Fc_detect, ...
                                head.Fc_ss, head.amp_cut, head.dur_cut, direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, head.showplot);
    head_yaw = stimSaccade(head_yaw, Stim(:,I.vel(kk)), false); % with approximate pattern position
    SACCADE.head_yaw(kk) = {head_yaw};
    SACCADE.head_roll(kk) = {head.roll};
    
    figure (13) ; cla ; hold on
        plot(tintrp, head.yaw_filt, 'Color', [0 0.1 0.6], 'LineWidth', 1)
        plot(tintrp, head.roll_filt, 'Color', [0.6 0.0 0.1], 'LineWidth', 1)
        pause
    
    if head_yaw.count == 0
        rep = 1;
    else
        HEAD_DATA{I.fly(kk),I.vel(kk)}(end+1,1) = head_yaw;
        rep = head_yaw.count;
    end
    VTable = table(D.vel(kk),'VariableNames',{'Vel'});
    ITable = [I(kk,:),VTable];
    ITable = repmat(ITable,rep,1);
    HEAD_SACCADE_STATS = [HEAD_SACCADE_STATS ; [ITable , head_yaw.SACD]];
    
    if head.showplot
        figure (1)
        pause
        close all
    end
    
    if head_yaw.count > 0
        % Yaw
    	[yaw_scd_pos,~,~,~] = getSaccade(head_yaw, head_yaw.position, scd_win, scd_win, true);
        SACCADE.yaw_scd_pos(kk) = {cat(2,yaw_scd_pos{:})};
        [yaw_scd_vel,~,~,~] = getSaccade(head_yaw, head_yaw.velocity, scd_win, scd_win, true);
        SACCADE.yaw_scd_vel(kk) = {cat(2,yaw_scd_vel{:})};
        
        % Roll
    	[roll_scd_pos,~,~,~] = getSaccade(head_yaw, head.roll_filt, scd_win, scd_win, true);
        SACCADE.roll_scd_pos(kk) = {cat(2,roll_scd_pos{:})};
        [roll_scd_vel,~,scd_time,~] = getSaccade(head_yaw, head.roll_vel_filt, scd_win, scd_win, true);
        SACCADE.roll_scd_vel(kk) = {cat(2,roll_scd_vel{:})};
        
        % Time
        SACCADE.scd_time(kk) = scd_time(1);
    end
end

% Fill in empty saccade trials
empty_idx = cellfun(@(x) isempty(x), HEAD_DATA);
HEAD_DATA(empty_idx) = {saccade_all(0*head.yaw, tintrp, head.thresh, head.true_thresh, head.Fc_detect, ...
                                head.Fc_ss, head.amp_cut, head.dur_cut, direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, false)};

%% Extract & group saccades & intervals by speed & by fly
fields = {'normpeak_saccade','norm_interval','normstart_interval','normend_interval',...
    'normstart_stimulus','error','int_error'};
nfield = length(fields);
norm_fields = {'time','position','velocity'};
norm_fields_stats = string(norm_fields) + "_stats";

center = 0; % normalization center for saccades & inter-saccade intervals
dim = 1; % dimension to center
dim_stats = 2; % dimension for statistics

FLY = []; % all trials per speed per fly
GRAND = []; % all trials per speed
for kk = 1:nfield % for each field in saccade structure
    if kk < 2 % even padding around center only for saccades (not intervals)
        even = true;
    else
        even = false;
    end
    
    FLY.(fields{kk}) = cellfun(@(x) struct_center([x.(fields{kk})], center, even, dim, norm_fields), ...
                            HEAD_DATA, 'UniformOutput', true);
    for jj = 1:N.vel
        GRAND.(fields{kk})(jj) = struct_center(FLY.(fields{kk})(:,jj), center, even , dim, norm_fields);
        
        % Calculate stats by fly
        for ff = 1:length(norm_fields)
            for ww = 1:N.fly
                FLY.(fields{kk})(ww,jj).(norm_fields_stats(ff)) = basic_stats( ...
                    FLY.(fields{kk})(ww,jj).(norm_fields{ff}), dim_stats);
            end
        end
    end
end

% Calculate stats for all
for kk = 1:nfield % for each field in saccade structure
    for jj = 1:N.vel
        for ff = 1:length(norm_fields)
            GRAND.(fields{kk})(jj).(norm_fields_stats(ff)) = basic_stats( ...
                GRAND.(fields{kk})(jj).(norm_fields{ff}), dim_stats);
        end
    end
end

%% SAVE %%
disp('Saving...')
save(['H:\DATA\Rigid_Data\Saccade\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','HEAD_DATA','SACCADE','HEAD_SACCADE_STATS','FLY','GRAND',...
      'Stim','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end