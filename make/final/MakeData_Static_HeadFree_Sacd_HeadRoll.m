function [] = MakeData_Static_HeadFree_Sacd_HeadRoll()
%% MakeData_Static_HeadFree_Sacd_HeadRoll:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%

% Data location
rootdir = 'H:\EXPERIMENTS\RIGID\Experiment_Static_forRoll';

% Output file name
% filename = ['Ramp_HeadFree_SACCD_Anti_HeadWingALL_Wave=' num2str(wave)];

% Setup Directories 
PATH.daq  = rootdir; % DAQ data location
PATH.head = fullfile(PATH.daq,'tracked_head_edge'); % tracked kinematic data location
% PATH.wing = fullfile(PATH.vid,'wing_filt', 'tracked_head_wing'); % tracked kinematic data location

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.head,'*.mat',false,'fly','trial','vel');
% [D,I,N,U,T,~,~,basename] = GetFileData(PATH.wing,'*.csv',false,'fly','trial','vel');

%% Get Data %%
clc
disp('Loading...')
Fs = 200; % sampling frequency [s]
tintrp = (0:(1/Fs):(10))'; % time vector for interpolation

% HEAD saccade detection parameters
head.showplot = false;
head.Fc_detect = [10 nan];
head.Fc_ss = [40 nan];
head.amp_cut = 4;
head.thresh = [80];
head.true_thresh = 250;
head.sacd_length = nan;
head.pks = [];
head.min_pkdist = 0.1;
head.min_pkwidth = 0.03;
head.min_pkprom = 75;
head.min_pkthresh = 0;
head.boundThresh = [0.25 50];
head_carry.Fc = 40;
[head_carry.b, head_carry.a] = butter(3, head_carry.Fc / (Fs/2) ,'low');

SACCADE = [I , splitvars(table(num2cell(zeros(N.file,4))))];
SACCADE.Properties.VariableNames(4:7) = {'head_yaw','head_roll','yaw2roll','roll2yaw'};
HEAD_DATA = cell(N.fly,N.vel);
HEAD_SACCADE_STATS = [];
Yaw_SCD = cell(N.fly,N.vel);
Roll_SCD = cell(N.fly,N.vel);
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p');
    load(fullfile(PATH.head, [basename{kk} '.mat']),'yaw','roll_idx');
 
  	% Sync frames and get pattern data
	[trig,pat] = sync_pattern_trigger(t_p, data(:,2), 10, data(:,1), true, round(0.5*Fs), false, false);
    
    % Get head data
    dx = 3;
    head.yaw = hampel(trig.time_sync, yaw, dx);
    head.roll = hampel(trig.time_sync, roll_idx, dx);
    
    head.yaw = interp1(trig.time_sync, head.yaw , tintrp, 'pchip');
    head.roll = interp1(trig.time_sync, 70*head.roll, tintrp, 'pchip');
    
    head.yaw_vel = diff(head.yaw) * Fs; head.yaw_vel = [head.yaw_vel(1) ; head.yaw_vel];
    head.roll_vel = diff(head.roll) * Fs; head.roll_vel = [head.roll_vel(1) ; head.roll_vel];
    
	% Head with filter
 	head.yaw_filt = filtfilt(head_carry.b, head_carry.a, head.yaw);
    head.roll_filt = filtfilt(head_carry.b, head_carry.a, head.roll);
    
 	head.yaw_vel_filt = diff(head.yaw_filt) * Fs;
    head.yaw_vel_filt = [head.yaw_vel_filt(1) ; head.yaw_vel_filt];
    head.yaw_vel_filt = filtfilt(head_carry.b, head_carry.a, head.yaw_vel_filt);
    
 	head.roll_vel_filt = diff(head.roll_filt) * Fs;
    head.roll_vel_filt = [head.roll_vel_filt(1) ; head.roll_vel_filt];
    head.roll_vel_filt = filtfilt(head_carry.b, head_carry.a, head.roll_vel_filt);
        
%     figure (13) ; cla ; hold on
%         plot(tintrp, head.yaw_filt, 'Color', [0 0.1 0.6], 'LineWidth', 1)
%         plot(tintrp, head.roll_filt, 'Color', [0.6 0.0 0.1], 'LineWidth', 1)
%         pause
       
    % Extract head yaw saccades
    direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    head_yaw = saccade_all(head.yaw_filt, tintrp, head.thresh, head.true_thresh, head.Fc_detect, ...
                                head.Fc_ss, head.amp_cut, direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, head.showplot);
    SACCADE.head_yaw(kk) = {head_yaw};
    SACCADE.head_roll(kk) = {head.roll};
    
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
        [yaw_scds,yaw_ints,~,~] = getSaccade(head_yaw, head.yaw_filt, 0.5, 0.5, true);
        [roll_scds,roll_ints,scd_time,int_time] = getSaccade(head_yaw, head.roll_filt, 0.5, 0.5, true);
        
        Yaw_SCD{I.fly(kk,:),I.vel(kk,:)} = cat(1, Yaw_SCD{I.fly(kk,:),I.vel(kk,:)}, yaw_scds);
        Roll_SCD{I.fly(kk,:),I.vel(kk,:)} = cat(1, Roll_SCD{I.fly(kk,:),I.vel(kk,:)}, roll_scds);
        %yaw2roll = saccade_interact(head_yaw, head_roll, 0.5, 0.1, true);
        %SACCADE{kk,6} = {yaw2roll};
        %SACCADE{kk,7} = {roll2yaw};
        %pause
    end
end

% Fill in empty saccade trials
empty_idx = cellfun(@(x) isempty(x), HEAD_DATA);
HEAD_DATA(empty_idx) = {saccade_all(0*head.yaw, tintrp, head.thresh, head.true_thresh, head.Fc_detect, ...
                                head.Fc_ss, head.amp_cut, direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, false)};
                            
%% Yaw vs Roll
% tt = int_time{1};
tt = scd_time{1};
Yaw.All = cellfun(@(x) [x{:}], Yaw_SCD, 'UniformOutput', false);
Roll.All = cellfun(@(x) [x{:}], Roll_SCD, 'UniformOutput', false);
Yaw.All = cellfun(@(x) x(:,~isnan(x(1,:))), Yaw.All, 'UniformOutput', false);
Roll.All = cellfun(@(x) x(:,~isnan(x(1,:))), Roll.All, 'UniformOutput', false);
for v = 1:N.vel
    for f = 1:N.fly
       for c = 1:size(Roll.All{f,v},2)
           Yaw.All{f,v}(:,c) = Yaw.All{f,v}(:,c) - 0*nanmean(Yaw.All{f,v}([1:20],c));
           Roll.All{f,v}(:,c) = Roll.All{f,v}(:,c) - 0*nanmean(Roll.All{f,v}([1:50],c));
       end
    end
end

Yaw.Fly = cellfun(@(x) nanmean(x,2), Yaw.All, 'UniformOutput', false);
Roll.Fly = cellfun(@(x) nanmean(x,2), Roll.All, 'UniformOutput', false);
Yaw.Fly = cat(2, Yaw.Fly{:});
Roll.Fly = cat(2, Roll.Fly{:});

Yaw.Fly_stats = basic_stats(Yaw.Fly, 2);
Roll.Fly_stats = basic_stats(Roll.Fly, 2);

fig = figure (102);
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 4 3])
movegui(fig, 'center')
clear ax

yaw_color = [0 0 1];
roll_color = [0.1 0.7 0.3];
ax(1) = subplot(1,1,1); cla ; hold on
%     plot(tt, Yaw.Fly, 'Color', 0.7*yaw_color, 'LineWidth', 0.5)
%     plot(tt, Roll.Fly, 'Color', 0.7*roll_color, 'LineWidth', 0.5)
    yline(0);
    xline(0);
    
    [~,h.mean(1)] = PlotPatch(Roll.Fly_stats(1).mean, Roll.Fly_stats(1).std, ...
        tt, 1, 1, roll_color, 0.7*roll_color, 0.3, 1);
    [~,h.mean(2)] = PlotPatch(Yaw.Fly_stats(1).mean, Yaw.Fly_stats(1).std, ...
        tt, 1, 1, yaw_color, 0.7*yaw_color, 0.3, 1);
    
    uistack(h.mean, 'top')
    
    set(ax, 'LineWidth', 1, 'Box', 'off')
    set(ax, 'XLim', 0.5*[-1 1], 'XTick', -0.5:0.1:0.5)
    xlabel('Time (s)')
    ylabel('Head Angle (°)')
%     ylim([-400 600])
    ylim(10*[-1 1])
%     xlim([0 0.5])
    
%%
clc
start_time_I = find(tt == 0);
end_time_I = find(tt == 0.1);

yaw_int_all = cat(2, Yaw.All_speed{:});
yaw_int_all = yaw_int_all(start_time_I:end_time_I,:);
yaw_int_all = yaw_int_all(~isnan(yaw_int_all));

roll_int_all = cat(2, Roll.All_speed{:});
roll_int_all = roll_int_all(start_time_I:end_time_I,:);
roll_int_all = roll_int_all(~isnan(roll_int_all));

[rho,pval] = corr(yaw_int_all, roll_int_all)

[rho,pval] = corr(Yaw.Vel_stats(1).mean(start_time_I:end_time_I), ...
    Roll.Vel_stats(1).mean(start_time_I:end_time_I))

%% Example Trial: fly 2 trial 2 I=735

idx = 2;

tt = SACCADE.head_yaw{idx}.time;
yaw = SACCADE.head_yaw{idx}.position;
roll = SACCADE.head_roll{idx};

[b,a] = butter(3, 20 / (Fs/2), 'low');
yaw = filtfilt(b, a, yaw);
roll = filtfilt(b, a, roll);

scd_frame = 735;
scd_time = tt(scd_frame);
scd_yaw = yaw(scd_frame);
scd_roll = roll(scd_frame);

fig = figure (102);
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 2.5])
movegui(fig, 'center')
clear ax

yaw_color = [0 0 1];
roll_color = [0.1 0.7 0.3];
ax(1) = subplot(1,1,1); cla ; hold on
    plot(tt, yaw, 'Color', yaw_color, 'LineWidth', 1)
    plot(tt, roll, 'Color', roll_color, 'LineWidth', 1)
    
    plot(scd_time, scd_yaw, 'r.', 'MarkerSize', 15, 'MarkerFaceColor', 'none')
    plot(scd_time, scd_roll, 'r.', 'MarkerSize', 15, 'MarkerFaceColor', 'none')
    
    xlim([0 10])
    ylim(30*[-1 1])
    set(ax, 'LineWidth', 1, 'Box', 'off')
    xlabel('Time (s)')
    ylabel('Head (°)')

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
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','HEAD_DATA','SACCADE','HEAD_SACCADE_STATS','WING_SACCADE_STATS','FLY','GRAND',...
      'Stim','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end