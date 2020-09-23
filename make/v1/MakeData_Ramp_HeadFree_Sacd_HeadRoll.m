function [] = MakeData_Ramp_HeadFree_Sacd_HeadRoll()
%% MakeData_Ramp_HeadFree_Sacd_HeadRoll:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%

% Data location
rootdir = 'H:\EXPERIMENTS\RIGID\Experiment_Ramp_forRoll';

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
head_carry.Fc = 60;
[head_carry.b, head_carry.a] = butter(3, head_carry.Fc / (Fs/2) ,'low');

Vel = U.vel{1}; % velocities
Stim = (Vel*tintrp')'; % stimuli
SACCADE = [I , splitvars(table(num2cell(zeros(N.file,4))))];
SACCADE.Properties.VariableNames(4:7) = {'head_yaw','head_roll','yaw2roll','roll2yaw'};
HEAD_DATA = cell(N.fly,N.vel);
HEAD_SACCADE_STATS = [];
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p');
    load(fullfile(PATH.head, [basename{kk} '.mat']),'yaw','roll_idx');
 
  	% Sync frames and get pattern data
	[trig,pat] = sync_pattern_trigger(t_p, data(:,2), 10, data(:,1), true, [], false, false);
    
    % Get head data
    head.yaw = hampel(trig.time_sync, yaw);
    head.roll = hampel(trig.time_sync, roll_idx);
    
    head.yaw = interp1(trig.time_sync, head.yaw , tintrp, 'pchip');
    head.roll = interp1(trig.time_sync, 36.33*head.roll, tintrp, 'pchip');
    
	% Head with filter
 	head.yaw_filt = filtfilt(head_carry.b, head_carry.a, head.yaw);
    head.roll_filt = filtfilt(head_carry.b, head_carry.a, head.roll);
       
    % Extract head yaw saccades
    direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    head_yaw = saccade_all(head.yaw_filt, tintrp, head.thresh, head.true_thresh, head.Fc_detect, ...
                                head.Fc_ss, head.amp_cut, direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, head.showplot);
    head_yaw = stimSaccade(head_yaw, Stim(:,I.vel(kk)), false); % with approximate pattern position
    SACCADE{kk,4} = {head_yaw}; % store data in cell
    
    
    % Extract head roll saccades
    direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    head_roll = saccade_all(head.roll_filt, tintrp, head.thresh, head.true_thresh, head.Fc_detect, ...
                                head.Fc_ss, head.amp_cut, direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, true);
    head_roll = stimSaccade(head_roll, Stim(:,I.vel(kk)), false); % with approximate pattern position
    SACCADE{kk,5} = {head_roll}; % store data in cell
    pause
    close all
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
        yaw2roll = saccade_interact(head_yaw, head_roll, 0.5, 0.1, true);
        SACCADE{kk,6} = {yaw2roll};
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
                            
%% Interact
name = 'yaw2roll';
keepI = cellfun(@(x) isstruct(x) | isobject(x), SACCADE.(name));
Saccade = SACCADE(keepI,:);
Saccade = Saccade(:,:);

velI = Saccade.vel;
cw = velI <= 1;
clear Yaw Roll Time

Time.align = cellfun(@(x) x.interval.time_align, Saccade.(name), 'UniformOutput', false);
Time.align = nanmean(cat(2,Time.align{:}),2);

Yaw.pos = cellfun(@(x) x.interval.in.pos , Saccade.(name), 'UniformOutput', false);
Yaw.vel = cellfun(@(x) x.interval.in.vel , Saccade.(name), 'UniformOutput', false);
Yaw.pos(cw) = cellfun(@(x) -x , Yaw.pos(cw), 'UniformOutput', false);
Yaw.vel(cw) = cellfun(@(x) -x , Yaw.vel(cw), 'UniformOutput', false);
Yaw = structfun(@(x) cat(2,x{:}), Yaw, 'UniformOutput', false);
Yaw.stats = structfun(@(x) basic_stats(x,2), Yaw, 'UniformOutput', false);

Roll.pos = cellfun(@(x) x.interval_rel.out.pos , Saccade.(name), 'UniformOutput', false);
Roll.vel = cellfun(@(x) x.interval_rel.out.vel , Saccade.(name), 'UniformOutput', false);
Roll.pos(cw) = cellfun(@(x) -x , Roll.pos(cw), 'UniformOutput', false);
Roll.vel(cw) = cellfun(@(x) -x , Roll.vel(cw), 'UniformOutput', false);
Roll = structfun(@(x) cat(2,x{:}), Roll, 'UniformOutput', false);
Roll.stats = structfun(@(x) basic_stats(x,2), Roll, 'UniformOutput', false);

fig = figure (102); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 5])
movegui(fig, 'center')
clear ax
yaw_color = [0.7 0 0];
roll_color = [0.1 1 0.5];
alpha = 0.3;
ax(1) = subplot(2,1,1); cla ; hold on ; ylabel('(°)')
    plot(Time.align, Roll.pos, 'Color', [0.7*roll_color alpha], 'LineWidth', 0.25)
    plot(Time.align, Yaw.pos, 'Color', [0.7*yaw_color alpha], 'LineWidth', 0.25)
    [~,h.mean(1)] = PlotPatch(Roll.stats.pos.mean, Roll.stats.pos.std, ...
        Time.align, 1, 1, roll_color, 0.7*roll_color, 0.3, 1);
    [~,h.mean(2)] = PlotPatch(Yaw.stats.pos.mean, Yaw.stats.pos.std, ...
        Time.align, 1, 1, yaw_color, 0.7*yaw_color, 0.3, 1);
    
ax(2) = subplot(2,1,2); cla ; hold on ; ylabel('(°/s)')
    %plot(Time.align, Roll.vel, 'Color', [0.7*roll_color alpha], 'LineWidth', 0.25)
    %plot(Time.align, Yaw.vel, 'Color', [0.7*yaw_color alpha], 'LineWidth', 0.25)
    [~,~] = PlotPatch(Roll.stats.vel.mean, Roll.stats.vel.std, ...
        Time.align, 1, 1, roll_color, 0.7*roll_color, 0.3, 1);
    [~,~] = PlotPatch(Yaw.stats.vel.mean, Yaw.stats.vel.std, ...
        Time.align, 1, 1, yaw_color, 0.7*yaw_color, 0.3, 1);
    
    xlabel('Time (s)')
    legend(h.mean, 'Roll', 'Yaw', 'Box', 'off')
    
set(ax, 'LineWidth', 1, 'XLim', 0.2*[-1 1], 'XTick', -0.2:0.05:0.2, 'Box', 'on')
set(ax(1), 'YLim', 40*[-1 1])
set(ax(2), 'YLim', [-600 1600])
linkaxes(ax,'x')

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