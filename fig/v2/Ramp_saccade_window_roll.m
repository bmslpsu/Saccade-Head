function [] = Ramp_saccade_window_roll()
%% Ramp_saccade_window_roll:
root = 'H:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'U','N','SACCADE')

%% Get windows around saccades
clc
clearvars -except SACCADE U N
keepI = cellfun(@(x) length(x) > 1, SACCADE.yaw_scd_pos);
Saccade = SACCADE(keepI,:);
if N.vel > 1
    n_speed = N.vel / 2;
else
    n_speed = N.vel;
end
velI = Saccade.vel;
velI(velI > n_speed) = velI(velI > n_speed) - n_speed;

[vel_fly_group, vel_group, fly_group] = findgroups(velI, Saccade.fly);

Head.time       = splitapply(@(x) {cat(2,x{:})}, Saccade.scd_time, vel_fly_group);
Head.yaw_pos    = splitapply(@(x) {cat(2,x{:})}, Saccade.yaw_scd_pos, vel_fly_group);
Head.yaw_vel    = splitapply(@(x) {cat(2,x{:})}, Saccade.yaw_scd_vel, vel_fly_group);
Head.roll_pos 	= splitapply(@(x) {cat(2,x{:})}, Saccade.roll_scd_pos, vel_fly_group);
Head.roll_vel 	= splitapply(@(x) {cat(2,x{:})}, Saccade.roll_scd_vel, vel_fly_group);

Head = structfun(@(x) splitapply(@(y) {y}, x, vel_group), Head, 'UniformOutput', false);
Head = structfun(@(x) cat(2,x{:}), Head, 'UniformOutput', false);
% Head.yaw_pos = cellfun(@(x) x - mean(x), Head.yaw_pos, 'UniformOutput', false);
Head.roll_pos = cellfun(@(x) x - mean(x), Head.roll_pos, 'UniformOutput', false);

Head.fly_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    Head, 'UniformOutput', false);

fnames = string(fieldnames(Head.fly_stats));
n_field = length(fnames);
for f = 1:n_field
    for v = 1:n_speed
        Head.fly_mean.(fnames(f)){v} = cat(2, Head.fly_stats.(fnames(f))(:,v).mean);
    end
end

Head.vel_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    Head.fly_mean, 'UniformOutput', false);

%% Plot Yaw vs Roll
clc
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', 2*[2 2 6/5 (2/3)*3.5])
movegui(fig, 'center')
clear ax H
yaw_color = [0 0 1];
roll_color = [0.3 0.2 0.6];
ax = gobjects(2,n_speed);
fly_lw = 0.5;
Lw = 1;
n_std = 1;
for v = 1:n_speed
    rowI = v + (0:1)*n_speed;
    ax(1,v) = subplot(2,n_speed,rowI(1)); hold on ; cla
%         plot(Head.fly_mean.time{v}, Head.fly_mean.yaw_pos{v}, ...
%             'Color', [0.7*yaw_cc 0.5], 'LineWidth', fly_lw)
        [~] = PlotPatch(Head.vel_stats.yaw_pos(v).mean, Head.vel_stats.yaw_pos(v).std, ...
            Head.vel_stats.time(v).mean, n_std, 1, yaw_color, yaw_color, 0.3, Lw);
        
%         plot(Head.fly_mean.time{v}, Head.fly_mean.roll_pos{v}, ...
%             'Color', [0.7*roll_cc 0.5], 'LineWidth', fly_lw)
        [~] = PlotPatch(Head.vel_stats.roll_pos(v).mean, Head.vel_stats.roll_pos(v).std, ...
            Head.vel_stats.time(v).mean, n_std, 1, roll_color, roll_color, 0.3, Lw);
        
    ax(2,v) = subplot(2,n_speed,rowI(2)); hold on ; cla
%         plot(Head.fly_mean.time{v}, Head.fly_mean.yaw_vel{v}, ...
%             'Color', [0.7*yaw_cc 0.5], 'LineWidth', fly_lw)
        [~] = PlotPatch(Head.vel_stats.yaw_vel(v).mean, Head.vel_stats.yaw_vel(v).std, ...
            Head.vel_stats.time(v).mean, n_std, 1, yaw_color, yaw_color, 0.3, Lw);
        
%         plot(Head.fly_mean.time{v}, Head.fly_mean.roll_vel{v}, ...
%             'Color', [0.7*roll_cc 0.5], 'LineWidth', fly_lw)
        [~] = PlotPatch(Head.vel_stats.roll_vel(v).mean, Head.vel_stats.roll_vel(v).std, ...
            Head.vel_stats.time(v).mean, n_std, 1, roll_color, roll_color, 0.3, Lw);
end
linkaxes(ax, 'x')
linkaxes(ax(1,:), 'y')
linkaxes(ax(2,:), 'y')
set(ax, 'XLim', 0.2*[-1 1])
set(ax(1,:), 'YLim', 10*[-1 1])
set(ax(2,:), 'YLim', [-400 600])
% set(ax(1:2,:), 'XTickLabel', [])
% set(ax(1:2,:), 'XColor', 'none')
% set(ax(:,2:n_speed), 'YTickLabel', [])
% set(ax(:,2:n_speed), 'YColor', 'none')
set(ax, 'LineWidth', 1)

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Position (°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Velocity (°/s)')
XLabelHC = get(ax(2,1), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')

%% Example Trial
idx = 2; % Ramp: fly 2 trial 2 I=735
% idx = 90; % Static: fly 

tt = SACCADE.head_yaw{idx}.time;
yaw = SACCADE.head_yaw{idx}.position;
roll = SACCADE.head_roll{idx};

[b,a] = butter(3, 20 / (SACCADE.head_yaw{idx}.Fs/2), 'low');
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

end