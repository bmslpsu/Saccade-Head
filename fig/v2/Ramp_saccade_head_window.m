function [] = Ramp_saccade_head_window()
%% Ramp_saccade_head_window:
root = 'H:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'U','N','SACCADE')

%% Get windows around saccades
clc
clearvars -except SACCADE U N
keepI = cellfun(@(x) length(x) > 1, SACCADE.head_scd_pos);
Saccade = SACCADE(keepI,:);
n_speed = N.vel / 2;
velI = Saccade.vel;
velI(velI > n_speed) = velI(velI > n_speed) - n_speed;
[vel_fly_group, vel_group, fly_group] = findgroups(velI, Saccade.fly);

Head.time   = splitapply(@(x) {cat(2,x{:})}, Saccade.scd_time, vel_fly_group);
Head.pos    = splitapply(@(x) {cat(2,x{:})}, Saccade.head_scd_pos, vel_fly_group);
Head.vel    = splitapply(@(x) {cat(2,x{:})}, Saccade.head_scd_vel, vel_fly_group);
Head.accel  = splitapply(@(x) {cat(2,x{:})}, Saccade.head_scd_accel, vel_fly_group);

Head = structfun(@(x) splitapply(@(y) {y}, x, vel_group), Head, 'UniformOutput', false);
Head = structfun(@(x) cat(2,x{:}), Head, 'UniformOutput', false);
% Head.pos = cellfun(@(x) x - mean(x), Head.pos, 'UniformOutput', false);

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

%% Plot by speed
clc
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 6 3.5])
movegui(fig, 'center')
clear ax H
cc_vel = hsv(n_speed);
ax = gobjects(3,n_speed);
fly_lw = 0.5;
vel_lw = 1;
for v = 1:n_speed
    rowI = v + (0:2)*n_speed;
    ax(1,v) = subplot(3,n_speed,rowI(1)); hold on ; cla
        plot(Head.fly_mean.time{v}, Head.fly_mean.pos{v}, ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', fly_lw)
        plot(Head.vel_stats.time(v).mean, Head.vel_stats.pos(v).mean, ...
            'Color', cc_vel(v,:), 'LineWidth', vel_lw)
        
    ax(2,v) = subplot(3,n_speed,rowI(2)); hold on ; cla
        plot(Head.fly_mean.time{v}, Head.fly_mean.vel{v}, ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', fly_lw)
        plot(Head.vel_stats.time(v).mean, Head.vel_stats.vel(v).mean, ...
            'Color', cc_vel(v,:), 'LineWidth', vel_lw)
        
    ax(3,v) = subplot(3,n_speed,rowI(3)); hold on ; cla
        plot(Head.fly_mean.time{v}, Head.fly_mean.accel{v}, ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', fly_lw)
        plot(Head.vel_stats.time(v).mean, Head.vel_stats.accel(v).mean, ...
            'Color', cc_vel(v,:), 'LineWidth', vel_lw)     
end
linkaxes(ax, 'x')
linkaxes(ax(1,:), 'y')
linkaxes(ax(2,:), 'y')
linkaxes(ax(3,:), 'y')
set(ax, 'XLim', 0.2*[-1 1])
set(ax(1,:), 'YLim', 13*[-1 1])
set(ax(2,:), 'YLim', [-200 700])
set(ax(3,:), 'YLim', 40000*[-1 1])
set(ax(1:2,:), 'XTickLabel', [])
set(ax(1:2,:), 'XColor', 'none')
set(ax(:,2:n_speed), 'YTickLabel', [])
set(ax(:,2:n_speed), 'YColor', 'none')
set(ax(3,2:n_speed), 'XColor', 'none')
set(ax, 'LineWidth', 1)

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Position (°)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Velocity (°/s)')
YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Accleration (°/s^{2})')
XLabelHC = get(ax(3,1), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')

%% Mean
clc
max_accel = nan(n_speed,2);
max_accel_std = nan(n_speed,2);
for v = 1:n_speed
   [max_accel(v,1),maxI] = max(Head.vel_stats.vel(v).mean);
   [max_accel(v,2),minI] = min(Head.vel_stats.vel(v).mean);
   max_accel_std(v,1) = Head.vel_stats.vel(v).std(maxI);
   max_accel_std(v,2) = Head.vel_stats.vel(v).std(minI);
end
disp(max_accel)
disp(['Mean Acceleration: ' num2str(mean(max_accel,1))])
disp(['STD Acceleration: ' num2str(mean(max_accel_std,1))])

%% Save
savedir = 'H:\DATA\Rigid_Data\Saccade\processed';
filename = ['Ramp_saccade_window_wave=' num2str(U.wave) '.mat'];
save(fullfile(savedir, filename), 'Head','U','N','-v7.3')

end