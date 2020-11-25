function [] = Ramp_saccade_head_wing_window()
%% Ramp_saccade_head_wing_window:
root = 'H:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'U','N','SACCADE')

%% Get windows around saccades
clc
clearvars -except SACCADE U N FILE
keepI = cellfun(@(x) isobject(x) || isstruct(x), SACCADE.head2wing);
Saccade = SACCADE(keepI,:);
if N.vel > 1
    n_speed = N.vel / 2;
else
    n_speed = N.vel;
end
velI = Saccade.vel;
velI(velI > n_speed) = velI(velI > n_speed) - n_speed;
[vel_fly_group, vel_group, fly_group] = findgroups(velI, Saccade.fly);
n_speed = length(unique(vel_group));

Scd.time        = splitapply(@(x) {cat(2,x{:})}, Saccade.scd_time, vel_fly_group);
Scd.head_pos    = splitapply(@(x) {cat(2,x{:})}, Saccade.head_scd_pos, vel_fly_group);
Scd.head_vel    = splitapply(@(x) {cat(2,x{:})}, Saccade.head_scd_vel, vel_fly_group);
Scd.head_accel  = splitapply(@(x) {cat(2,x{:})}, Saccade.head_scd_accel, vel_fly_group);
Scd.wing_pos    = splitapply(@(x) {cat(2,x{:})}, Saccade.wing_scd_pos, vel_fly_group);
Scd.wing_vel    = splitapply(@(x) {cat(2,x{:})}, Saccade.wing_scd_vel, vel_fly_group);

% head_pos_sync = cellfun(@(x) x.hsacd.pos, Saccade.head2wing, 'UniformOutput', false);
% head_vel_sync = cellfun(@(x) x.hsacd.vel, Saccade.head2wing, 'UniformOutput', false);
% head_pos_desync = cellfun(@(x) x.hsacd.pos_desync, Saccade.head2wing, 'UniformOutput', false);
% head_vel_desync = cellfun(@(x) x.hsacd.vel_desync, Saccade.head2wing, 'UniformOutput', false);
% 
% wing_pos_sync = cellfun(@(x) x.wsacd.pos, Saccade.head2wing, 'UniformOutput', false);
% wing_vel_sync = cellfun(@(x) x.wsacd.vel, Saccade.head2wing, 'UniformOutput', false);
% wing_pos_desync = cellfun(@(x) x.wsacd.pos_desync, Saccade.head2wing, 'UniformOutput', false);
% wing_vel_desync = cellfun(@(x) x.wsacd.vel_desync, Saccade.head2wing, 'UniformOutput', false);
% 
% Scd.head_pos_sync = splitapply(@(x) {cat(2,x{:})}, head_pos_sync, vel_fly_group);
% Scd.head_vel_sync = splitapply(@(x) {cat(2,x{:})}, head_vel_sync, vel_fly_group);
% Scd.head_pos_desync = splitapply(@(x) {cat(2,x{:})}, head_pos_desync, vel_fly_group);
% Scd.head_vel_desync = splitapply(@(x) {cat(2,x{:})}, head_vel_desync, vel_fly_group);
% 
% Scd.wing_pos_sync = splitapply(@(x) {cat(2,x{:})}, wing_pos_sync, vel_fly_group);
% Scd.wing_vel_sync = splitapply(@(x) {cat(2,x{:})}, wing_vel_sync, vel_fly_group);
% Scd.wing_pos_desync = splitapply(@(x) {cat(2,x{:})}, wing_pos_desync, vel_fly_group);
% Scd.wing_vel_desync = splitapply(@(x) {cat(2,x{:})}, wing_vel_desync, vel_fly_group);

Scd = structfun(@(x) splitapply(@(y) {y}, x, vel_group), Scd, 'UniformOutput', false);
Scd = structfun(@(x) cat(2,x{:}), Scd, 'UniformOutput', false);

% Normalize position
Fs = round(Saccade.head_saccade{1}.Fs);
winI = 1 * Fs;
% Scd.wing_pos = cellfun(@(x) x - mean(x([1:winI, length(x)-winI:length(x)],:),1), ...
%     Scd.wing_pos, 'UniformOutput', false);
Scd.wing_pos = cellfun(@(x) x - mean(x(1:winI,:),1), ...
    Scd.wing_pos, 'UniformOutput', false);

Scd.fly_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    Scd, 'UniformOutput', false);

fnames = string(fieldnames(Scd.fly_stats));
n_field = length(fnames);
for f = 1:n_field
    for v = 1:n_speed
        Scd.fly_mean.(fnames(f)){v} = cat(2, Scd.fly_stats.(fnames(f))(:,v).mean);
    end
end

Scd.vel_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    Scd.fly_mean, 'UniformOutput', false);

%% Plot ALL by speed
clc
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 3 4.5])
movegui(fig, 'center')
clear ax H
ax = gobjects(3,n_speed);
fly_lw = 0.5;
vel_lw = 1;
head_color = [0 0 1];
wing_color = [1 0 0];
n_std = 1;
for v = 1:n_speed
    rowI = v + (0:2)*n_speed;
    ax(1,v) = subplot(3,n_speed,rowI(1)); hold on ; cla
%         plot(Scd.fly_mean.time{v}, Scd.fly_mean.head_pos{v}, ...
%             'Color', [0.7*head_color 0.5], 'LineWidth', fly_lw)
        
%         plot(Scd.fly_mean.time{v}, Scd.fly_mean.wing_pos{v}, ...
%             'Color', [0.7*wing_color 0.5], 'LineWidth', fly_lw)

        [~] = PlotPatch(Scd.vel_stats.head_pos(v).mean, Scd.vel_stats.head_pos(v).std, ...
            Scd.vel_stats.time(v).mean, n_std, 1, head_color, head_color, 0.3, vel_lw);
        
%         [~] = PlotPatch(Scd.vel_stats.wing_pos(v).mean, Scd.vel_stats.wing_pos(v).std, ...
%             Scd.vel_stats.time(v).mean, n_std, 1, wing_color, wing_color, 0.3, vel_lw);
        
    ax(2,v) = subplot(3,n_speed,rowI(2)); hold on ; cla
%         plot(Scd.fly_mean.time{v}, Scd.fly_mean.head_vel{v}, ...
%             'Color', [0.7*head_color 0.5], 'LineWidth', fly_lw)
%         
%         plot(Scd.fly_mean.time{v}, Scd.fly_mean.wing_vel{v}, ...
%             'Color', [0.7*wing_color 0.5], 'LineWidth', fly_lw)
        
        [~] = PlotPatch(Scd.vel_stats.head_vel(v).mean, Scd.vel_stats.head_vel(v).std, ...
            Scd.vel_stats.time(v).mean, n_std, 1, head_color, head_color, 0.3, vel_lw);
        
%         [~] = PlotPatch(Scd.vel_stats.wing_vel(v).mean, Scd.vel_stats.wing_vel(v).std, ...
%             Scd.vel_stats.time(v).mean, n_std, 1, wing_color, wing_color, 0.3, vel_lw);
 
    ax(3,v) = subplot(3,n_speed,rowI(3)); hold on ; cla
        plot(Scd.fly_mean.time{v}, Scd.fly_mean.head_accel{v}, ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', fly_lw)       
        
        [~] = PlotPatch(Scd.vel_stats.head_accel(v).mean, Scd.vel_stats.head_accel(v).std, ...
            Scd.vel_stats.time(v).mean, n_std, 1, head_color, head_color, 0.3, vel_lw);
end
linkaxes(ax, 'x')
linkaxes(ax(1,:), 'y')
linkaxes(ax(2,:), 'y')
linkaxes(ax(3,:), 'y')
set(ax, 'XLim', 0.2*[-1 1])
set(ax(1,:), 'YLim', 15*[-1 1])
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
set([YLabelHC], 'String', 'Accleratin (°/s^{2})')
XLabelHC = get(ax(3,1), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')

%% Save
savedir = 'H:\DATA\Rigid_Data\Saccade\processed';
filedata = textscan(char(FILE), '%s', 'delimiter', '_');
fclass = filedata{1}{1};
filename = [fclass '_saccade_window_head_wing.mat'];
save(fullfile(savedir, filename), 'Scd','U','N','-v7.3')

end