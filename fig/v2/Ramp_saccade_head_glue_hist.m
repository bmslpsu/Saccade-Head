function [] = Ramp_saccade_head_glue_hist()
%% Ramp_saccade_head_window:
root = 'E:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'U','N','SACCADE')

%% Histogram
clc
clearvars -except SACCADE U N
keepI = cellfun(@(x) length(x) > 1, SACCADE.head_scd_pos);
Saccade = SACCADE(keepI,:);
G = sign(Saccade.vel);
T = table(cellfun(@(x) x.position, Saccade.head_saccade, 'UniformOutput', false), ...
    'VariableNames', {'head_all'});
T = [Saccade , T];
head_all_glue = table_fly_stats(T, G, 15);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 4 2])
movegui(fig, 'center')
clear ax H
% hbins = -21.25:0.5:22.15;
hbins = -21.5:1:21.5;
H = gobjects(1,2);
cc = [1 0 0; 0 0 1];

ax = subplot(1,1,1) ; cla ; hold on
ylabel('Probability')
xlabel('Head (°)')
for k = 1:2
    Y = head_all_glue.comb.val_all.head_all{1,k}(:);
    Y = Y(~isnan(Y));
    H(k) = histogram(Y, hbins, 'Normalization', 'probability', 'EdgeColor', 'none', ...
        'FaceColor', cc(k,:));
end
xline(0, '--', 'Color', -0.5+[0.5 0.5 0.5])
ax.YLim(1) = -0.05*ax.YLim(2);

xx = [5 5 15 15];
yy = [ax.YLim(1) ax.YLim(2) ax.YLim(2) ax.YLim(1)];
patch(xx, yy, 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.2)

%% Histogram norm
clc
clearvars -except SACCADE U N
keepI = cellfun(@(x) length(x) > 1, SACCADE.head_scd_pos);
speedI = any(SACCADE.vel == [1 2 6 7],2);
Saccade = SACCADE(keepI & speedI,:);
% G = sign(Saccade.vel);
direction = Saccade.vel;
direction(any(direction==(1:5),2)) = 1;
direction(any(direction==(6:10),2)) = -1;
% T = table(cellfun(@(x) x.position, Saccade.head_saccade, 'UniformOutput', false), ...
%     'VariableNames', {'head_all'});
T = table(cellfun(@(x,y) y*x.position, Saccade.head_saccade, num2cell(direction), ...
    'UniformOutput', false), 'VariableNames', {'head_all'});
T = [Saccade , T];
G = ones(size(direction));
head_all_glue = table_fly_stats(T, G, 15);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 4 2])
movegui(fig, 'center')
clear ax H
% hbins = -21.25:0.5:22.15;
hbins = -23.5:1:23.5;
H = gobjects(1,2);
cc = [0 0 0];

ax = subplot(1,1,1) ; cla ; hold on
ylabel('Probability')
xlabel('Head (°)')
for k = 1
    Y = head_all_glue.comb.val_all.head_all{1,k}(:);
    Y = Y(~isnan(Y));
    H(k) = histogram(Y, hbins, 'Normalization', 'probability', 'EdgeColor', 'none', ...
        'FaceColor', cc(k,:));
end
xline(0, '--', 'Color', -0.5+[0.5 0.5 0.5])
ylim([0 0.1])
ax.YLim(1) = -0.05*ax.YLim(2);

xx = [8 8 20 20];
yy = [ax.YLim(1) ax.YLim(2) ax.YLim(2) ax.YLim(1)];
patch(xx, yy, 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.2)


%% Get windows around saccades
clc
clearvars -except SACCADE U N
keepI = cellfun(@(x) length(x) > 1, SACCADE.head_scd_pos);
Saccade = SACCADE(keepI,:);
n_speed = N.vel / 2;
velI = findgroups(Saccade.vel);
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


%% Save
% savedir = 'E:\DATA\Rigid_Data\Saccade\processed';
% filename = ['Ramp_saccade_window_wave=' num2str(U.wave) '.mat'];
% save(fullfile(savedir, filename), 'Head','U','N','-v7.3')

end