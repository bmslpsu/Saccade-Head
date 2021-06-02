function [] = Ramp_saccade_head_wing_window()
%% Ramp_saccade_head_wing_window:
root = 'E:\DATA\Rigid_Data\Saccade';
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
velI = ones(size(velI));
velI(velI > n_speed) = velI(velI > n_speed) - n_speed;
[vel_fly_group, vel_group, fly_group] = findgroups(velI, Saccade.fly);
n_speed = length(unique(vel_group));

Scd.time        = splitapply(@(x) {cat(2,x{:})}, Saccade.scd_time, vel_fly_group);
Scd.head_pos    = splitapply(@(x) {cat(2,x{:})}, Saccade.head_scd_pos, vel_fly_group);
Scd.head_vel    = splitapply(@(x) {cat(2,x{:})}, Saccade.head_scd_vel, vel_fly_group);
Scd.head_accel  = splitapply(@(x) {cat(2,x{:})}, Saccade.head_scd_accel, vel_fly_group);
Scd.wing_pos    = splitapply(@(x) {cat(2,x{:})}, Saccade.wing_scd_pos, vel_fly_group);
Scd.wing_vel    = splitapply(@(x) {cat(2,x{:})}, Saccade.wing_scd_vel, vel_fly_group);

Scd = structfun(@(x) splitapply(@(y) {y}, x, vel_group), Scd, 'UniformOutput', false);
Scd = structfun(@(x) cat(2,x{:}), Scd, 'UniformOutput', false);

% Normalize position
Fs = round(Saccade.head_saccade{1}.Fs);
winI = (0.45*Fs):(0.47*Fs);
Scd.wing_pos = cellfun(@(x) x - mean(x(winI,:),1), ...
    Scd.wing_pos, 'UniformOutput', false);
Scd.head_pos = cellfun(@(x) x - mean(x(winI,:),1), ...
    Scd.head_pos, 'UniformOutput', false);
% winI = 0.4 * Fs;
% % Scd.wing_pos = cellfun(@(x) x - mean(x([1:winI],:),1), ...
% %     Scd.wing_pos, 'UniformOutput', false);
% Scd.wing_pos = cellfun(@(x) x - mean(x([1:winI, length(x)-winI:length(x)],:),1), ...
%     Scd.wing_pos, 'UniformOutput', false);

% Head-wing cross correlation
time = Scd.time{1}(:,1);
time_span = [-0.05 0.05];
[~,sI] = min(abs(time - time_span(1)));
[~,eI] = min(abs(time - time_span(2)));
span = sI:eI;
Scd.head_wing_cc = cell(N.fly,n_speed);
Scd.head_wing_max_cc = cell(N.fly,n_speed);
Scd.head_wing_time_diff = cell(N.fly,n_speed);
Scd.head_unsync = cell(N.fly,n_speed);
Scd.wing_unsync = cell(N.fly,n_speed);
Scd.cc_unsync = cell(N.fly,n_speed);
for v = 1:n_speed
    for n = 1:N.fly
       n_scd = size(Scd.head_pos{n,v},2);
       pp = 1;
       for s = 1:n_scd
            hpos = Scd.head_pos{n,v}(span,s);
            wpos = Scd.wing_pos{n,v}(span,s);
            [Scd.head_wing_cc{n,v}(:,s),lags] = xcorr(hpos, wpos, 'Normalized');
            if all(~isnan(Scd.head_wing_cc{n,v}(:,s)))
                [~,maxI] = max(abs(Scd.head_wing_cc{n,v}(:,s)));
                time_lag = lags/Fs;
                Scd.head_wing_max_cc{n,v}(1,pp) = Scd.head_wing_cc{n,v}(maxI,s);
                Scd.head_wing_time_diff{n,v}(1,pp) = time_lag(maxI);
                
                if Scd.head_wing_max_cc{n,v}(1,pp) < 0.5 %Scd.head_wing_time_diff{n,v}(1,pp) < 0
                    Scd.head_unsync{n,v}(:,end+1) = Scd.head_pos{n,v}(:,s);
                    Scd.wing_unsync{n,v}(:,end+1) = Scd.wing_pos{n,v}(:,s);
                    Scd.cc_unsync{n,v}(:,end+1) = Scd.head_wing_cc{n,v}(:,s);
                end             
                
                pp = pp + 1;
            end
       end
       if isempty(Scd.head_unsync{n,v})
          Scd.head_unsync{n,v} = nan(size(Scd.head_pos{n,v}(:,s)));
          Scd.wing_unsync{n,v} = nan(size(Scd.wing_pos{n,v}(:,s)));
          Scd.cc_unsync{n,v} = nan(size(Scd.head_wing_cc{n,v}(:,s)));
       end
    end
end

% Fly stats
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

%% Plot ALL by speed small window
clc
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 n_speed*2.25 4.5])
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
        plot(Scd.fly_mean.time{v}, Scd.fly_mean.head_pos{v}, ...
            'Color', [0.7*head_color 0.5], 'LineWidth', fly_lw)
        
        plot(Scd.fly_mean.time{v}, Scd.fly_mean.wing_pos{v}, ...
            'Color', [0.7*wing_color 0.3], 'LineWidth', fly_lw)

        [~] = PlotPatch(Scd.vel_stats.head_pos(v).mean, Scd.vel_stats.head_pos(v).std, ...
            Scd.vel_stats.time(v).mean, 1, 1, head_color, head_color, 0.3, vel_lw);
        
        [~] = PlotPatch(Scd.vel_stats.wing_pos(v).mean, Scd.vel_stats.wing_pos(v).std, ...
            Scd.vel_stats.time(v).mean, 1, 1, wing_color, wing_color, 0.3, vel_lw);
        
    ax(2,v) = subplot(3,n_speed,rowI(2)); hold on ; cla
        plot(Scd.fly_mean.time{v}, Scd.fly_mean.head_vel{v}, ...
            'Color', [0.7*head_color 0.5], 'LineWidth', fly_lw)
        
        plot(Scd.fly_mean.time{v}, Scd.fly_mean.wing_vel{v}, ...
            'Color', [0.7*wing_color 0.5], 'LineWidth', fly_lw)
        
        [~] = PlotPatch(Scd.vel_stats.head_vel(v).mean, Scd.vel_stats.head_vel(v).std, ...
            Scd.vel_stats.time(v).mean, n_std, 1, head_color, head_color, 0.3, vel_lw);
        
        [~] = PlotPatch(Scd.vel_stats.wing_vel(v).mean, Scd.vel_stats.wing_vel(v).std, ...
            Scd.vel_stats.time(v).mean, n_std, 1, wing_color, wing_color, 0.3, vel_lw);
 
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
set(ax, 'XLim', 0.05*[-1 1])
set(ax(1,:), 'YLim', 20*[-0.2 1])
set(ax(2,:), 'YLim', [-200 600])
set(ax(3,:), 'YLim', 40000*[-1 1])
set(ax(1:2,:), 'XTickLabel', [])
set(ax, 'XTick', -0.05:0.01:0.05)
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

%% Cross corelation
clc
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 n_speed*2 3*1.5])
movegui(fig, 'center')
clear ax H
ax = gobjects(1,n_speed);
fly_lw = 0.5;
vel_lw = 1.5;
cc = [0.4 0.2 0.8];
n_std = 0;
for v = 1:n_speed
    rowI = v + (0:2)*n_speed;
    ax(1,v) = subplot(3,n_speed,rowI(1)); hold on ; cla
        yline(0, '--k', 'LineWidth', 0.5)
        xline(0, '--k', 'LineWidth', 0.5)
        
        plot(1000*time_lag, cat(2,Scd.head_wing_cc{3,v}), ...
            'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.25)
        
        plot(1000*time_lag, Scd.fly_mean.head_wing_cc{v}, ...
            'Color', [0.7*cc 0.5], 'LineWidth', fly_lw)
        
        [~] = PlotPatch(Scd.vel_stats.head_wing_cc(v).mean, Scd.vel_stats.head_wing_cc(v).std, ...
            1000*time_lag, n_std, 1, cc, 0.7*cc, 0.3, vel_lw);
        
    ax(2,v) = subplot(3,n_speed,rowI(2)); hold on ; cla
        histogram(cat(2,Scd.head_wing_max_cc{:,v}), -1:0.05:1, ...
            'Normalization', 'Probability', 'FaceColor', 'k', 'EdgeColor', 'none')
        ax(2,v).YLim(1) = -0.01;
        
    ax(3,v) = subplot(3,n_speed,rowI(3)); hold on ; cla
        hbins = 1000*(-0.10125:0.0025:0.10125);
        histogram(1000*cat(2,Scd.head_wing_time_diff{:,v}), hbins, ...
            'Normalization', 'Probability', 'FaceColor', 'k', 'EdgeColor', 'none')
        ax(3,v).YLim(1) = -0.02;
        
end
set(ax, 'Color', 'none', 'LineWidth', 1, 'FontSize', 8)
set(ax(1,:), 'XLim', 1000*0.1*[-1 1], 'YLim', [-1 1])
for a = 1:size(ax,1)
   linkaxes(ax(a,:),'xy')
end
set(ax(2,:), 'XLim', [-1 1])
set(ax(3,:), 'XLim', 1000*0.1*[-1 1])
set(ax([1,3],:), 'XTick', -100:10:100)

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Cross correlation')
XLabelHC = get(ax(1,1), 'XLabel');
set([XLabelHC], 'String', 'Time lag (ms)')

YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Probability')
XLabelHC = get(ax(2,1), 'XLabel');
set([XLabelHC], 'String', 'Max correlation')

YLabelHC = get(ax(3,1), 'YLabel');
set([YLabelHC], 'String', 'Probability')
XLabelHC = get(ax(3,1), 'XLabel');
set([XLabelHC], 'String', 'Time lag (ms)')

%% Desync
fig = figure (3) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 n_speed*2.25 2*4.5/3])
movegui(fig, 'center')
clear ax H
ax = gobjects(1,n_speed);
fly_lw = 0.25;
head_color = [0 0 1];
wing_color = [1 0 0];
time_all = Scd.time{1}(:,1);
for v = 1:n_speed
    rowI = v + (0:1)*n_speed;
    ax(1,v) = subplot(1,n_speed,rowI(1)); hold on ; cla
%         plot(time_all, cat(2, Scd.head_unsync{:,v}), ...
%             'Color', [0.7*head_color 0.5], 'LineWidth', fly_lw)
%         
%         plot(time_all, cat(2, Scd.wing_unsync{:,v}), ...
%             'Color', [0.7*wing_color 0.3], 'LineWidth', fly_lw)
        
        plot(time_all, Scd.fly_mean.head_unsync{v}, ...
            'Color', [0.7*head_color 0.5], 'LineWidth', fly_lw)
        
        plot(time_all, Scd.fly_mean.wing_unsync{v}, ...
            'Color', [0.7*wing_color 0.3], 'LineWidth', fly_lw)
        
        [~] = PlotPatch(Scd.vel_stats.head_unsync(v).mean, Scd.vel_stats.head_unsync(v).std, ...
            time_all, 0, 1, head_color, 0.7*head_color, 0.3, 1);
        
        [~] = PlotPatch(Scd.vel_stats.wing_unsync(v).mean, Scd.vel_stats.wing_unsync(v).std, ...
            time_all, 0, 1, wing_color, 0.7*wing_color, 0.3, 1);
end
linkaxes(ax, 'x')
linkaxes(ax(1,:), 'y')
set(ax, 'XLim', 0.05*[-1 1])
set(ax(1,:), 'YLim', 20*[-1 1])
% set(ax(1:2,:), 'XTickLabel', [])
set(ax, 'XTick', -0.1:0.025:0.1)
% set(ax(1:2,:), 'XColor', 'none')
% set(ax(:,2:n_speed), 'YTickLabel', [])
% set(ax(:,2:n_speed), 'YColor', 'none')
set(ax, 'LineWidth', 1)

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Position (°)')
XLabelHC = get(ax(end,1), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')

%% Save
Scd.time_lag = time_lag;
savedir = 'E:\DATA\Rigid_Data\Saccade\processed';
filedata = textscan(char(FILE), '%s', 'delimiter', '_');
fclass = filedata{1}{1};
filename = [fclass '_saccade_window_head_wing.mat'];
save(fullfile(savedir, filename), 'Scd','U','N','-v7.3')
end