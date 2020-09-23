function [] = Ramp_Static_stats_heatmap()
%% Ramp_Static_stats_heatmap:
Ramp = load('H:\DATA\Rigid_Data\Saccade\combined\Ramp_All_Stats.mat');
Static = load('H:\DATA\Rigid_Data\Saccade\combined\Static_All_Stats.mat');

%% Saccade Statistics %%
CC = [0.3 0.1 0.6 ; 0.2 0.4 0.9];
% ramp.fly_group_all = Ramp.All_Stats.fly;

static_wave = [2:4];
Static.All_Stats = Static.All_Stats(any(Static.All_Stats.wave==static_wave,2),:);
Static.Count_Stats = Static.Count_Stats(any(Static.Count_Stats.wave==static_wave,2),:);
% static.fly_group_all = Static.All_Stats.fly;

stat_names = ["PeakVel", "Amplitude", "Duration","count"];
ylabel_names = ["Peak Velocity (°/s)", "Amplitude (°)", "Duration (s)","Frequency (Hz)"];
ylim_list = {[0 1500],[0 35],[0 0.1],[-0.1 3]};
n_plot = length(stat_names);

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 n_plot*2 1.5];
ax = gobjects(n_plot,1);
for ww = 1:n_plot
    ax(ww) = subplot(1,n_plot,ww);  axis tight
    
    if ww < 4
        ramp_data = abs(Ramp.All_Stats.(stat_names(ww)));
        static_data = abs(Static.All_Stats.(stat_names(ww)));
       	%ramp_data = splitapply(@(x) nanmean(x,1), ramp_data, ramp.fly_group_all);
        %static_data = splitapply(@(x) nanmean(x,1), static_data, static.fly_group_all);
    else
        ramp_data = abs(Ramp.Count_Stats.(stat_names(ww))) / 10;
        static_data = abs(Static.Count_Stats.(stat_names(ww))) / 10;
    end
    
    all_data = [ramp_data ; static_data];
    G = [1*ones(size(ramp_data)) ; 2*ones(size(static_data))];
    
    bx = boxplot(all_data, G, 'Labels', {'Ramp','Static'}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    ylabel(ylabel_names(ww))

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2}, CC(kk,:));
    end

    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);
    ylim(ylim_list{ww})
end

set(ax ,'LineWidth', 1, 'Box', 'off')

%% Heat Maps
FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = 2*[2 2 3 3];
movegui(FIG,'center')
colormap(bone)

edges.dur = (0:5:100)';
edges.amp = (4:0.5:30)';
edges.pkv = (250:10:1000)';

h = gobjects(6,1);
ax(1) = subplot(3,2,1);
    X = abs(Ramp.All_Stats.Amplitude);
    Y = abs(Ramp.All_Stats.PeakVel);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    [rho_move.amp_pkv,pval_move.amp_pkv] = corr(X,Y);
    h(1) = histogram2(X, Y, edges.amp, edges.pkv, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'ShowEmptyBins','on', 'Normalization','Probability','DisplayStyle','tile');
    axis tight
    grid off
    xlabel('Amplitude (°)')
    ylabel('Peak Velocity (°/s)')
    zlabel('Probability')
    cbar = colorbar('Location', 'eastoutside');
    cbar.Ticks = [0 max(h(1).Values,[],'all')];
	title('Moving')

ax(2) = subplot(3,2,2);
    X = abs(Static.All_Stats.Amplitude);
    Y = abs(Static.All_Stats.PeakVel);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    [rho_static.amp_pkv,pval_static.amp_pkv] = corr(X,Y);
    h(2) = histogram2(X, Y, edges.amp, edges.pkv, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'ShowEmptyBins','on', 'Normalization','Probability','DisplayStyle','tile');
    grid off
    %xlabel('Amplitude (°)')
    %ylabel('Peak Velocity (°/s)')
    zlabel('Probability')
    colorbar('Location', 'eastoutside')
    cbar = colorbar('Location', 'eastoutside');
    cbar.Ticks = [0 max(h(2).Values,[],'all')];
    title('Static')
   
ax(3) = subplot(3,2,3);
    X = abs(Ramp.All_Stats.Amplitude);
    Y = 1000*abs(Ramp.All_Stats.Duration);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    [rho_move.amp_dur,pval_move.amp_dur] = corr(X,Y);
    h(3) = histogram2(X, Y, edges.amp, edges.dur, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'ShowEmptyBins','on', 'Normalization','Probability','DisplayStyle','tile');
    grid off
    xlabel('Amplitude (°)')
    ylabel('Duration (s)')
    zlabel('Probability')
    cbar = colorbar('Location', 'eastoutside');
    cbar.Ticks = [0 max(h(3).Values,[],'all')];
    
ax(4) = subplot(3,2,4);
    X = abs(Static.All_Stats.Amplitude);
    Y = 1000*abs(Static.All_Stats.Duration);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    [rho_static.amp_dur,pval_static.amp_dur] = corr(X,Y);
    h(4) = histogram2(X, Y, edges.amp, edges.dur, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'ShowEmptyBins','on', 'Normalization','Probability','DisplayStyle','tile');
    grid off
    %xlabel('Amplitude (°)')
    %ylabel('Duration (s)')
    zlabel('Probability')
    cbar = colorbar('Location', 'eastoutside');
    cbar.Ticks = [0 max(h(4).Values,[],'all')];
    
ax(5) = subplot(3,2,5);
    X = abs(Ramp.All_Stats.PeakVel);
    Y = 1000*abs(Ramp.All_Stats.Duration);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    [rho_move.pkv_dur,pval_move.pkv_dur] = corr(X,Y);
    h(5) = histogram2(X, Y, edges.pkv, edges.dur, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'ShowEmptyBins','on', 'Normalization','Probability','DisplayStyle','tile');
    grid off
    xlabel('Peak Velocity (°/s)')
    ylabel('Duration (s)')
    zlabel('Probability')
    cbar = colorbar('Location', 'eastoutside');
    cbar.Ticks = [0 max(h(5).Values,[],'all')];
    
ax(6) = subplot(3,2,6);
    X = abs(Static.All_Stats.PeakVel);
    Y = 1000*abs(Static.All_Stats.Duration);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    [rho_static.pkv_dur,pval_static.pkv_dur] = corr(X,Y);
    h(6) = histogram2(X, Y, edges.pkv, edges.dur, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'ShowEmptyBins','on', 'Normalization','Probability','DisplayStyle','tile');
    grid off
    %xlabel('Peak Velocity (°/s)')
    %ylabel('Duration (s)')
    zlabel('Probability')
    cbar = colorbar('Location', 'eastoutside');
    cbar.Ticks = [0 max(h(6).Values,[],'all')];
    
set(ax, 'LineWidth', 1, 'Visible', 'on', 'Color', 'none', 'XColor', 'k', 'YColor', 'k')
% set(h, 'ShowEmptyBins', 'off')
% set(h,'EdgeColor','none')

% hp = get(subplot(3,2,4),'Position');
% cbar = colorbar('Position', [hp(1)+hp(3)+0.02  0.15+hp(2)  0.03  hp(2)+hp(3)*1]);
% cbar.Label.String = 'Probability';

% hp = get(subplot(3,2,6),'Position');
% cbar = colorbar('Position', [hp(1)+hp(3)+0.02  0.15+hp(2)  0.03  hp(2)+hp(3)*1]);
% cbar.Label.String = 'Probability';

%% Static Saccade Statistics by fly
clc

static_wave = [2:4];
Static.All_Stats = Static.All_Stats(any(Static.All_Stats.wave==static_wave,2),:);
Static.Count_Stats = Static.Count_Stats(any(Static.Count_Stats.wave==static_wave,2),:);

wave_group_all = Static.All_Stats.wave;
fly_group_all = Static.All_Stats.fly;

Wave = Static.U.wave{1};
n_wave = length(Wave);
CC = repmat(hsv(n_wave),2,1);

% Combine wavelengths
wave_group_all = wave_group_all ./ wave_group_all;
Wave = 1;

[fly_wave_group, wave_group, fly_group] = findgroups(wave_group_all, fly_group_all);

stat_names = ["PeakVel", "Amplitude", "Duration"];
ylabel_names = ["Peak Velocity (°/s)", "Amplitude (°)", "Duration (s)"];
ylim_list = {[0 1500],[0 35],[0 0.1],[-0.1 3]};
n_plot = length(stat_names);

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 (n_plot+1)*2 1.5];
ax = gobjects(n_plot+1,1);
for ww = 1:n_plot
    ax(ww) = subplot(1,n_plot+1,ww);  axis tight
    if ww == 3
        flip = 1;
    else
        flip = Static.All_Stats.Direction;
    end
    data = flip .* Static.All_Stats.(stat_names(ww));
    fly_mean = splitapply(@(x) nanmean(x,1), data, fly_wave_group);
    
    %bx = boxplot(data, wave_group_all, 'Labels', {Wave}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    bx = boxplot(fly_mean, wave_group, 'Labels', {Wave}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    ylabel(ylabel_names(ww))

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2}, CC(kk,:));
    end

    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);
    ylim(ylim_list{ww})
end

% Saccade Count/Rate
fly_group = Static.Count_Stats.fly;
wave_group_all = Static.Count_Stats.vel;
[fly_wave_group,fly_group,wave_group]  = findgroups(fly_group, wave_group_all);

count = Static.Count_Stats.count;
count_wave_fly_mean = splitapply(@(x) nanmean(x,1), count, fly_wave_group);

ww = 4;
ax(ww) = subplot(1,4,4); hold on

% bx = boxplot(count./10, wave_group_all, 'Labels', {Wave}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
bx = boxplot(count_wave_fly_mean./10, wave_group, 'Labels', {Wave}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
ylabel('Frequency (Hz)')

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1},h{kk,2},CC(kk,:));
end

set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(ww).Children = ax(ww).Children([end 1:end-1]);
ylim(ylim_list{ww})

set(ax, 'LineWidth', 1, 'Box', 'off')

end