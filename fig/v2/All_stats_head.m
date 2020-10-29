function [] = All_stats_head()
%% All_stats_head:
Ramp = load('H:\DATA\Rigid_Data\Saccade\combined\Ramp_All_Stats.mat');
Static = load('H:\DATA\Rigid_Data\Saccade\combined\Static_All_Stats.mat');

%% Head saccade stats
clearvars -except Ramp Static
clc
n_speed = 5;
vel_group_all = Ramp.All_Stats.vel;
vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
vel_group_all = [0*Static.All_Stats.vel ; vel_group_all];

stats = ["Amplitude", "PeakVel", "Duration", "Skew", "StartPos", "EndPos"];
yL = [0 35; 0 1200; 0 0.14; 0 1; -25 15; -15 25];
n_stat = length(stats);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 n_stat*2 2])
movegui(fig, 'center')
ax = gobjects(n_stat,1);
Y = [];
for n = 1:n_stat
    ax(n) = subplot(1,n_stat,n); hold on ; 
        ylim(yL(n,:))
        ylabel(stats(n))
        y = [Static.All_Stats.(stats(n)) ; Ramp.All_Stats.(stats(n))];
        if any(y < 0)
           y = y .* [Static.All_Stats.Direction ; Ramp.All_Stats.Direction];
        end
        vel_group_temp = vel_group_all(~isnan(y));
        y = y(~isnan(y));
        Y.(stats(n)) = y;
        G = vel_group_temp;
        boxchart(y, 'GroupByColor', G, ...
            'LineWidth', 0.5, 'MarkerStyle', '.', 'JitterOutliers', 'on', 'Notch', 'off');
end

set(ax, 'LineWidth', 1, 'Color', 'none', 'XColor', 'none')

%% Stats
stat_name = 'EndPos';
[p,tbl,stats] = anova1(Y.(stat_name), G);
% [p,tbl,stats] = kruskalwallis(Y.(stat_name), G);
c = multcompare(stats, 'alpha', 0.001);

%% Saccade frequency
vel_group_all = Ramp.Count_Stats.vel;
vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
vel_group_all = [0*Static.Count_Stats.vel ; vel_group_all];

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 2 2])
movegui(fig, 'center')
clear ax
ax(1) = subplot(1,1,1); hold on ; 
    ylim([0 3])
    ylabel('Frequency (Hz)')
    Count = [Static.Count_Stats.count ; Ramp.Count_Stats.count] ./ 10;
    G_count = vel_group_all;
    boxchart(Count, 'GroupByColor', G_count, ...
        'LineWidth', 0.5, 'MarkerStyle', '.', 'JitterOutliers', 'on', 'Notch', 'off');

set(ax, 'LineWidth', 1, 'Color', 'none', 'XColor', 'none')

%% Stats
[p,tbl,stats] = anova1(Count, G_count);
% [p,tbl,stats] = kruskalwallis(Count, G_count);
c = multcompare(stats, 'alpha', 0.001);

%% Save
fname = 'Head_all_stats';
root = 'H:\DATA\Rigid_Data\Saccade';
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'Y', 'G', 'Count', 'G_count');

end