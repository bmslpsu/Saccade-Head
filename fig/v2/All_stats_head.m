function [] = All_stats_head()
%% All_stats_head:
Ramp = load('E:\DATA\Rigid_Data\Saccade\processed\Ramp_All_Stats.mat');
Static = load('E:\DATA\Rigid_Data\Saccade\processed\Static_All_Stats.mat');

%% Head saccade stats
clearvars -except Ramp Static
clc

waveI = any(Static.All_Stats.Wave == [22.5 30 60],2);
Static.All_stats = Static.All_Stats(waveI,:);

n_speed = 5;
ramp_vel_group_all = Ramp.All_Stats.vel;
ramp_vel_group_all(ramp_vel_group_all > n_speed) = ...
    ramp_vel_group_all(ramp_vel_group_all > n_speed) - n_speed;
% ramp_vel_group_all = ones(size(ramp_vel_group_all));
vel_group_all = [0*Static.All_Stats.vel ; ramp_vel_group_all];

[ramp_fly_vel_group,ramp_vel_group,ramp_fly_group] = ...
    findgroups(ramp_vel_group_all, Ramp.All_Stats.fly);
[static_fly_vel_group,static_vel_group,static_fly_group] = ...
    findgroups(Static.All_Stats.vel, Static.All_Stats.fly);
fly_vel_group_all = [0*static_vel_group ; ramp_vel_group];

stats = ["Amplitude", "PeakVel", "Duration", "Skew", "StartPos", "EndPos"];
yL = [0 35; 0 1200; 0 0.14; 0 1; -25 15; -15 25];
n_stat = length(stats);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 n_stat*2 2])
movegui(fig, 'center')
ax = gobjects(n_stat,1);
Y = [];
for n = 1:n_stat
        y = [Static.All_Stats.(stats(n)) ; Ramp.All_Stats.(stats(n))];
        if any(y < 0)
            y = y .* [Static.All_Stats.Direction ; Ramp.All_Stats.Direction];
            y_fly_ramp = splitapply(@(x) nanmean(x,1), ...
                Ramp.All_Stats.Direction.*Ramp.All_Stats.(stats(n)), ramp_fly_vel_group);
            y_fly_static = splitapply(@(x) nanmean(x,1), ...
                Static.All_Stats.Direction.*Static.All_Stats.(stats(n)), static_fly_vel_group);
        else
            y_fly_ramp = splitapply(@(x) nanmean(x,1), ...
                Ramp.All_Stats.(stats(n)), ramp_fly_vel_group);
            y_fly_static = splitapply(@(x) nanmean(x,1), ...
                Static.All_Stats.(stats(n)), static_fly_vel_group);
        end
        y_fly = [y_fly_static ; y_fly_ramp];

        vel_group_temp = vel_group_all(~isnan(y));
        y = y(~isnan(y));
        Y.(stats(n)) = y;
        
        fly_vel_group_temp = fly_vel_group_all(~isnan(y_fly));
        Y.fly.(stats(n)) = y_fly;
        G_fly = fly_vel_group_temp;
    
    ax(n) = subplot(1,n_stat,n); hold on ; 
        ylim(yL(n,:))
        ylabel(stats(n))
        G = vel_group_temp;
        boxchart(y, 'GroupByColor', G, ...
            'LineWidth', 0.5, 'MarkerStyle', '.', 'JitterOutliers', 'on', 'Notch', 'off');
%         boxchart(y_fly, 'GroupByColor', G_fly, ...
%             'LineWidth', 0.5, 'MarkerStyle', '.', 'JitterOutliers', 'on', 'Notch', 'off');
end

set(ax, 'LineWidth', 1, 'Color', 'none', 'XColor', 'none')

%% Stats
clc
stat_name = 'Duration';
speedI = [1:5];
keepI = any(G==speedI,2);
G_speedI = G(keepI);
temp = Y.(stat_name)(keepI);
% [p,tbl,stats] = anova1(temp, G_speedI);
[p,tbl,stats] = kruskalwallis(temp, G_speedI);
c = multcompare(stats, 'alpha', 0.01);
M = splitapply(@(x) mean(x), temp, findgroups(G_speedI));

%% Fly stats
clc
stat_name = 'Duration';
speedI = [0:5];
keepI = any(G_fly==speedI,2);
G_speedI = G_fly(keepI);
temp = Y.fly.(stat_name)(keepI);
% [p,tbl,stats] = anova1(temp, G_speedI);
[p,tbl,stats] = kruskalwallis(temp, G_speedI);
c = multcompare(stats, 'alpha', 0.001);
M = splitapply(@(x) mean(x), temp, findgroups(G_speedI));

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
G_no_static = G_count(G_count~=0);
temp = Count(G_count~=0);
% [p,tbl,stats] = anova1(Count, G_count);
[p,tbl,stats] = kruskalwallis(temp, G_no_static);
c = multcompare(stats, 'alpha', 0.01);

%% Save
fname = 'Head_all_stats';
root = 'H:\DATA\Rigid_Data\Saccade';
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'Y', 'G', 'Count', 'G_count');

end