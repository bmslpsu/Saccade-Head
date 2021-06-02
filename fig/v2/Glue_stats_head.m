function [] = Glue_stats_head()
%% Glue_stats_head:
Ramp = load('E:\DATA\Rigid_Data\Saccade\processed\Ramp_All_Stats.mat');
Static = load('E:\DATA\Rigid_Data\Saccade\processed\Static_All_Stats.mat');
Glue = load('E:\DATA\Rigid_Data\Saccade\processed\Ramp_Glue_All_Stats.mat');

%% Head saccade stats
clearvars -except Ramp Static Glue
clc

n_speed = 5;
vel_group_all = Ramp.All_Stats.Direction + 2;
vel_group_all = ones(size(vel_group_all));
% vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
vel_group_all = [0*Static.All_Stats.vel ; vel_group_all ; Glue.All_Stats.Direction+5];

stats = ["Amplitude", "PeakVel", "Duration", "Skew", "StartPos", "EndPos"];
yL = [0 35; 0 1200; 0 0.14; 0 1; -25 25; -15 25];
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
        y = [Static.All_Stats.(stats(n)) ; Ramp.All_Stats.(stats(n)) ; Glue.All_Stats.(stats(n))];
        if strcmp(stats(n), 'StartPos')
            sense = -1;
        else
            sense = 1;
        end
        
        if any(y < 0)
           y = sense*y .* [Static.All_Stats.Direction ; Ramp.All_Stats.Direction ; Glue.All_Stats.Direction];
        end
        vel_group_temp = vel_group_all(~isnan(y));
        y = y(~isnan(y));
        Y.(stats(n)) = y;
        G = vel_group_temp;
        boxchart(y, 'GroupByColor', G, ...
            'LineWidth', 0.5, 'MarkerStyle', 'none', 'JitterOutliers', 'on', 'Notch', 'off');
end

set(ax, 'LineWidth', 1, 'Color', 'none', 'XColor', 'none')

%% Stats
stat_name = 'StartPos';
[p,tbl,stats] = anova1(Y.(stat_name), G);
% [p,tbl,stats] = kruskalwallis(Y.(stat_name), G);
figure
c = multcompare(stats, 'alpha', 0.001);

%% Saccade frequency
% vel_group_all = Ramp.Count_Stats.vel;
% vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
% vel_group_all = [0*Static.Count_Stats.vel ; vel_group_all];
velI = [1 2 6 7];
velI_keep = any(Ramp.Count_Stats.vel == velI, 2);
Ramp_temp = Ramp.Count_Stats(velI_keep,:);
vel_group_all = Ramp_temp.vel;
n_speed = 5;
vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
vel_group_all = [0*Static.Count_Stats.vel ; vel_group_all ; findgroups(Glue.Count_Stats.vel) + 10];

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 2 2])
movegui(fig, 'center')
clear ax
ax(1) = subplot(1,1,1); hold on ; 
    ylim([0 2])
    ylabel('Frequency (Hz)')
    Count = [Static.Count_Stats.count ; Ramp_temp.count ; Glue.Count_Stats.count] ./ 10;
    G_count = vel_group_all;
    boxchart(Count, 'GroupByColor', G_count, ...
        'LineWidth', 0.5, 'MarkerStyle', 'none', 'JitterOutliers', 'on', 'Notch', 'off');

set(ax, 'LineWidth', 1, 'Color', 'none', 'XColor', 'none')

%% Stats
[p,tbl,stats] = anova1(Count, G_count);
% [p,tbl,stats] = kruskalwallis(Count, G_count);
figure
c = multcompare(stats, 'alpha', 0.473);

%% T-test
% close all
clc
g1 = 0;
g2 = [11 12];

G = G_count;
% G = findgroups(G_count);

G1 = any(G == g1,2);
Y1 = Count(G1);
G1 = G(G1);
G2 = any(G == g2,2);
Y2 = Count(G2);
G2 = G(G2);

figure (100)
boxplot([Y1;Y2], [1*ones(size(G1)); 2*ones(size(G2))])

[~,p] = ttest2(Y1,Y2)

% anova1([Y1;Y2], [G1;G2])

%% Save
fname = 'Head_all_stats';
root = 'H:\DATA\Rigid_Data\Saccade';
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'Y', 'G', 'Count', 'G_count');

end