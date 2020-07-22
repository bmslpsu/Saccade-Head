function [] = Saccade_Stats()
%% Saccade_Stats:

root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE','HEAD_SACCADE_STATS','U','N')

% clearvars -except clms CC Vel PATH COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

n_speed = N.vel/2;
CC = repmat(hsv(n_speed),2,1);
Vel = U.vel{1};
Speed = Vel(1:n_speed);

%% Saccade Statistics %%
vel_group_all = HEAD_SACCADE_STATS.vel;
fly_group_all = HEAD_SACCADE_STATS.fly;
vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
[fly_vel_group, vel_group, fly_group] = findgroups(vel_group_all, fly_group_all);
% vel_group(vel_group > n_speed) = vel_group(vel_group > n_speed) - n_speed;

stat_names = ["PeakVel", "Amplitude", "Duration"];
ylabel_names = ["Peak Velocity (°/s)", "Amplitude (°)", "Duration (s)"];
n_plot = length(stat_names);

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 (n_plot+1)*2 1.5];
ax = gobjects(n_plot+1,1);
for ww = 1:n_plot
    ax(ww) = subplot(1,n_plot+1,ww); axis tight
    
    data = abs(HEAD_SACCADE_STATS.(stat_names(ww)));
    fly_mean = splitapply(@(x) nanmean(x,1), data, fly_vel_group);
    
    bx = boxplot(data, vel_group_all, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    xlabel('Stimulus Speed (°/s)')
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
    ax(ww).YLim(1) = 0;
end

% Saccade Count/Rate
Saccade = SACCADE(:,:);
fly_group = Saccade.fly;
vel_group_all = Saccade.vel;
vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
[fly_vel_group,fly_group,vel_group]  = findgroups(fly_group, vel_group_all);

count = cellfun(@(x) x.count, Saccade.head_saccade, 'UniformOutput', true);
count_vel_fly_mean = splitapply(@(x) nanmean(x,1), count, fly_vel_group);

ww = 4;
ax(ww) = subplot(1,4,4); hold on

bx = boxplot(count./10, vel_group_all, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
xlabel('Stimulus Speed (°/s)')
ylabel('Frequency (Hz)')
ylim([-0.1 3])

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1},h{kk,2},CC(kk,:));
end

set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(ww).Children = ax(ww).Children([end 1:end-1]);

set(ax, 'LineWidth', 1)

%% ANOVA
vel_group_all = HEAD_SACCADE_STATS.vel;
fly_group_all = HEAD_SACCADE_STATS.fly;
vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
[fly_vel_group, vel_group, fly_group] = findgroups(vel_group_all, fly_group_all);

%% Amplitude
data = abs(HEAD_SACCADE_STATS.Amplitude);
[p,tbl,stats,terms] = anovan(data, {vel_group_all,fly_group_all});

data = splitapply(@(x) nanmean(x,1), data, fly_vel_group);
[p,tbl,stats,terms] = anovan(data, {vel_group,fly_group});

%% Peak Velocity
data = abs(HEAD_SACCADE_STATS.PeakVel);
[p,tbl,stats,terms] = anovan(data, {vel_group_all,fly_group_all});

data = splitapply(@(x) nanmean(x,1), data, fly_vel_group);
[p,tbl,stats,terms] = anovan(data, {vel_group_all,fly_group});

%% Duration
data = abs(HEAD_SACCADE_STATS.Duration);
[p,tbl,stats,terms] = anovan(data, {vel_group_all,fly_group_all});

data = splitapply(@(x) nanmean(x,1), data, fly_vel_group);
[p,tbl,stats,terms] = anovan(data, {vel_group_all,fly_group});

%% Rate
Saccade = SACCADE(:,:);
fly_group_all = Saccade.fly;
vel_group_all = Saccade.vel;
vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
[fly_vel_group,fly_group,vel_group_all]  = findgroups(fly_group, vel_group_all);

count = cellfun(@(x) x.count, Saccade.head_saccade, 'UniformOutput', true);
count_vel_fly_mean = splitapply(@(x) nanmean(x,1), count, fly_vel_group);

[p,tbl,stats,terms] = anovan(count, {vel_group_all, fly_group_all});
[p,tbl,stats,terms] = anovan(count_vel_fly_mean, {vel_group_all, fly_group});


end