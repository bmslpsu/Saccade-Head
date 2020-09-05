function [] = Saccade_Stats_new()
%% Saccade_Stats_new:

load('C:\Users\BC\Box\Research\Manuscripts\Head Saccade\Data\Ramp_All_Stats.mat')

%% Saccade Statistics %%
vel_group_all = All_Stats.vel;
fly_group_all = All_Stats.fly;

Vel = unique(All_Stats.Vel);
n_speed = length(Vel)/2;
Speed = Vel(n_speed+1:end);
CC = repmat(hsv(n_speed),2,1);

vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
[fly_vel_group, vel_group, fly_group] = findgroups(vel_group_all, fly_group_all);
% vel_group(vel_group > n_speed) = vel_group(vel_group > n_speed) - n_speed;

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
    data = abs(All_Stats.(stat_names(ww)));
    fly_mean = splitapply(@(x) nanmean(x,1), data, fly_vel_group);
    
    %bx = boxplot(data, vel_group_all, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    bx = boxplot(fly_mean, vel_group, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    xlabel('Stimulus Speed (°/s)')
    ylabel(ylabel_names(ww))

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1}, h{kk,2}, CC(kk,:),  'EdgeColor', 'none'););
    end

    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);
    ylim(ylim_list{ww})
end

% Saccade Count/Rate
fly_group = Count_Stats.fly;
vel_group_all = Count_Stats.vel;
vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
[fly_vel_group,fly_group,vel_group]  = findgroups(fly_group, vel_group_all);

count = Count_Stats.count;
count_vel_fly_mean = splitapply(@(x) nanmean(x,1), count, fly_vel_group);

ww = 4;
ax(ww) = subplot(1,4,4); hold on

% bx = boxplot(count./10, vel_group_all, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
bx = boxplot(count_vel_fly_mean./10, vel_group, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
xlabel('Stimulus Speed (°/s)')
ylabel('Frequency (Hz)')

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1}, h{kk,2}, CC(kk,:), 'EdgeColor', 'none'););
end

set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(ww).Children = ax(ww).Children([end 1:end-1]);
ylim(ylim_list{ww})

set(ax, 'LineWidth', 1, 'Box', 'off')

%% Saccade Trigger Polar Plot, normalized to stimulus direction
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 3 4];
FIG.Name = 'Normalized Saccade Trigger';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h

stim_dir = sign(All_Stats.Vel);
CW  = stim_dir ==  1;
CCW = stim_dir == -1;
edges = deg2rad(-25:1:25);

ax(1) = subplot(1,1,1,polaraxes); grid off ; axis tight
    Trig_Start = [ -All_Stats.StartPos(CW); All_Stats.StartPos(CCW)];
    Trig_End = [ -All_Stats.EndPos(CW);  All_Stats.EndPos(CCW)];
    h(1) = polarhistogram(deg2rad(Trig_Start),edges,...
        'FaceColor','g','FaceAlpha',0.5,'Normalization','Probability'); hold on
    h(2) = polarhistogram(deg2rad(Trig_End),'BinEdges', edges, ...
        'FaceColor','r','FaceAlpha',0.5,'Normalization','Probability');
    title('CW Stimulus Direction   ===>')
    ax(1).ThetaAxis.Label.String = 'Head Position (°)';
    
set(h,'EdgeColor','none')
set(ax,'FontSize',8);
set(ax,'Color','w');
set(ax,'ThetaLim',[-20 20]);
set(ax,'RLim',[0 0.1]);
set(ax,'ThetaDir','clockwise')
set(ax,'ThetaTick',-20:10:20);
set(ax,'ThetaZeroLocation','top');
% leg = legend('Start','End','Location','North');
% leg.Location = 'northwest';
% leg.Box = 'off';

set(ax, 'LineWidth', 1)

%% Error
vel_group_all = All_Stats.vel;
fly_group_all = All_Stats.fly;

Vel = unique(All_Stats.Vel);
n_speed = length(Vel)/2;
Speed = Vel(n_speed+1:end);
CC = repmat(hsv(n_speed),2,1);

% Flip_Stats = All_Stats;
% Flip_Stats{vel_group_all > n_speed,15:25} = -Flip_Stats{vel_group_all > n_speed,15:25};

vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
[fly_vel_group, vel_group, fly_group] = findgroups(vel_group_all, fly_group_all);
% vel_group(vel_group > n_speed) = vel_group(vel_group > n_speed) - n_speed;

stat_names = ["ErrorPos", "IntErrorPos", "ErrorVel", "IntErrorVel", "StartPos", "EndPos"];
ylabel_names = ["Position Error (°)", "Integrated Position Error (°s)", ...
    "Velocity Error (°/s)","Integrated Velocity Error (°)", ...
    "Start Pos (°)", "End Pos (°)"];
ylim_list = {[0 1500],[0 35],[0 0.1],[-0.1 3]};
n_plot = length(stat_names);

FIG = figure (3) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 (n_plot+1)*2 1.5];
ax = gobjects(n_plot+1,1);
for ww = 1:n_plot
    ax(ww) = subplot(1,n_plot+1,ww);  axis tight
    %data = All_Stats.Direction .* (All_Stats.(stat_names(ww)));
    data = -All_Stats.Direction .* (All_Stats.(stat_names(ww)));
    fly_mean = splitapply(@(x) nanmean(x,1), data, fly_vel_group);
    
  	bx = boxplot(fly_mean, vel_group, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    %bx = boxplot(data, vel_group_all, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
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
    %ylim(ylim_list{ww})
end

%% ANOVA
vel_group_all = All_Stats.vel;
fly_group_all = All_Stats.fly;
vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
[fly_vel_group, vel_group, fly_group] = findgroups(vel_group_all, fly_group_all);

%% Amplitude
data = abs(All_Stats.Amplitude);
[p,tbl,stats,terms] = anovan(data, {vel_group_all,fly_group_all});

data = splitapply(@(x) nanmean(x,1), data, fly_vel_group);
[p,tbl,stats,terms] = anovan(data, {vel_group,fly_group});

%% Peak Velocity
data = abs(All_Stats.PeakVel);
[p,tbl,stats,terms] = anovan(data, {vel_group_all,fly_group_all});

data = splitapply(@(x) nanmean(x,1), data, fly_vel_group);
[p,tbl,stats,terms] = anovan(data, {vel_group,fly_group});

%% Duration
data = abs(All_Stats.Duration);
[p,tbl,stats,terms] = anovan(data, {vel_group_all,fly_group_all});

data = splitapply(@(x) nanmean(x,1), data, fly_vel_group);
[p,tbl,stats,terms] = anovan(data, {vel_group,fly_group});

%% Rate
fly_group_all = Count_Stats.fly;
vel_group_all = Count_Stats.vel;
vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
[fly_vel_group,fly_group,vel_group]  = findgroups(fly_group_all, vel_group_all);

count = Count_Stats.count;
count_vel_fly_mean = splitapply(@(x) nanmean(x,1), count, fly_vel_group);

[p,tbl,stats,terms] = anovan(count, {vel_group_all, fly_group_all});
[p,tbl,stats,terms] = anovan(count_vel_fly_mean, {vel_group, fly_group});


end