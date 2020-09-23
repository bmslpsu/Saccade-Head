function [] = Ramp_saccade_stats()
%% Ramp_saccade_stats:
load('H:\DATA\Rigid_Data\Saccade\combined\Ramp_All_Stats.mat')

%% Saccade Statistics by velocity
wave_group_all = All_Stats.wave;
wave = unique(wave_group_all);
% All_Stats = All_Stats(wave_group_all==1,:);

vel_group_all = All_Stats.vel;
fly_group_all = All_Stats.fly;

Vel = unique(All_Stats.Vel);
n_speed = length(Vel)/2;
Speed = Vel(n_speed+1:end);
CC = repmat(hsv(n_speed),2,1);

vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;

[fly_vel_group, vel_group, fly_group] = findgroups(vel_group_all, fly_group_all);
[wave_vel_group, ~, wave_group] = findgroups(vel_group_all, wave_group_all);
% vel_group(vel_group > n_speed) = vel_group(vel_group > n_speed) - n_speed;

stat_names = ["PeakVel", "Amplitude", "Duration", "Skew", "StartPos", "EndPos"];
ylabel_names = ["Peak Velocity (°/s)", "Amplitude (°)", "Duration (s)", "Skew", "Start (°)", "End (°)"];
ylim_list = {[0 1500],[0 35],[0 0.1],[0 1],[-30 30], [-30 30], [-0.1 3],};
n_plot = length(stat_names);

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 (n_plot+1)*2 2];
ax = gobjects(n_plot+1,1);
for ww = 1:n_plot
    ax(ww) = subplot(1,n_plot+1,ww);  axis tight
    data = All_Stats.(stat_names(ww));
    if any(data < 0)
        data = data .* All_Stats.Direction;
    end
    fly_mean = splitapply(@(x) nanmean(x,1), data, fly_vel_group);
    
    bx = boxplot(data, vel_group_all, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    %bx = boxplot(fly_mean, vel_group, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    xlabel('Stimulus Speed (°/s)')
    ylabel(ylabel_names(ww))

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1}, h{kk,2}, CC(kk,:),  'EdgeColor', 'none');
    end

    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',0.75);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);
    ylim(ylim_list{ww})
end

%% Saccade Count/Rate
fly_group = Count_Stats.fly;
vel_group_all = Count_Stats.vel;
vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
[fly_vel_group,fly_group,vel_group]  = findgroups(fly_group, vel_group_all);

count = Count_Stats.count;
count_vel_fly_mean = splitapply(@(x) nanmean(x,1), count, fly_vel_group);

ww = length(ax);
ax(ww) = subplot(1,n_plot+1,ww); hold on

bx = boxplot(count./10, vel_group_all, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
% bx = boxplot(count_vel_fly_mean./10, vel_group, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
xlabel('Stimulus Speed (°/s)')
ylabel('Frequency (Hz)')

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1}, h{kk,2}, CC(kk,:), 'EdgeColor', 'none');
end

set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(ww).Children = ax(ww).Children([end 1:end-1]);
ylim(ylim_list{ww})

set(ax, 'LineWidth', 0.75, 'Box', 'off')

%% By Wavelength
cc_wave = parula(length(wave));
FIG = figure (11) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 (n_plot+1)*2 2*1.5];
ax = gobjects(n_plot+1,1);
for ww = 1:n_plot
    ax(ww) = subplot(2,ceil(1 + n_plot/2),ww);  axis tight
    data = All_Stats.(stat_names(ww));
    if any(data < 0)
        data = data .* All_Stats.Direction;
    end
    
    bx = boxplot(data, wave_group_all, 'Labels', {wave}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
%     bx = boxplot(data, wave_vel_group, 'Labels', {wave_group}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    xlabel('Stimulus Speed (°/s)')
    ylabel(ylabel_names(ww))

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1}, h{kk,2}, cc_wave(kk,:),  'EdgeColor', 'none');
    end

    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);
    ylim(ylim_list{ww})
end

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

end