function [] = Static_Saccade_Stats()
%% Static_Saccade_Stats:

load('C:\Users\BC\Box\Research\Manuscripts\Head Saccade\Data\Static_All_Stats.mat')

%% Saccade Statistics %%
All_Stats = All_Stats(All_Stats.Duration < 0.08,:);
wave_group_all = All_Stats.wave;
fly_group_all = All_Stats.fly;

Wave = U.wave{1};
n_wave = length(Wave);
CC = repmat(jet(n_wave),2,1);

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
    data = abs(All_Stats.(stat_names(ww)));
    fly_mean = splitapply(@(x) nanmean(x,1), data, fly_wave_group);
    
    bx = boxplot(data, wave_group_all, 'Labels', {Wave}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    xlabel('Wavelength (°)')
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
Count_Stats = Count_Stats(Count_Stats.count < 10,:);

fly_group = Count_Stats.fly;
wave_group_all = Count_Stats.wave;
[fly_wave_group,fly_group,wave_group] = findgroups(fly_group, wave_group_all);

count = Count_Stats.count;
count_vel_fly_mean = splitapply(@(x) nanmean(x,1), count, fly_wave_group);

ww = 4;
ax(ww) = subplot(1,4,4); hold on

bx = boxplot(count./10, wave_group_all, 'Labels', {Wave}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
xlabel('Wavelength (°)')
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

set(ax, 'LineWidth', 1)

%% Saccade Trigger Polar Plot, normalized to stimulus direction
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 3 4];
FIG.Name = 'Normalized Saccade Trigger';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h

stim_dir = sign(All_Stats.Direction);
CW  = stim_dir ==  1;
CCW = stim_dir == -1;
edges = deg2rad(-20:1:20);

ax(1) = subplot(1,1,1,polaraxes); grid off ; axis tight
    Trig_Positive = [ All_Stats.StartPos(CW); -All_Stats.StartPos(CCW)];
    Trig_Negative = [ All_Stats.EndPos(CW);  -All_Stats.EndPos(CCW)];
    h(1) = polarhistogram(deg2rad(Trig_Positive),edges,...
        'FaceColor','g','FaceAlpha',0.5,'Normalization','Probability'); hold on
    h(2) = polarhistogram(deg2rad(Trig_Negative),'BinEdges', edges, ...
        'FaceColor','r','FaceAlpha',0.5,'Normalization','Probability');
    title('CW Stimulus Direction   ===>')
    ax(1).ThetaAxis.Label.String = 'Head Position (°)';
    
set(h,'EdgeColor','none')
set(ax,'FontSize',8);
set(ax,'Color','w');
set(ax,'ThetaLim',[-20 20]);
set(ax,'RLim',[0 0.105]);
set(ax,'RTick',[0:0.02:0.1]);
set(ax,'ThetaDir','clockwise')
set(ax,'ThetaTick',-20:10:20);
set(ax,'ThetaZeroLocation','top');

set(ax, 'LineWidth', 1)

%% ANOVA
% [p,tbl,stats,terms] = anovan(count, {wave_group_all, fly_group_all});
% [p,tbl,stats,terms] = anovan(count_vel_fly_mean, {wave_group, fly_group});

end