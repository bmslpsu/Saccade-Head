function [] = Saccade_Stats_SS()
%% Saccade_Stats_SS:

load('C:\Users\BC\Box\Research\Manuscripts\Head Saccade\Data\SS_All_Stats.mat')

%% Saccade Rate by amplitude
freqI = 3;
All_Stats = All_Stats(All_Stats.Duration < 0.075,:);
All_Stats = All_Stats(All_Stats.freq==freqI,:);
[amp_group,Amp] = findgroups(All_Stats.amp);

n_amp = length(Amp);
CC = repmat(prism(n_amp),2,1);

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
    
    bx = boxplot(data, amp_group, 'Labels', {Amp}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    xlabel('Amplitude (Hz)')
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
Count_Stats = Count_Stats(Count_Stats.freq==freqI,:);
[amp_group,Amp] = findgroups(Count_Stats.amp);

count = Count_Stats.count;

ww = 4;
ax(ww) = subplot(1,4,4); hold on

bx = boxplot(count./10, amp_group, 'Labels', {Amp}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
xlabel('Amplitude (Hz)')
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


%% Saccade Statistics
amp = 15;
All_Stats = All_Stats(All_Stats.amp==amp,:);

All_Stats = All_Stats(All_Stats.Duration < 0.075,:);


freq_group_all = All_Stats.freq;
fly_group_all = All_Stats.fly;

Freq = unique(All_Stats.Freq);
n_freq = length(Freq);
CC = repmat(hsv(n_freq),2,1);

[fly_vel_group, freq_group, fly_group] = findgroups(freq_group_all, fly_group_all);

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
    
    bx = boxplot(data, freq_group_all, 'Labels', {Freq}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    xlabel('Amplitude (°)')
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
Count_Stats = Count_Stats(Count_Stats.amp==amp,:);
% Count_Stats = Count_Stats(Count_Stats.count < 14,:);
fly_group = Count_Stats.fly;
freq_group_all = Count_Stats.freq;
[fly_vel_group,fly_group,freq_group]  = findgroups(fly_group, freq_group_all);

count = Count_Stats.count;
count_freq_fly_mean = splitapply(@(x) nanmean(x,1), count, fly_vel_group);

ww = 4;
ax(ww) = subplot(1,4,4); hold on

bx = boxplot(count./10, freq_group_all, 'Labels', {Freq}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
xlabel('Frequency (Hz)')
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

%% Saccade Trigger Polar Plot, normalized to stimulus direction
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 3 4];
FIG.Name = 'Normalized Saccade Trigger';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h

amp = 18.75;
All_Stats = All_Stats(All_Stats.amp==amp,:);

% freq = 3;
% All_Stats = All_Stats(All_Stats.freq==freq,:);

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
set(ax,'RLim',[0 0.1]);
set(ax,'ThetaDir','clockwise')
set(ax,'ThetaTick',-20:10:20);
set(ax,'ThetaZeroLocation','top');
% leg = legend('Start','End','Location','North');
% leg.Location = 'northwest';
% leg.Box = 'off';

set(ax, 'LineWidth', 1)

%% ANOVA

end