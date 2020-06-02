function [] = MovingvsStatic_Saccade_Stats()
%% MovingvsStatic_Saccade_Stats:

root = 'H:\DATA\Rigid_Data\';

[Move,MovePATH] = uigetfile({'*.mat'},'Select moving data', root, 'MultiSelect','off');
[Static,StaticPATH] = uigetfile({'*.mat'},'Select static data', root, 'MultiSelect','off');

Move = load(fullfile(MovePATH,Move),'PATH','COUNT','SACCADE_STATS','D','I','U','N');
Static = load(fullfile(StaticPATH,Static),'PATH','COUNT','SACCADE_STATS','D','I','U','N');

%% Boxplots
CC = [0.8 0 0 ; 0 0 0.8];
stats_list = {'Amplitude','PeakVel','Duration'};
ylim_list = [35 1200 100];

YY = [  "Amplitude (°)",...
        "Peak Velocity (°/s)",...
        "Duration (ms)"];

n_plot = length(YY);
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches');
fig.Position = [2 4 4 4];
ax = gobjects(n_plot+1,1);
for ww = 1:n_plot
    ax(ww) = subplot(2,2,ww); axis tight
    move = abs(Move.SACCADE_STATS.(stats_list{ww}));
    static = abs(Static.SACCADE_STATS.(stats_list{ww}));
    G = [ones(length(move),1) ; 2*ones(length(static),1)];
    
    data = [move ; static];
    if ww==3
        data = 1000*data;
    end
    
    bx = boxplot(data,G,'Labels',{'Moving','Static'},'Width',0.5,'Symbol','.','Whisker',2);
    % xlabel('Stimulus Speed (°/s)')
    ylabel(YY(ww))
    box off

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
    ax(ww).YLim(2) = ylim_list(ww);
end

% Saccade Count/Rate
move_count = cat(1, Move.COUNT {:});
static_count = cat(1, Static.COUNT{:});
count_all = [move_count ; static_count];
G = [ones(length(move_count),1) ; 2*ones(length(static_count),1)];

ax(ww+1) = subplot(2,2,ww+1); axis tight
bx = boxplot(count_all./10, G, 'Labels', {'Moving','Static'}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
% xlabel('Stimulus Speed (°/s)')
ylabel('Rate (#/s)')
box off
ylim([0 3])

ax(2).YTick = [ax(2).YTick 1200];
ax(1).YTick = [ax(1).YTick 35];

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1},h{kk,2},CC(kk,:));
end

set(findobj(ax(ww+1),'tag','Median'), 'Color', 'w','LineWidth',1.5);
set(findobj(ax(ww+1),'tag','Box'), 'Color', 'none');
set(findobj(ax(ww+1),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(ww+1),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(ww+1).Children = ax(ww+1).Children([end 1:end-1]);

set(ax,'LineWidth',1,'FontWeight','bold')

end