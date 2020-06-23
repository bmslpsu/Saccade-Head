function [] = SSvsStatic_Saccade_Stats()
%% SSvsStatic_Saccade_Stats:

root = 'H:\DATA\Rigid_Data\';

[SS,SSPATH] = uigetfile({'*.mat'},'Select SS data', root, 'MultiSelect','off');
[Static,StaticPATH] = uigetfile({'*.mat'},'Select static data', root, 'MultiSelect','off');

SS = load(fullfile(SSPATH,SS),'PATH','COUNT','SACCADE_STATS','D','I','U','N');
Static = load(fullfile(StaticPATH,Static),'PATH','COUNT','SACCADE_STATS','D','I','U','N');

%% Boxplots
CC = [0 0.8 0 ; 0 0 0.8];

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches');
fig.Position = [2 2 2 2];

% Saccade Count/Rate
ss_count = cat(1, SS.COUNT {:});
static_count = cat(1, Static.COUNT{:});
count_all = [ss_count ; static_count];
G = [ones(length(ss_count),1) ; 2*ones(length(static_count),1)];

ax(1) = subplot(1,1,1); axis tight
bx = boxplot(count_all./10, G, 'Labels', {'Sine','Static'}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
ylabel('Rate (#/s)')
box off
ylim([-0.1 2])

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1},h{kk,2},CC(kk,:));
end

set(findobj(ax(1),'tag','Median'), 'Color', 'w','LineWidth',1.5);
set(findobj(ax(1),'tag','Box'), 'Color', 'none');
set(findobj(ax(1),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(1),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(1).Children = ax(1).Children([end 1:end-1]);
set(ax,'LineWidth',1,'FontWeight','bold')

end