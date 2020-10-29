function [] = Passive_vs_Saccade_time_constant()
%% Passive_vs_Saccade_time_constant:
root = 'H:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','off');
Passive = load(fullfile(PATH,FILE),'U','N','DATA');

Ramp = load("H:\DATA\Rigid_Data\Saccade\processed\Ramp_time_constant_wave=30.mat");

%% Make table
ramp_table = Ramp.Saccade_Table(Ramp.Saccade_Table.r2 > 0.8, 4:5);
A = [table(ones(size(Passive.DATA,1),1),'VariableNames', {'class'}) , Passive.DATA(:,[12,11])];
B = [table(2*ones(size(ramp_table,1),1),'VariableNames', {'class'}) , ramp_table];
C = [A ; B];

%% Anova
[p,tb,stats] = anova1(C.tau, C.class, 'varnames', {'passive', 'active'});
% [h, p] = ttest2(A.tau, B.tau)
[c,m] = multcompare(stats);

%% Time constant
fig = figure (2) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', 1*[2 2 2 2])
clear ax
ax(1) = subplot(1,1,1); hold on
ylabel('Time Constant (ms)')
bx = boxplot(C.tau, C.class, 'Labels', {'Passive', 'Saccade'}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);

h = get(bx(5,:),{'XData','YData'});
cc = {'k', 'r'};
for kk = 1:size(h,1)
   patch(h{kk,1}, h{kk,2}, cc{kk}, 'EdgeColor', 'none');
end
ww = 1;
set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth', 1);
set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(ww).Children = ax(ww).Children([end 1:end-1]);
ylim([0 60])
set(ax , 'LineWidth', 1, 'Box', 'off', 'XColor', 'none')

%% R^2
fig = figure (2) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', 1*[2 2 2 2])
clear ax
ax(1) = subplot(1,1,1); hold on
ylabel('Time Constant (ms)')
bx = boxplot(C.r2, C.class, 'Labels', {'Passive', 'Saccade'}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);

h = get(bx(5,:),{'XData','YData'});
cc = {'k', 'r'};
for kk = 1:size(h,1)
   patch(h{kk,1}, h{kk,2}, cc{kk}, 'EdgeColor', 'none');
end
ww = 1;
set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth', 1);
set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(ww).Children = ax(ww).Children([end 1:end-1]);
ylim([0 1])
set(ax , 'LineWidth', 1, 'Box', 'off', 'XColor', 'none')

end