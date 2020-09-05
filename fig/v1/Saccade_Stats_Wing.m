function [] = Saccade_Stats_Wing()
%% Saccade_Stats_Wing:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE','WING_SACCADE_STATS','U','N')

%% Saccade Count/Rate
keepI = cellfun(@(x) isstruct(x) | isobject(x), SACCADE.wing_saccade);
Saccade = SACCADE(keepI,:);
n_speed = N.vel/2;

fly_group = Saccade.fly;
vel_group_all = Saccade.vel;

vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
[fly_vel_group,fly_group,vel_group] = findgroups(fly_group, vel_group_all);

count = cellfun(@(x) x.count, Saccade.wing_saccade, 'UniformOutput', true);
count_vel_fly_mean = splitapply(@(x) nanmean(x,1), count, fly_vel_group);

fig = figure (13) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2 2])
movegui(fig, 'center')

ww = 1;
clear ax
ax(1) = subplot(1,1,1); hold on
    bx = boxplot(count./10, vel_group_all, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    ylabel('Frequency (Hz)')
    ylim([-0.1 3])

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2},'k');
    end

    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);

set(ax, 'LineWidth', 1, 'Box', 'off')

end