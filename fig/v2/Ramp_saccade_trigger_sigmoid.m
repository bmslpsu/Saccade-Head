function [] = Ramp_saccade_trigger_sigmoid()
%% Ramp_saccade_stats:
load('H:\DATA\Rigid_Data\Saccade\combined\Ramp_All_Stats.mat')

%% Saccade Statistics by velocity
clc
All_Stats = All_Stats(~isnan(All_Stats.StartPos),:);

all_pos = [ cat(2, Data(1).pos{:}) , cat(2, Data(3).pos{:}) , cat(2, Data(3).pos{:}) ];

Vel = unique(All_Stats.Vel);
n_speed = length(Vel)/2;
Speed = Vel(n_speed+1:end);
CC = repmat(hsv(n_speed),2,1);

vel_group_all = All_Stats.vel;
vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;

clear ax h H
fig = figure (1);
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 3.5])
h = gobjects(n_speed,1);
ax(1) = subplot(2,1,1); cla; hold on
    edges = -25:25;
    starts = -1 * All_Stats.StartPos .* All_Stats.Direction;
%     H = histogram(all_pos(:), edges, 'Normalization', 'pdf', ...
%         'FaceColor', 'r', 'FaceAlpha', 1, 'EdgeColor', 'none');
    H = histogram(starts, edges, 'Normalization', 'pdf', ...
        'FaceColor', 'g', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
    xline(median(starts), 'k', 'LineWidth', 1)
    xline(0, 'k', 'LineWidth', 1)
    ylabel('Probability')
    
%     for v = 1:n_speed
%         h(v) = histogram(starts(vel_group_all==v), edges, 'Normalization', 'pdf', ...
%             'FaceColor', CC(v,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%         xline(median(starts(vel_group_all==v)), 'Color', CC(v,:), 'LineWidth', 1)
%     end
    
ax(2) = subplot(2,1,2); cla; hold on
    [h_cdf,stats] = cdfplot(starts);
    set(h_cdf, 'Color', 'k', 'LineWidth', 1)
    title('')
    xlabel('Trigger Position (°)')
    ylabel('Probability')
    xline(0, 'k', 'LineWidth', 1)
    
    [fitresult, gof] = sigmoid_fit(h_cdf.XData, h_cdf.YData);
    
 	xfit = -11:0.1:25;
    yfit = 1 ./ (1 + fitresult.a*exp(fitresult.b*xfit));
    h_fit = plot(xfit, yfit, 'r', 'LineWidth', 1);
    
    legend([h_cdf h_fit], 'CDF', 'Sigmoid Fit', 'Box', 'off')
    
set(ax, 'LineWidth', 1, 'Box', 'off', 'XGrid', 'off', 'YGrid', 'off')
linkaxes(ax, 'x')
set(ax, 'XTick', -25:5:25)

end