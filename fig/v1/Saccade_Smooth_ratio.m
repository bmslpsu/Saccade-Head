function [] = Saccade_Smooth_ratio()
%% Saccade_Smooth_ratio:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE','D','I','U','N')

%% Saccade Percent
clearvars -except SACCADE HEAD_SACCADE_STATSD I U N
clc

n_speed = N.vel/2;
Speed = U.vel{1}(1:n_speed);
CC = repmat(hsv(n_speed),2,1);

speedI = I.vel;
speedI(speedI > n_speed) = speedI(speedI > n_speed) - n_speed;
Scd_trial = cell(N.fly,n_speed);
%Scd_fly = cell(N.fly,n_speed);
for n = 1:N.file
    trial = SACCADE.head_saccade{n};
    if trial.count > 0
        scd_points = size(trial.saccades_all,1);
    else
       scd_points = 0; 
    end

    scd_percent = 100 * scd_points / trial.n;
    Scd_trial{I.fly(n),speedI(n)}(end+1,1) = scd_percent;
    %Scd_fly{I.fly(n),speedI(n)}(end+1,1) = scd_points;
end

Scd_trial_mean = cellfun(@(x) mean(x,1), Scd_trial);

Scd_vel_all = cell(1,n_speed);
for v = 1:n_speed
    Scd_vel_all{v} = cat(1,Scd_trial{:,v});
end

g = num2cell(1:n_speed);
G = cellfun(@(x,y) y*ones(size(x)), Scd_vel_all, g, 'UniformOutput', false);
G = cat(1,G{:});
Scd_vel_all_g = cat(1,Scd_vel_all{:});


fig = figure(1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2 2])
ww = 1;
ax(ww) = subplot(1,1,1);
    bx = boxplot(Scd_vel_all_g, G, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    %bx = boxplot(Scd_trial_mean, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2)
    xlabel('Stimulus Speed (°/s)')
    ylabel('Saccade (%)')

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1}, h{kk,2}, CC(kk,:), 'EdgeColor', 'none');
    end

    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);

set(ax , 'LineWidth', 0.75, 'FontSize', 8, 'Color', 'none')

end