function [] = Ramp_saccade_head_leg_new()
%% Ramp_saccade_head_leg_new:
root = 'E:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'SACCADE','U','N')

%% Saccade during leg extension
clearvars -except U N SACCADE
clc
leg_thresh = 0.9;
n_speed = N.vel/2;
n_trial = size(SACCADE,1);
velI = SACCADE.vel;
velI(velI > n_speed) = velI(velI > n_speed) - n_speed;
[vel_fly_group, vel_group, fly_group] = findgroups(velI, SACCADE.fly);

DATA = splitvars(table(nan(n_trial,7)));
DATA.Properties.VariableNames = {'leg_percent', 'count', 'count_leg', ...
    'sacd_leg_percent', 'sacd_leg_norm', 'sacd_rate', 'sacd_leg_rate'};
DATA = [SACCADE(:,1:4), DATA];
for n = 1:n_trial    
   	leg_class = max(SACCADE.leg{n} > leg_thresh, [], 2); % classify time series into leg extension / stabilization
  	leg_class = medfilt2(leg_class,[200,1]); % filter classification data
    DATA.leg_percent(n) = 100 * sum(leg_class) / SACCADE.head_saccade{n}.n; % percent of time fly is extending legs
    
   	n_sacd = SACCADE.head_saccade{n}.count;
    DATA.count(n) = n_sacd; % total # of saccades
    if n_sacd > 0
        sacd_idx = SACCADE.head_saccade{n}.SACD.PeakIdx; % saccade indicies
        %sacd_idx = SACCADE.head_saccade{n}.saccades_all.Index;
        sacd_locs = false(SACCADE.head_saccade{n}.n, 1); % saccade class vector
        sacd_locs(sacd_idx) = true; % classify saccades
        sacd_class = sacd_locs & leg_class; % saccades during leg extension
        DATA.count_leg(n) = sum(sacd_class);
        DATA.sacd_leg_percent(n) = 100 * DATA.count_leg(n) / DATA.count(n);
    else % no saccades
        DATA.count_leg(n) = 0;
    end
    
    DATA.sacd_rate(n) = (DATA.count(n) - DATA.count_leg(n)) / (10*(1-DATA.leg_percent(n)/100));
    DATA.sacd_leg_rate(n) = DATA.count_leg(n) / (10*DATA.leg_percent(n)/100);
    
    if DATA.leg_percent(n) == 0
        DATA.sacd_leg_norm(n) = nan;
    else
        DATA.sacd_leg_norm(n) = DATA.sacd_leg_percent(n) / DATA.leg_percent(n);
    end
    
end

%%
clear ALL
ALL.leg_percent = splitapply(@(x) {cat(2,x)}, DATA.leg_percent, vel_fly_group);
ALL.sacd_leg_percent = splitapply(@(x) {cat(2,x)}, DATA.sacd_leg_percent, vel_fly_group);
ALL.sacd_rate = splitapply(@(x) {cat(2,x)}, DATA.sacd_rate, vel_fly_group);
ALL.sacd_leg_rate = splitapply(@(x) {cat(2,x)}, DATA.sacd_leg_rate, vel_fly_group);

ALL = structfun(@(x) splitapply(@(y) {y}, x, vel_group), ALL, 'UniformOutput', false);
ALL = structfun(@(x) cat(2,x{:}), ALL, 'UniformOutput', false);

ALL.fly_stats = structfun(@(x) cellfun(@(y) basic_stats(y,1), x, 'UniformOutput', true), ...
    ALL, 'UniformOutput', false);

fnames = string(fieldnames(ALL.fly_stats));
n_field = length(fnames);
for f = 1:n_field
    for v = 1:n_speed
        ALL.fly_mean.(fnames(f)){v} = cat(1, ALL.fly_stats.(fnames(f))(:,v).mean);
        ALL.vel_trial.(fnames(f)){v} = cat(1, ALL.(fnames(f)){:,v});
    end
    ALL.fly_mean_all.(fnames(f)) = cat(2, ALL.fly_mean.(fnames(f)){:});
    ALL.vel_trial_all.(fnames(f)) = cat(1, ALL.vel_trial.(fnames(f)){:});
%     ALL.vel_trial_group.(fnames(f)) = cellfun(@(x,y) y*ones(size(x)), ALL.vel_trial.(fnames(f)), ...
%         num2cell(1:n_speed), 'UniformOutput', false);
%     ALL.vel_trial_group.(fnames(f)) = cat(1, ALL.vel_trial_group.(fnames(f)){:});
end

ALL.vel_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    ALL.fly_mean, 'UniformOutput', false);

%% Percent Leg Extension by trial
fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 2 2])
ax(1) = subplot(1,1,1); cla ; hold on ; axis tight
    speeds = U.vel{1}(1:n_speed);
    CC = hsv(n_speed);
    
    bx = boxplot(DATA.leg_percent, velI, 'Labels', {speeds}, ...
        'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    
    xlabel('Stimulus Speed (°/s)')
    ylabel('Percentage Leg Extension')

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2},CC(kk,:), 'EdgeColor', 'none');
    end

    set(findobj(ax(1),'tag','Median'), 'Color', 'k', 'LineWidth', 1.5);
    set(findobj(ax(1),'tag','Box'), 'Color', 'none');
    set(findobj(ax(1),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(1),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(1).Children = ax(1).Children([end 1:end-1]);
    ylim([-5 100])

    set(ax , 'LineWidth', 1, 'Box', 'off')

%% Percent Leg Extension by fly
fig = figure (2) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 3 3])
ax(1) = subplot(1,1,1); cla ; hold on ; axis tight
    speeds = U.vel{1}(1:n_speed);
    CC = hsv(n_speed);
    
    set(ax , 'LineWidth', 1, 'YLim', [0 100])
    xlabel('Stimulus Speed (°/s)')
    ylabel('Leg Extension Percent')
    xticks(speeds)
    
    bx = boxplot(ALL.fly_mean_all.leg_percent, 'Labels', {speeds}, ...
        'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    xlabel('Stimulus Speed (°/s)')
    ylabel('Percentage Leg Extension')

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2},CC(kk,:), 'EdgeColor', 'none');
    end

    set(findobj(ax(1),'tag','Median'), 'Color', 'k', 'LineWidth', 1.5);
    set(findobj(ax(1),'tag','Box'), 'Color', 'none');
    set(findobj(ax(1),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(1),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(1).Children = ax(1).Children([end 1:end-1]);
    ylim([-5 100])
    
set(ax , 'LineWidth', 1, 'Box', 'off')

%% Saccade rate land vs fly
fig = figure (3) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 2 2])
ax(1) = subplot(1,1,1); cla ; hold on ; axis tight   
    %bx = boxplot(ALL.fly_mean_all.sacd_rate, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    %bx = boxplot(ALL.fly_mean_all.sacd_leg_rate, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    data = [ALL.fly_mean_all.sacd_rate(:), ALL.fly_mean_all.sacd_leg_rate(:)];
    bx = boxplot(data, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    
    ylabel('Saccades frequency (Hz)')

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2}, 'k', 'EdgeColor', 'none');
    end

    set(findobj(ax(1),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(1),'tag','Box'), 'Color', 'none');
    set(findobj(ax(1),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(1),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(1).Children = ax(1).Children([end 1:end-1]);
    ylim([-0.1 2])

set(ax , 'LineWidth', 1)

end