function [] = Head_Leg_Saccade_v2()
%% Head_Leg_Saccade_v2:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N')

clearvars -except clms CC Vel PATH COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

%% Saccade distribution
leg_thresh = 0.9;
n_speed = N.vel/2;
n_trial = size(SACCADE,1);
sacd_leg_ratio = cell(N.fly,n_speed);
vel_idx = SACCADE.vel;
vel_idx(vel_idx > n_speed) = vel_idx(vel_idx > n_speed) - n_speed; % speeds
for n = 1:n_trial
	n_sacd = SACCADE.head_saccade{n}.count;
    if n_sacd > 0
        leg_class = max(SACCADE.leg{n} > leg_thresh,[],2); % classify time series into leg extension / stabilization
        leg_class = medfilt2(leg_class,[200,1]); % filter classification data
        all_times = SACCADE.head_saccade{n}.time; % time vetor
        % sacd_idx = SACCADE.head_saccade{n}.SACD.PeakIdx; % saccade indicies
        sacd_idx = SACCADE.head_saccade{n}.saccades_all.Index;
        sacd_locs = false(length(all_times),1); % saccade class vector
        sacd_locs(sacd_idx) = true; % classify saccades
        sacd_class = sacd_locs & leg_class; % saccades during leg extension
        sacd_leg_ratio{SACCADE.fly(n),vel_idx(n)}(end+1,1) = length(sacd_idx); % # of saccades in trial
        
        pp = size(sacd_leg_ratio{SACCADE.fly(n),vel_idx(n)},1);
        sacd_leg_ratio{SACCADE.fly(n),vel_idx(n)}(pp,2) = sum(sacd_class); % # of saccades occuring with leg extension
        sacd_leg_ratio{SACCADE.fly(n),vel_idx(n)}(pp,3) = sum(leg_class); % leg extension data points
        sacd_leg_ratio{SACCADE.fly(n),vel_idx(n)}(pp,4) = length(leg_class); % total data points
    else
        % no saccades
    end
end

comb_vel_sacd_leg_percent = cell(N.fly,1);
for f = 1:N.fly
   comb_vel_sacd_leg_percent{f} = cat(1,sacd_leg_ratio{f,:}); % combine speeds 
end

percent_leg_ext = cellfun(@(x) 100*sum(x(:,3))/sum(x(:,4)), sacd_leg_ratio, 'UniformOutput', true);
percent_leg_ext_stats = basic_stats(percent_leg_ext,1);

percet_sacd_leg_ext = cellfun(@(x) 100 * sum(x(:,2)) / sum(x(:,1)), comb_vel_sacd_leg_percent, 'UniformOutput', true);
percet_sacd_leg_ext_stats = basic_stats(percet_sacd_leg_ext,1);

%% Percent Leg Extension
fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 3 3])
ax(1) = subplot(1,1,1); cla ; hold on ; axis tight
    speeds = U.vel{1}(1:n_speed);
    CC = hsv(n_speed);
    
    [h1,h2] = PlotPatch(percent_leg_ext_stats.mean, percent_leg_ext_stats.std, speeds, 1, 1, ...
        'k', [0.5 0.5 0.5], 0.5, 2);
    scatter(speeds, percent_leg_ext_stats.mean, 200, CC, '.');    
    
    set(ax , 'LineWidth', 1, 'YLim', [0 100])
    xlabel('Stimulus Speed (�/s)')
    ylabel('Leg Extension Percent')
    xticks(speeds)
    
%     bx = boxplot(percent_leg_ext, 'Labels', {speeds}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
%     xlabel('Stimulus Speed (�/s)')
%     ylabel('Percentage Leg Extension')
% 
%     h = get(bx(5,:),{'XData','YData'});
%     for kk = 1:size(h,1)
%        patch(h{kk,1},h{kk,2},CC(kk,:));
%     end
% 
%     set(findobj(ax(1),'tag','Median'), 'Color', 'w','LineWidth',1.5);
%     set(findobj(ax(1),'tag','Box'), 'Color', 'none');
%     set(findobj(ax(1),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
%     set(findobj(ax(1),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
%     ax(1).Children = ax(1).Children([end 1:end-1]);
%     %ax(1).YLim(1) = -0.1;
%     %ax(1).YLim(2) = ylim_list(1);


%% Percent Leg Extension Saccades
fig = figure (2) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 2 2])
ax(1) = subplot(1,1,1); cla ; hold on ; axis tight   
    bx = boxplot(percet_sacd_leg_ext, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    ylabel('Saccades During Leg Extension (%)')

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2}, [0.5 0.5 0.5]);
    end

    set(findobj(ax(1),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(1),'tag','Box'), 'Color', 'none');
    set(findobj(ax(1),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(1),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(1).Children = ax(1).Children([end 1:end-1]);

set(ax , 'LineWidth', 1)

end