function [] = Head_Leg_Saccade()
%% Head_Leg:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N')

clearvars -except clms CC Vel PATH COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

%% Saccade distribution
leg_thresh = 0.9;
n_speed = N.vel/2;
n_trial = size(SACCADE,1);
sacd_leg_ratio = nan(n_trial,3);
for n = 1:n_trial
	n_sacd = SACCADE.head_saccade{n}.count;
    vel_idx = SACCADE.vel(n);
    vel_idx(vel_idx > n_speed) = vel_idx(vel_idx > n_speed) - n_speed; % speeds
    if n_sacd > 0
        leg_class = max(SACCADE.leg{n} > leg_thresh,[],2);
        leg_class = medfilt2(leg_class,[200,1]);
        all_times = SACCADE.head_saccade{n}.time;
        sacd_idx = SACCADE.head_saccade{n}.SACD.PeakIdx;
        sacd_locs = false(length(all_times),1);
        sacd_locs(sacd_idx) = true;
        sacd_class = sacd_locs & leg_class;
        sacd_leg_ratio(n,1) = n_sacd; % # of saccades in trial
        sacd_leg_ratio(n,2) = sum(sacd_class); % # of saccades occuring with leg extension
        sacd_leg_ratio(n,3) = sum(leg_class) / length(leg_class); % percent leg movement
        sacd_leg_ratio(n,4) = vel_idx; % speed index
    else
        sacd_leg_ratio(n,1:4) = nan;
    end
   
end
sacd_leg_ratio = sacd_leg_ratio(~isnan(sum(sacd_leg_ratio,2)),:);

% For each speed
sacd_leg_ratio_vel = cell(n_speed,1);
for v = 1:n_speed 
    sacd_leg_ratio_vel{v} = sacd_leg_ratio(sacd_leg_ratio(:,4)==v,:);
end

percet_sacd_leg_ext = 100 * sum(sacd_leg_ratio(:,2)) / sum(sacd_leg_ratio(:,1));

percent_leg_ext = cellfun(@(x) basic_stats(100*x(:,3),1), sacd_leg_ratio_vel, 'UniformOutput', true);
percet_sacd_leg_ext = cellfun(@(x) 100 * sum(x(:,2)) / sum(x(:,1)), sacd_leg_ratio_vel, 'UniformOutput', true);

percent_sacd_trial = sacd_leg_ratio(:,2) ./ sacd_leg_ratio(:,1);

%%
fig = figure (3) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 3 3])
ax(1) = subplot(1,1,1); cla ; hold on ; axis tight
    histogram(100*percent_sacd_trial, 30, 'Normalization', 'Probability', ...
        'FaceColor', 'k', 'FaceAlpha', 1, 'EdgeColor', 'none')
    xlim([-10 100])
    ylim([-0.1 1])
    xlabel('Saccade during Leg Extension Rate (%)')


end