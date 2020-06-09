function [] = Head_Leg()
%% Head_Leg:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N')

clearvars -except clms CC Vel PATH COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

%% Land %%
Vel = U.vel{1};
n_vel = length(Vel);
CC = repmat(hsv(n_vel/2),2,1);

vel_idx = SACCADE.vel;
% head_data = cellfun(@(x) x.position, SACCADE.head_saccade, 'UniformOutput', false);
head_data = cellfun(@(x) x.velocity, SACCADE.head_saccade, 'UniformOutput', false);
% head_data = SACCADE.dWBA;
head_data = cat(2,head_data{:});
leg_data = cat(2,SACCADE.leg{:});

stable_head = head_data;
land_head = head_data;
stable_head(leg_data) = nan;
land_head(~leg_data) = nan;

stable_head_vel = cell(n_vel,1);
land_head_vel = cell(n_vel,1);
for n = 1:n_vel
    if n > n_vel/2
        flip = -1;
    else
        flip = 1;
    end
    stable_head_vel{n} = flip*stable_head(:,vel_idx == n);
    land_head_vel{n} = flip*land_head(:,vel_idx == n);
end
stable_head_vel = reshape(stable_head_vel, n_vel/2, 2);
land_head_vel = reshape(land_head_vel, n_vel/2, 2);

stable_head_vel = cellfun(@(x,y) [x y], stable_head_vel(:,1), stable_head_vel(:,2), 'UniformOutput', false);
land_head_vel = cellfun(@(x,y) [x y], land_head_vel(:,1), land_head_vel(:,2), 'UniformOutput', false);

stable_head_vel = cellfun(@(x) x(~isnan(x)), stable_head_vel, 'UniformOutput', false);
land_head_vel = cellfun(@(x) x(~isnan(x)), land_head_vel, 'UniformOutput', false);

land_ratio = cellfun(@(x,y) length(y) / (length(x) + length(y)), stable_head_vel, ...
                            land_head_vel, 'UniformOutput', true);
                        
stable_head_all = cat(1,land_head_vel{:});
land_head_all = cat(1,stable_head_vel{:});

%% Speeds
fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 3 8])
ax = gobjects(N.vel/2,1);
bins = 0:5:300;
% bins = -25:0.5:25;
pp = 1;
for v = 1:n_vel/2
    ax(v) = subplot(N.vel/2,1,pp); hold on ; title([ num2str(Vel(v)) '(�/s)'])
        h_stable = histogram(stable_head_vel{v}, bins, 'Normalization', 'probability', ...
            'FaceColor', CC(v,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        h_land = histogram(land_head_vel{v}, bins, 'Normalization', 'probability', ...
            'FaceColor', 'k', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        
        % xlabel('Position (�)')
        xlabel('Velocity (�/s)')
        
        ylabel('Probability')
        axis tight
	pp = pp + 1;
end
linkaxes(ax,'xy')
set(ax, 'LineWidth', 1.5, 'Box', 'on', 'YLim', [-0.005 ax(1).YLim(2)])

%% All
fig = figure (2) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 3 3])
bins = 0:5:300;
% bins = 0:0.5:25;
ax = subplot(1,1,1); hold on
    h_stable = histogram(abs(stable_head_all), bins, 'Normalization', 'probability', ...
        'FaceColor', 'k', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    h_land = histogram(abs(land_head_all), bins, 'Normalization', 'probability', ...
        'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none');

    xlabel('Position (�)')
    xlabel('Velocity (�/s)')
    
    ylabel('Probability')
    axis tight
    
    legend('Landing State', 'Stabilizing State')

set(ax, 'LineWidth', 1.5, 'Box', 'on', 'YLim', [-0.005 1.1*ax(1).YLim(2)])



end