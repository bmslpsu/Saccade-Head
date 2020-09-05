function [] = Head_Leg_sacd()
%% Head_Leg_sacd:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE','U','N')

%% Saccade distribution
clearvars -except U N SACCADE

leg_thresh = 0.9;
n_speed = N.vel/2;
n_trial = size(SACCADE,1);
sacd_leg = cell(N.fly,n_speed);
vel_idx = SACCADE.vel;
vel_idx(vel_idx > n_speed) = vel_idx(vel_idx > n_speed) - n_speed; % speeds
for n = 1:n_trial
	n_sacd = SACCADE.head_saccade{n}.count;
    if n_sacd > 0
        leg_class = max(SACCADE.leg{n} > leg_thresh,[],2); % classify time series into leg extension / stabilization
        leg_class = medfilt2(leg_class,[100,1]); % filter classification data
        all_times = SACCADE.head_saccade{n}.time; % time vetor
        sacd_idx = SACCADE.head_saccade{n}.SACD.PeakIdx; % saccade indicies
        %sacd_idx = SACCADE.head_saccade{n}.saccades_all.Index;
        sacd_locs = false(length(all_times),1); % saccade class vector
        sacd_locs(sacd_idx) = true; % classify saccades
        sacd_class = sacd_locs & leg_class; % saccades during leg extension
        
        sacd_leg{SACCADE.fly(n),vel_idx(n)}(end+1,1) = length(sacd_idx); % # of saccade points in trial
        pp = size(sacd_leg{SACCADE.fly(n),vel_idx(n)},1);
        
        sacd_leg{SACCADE.fly(n),vel_idx(n)}(pp,2) = sum(sacd_class); % # of saccades occuring with leg extension
        sacd_leg{SACCADE.fly(n),vel_idx(n)}(pp,3) = sum(leg_class); % leg extension data points
        sacd_leg{SACCADE.fly(n),vel_idx(n)}(pp,4) = length(leg_class); % total data points
        sacd_leg{SACCADE.fly(n),vel_idx(n)}(pp,5) = sum(sacd_class) / length(sacd_idx); % percent sacd during leg extension
        sacd_leg{SACCADE.fly(n),vel_idx(n)}(pp,6) = sum(leg_class) / length(leg_class); % percent leg extension per trial
    else
        % no saccades
    end
end
%%
% comb_vel_sacd_leg = cell(N.fly,1);
% for f = 1:N.fly
%    comb_vel_sacd_leg{f} = cat(1,sacd_leg{f,:}); % combine speeds 
% end

comb_fly_sacd_leg = cell(n_speed,1);
for v = 1:n_speed
   comb_fly_sacd_leg{v} = cat(1,sacd_leg{:,v}); % combine speeds 
end

percent_leg_ext = cellfun(@(x) 100*sum(x(:,3))/sum(x(:,4)), sacd_leg, 'UniformOutput', true);
percent_leg_ext_stats = basic_stats(percent_leg_ext,1);
percent_leg_ext_trial = cat(1,comb_fly_sacd_leg{:});

% percent_sacd_leg_ext = cellfun(@(x) 100 * sum(x(:,2)) / sum(x(:,1)), comb_vel_sacd_leg, 'UniformOutput', true);
percent_sacd_leg_ext = cellfun(@(x) 100 * sum(x(:,2)) / sum(x(:,1)), sacd_leg, 'UniformOutput', true);
percent_sacd_leg_ext_stats = basic_stats(percent_sacd_leg_ext,1);

% percent_sacd_leg_ext = cellfun(@(x) 100 * sum(x(:,2)) / sum(x(:,1)), comb_fly_sacd_leg, 'UniformOutput', true);

vel_group = cellfun(@(x,y) y*ones(size(x,1),1), comb_fly_sacd_leg, num2cell(1:n_speed)', ...
    'UniformOutput', false);
vel_group = cat(1, vel_group{:});

%% Percent Leg Extension
fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 3 3])
ax(1) = subplot(1,1,1); cla ; hold on ; axis tight
    speeds = U.vel{1}(1:n_speed);
    CC = hsv(n_speed);
    
%     bx = boxplot(100*percent_leg_ext_trial(:,6), vel_group, 'Labels', {speeds}, ...
%                 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    bx = boxplot(percent_leg_ext, 'Labels', {speeds}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2},CC(kk,:));
    end

    set(findobj(ax(1),'tag','Median'), 'Color', 'k','LineWidth',1.5);
    set(findobj(ax(1),'tag','Box'), 'Color', 'none');
    set(findobj(ax(1),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(1),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(1).Children = ax(1).Children([end 1:end-1]);
    ylim([-5 100])

 	set(ax , 'LineWidth', 1)
    xlabel('Stimulus Speed (°/s)')
    ylabel('Leg Extension (%)')
    xticks(speeds)

%% Percent Leg Extension Saccades
fig = figure (2) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 2 2])
ax(1) = subplot(1,1,1); cla ; hold on ; axis tight   
    bx = boxplot(percent_sacd_leg_ext, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    %bx = boxplot(100*percent_leg_ext_trial(:,5), vel_group, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    
    ylabel('Saccades During Leg Extension (%)')

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2}, [0.5 0.5 0.5]);
    end

    set(findobj(ax(1),'tag','Median'), 'Color', 'k','LineWidth',1.5);
    set(findobj(ax(1),'tag','Box'), 'Color', 'none');
    set(findobj(ax(1),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(1),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(1).Children = ax(1).Children([end 1:end-1]);
    ylim([-5 100])

set(ax , 'LineWidth', 1)

end