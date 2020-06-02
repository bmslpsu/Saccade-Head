function [] = Interval_KMeans()
%% Interval_KMeans:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N')

clms = N.vel/2;
CC = repmat(hsv(clms),2,1);

Vel = U.vel{1};
Speed = Vel(1:N.vel/2);

clearvars -except clms CC Vel PATH COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

%%
Ts = SACCADE.saccade{1}.Ts; % sampling time [s]

% Get rid of intervals before 1st saccade
valid_time = cell(N.vel,1);
length_idx = cell(N.vel,1);
int_times = cell(N.vel,1);
int_stats = cell(N.vel,1);
norm_interval = GRAND.norm_interval;
normstart_interval = GRAND.normstart_interval;
error = GRAND.error;
int_error = GRAND.int_error;
for jj = 1:N.vel
    allnan_idx = isnan(norm_interval(jj).time(1,:));
    field_name = string(fieldnames(norm_interval));
    for f = 1:length(field_name)
        norm_interval(jj).(field_name(f)) = norm_interval(jj).(field_name(f))(:,~allnan_idx);
        normstart_interval(jj).(field_name(f)) = normstart_interval(jj).(field_name(f))(:,~allnan_idx);
        error(jj).(field_name(f)) = error(jj).(field_name(f))(:,~allnan_idx);
        int_error(jj).(field_name(f)) = int_error(jj).(field_name(f))(:,~allnan_idx);
    end
    
    % Sort intervals
    valid_time{jj} = ~isnan(norm_interval(jj).time); % valid times in intervals
    length_idx{jj} = sum(valid_time{jj},1); % how long intervals are in samples
    [length_idx{jj},order] = sort(length_idx{jj},'ascend'); % sort by lengths by length
    max_length = max(length_idx{jj});
    valid_time{jj} = valid_time{jj}(1:max_length,order); % valid times by length
    
 	% Sort intervals by length & get rid of trailing nan's
    for f = 1:length(field_name)
        norm_interval(jj).(field_name(f)) = ...
            norm_interval(jj).(field_name(f))(1:max_length,order);
        normstart_interval(jj).(field_name(f)) = ...
            normstart_interval(jj).(field_name(f))(1:max_length,order);
      	error(jj).(field_name(f)) = ...
            error(jj).(field_name(f))(1:max_length,order);
        int_error(jj).(field_name(f)) = ...
            int_error(jj).(field_name(f))(1:max_length,order);
    end

    int_times{jj} = length_idx{jj} * Ts;
    int_stats{jj} = basic_stats(int_times{jj},2);
end

% Cut off based on times
head_amp = 15;
head_gain = 0.4;
cut_time = head_amp ./ (head_gain*abs(Vel));

norm_interval_cut = norm_interval;
normstart_interval_cut = normstart_interval;
error_cut = normstart_interval;
int_error_cut = int_error;
for jj = 1:N.vel
    cut_idx = find(norm_interval(jj).time(:,end) <= cut_time(jj), 1, 'last');
    rmv_idx = sum(valid_time{jj},1) > cut_idx;
    
 	% Remove intervals longer than cut-time
    for f = 1:length(field_name)
        norm_interval_cut(jj).(field_name(f)) = ...
            norm_interval(jj).(field_name(f))(:,~rmv_idx);
        normstart_interval_cut(jj).(field_name(f)) = ...
            normstart_interval(jj).(field_name(f))(:,~rmv_idx);
      	error_cut(jj).(field_name(f)) = ...
            error(jj).(field_name(f))(:,~rmv_idx);
        int_error_cut(jj).(field_name(f)) = ...
            int_error(jj).(field_name(f))(:,~rmv_idx);
    end
    
    valid_time{jj} = ~isnan(norm_interval_cut(jj).time);
 	max_length = max(sum(valid_time{jj},1));
    
  	% Get rid of trailing nan's
    for f = 1:length(field_name)
        norm_interval_cut(jj).(field_name(f)) = ...
            norm_interval_cut(jj).(field_name(f))(1:max_length,:);
        normstart_interval_cut(jj).(field_name(f)) = ...
            normstart_interval_cut(jj).(field_name(f))(1:max_length,:);
      	error_cut(jj).(field_name(f)) = ...
            error_cut(jj).(field_name(f))(1:max_length,:);
        int_error_cut(jj).(field_name(f)) = ...
            int_error_cut(jj).(field_name(f))(1:max_length,:);
    end
end

%%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = 2.0*[2 2 clms*(4/3) 3];
FIG.Name = 'Interval Time Visualization';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel,1);
clear h

for jj = 1:N.vel/2
    ax(jj) = subplot(2,clms,jj) ; hold on
    title([num2str(Vel(jj)) ' (°/s)']);
        time = normstart_interval(jj).time;
        pos = normstart_interval(jj).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);

        time = normstart_interval(jj + N.vel/2).time;
        pos = normstart_interval(jj + N.vel/2).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);

        plot([cut_time(jj) cut_time(jj)],40*[-1 1],'-','Color',[0.5 0.5 0.5],'LineWidth', 1);
    
    ax(jj + N.vel/2) = subplot(2,clms,jj + N.vel/2) ; hold on
        time = normstart_interval_cut(jj).time;
        pos = normstart_interval_cut(jj).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);

        time = normstart_interval_cut(jj + N.vel/2).time;
        pos = normstart_interval_cut(jj + N.vel/2).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
        ax(jj + N.vel/2).XLim = [0 cut_time(jj)];
end
linkaxes(ax,'y')
set(ax,'LineWidth',1.5,'YLim',30*[-1 1])
set(ax(1:clms),'XLim',[0 1.7])
%% Visualize interval lengths by velocity
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = 1.0*[2 2 clms*(4/3) 3];
FIG.Name = 'Interval Time Visualization: Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel,1);
clear h
colormap(jet)

Groups = cell(N.vel,1);
K = 2;
for jj = 1:N.vel
    ax(jj) = subplot(ceil(N.vel/clms),clms,jj); 
        imagesc(valid_time{jj}) ; hold on ; title([num2str(Vel(jj)) ' (°/s)'])
        n_trial = size(valid_time{jj},2);
        % cut_time = int_stats{jj}.median + 1*int_stats{jj}.std;
        
        % Perform k-means clustering to seperate groups
        train_data = [int_times{jj}]';
        [Groups{jj}] = kmeans(train_data, 2, 'Distance', 'cityblock');
        
        % Ensure cluster labels start with "1"
        if Groups{jj}(1)==2
            k1 = Groups{jj}==1;
            k2 = Groups{jj}==2;
            Groups{jj}(k1) = 2;
            Groups{jj}(k2) = 1;
        end
        
     	% Find the transition cutoff
        cut_trial = find(Groups{jj} == K, 1, 'first');
        cut_time = int_times{jj}(cut_trial);

        cut_line1 = repmat(cut_time / Ts, 1, n_trial);
        cut_line2 = repmat(cut_trial, 2, 1);
        
        plot(1:n_trial, cut_line1, 'g', 'LineWidth', 1)
        plot(cut_line2, [1 , max(length_idx{jj})], 'y', 'LineWidth', 1)
        
        ax(jj).YTick = 1:200:2000;
        ax(jj).YTickLabels = string((ax(jj).YTick -1 )*Ts); % convert lables to time
end
linkaxes(ax,'y')
set(ax,'XTick',[1, 100:100:500])
YLabelHC = get(ax([1,clms+1]), 'YLabel');
set([YLabelHC{:}], 'String', 'Time (s)')
XLabelHC = get(ax(clms+1:2*clms), 'XLabel');
set([XLabelHC{:}], 'String', 'Trial #')

set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','k',...
    'YColor','k','XColor','k')

%% Visualize interval lengths by speed
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = 1.5*[2 2 clms*(4/3) 3/2];
FIG.Name = 'Interval Time Visualization: Speed';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(clms,1);
clear h
colormap(jet)

valid_time_speed = cell(N.vel/2,1);
length_idx_speed = cell(N.vel/2,1);
int_times_speed = cell(N.vel/2,1);
Groups = cell(N.vel/2,1);
K = 2;
for jj = 1:N.vel/2
    valid_time_speed{jj} = nancat(2,valid_time{jj}, valid_time{jj + N.vel/2});  % both directions
   	length_idx_speed{jj} = nansum(valid_time_speed{jj},1); % how long intervals are in samples
    [length_idx_speed{jj},order] = sort(length_idx_speed{jj},'ascend'); % sort by lengths by length
 	max_length = max(length_idx{jj});
    valid_time_speed{jj} = valid_time_speed{jj}(1:max_length,order); % valid times by length
    int_times_speed{jj} = [int_times{jj}, int_times{jj + N.vel/2}];  % both directions
    int_times_speed{jj} = int_times_speed{jj}(order);
    
    ax(jj) = subplot(1,clms,jj); 
        imagesc(valid_time_speed{jj}) ; hold on ; title([num2str(Vel(jj)) ' (°/s)'])
        n_trial = size(valid_time_speed{jj},2);
        % cut_time = int_stats{jj}.median + 1*int_stats{jj}.std;
        
        % Perform k-means clustering to seperate groups
        train_data = int_times_speed{jj}';
        [Groups{jj}] = kmeans(train_data, 2, 'Distance', 'cityblock');
        
        % Ensure cluster labels start with "1"
        if Groups{jj}(1)==2
            k1 = Groups{jj}==1;
            k2 = Groups{jj}==2;
            Groups{jj}(k1) = 2;
            Groups{jj}(k2) = 1;
        end
        
     	% Find the transition cutoff
        cut_trial = find(Groups{jj} == K, 1, 'first');
        cut_time = int_times_speed{jj}(cut_trial);
        
        cut_line1 = repmat(cut_time / Ts, 1, n_trial);
        cut_line2 = repmat(cut_trial, 2, 1);
        
        plot(1:n_trial, cut_line1, 'g', 'LineWidth', 1)
        plot(cut_line2, [1 , max(length_idx{jj})], 'y', 'LineWidth', 1)
        
        ax(jj).XTick = [1 n_trial];
        ax(jj).YTick = 1:200:2000;
        ax(jj).YTickLabels = string((ax(jj).YTick -1 )*Ts); % convert lables to time
end
linkaxes(ax,'y')
YLabelHC = get(ax([1,2]), 'YLabel');
set([YLabelHC{:}], 'String', 'Time (s)')
% XLabelHC = get(ax(clms+1:2*clms), 'XLabel');
% set([XLabelHC{:}], 'String', 'Trial #')

set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','k',...
    'YColor','k','XColor','k')


%% Visualize interval lengths all
FIG = figure (3) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 3 7];
FIG.Name = 'Interval Time Visualization: All';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear h ax
colormap(jet)

speed_train = 1:5;

valid_time_all = nancat(2,valid_time{speed_train});  % both directions
length_idx_all= nansum(valid_time_all,1); % how long intervals are in samples
[length_idx_all,order] = sort(length_idx_all,'ascend'); % sort by lengths by length
max_length = max(length_idx_all);
valid_time_all = valid_time_all(1:max_length,order); % valid times by length

int_times_all = cat(2,int_times{speed_train});  % both directions
int_times_all = int_times_all(order);
n_trial_speed = cellfun(@(x) length(x), int_times_speed);

Groups = [];
K = 2;
ax(1) = subplot(3,1,1);
    imagesc(valid_time_all) ; hold on
    n_trial = size(valid_time_all,2);

    % Perform k-means clustering to seperate groups
    train_data = int_times_all';
    [Groups] = kmeans(train_data, 2, 'Distance', 'cityblock');

    % Ensure cluster labels start with "1"
    if Groups(1)==2
        k1 = Groups==1;
        k2 = Groups==2;
        Groups(k1) = 2;
        Groups(k2) = 1;
    end

    % Find the transition cutoff
    cut_trial = find(Groups == K, 1, 'first');
    cut_time = int_times_all(cut_trial);

    cut_line1 = repmat(cut_time / Ts, 1, n_trial);
    cut_line2 = repmat(cut_trial, 2, 1);

    plot(1:n_trial, cut_line1, 'g', 'LineWidth', 1)
    plot(cut_line2, [1 , max(length_idx_all)], 'y', 'LineWidth', 1)

    ax(1).YTick = 1:200:2000;
    ax(1).YTickLabels = string((ax(1).YTick -1 )*Ts); % convert lables to time
    ylabel('Time (s)')
    
ax(2) = subplot(3,1,2) ; hold on ; axis tight ; title('All speeds')
    for jj = 1:N.vel/2
        trial_norm = linspace(0,100,n_trial_speed(jj));
        plot(trial_norm,int_times_speed{jj}, '.', 'Color', CC(jj,:), 'LineWidth', 1)
    end
    plot([0 ax(2).XLim(2)], [cut_time cut_time], '--', 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
    ylabel('Time (s)')
    xlabel('Percent of Saccades')

ax(3) = subplot(3,1,3) ; hold on ; axis tight ; title('Combined speeds')
    n_trial_all = length(int_times_all);
    tria_norm  = linspace(0,100,n_trial_all);
 	plot(tria_norm, int_times_all, '.', 'Color', 'k', 'LineWidth', 2)
    plot([0 ax(3).XLim(2)], [cut_time cut_time], '--', 'LineWidth', 1, 'Color', [0.5 0.5 0.5])
    ylabel('Time (s)')
    xlabel('Percent of Saccades')
    
set(ax,'LineWidth',1,'FontWeight','normal','FontSize',8,'Color','w',...
    'YColor','k','XColor','k')

set(ax(2:3),'YTick',0:10)
linkaxes(ax(2:3),'xy')

%% Interval Position %%
FIG = figure (10) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 clms*(4/3) 3];
FIG.Name = 'Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N{1,3},1);
clear ax h

med_span = 5;
for jj = 1:N.vel
    ax(jj) = subplot(ceil(N{1,3}/clms),clms,jj) ; hold on
    ax(jj).Title.String = [num2str(Vel(jj)) ' (°/s)'];
    ax(jj).Title.Color = 'k';
    
    med_time = GRAND.norm_interval(jj).time;
    med_time = nanmedian(med_time,2);
    vel = GRAND.norm_interval(jj).position;
    med_vel = nanmedian(vel,2);
    std_vel = nanstd(vel,[],2);
    
    
    h.trial = plot(med_time, vel, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
                                        
%     h.patch = PlotPatch(med_vel, std_vel, med_time, 1, 1, ...
%         CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
%     h.patch.EdgeColor = CC(jj,:);
end
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',[0 10],'YLim',20*[-1 1])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax([1,clms+1]), 'YLabel');
set([YLabelHC{:}], 'String', ['Head Position (' char(176) ')'])
linkaxes(ax(1:clms),'y')
linkaxes(ax((clms+1):N{1,3}),'y')
linkaxes(ax,'x')

end