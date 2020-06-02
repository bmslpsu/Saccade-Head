function [] = Interval_Cluster()
%% Interval_Cluster:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N')

clms = N.vel/2;
CC = repmat(hsv(clms),2,1);

Vel = U.vel{1};
Speed = Vel(1:N.vel/2);

%% Cluster saccades
clearvars -except clms CC Speed Vel PATH COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

Ts = SACCADE.saccade{1}.Ts; % sampling time [s]

IntAll = GRAND;
IntSort = IntAll; % grouped saccade & intervals
main_fields = string(fieldnames(IntSort)); % all fields
% IntSort = rmfield(IntSort,main_fields(1)); % only intervals
% main_fields = main_fields(2:end); % only intervals
n_main_fields = length(main_fields); % # of fields
field_time = main_fields(1);

% Get rid of intervals before 1st saccade
valid_time = cell(N.vel,1);
length_idx = cell(N.vel,1);
int_times = cell(N.vel,1);
int_stats = cell(N.vel,1);
for jj = 1:N.vel
    for kk = 1:n_main_fields
        time_field = IntAll.(field_time)(jj).time;
        allnan_idx = isnan(time_field(1,:)); % use time of first field to find nan's
    
        field_names = string(fieldnames(IntSort.(main_fields(kk))));
        for f = 1:length(field_names)
            IntSort.(main_fields(kk))(jj).(field_names(f)) = ...
                IntSort.(main_fields(kk))(jj).(field_names(f))(:,~allnan_idx);
        end

        % Sort intervals
        time_field = IntSort.(main_fields(kk))(jj).time;
        valid_time{jj} = ~isnan(time_field); % valid times in intervals
        length_idx{jj} = sum(valid_time{jj},1); % how long intervals are in samples
        [length_idx{jj},order] = sort(length_idx{jj},'ascend'); % sort by lengths by length
        max_length = max(length_idx{jj});
        valid_time{jj} = valid_time{jj}(1:max_length,order); % valid times by length

        % Sort intervals by length & get rid of trailing nan's
        for f = 1:length(field_names)
            IntSort.(main_fields(kk))(jj).(field_names(f)) = ...
                IntSort.(main_fields(kk))(jj).(field_names(f))(1:max_length,order);
        end

        int_times{jj} = length_idx{jj} * Ts;
        int_stats{jj} = basic_stats(int_times{jj},2);
    end
end

% Filtering
Fs = 1/Ts;
Fc = 2;
n = 2;
[b,a] = butter(n, Fc/(Fs/2), 'low');
for jj = 1:N.vel
    for kk = 1:n_main_fields
        field_names = string(fieldnames(IntSort.(main_fields(kk))));
        field_names = field_names(2:3); % don't filter "time"
        for f = 1:length(field_names)
        	all = IntSort.(main_fields(kk))(jj).(field_names(f));
            for ii = 1:size(all,2)
              	vals = all(:,ii);
                filtidx = ~isnan(vals);
                vals(filtidx) = filtfilt(b,a,vals(filtidx));
                IntSort.(main_fields(kk))(jj).(field_names(f))(filtidx,ii) = vals(filtidx);
            end
        end
    end
end

% Cut off based on times
head_amp = 20;
head_gain = 0.3;
% head_gain = repmat( 20 ./ [30 60 60 60 60]' , 2, 1);
cut_time = head_amp ./ (head_gain.*abs(Vel));

IntBelow = IntSort;
IntAbove = IntSort;
int_time_below = int_times;
int_time_above = int_times;
below_valid_time = valid_time;
for jj = 1:N.vel
    for kk = 1:n_main_fields
    	time_field = IntSort.(field_time)(jj).time;
        cut_idx = find(time_field(:,end) <= cut_time(jj), 1, 'last');
        rmv_idx = sum(valid_time{jj},1) > cut_idx;

        % Seperate intervals shorter & longer than cut-time
        field_names = string(fieldnames(IntSort.(main_fields(kk))));
        for f = 1:length(field_names)
            IntBelow.(main_fields(kk))(jj).(field_names(f)) = ...
                IntSort.(main_fields(kk))(jj).(field_names(f))(:,~rmv_idx);
            IntAbove.(main_fields(kk))(jj).(field_names(f)) = ...
                IntAbove.(main_fields(kk))(jj).(field_names(f))(:,rmv_idx);
        end
        int_time_below{jj} = Ts * sum(~isnan(IntBelow.norm_interval(jj).time),1);
        int_time_above{jj} = Ts * sum(~isnan(IntAbove.norm_interval(jj).time),1);

        below_valid_time{jj} = ~isnan(IntBelow.(main_fields(kk))(jj).time);
        max_length = max(sum(below_valid_time{jj},1));

        % Get rid of trailing nan's in shorter intervals
        for f = 1:length(field_names)
            IntBelow.(main_fields(kk))(jj).(field_names(f)) = ...
                IntBelow.(main_fields(kk))(jj).(field_names(f))(1:max_length,:);
        end
    end
end

%% Split intervals
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = 1.0*[2 2 clms*(4/3) 3*(3/2)];
FIG.Name = 'Interval Time Visualization';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(3*(N.vel/2),1);
clear h

for jj = 1:N.vel/2
    ax(jj) = subplot(3,clms,jj) ; hold on
    title([num2str(Vel(jj)) ' (°/s)']);
        time = IntSort.normstart_interval(jj).time;
        pos = IntSort.normstart_interval(jj).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);

        time = IntSort.normstart_interval(jj + N.vel/2).time;
        pos = IntSort.normstart_interval(jj + N.vel/2).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);

        plot([cut_time(jj) cut_time(jj)],40*[-1 1],'-','Color',[0.5 0.5 0.5],'LineWidth', 1);
    
    ax(jj + clms) = subplot(3,clms,jj + clms) ; hold on
        time = IntBelow.normstart_interval(jj).time;
        pos = IntBelow.normstart_interval(jj).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);

        time = IntBelow.normstart_interval(jj + N.vel/2).time;
        pos = IntBelow.normstart_interval(jj + N.vel/2).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
        ax(jj + clms).XLim = [0 cut_time(jj)];
        
    ax(jj + 2*clms) = subplot(3,clms,jj + 2*clms) ; hold on
        time = IntAbove.normstart_interval(jj).time;
        pos = IntAbove.normstart_interval(jj).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);

        time = IntAbove.normstart_interval(jj + N.vel/2).time;
        pos = IntAbove.normstart_interval(jj + N.vel/2).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
end
linkaxes(ax,'y')
set(ax,'LineWidth',1.5,'YLim',30*[-1 1])
% set(ax(1:clms),'XLim',[0 1.7])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax([1,clms+1,2*clms+1]), 'YLabel');
set([YLabelHC{:}], 'String', 'Head Position (°)')

%% Visualize interval position/time
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 clms*(4/3) 3];
FIG.Name = 'Interval Time Visualization: Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel,1);
clear h

for jj = 1:N.vel
    ax(jj) = subplot(ceil(N.vel/clms),clms,jj);
        img = IntSort.normstart_interval(jj).position;
        if U.vel{1}(jj) < 0
            img = -img;
        end
        blank = isnan(img);
        h = imagesc(img,'AlphaData',~blank,[0 15]); hold on

        title([num2str(Vel(jj)) ' (°/s)'])
        
        n_trial = size(valid_time{jj},2);
        
     	% Find the transition cutoff
        cut_trial = find(int_times{jj} >= cut_time(jj), 1, 'first');
        if isempty(cut_trial)
            cut_trial = size(int_times{jj},2);
        end

        cut_line1 = repmat(cut_time(jj) / Ts, 1, n_trial);
        cut_line2 = repmat(cut_trial, 2, 1);
        
        %plot(1:n_trial, cut_line1, 'g', 'LineWidth', 1)
        %plot(cut_line2, [1 , max(length_idx{jj})], 'g', 'LineWidth', 1)
        
        ax(jj).XTick = [1 n_trial];
        ax(jj).YTick = 1:200:2000;
        ax(jj).YTickLabels = string((ax(jj).YTick -1 )*Ts); % convert lables to time
end
linkaxes(ax,'y')
YLabelHC = get(ax([1,clms+1]), 'YLabel');
set([YLabelHC{:}], 'String', 'Time (s)')
XLabelHC = get(ax(clms+1:2*clms), 'XLabel');
set([XLabelHC{:}], 'String', 'Interval #')

set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','k',...
    'YColor','k','XColor','k')

colormap(jet)

hp = get(subplot(2,5,10),'Position');
cbar = colorbar('Position', [hp(1)+hp(3)+0.02  0.05+hp(2)  0.01  hp(2)+hp(3)*1],'Ticks',[-15 15]);
cbar.Label.String = 'Head Position (°)';

%% Interval Position %%
FIG = figure (3) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 4 3];
FIG.Name = 'Interval Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N{1,3},1);
clear ax h
ax(1) = subplot(1,1,1) ; hold on

pp = 1;
pLim = 0.9;
% test = [1:3,6:8];
for jj = 1:N.vel   
    time = IntBelow.normstart_interval(jj).time;
    pos = IntBelow.normstart_interval(jj).position;
    med_time = nanmedian(time,2);
    med_pos = nanmedian(pos,2);
    std_pos = nanstd(pos,[],2);
    
   	tLim = sum(isnan(time),2)./(size(time,2)-1);
    span = 1:length(tLim(tLim<pLim));
   	
    h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.4], 'LineWidth', 0.6);
	                              
    [h.patch(pp),h.med(pp)] = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
        CC(jj,:), 0.5*CC(jj,:), 0.2, 1.5);
    % h.patch.EdgeColor = CC(jj,:);
    
    pp = pp + 1;
end
xlabel('Time (s)')
ylabel('Head Position (°)')
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',[0 2])
set(ax,'YLim',30*[-1 1])
uistack(h.med','top')

%% Interval Position Error %%
FIG = figure (4) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 4 3];
FIG.Name = 'Interval Position Error';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N{1,3},1);
clear ax h
ax(1) = subplot(1,1,1) ; hold on

pp = 1;
pLim = 0.9;
% test = [1:3,6:8];
for jj = 1:N.vel   
    time = IntBelow.error(jj).time;
    pos = IntBelow.error(jj).position;
    med_time = nanmedian(time,2);
    med_pos = nanmedian(pos,2);
    std_pos = nanstd(pos,[],2);
    
   	tLim = sum(isnan(time),2)./(size(time,2)-1);
    span = 1:length(tLim(tLim<pLim));
   	
    h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.4], 'LineWidth', 0.6);
	                              
%     [h.patch(pp),h.med(pp)] = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
%         CC(jj,:), 0.5*CC(jj,:), 0.2, 1.5);
    % h.patch.EdgeColor = CC(jj,:);
    
    pp = pp + 1;
end
xlabel('Time (s)')
ylabel('Head Position Error (°)')
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',[0 2])
set(ax,'YLim',50*[-1 1])
% uistack(h.med','top')

%% Interval Velocity %%
FIG = figure (5) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 1.5*clms*(4/3) 3];
FIG.Name = 'Interval Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N{1,3},1);
clear ax h
ax(1) = subplot(1,1,1) ; hold on

pp = 1;
pLim = 0.9;
for jj = 1:N.vel
    ax(pp) = subplot(2,clms,pp) ; hold on
    title([num2str(Vel(jj)) ' (°/s)']);
    
    time = IntBelow.normstart_interval(jj).time;
    vel = IntBelow.normstart_interval(jj).velocity;
    med_time = nanmedian(time,2);
    med_vel = nanmedian(vel,2);
    std_vel = nanstd(vel,[],2);
    
   	tLim = sum(isnan(time),2)./(size(time,2)-1);
    span = 1:length(tLim(tLim<pLim));
   	
    h.trial = plot(time, vel, 'Color', [0.7*CC(jj,:) , 0.4], 'LineWidth', 0.6);
	                              
    [h.patch(1),h.med(1)] = PlotPatch(med_vel(span), std_vel(span), med_time(span), 1, 1, ...
        CC(jj,:), 0.5*CC(jj,:), 0.2, 1.5);
    h.patch(1).EdgeColor = CC(jj,:);
    
%     time = IntBelow.normstart_interval(jj + clms).time;
%     pos = IntBelow.normstart_interval(jj + clms).velocity;
%     med_time = nanmedian(time,2);
%     med_pos = nanmedian(pos,2);
%     std_pos = nanstd(pos,[],2);
%     
%    	tLim = sum(isnan(time),2)./(size(time,2)-1);
%     span = 1:length(tLim(tLim<pLim));
%    	
%     h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.4], 'LineWidth', 0.6);
% 	                              
%     [h.patch(2),h.med(2)] = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
%         CC(jj,:), [0.5 0.5 0.5], 0.2, 1.5);
%     h.patch(2).EdgeColor = CC(jj,:);
    
    uistack(h.patch,'top')
    uistack(h.med,'top')
    
    pp = pp + 1;
end
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (s)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Head Velocity (°/s)')

set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',[0 2])
set(ax,'YLim',100*[-1 1])
linkaxes(ax,'xy')

%% Interval Times %%
int_time_below_all = cat(2,int_time_below{:})';
% groups = (1:N.vel)';
groups = repmat((1:N.vel/2)',2,1);
G = cellfun(@(x,y) y*ones(size(x)), int_time_below, num2cell(groups), 'UniformOutput', false);
G = cat(2,G{:})';

FIG = figure (5) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 2 2];
FIG.Name = 'Interval Times';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel,1);
clear ax h
ax(1) = subplot(1,1,1) ; hold on
    bx = boxplot(int_time_below_all,G,'Labels',{Speed},'Width',0.5,'Symbol','.','Whisker',2);
    box off
    xlabel('Stimulus Speed (°/s)')
    ylabel('Interval Times (s)')
    
    h = get(bx(5,:),{'XData','YData'});
    for jj = 1:size(h,1)
       patch(h{jj,1},h{jj,2},CC(jj,:));
    end

    set(findobj(ax,'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax,'tag','Box'), 'Color', 'none');
    set(findobj(ax,'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax,'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax.Children = ax.Children([end 1:end-1]);
    
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8)

%% Interval Gains %%
FIG = figure (100) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 5 5];
FIG.Name = 'Velocity Response';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = subplot(1,1,1) ; hold on
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8)

PeakGain = cell(N.vel,1);
PeakTime = cell(N.vel,1);
Gain = cell(N.vel,1);
for jj = 1:N.vel
    for ii = 1:size(IntBelow.normstart_interval(jj).time,2)
        time = IntBelow.normstart_interval(jj).time(:,ii);
        vel = IntBelow.normstart_interval(jj).velocity(:,ii);
        Vel_temp = U.vel{1}(jj);
        if Vel_temp > 0
            [max_vel,max_idx] = max(vel);
        else
            [max_vel,max_idx] = min(vel);
        end
        
        peak_time = time(max_idx);
        
        PeakGain{jj}(ii,1) = max_vel./ Vel_temp;
        PeakTime{jj}(ii,1) = peak_time;     
        Gain{jj}(ii,1) = nanmedian(vel) / Vel_temp;
        
%         cla
%         plot(time,vel,'Color', CC(jj,:),'LineWidth', 1)
%         plot(peak_time, max_vel,'.k','MarkerSize',15)  
%         pause
    end
end
close(FIG)

PeakGain_all = cat(1,PeakGain{:});
PeakTime_all = cat(1,PeakTime{:});
Gain_all = cat(1,Gain{:});

% groups = (1:N.vel)';
groups = repmat((1:N.vel/2)',2,1);
G = cellfun(@(x,y) y*ones(size(x)), PeakGain, num2cell(groups), 'UniformOutput', false);
G = cat(1,G{:});

FIG = figure (6) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 2 2];
FIG.Name = 'Interval Times';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel,1);
clear ax h
ax(1) = subplot(1,1,1) ; hold on
    bx = boxplot(Gain_all,G,'Labels',{Speed},'Width',0.5,'Symbol','.','Whisker',2);
    box off
    xlabel('Stimulus Speed (°/s)')
    ylabel('Interval Peak Gain (°/s / °/s)')
    
    h = get(bx(5,:),{'XData','YData'});
    for jj = 1:size(h,1)
       patch(h{jj,1},h{jj,2},CC(jj,:));
    end

    set(findobj(ax,'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax,'tag','Box'), 'Color', 'none');
    set(findobj(ax,'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax,'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax.Children = ax.Children([end 1:end-1]);
    
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8)

end