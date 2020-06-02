function [] = Interval_Cluster_v2()
%% Interval_Cluster:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'GRAND','U','N','SACCADE')

clms = N.vel/2;
CC = repmat(hsv(clms),2,1);

Vel = U.vel{1};
Speed = Vel(1:N.vel/2);

%% Cluster saccades
clearvars -except clms CC Speed Vel COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

Ts = SACCADE.head_saccade{1}.Ts; % sampling time [s]

IntAll = GRAND;
IntSort = IntAll; % grouped saccade & intervals
main_fields = string(fieldnames(IntSort)); % all fields
n_main_fields = length(main_fields); % # of fields
field_time = main_fields(2);

% Get rid of intervals before 1st saccade
length_idx = struct('all',[],'below',[],'above',[],'speed',[]);
valid_time = struct('all',[],'below',[],'above',[],'speed',[]);
int_times  = struct('all',[],'below',[],'above',[],'speed',[]);
int_stats  = struct('all',[],'below',[],'above',[],'speed',[]);
for jj = 1:N.vel
    time_field = IntAll.(field_time)(jj).time;
	allnan_idx = isnan(time_field(1,:)); % use time of first field to find nan's
    for kk = 1:n_main_fields
        field_names = string(fieldnames(IntSort.(main_fields(kk))));
        for f = 1:length(field_names)
            IntSort.(main_fields(kk))(jj).(field_names(f)) = ...
                IntSort.(main_fields(kk))(jj).(field_names(f))(:,~allnan_idx);
        end

        % Sort intervals
        time_field = IntSort.(main_fields(kk))(jj).time;
        valid_time(jj).all = ~isnan(time_field); % valid times in intervals
        temp = sum(valid_time(jj).all, 1); % how long intervals are in samples
        [temp,order] = sort(temp,'ascend'); % sort by lengths by length
        length_idx(jj).all = temp;
        max_length = max(length_idx(jj).all);
        valid_time(jj).all =valid_time(jj).all(1:max_length,order); % valid times by length

        % Sort intervals by length & get rid of trailing nan's
        for f = 1:length(field_names)
            IntSort.(main_fields(kk))(jj).(field_names(f)) = ...
                IntSort.(main_fields(kk))(jj).(field_names(f))(1:max_length,order);
        end

        int_times(jj).all = length_idx(jj).all * Ts;
        int_stats(jj).all = basic_stats(int_times(jj).all,2);
    end
end

% Filtering
Fs = 1/Ts;
Fc = 2;
n = 2;
[b,a] = butter(n, Fc/(Fs/2), 'low');
for jj = 1:N.vel
    for kk = 1:n_main_fields
        if kk ~= 1 && kk~=4
            field_names = string(fieldnames(IntSort.(main_fields(kk))));
            field_names = field_names(2); % don't filter time or velocity (only position)
            for f = 1:length(field_names)
                all = IntSort.(main_fields(kk))(jj).(field_names(f));
                for ii = 1:size(all,2)
                    vals = all(:,ii);
                    filtidx = ~isnan(vals);
                    vals(filtidx) = filtfilt(b,a,vals(filtidx));
                    IntSort.(main_fields(kk))(jj).(field_names(f))(filtidx,ii) = vals(filtidx);
                    
                    dvals = diff(vals(filtidx),1)*Fs;
                    IntSort.(main_fields(kk))(jj).('velocity')(filtidx,ii) = [dvals(1) ; dvals];
                end
            end
        end
    end
end

% Cut off based on times
head_amp = 20;
% head_gain = [0.5 0.5 0.5 0.5 0.5]';
head_gain = [0.5 0.45 0.4 0.27 0.17]'; % from previous calculations
cut_time = head_amp ./ (head_gain.*Speed);
cut_time = [cut_time ; cut_time];
% cut_amp = 10;

IntBelow = IntSort;
IntAbove = IntSort;
Track = struct('tracking',[],'landing',[],'ratio',[]);
for jj = 1:N.vel
	time_field = IntSort.(field_time)(jj).time;
    % amp_field = IntSort.norm_interval(jj).position;
    % amp_end = indexbyarray(amp_field, length_idx(jj).all);
    cut_time_idx = find(time_field(:,end) <= cut_time(jj), 1, 'last');
    rmv_idx = sum(valid_time(jj).all ,1) > cut_time_idx;
    Track(jj).tracking = sum(~rmv_idx);
    Track(jj).landing = sum(rmv_idx);
    Track(jj).ratio = sum(~rmv_idx)/length(rmv_idx);
    for kk = 1:n_main_fields
        % Seperate intervals shorter & longer than cut-time
        field_names = string(fieldnames(IntSort.(main_fields(kk))));
        for f = 1:length(field_names)
            IntBelow.(main_fields(kk))(jj).(field_names(f)) = ...
                IntSort.(main_fields(kk))(jj).(field_names(f))(:,~rmv_idx);
            IntAbove.(main_fields(kk))(jj).(field_names(f)) = ...
                IntAbove.(main_fields(kk))(jj).(field_names(f))(:,rmv_idx);
        end
        int_times(jj).below = Ts * sum(~isnan(IntBelow.norm_interval(jj).time),1);
        int_times(jj).above = Ts * sum(~isnan(IntAbove.norm_interval(jj).time),1);
        int_stats(jj).below = basic_stats(int_times(jj).below,2);
        int_stats(jj).above = basic_stats(int_times(jj).above,2);
        
        valid_time(jj).below = ~isnan(IntBelow.(main_fields(kk))(jj).time);
        valid_time(jj).above = ~isnan(IntAbove.(main_fields(kk))(jj).time);
        max_length = max(sum(valid_time(jj).below,1));

        % Get rid of trailing nan's in shorter intervals
        for f = 1:length(field_names)
            IntBelow.(main_fields(kk))(jj).(field_names(f)) = ...
                IntBelow.(main_fields(kk))(jj).(field_names(f))(1:max_length,:);
        end
    end
    length_idx(jj).below = sum(~isnan(IntBelow.(field_time)(jj).time));
    length_idx(jj).above = sum(~isnan(IntAbove.(field_time)(jj).time));
end

% Group & sort by speed & for all
IntSpeed = [];
for jj = 1:N.vel/2
    for kk = 1:n_main_fields
        % Concatenate intervals with same speed
        field_names = string(fieldnames(IntSort.(main_fields(kk))));
        for f = 1:length(field_names)
            IntSpeed.(main_fields(kk))(jj).(field_names(f)) = nancat(2, ...
                IntSort.(main_fields(kk))(jj).(field_names(f)), ...
                 -IntSort.(main_fields(kk))(jj + clms).(field_names(f)));
        end
        
        % Sort intervals
        time_field = IntSpeed.(main_fields(kk))(jj).time;
        valid_time(jj).speed = ~isnan(time_field); % valid times in intervals
        length_idx(jj).speed = sum(valid_time(jj).speed,1); % how long intervals are in samples
        [length_idx(jj).speed,order] = sort(length_idx(jj).speed,'ascend'); % sort by lengths by length
        max_length = max(length_idx(jj).speed);
        valid_time(jj).speed = valid_time(jj).speed(1:max_length,order); % valid times by length

        % Sort intervals by length & get rid of trailing nan's
        for f = 1:length(field_names)
            IntSpeed.(main_fields(kk))(jj).(field_names(f)) = ...
                IntSpeed.(main_fields(kk))(jj).(field_names(f))(1:max_length,order);
        end

        int_times(jj).speed = length_idx(jj).speed * Ts;
        int_stats(jj).speed = basic_stats(int_times(jj).speed,2);
    end
end
int_times(1).all_speed = cat(2,int_times(:).speed);
int_times(1).all_speed = sort(int_times(1).all_speed);
int_stats(1).all_speed = basic_stats(int_times(jj).speed,2);

%% Visualize interval lengths
FIG = figure (30) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 2*(4/3) 4];
FIG.Name = 'Interval Time Visualization: Speed';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear h
h = gobjects(N.vel/2,1);
ax = subplot(1,1,1); hold on
for jj = 1:N.vel/2
    times = int_times(jj).speed;
    span = linspace(0,1,length(times));
    h(jj) = plot(span, times, '.', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', CC(jj,:));
end
span = linspace(0,1,length(int_times(1).all_speed)); 
h(end+1) = plot(span, int_times(1).all_speed, '.', 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'k');
leg = legend(h, [ string(Vel(1:clms)) ; "ALL"], 'Box', 'off');
leg.Title.String = 'Speed (°/s)';
xlabel('Normalized Trial #')
ylabel('Interval Time')
set(ax,'LineWidth',1.5)

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
        
    ax(jj + 1*clms) = subplot(3,clms,jj + 1*clms) ; hold on
        time = IntAbove.normstart_interval(jj).time;
        pos = IntAbove.normstart_interval(jj).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);

        time = IntAbove.normstart_interval(jj + N.vel/2).time;
        pos = IntAbove.normstart_interval(jj + N.vel/2).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
                
    ax(jj + 2*clms) = subplot(3,clms,jj + 2*clms) ; hold on
        time = IntBelow.normstart_interval(jj).time;
        pos = IntBelow.normstart_interval(jj).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);

        time = IntBelow.normstart_interval(jj + N.vel/2).time;
        pos = IntBelow.normstart_interval(jj + N.vel/2).position;
        h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
        ax(jj + 2*clms).XLim = [0 cut_time(jj)];
        ax(jj + 2*clms).XTick = [0 cut_time(jj)];
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
        
        n_trial = size(valid_time(jj).all,2);
        
     	% Find the transition cutoff
        cut_trial = find(int_times(jj).all >= cut_time(jj), 1, 'first');
        if isempty(cut_trial)
            cut_trial = size(int_times{jj},2);
        end

        cut_line1 = repmat(cut_time(jj) / Ts, 1, n_trial);
        cut_line2 = repmat(cut_trial, 2, 1);
        
        plot(1:n_trial, cut_line1, 'g', 'LineWidth', 1)
        plot(cut_line2, [1 , max(length_idx(jj).all)], 'g', 'LineWidth', 1)
        
        ax(jj).XTick = [1 n_trial];
        ax(jj).YTick = 1:200:2000;
        ax(jj).YTickLabels = string((ax(jj).YTick -1 )*Ts); % convert lables to time
end
set(ax,'YLim',[0 4*Fs])
linkaxes(ax,'y')
YLabelHC = get(ax([1,clms+1]), 'YLabel');
set([YLabelHC{:}], 'String', 'Time (s)')
XLabelHC = get(ax(clms+1:2*clms), 'XLabel');
set([XLabelHC{:}], 'String', 'Interval #')

set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','k',...
    'YColor','k','XColor','k')

colormap(jet)

hp = get(subplot(2,5,10),'Position');
cbar = colorbar('Position', [hp(1)+hp(3)+0.02  0.05+hp(2)  0.01  hp(2)+hp(3)*1],'Ticks',[0 15]);
cbar.Label.String = 'Absolute Value of Head Position (°)';

%% Visualize interval position/time by speed
FIG = figure (3) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 clms*(4/3) 1.5];
FIG.Name = 'Interval Time Visualization: Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel/2,1);
clear h

for jj = 1:N.vel/2
    ax(jj) = subplot(1,clms,jj);
        img = IntSpeed.normstart_interval(jj).position;
        blank = isnan(img);
        h = imagesc(img,'AlphaData',~blank,[0 15]); hold on

        title([num2str(Vel(jj)) ' (°/s)'])
        
        n_trial = size(valid_time(jj).speed,2);
        
     	% Find the transition cutoff
        cut_trial = find(int_times(jj).speed >= cut_time(jj), 1, 'first');
        if isempty(cut_trial)
            cut_trial = size(int_times(jj).speed,2);
        end

        cut_line1 = repmat(cut_time(jj) / Ts, 1, n_trial);
        cut_line2 = repmat(cut_trial, 2, 1);
        
        plot(1:n_trial, cut_line1, 'g', 'LineWidth', 1)
        plot(cut_line2, [1 , max(length_idx(jj).speed)], 'g', 'LineWidth', 1)
        
        ax(jj).XTick = [1 n_trial];
        ax(jj).YTick = 1:200:2000;
        ax(jj).YTickLabels = string((ax(jj).YTick -1 )*Ts); % convert lables to time
end
% set(ax,'YLim',[0 4.01*Fs])
linkaxes(ax,'y')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Time (s)')
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Interval #')
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','k',...
    'YColor','k','XColor','k')

colormap(jet)

hp = get(subplot(1,5,5),'Position');
cbar = colorbar('Position', [hp(1)+hp(3)+0.02  0.05+hp(2)  0.01  hp(2)+hp(3)*1],'Ticks',[0 15]);
cbar.Label.String = 'Absolute Value of Head Position (°)';

%% Interval Position %%
FIG = figure (4) ; clf
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
	                              
%     [h.patch(pp),h.med(pp)] = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
%         CC(jj,:), 0.5*CC(jj,:), 0.2, 1.5);
    % h.patch.EdgeColor = CC(jj,:);
    
    pp = pp + 1;
end
xlabel('Time (s)')
ylabel('Head Position (°)')
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',[0 1.3])
set(ax,'YLim',30*[-1 1])
% uistack(h.med','top')

%% Interval Position Error %%
FIG = figure (5) ; clf
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

%% Interval Integrated Position Error %%
FIG = figure (6) ; clf
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
    time = IntBelow.int_error(jj).time;
    pos = IntBelow.int_error(jj).position;
    med_time = nanmedian(time,2);
    med_pos = nanmedian(pos,2);
    std_pos = nanstd(pos,[],2);
    
   	tLim = sum(isnan(time),2)./(size(time,2)-1);
    span = 1:length(tLim(tLim<pLim));
   	
    h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.4], 'LineWidth', 0.6);
	                              
%     [h.patch(pp),h.med(pp)] = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
%         CC(jj,:), 0.5*CC(jj,:), 0.2, 1.5);
%     h.patch(pp).EdgeColor = CC(jj,:);
    
    pp = pp + 1;
end
xlabel('Time (s)')
ylabel('Head Integrated Position Error (°)')
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',[0 2])
set(ax,'YLim',30*[-1 1])
% uistack(h.med','top')

%% Interval Velocity %%
FIG = figure (5) ; clf
FIG.Units = 'inches';
FIG.Position = 1.0*[2 2 clms*(4/3) 2*(3/2)];
FIG.Name = 'Interval Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel,1);
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
        CC(jj,:), [0.4 0.4 0.6], 0.7, 1.5);
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
    
	plot([0 2],  Vel(jj)*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',1)
 	% plot([0 2], -Vel(jj)*[1 1],'--','Color',[0.5 0.5 0.5],'LineWidth',1)

	plot([0 2],0*[1 1],'-','Color','k','LineWidth',1)

    pp = pp + 1;
end
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (s)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Head Velocity (°/s)')

set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',[0 1])
set(ax(1:5), 'YLim', [-50 155])
set(ax(6:10),'YLim', [-155 50])
linkaxes(ax(1:5),'y')
linkaxes(ax(6:10),'y')
% set(ax,'YTick',sort([0,-150:60:150]))
% linkaxes(ax,'xy')
% set(ax(1:5), 'YLim', 155*[-1 1])

%% Interval Times %%
int_time_below_all = cat(2,int_times(:).below)';
% groups = (1:N.vel)';
% Labels = Vel;
groups = repmat((1:N.vel/2)',2,1);
Labels = Speed;
G = cellfun(@(x,y) y*ones(size(x)), {int_times(:).below}', num2cell(groups), ...
    'UniformOutput', false);
G = cat(2,G{:})';

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
    bx = boxplot(int_time_below_all,G,'Labels',{Labels},'Width',0.5,'Symbol','.','Whisker',2);
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
    
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'YLim',[0 1.4],'YTick',[0:0.2:1.4])

%% Interval Error %%
int_below_error = {IntBelow.error(:).position}';
int_below_error = cellfun(@(x,y) indexbyarray(x,y), int_below_error, {length_idx(:).below}', ...
    'UniformOutput', false);
% groups = (1:N.vel)';
% Labels = Vel;
groups = repmat((1:N.vel/2)',2,1);
Labels = Speed;
G = cellfun(@(x,y) y*ones(size(x)), int_below_error, num2cell(groups), ...
    'UniformOutput', false);
G = cat(2,G{:})';
int_below_error = abs(cat(2,int_below_error{:}));

FIG = figure (7) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 2 2];
FIG.Name = 'Interval Error';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel,1);
clear ax h
ax(1) = subplot(1,1,1) ; hold on
    bx = boxplot(int_below_error,G,'Labels',{Labels},'Width',0.5,'Symbol','.','Whisker',2);
    box off
    xlabel('Stimulus Speed (°/s)')
    ylabel('Interval Error (°)')
    
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

%% Interval Integrated Error %%
int_below_interror = {IntBelow.int_error(:).position}';
int_below_interror = cellfun(@(x,y) indexbyarray(x,y), int_below_interror, {length_idx(:).below}', ...
    'UniformOutput', false);
% groups = (1:N.vel)';
% Labels = Vel;
groups = repmat((1:N.vel/2)',2,1);
Labels = Speed;
G = cellfun(@(x,y) y*ones(size(x)), int_below_interror, num2cell(groups), ...
    'UniformOutput', false);
G = cat(2,G{:})';
int_below_interror = abs(cat(2,int_below_interror{:}));

FIG = figure (7) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 2 2];
FIG.Name = 'Interval Error';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel,1);
clear ax h
ax(1) = subplot(1,1,1) ; hold on
    bx = boxplot(int_below_interror,G,'Labels',{Labels},'Width',0.5,'Symbol','.','Whisker',2);
    box off
    xlabel('Stimulus Speed (°/s)')
    ylabel('Integrated Error (°*s)')
    
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
% set(ax,'YLim',[0 30])

%% Interval Gains %%
FIG = figure (9) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 2 2];
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

Gain_stats = cellfun(@(x) basic_stats(x,1), Gain,'UniformOutput', true);
Gain_med = mean(reshape([Gain_stats(:).median]',N.vel/2,2),2);
Gain_std = mean(reshape([Gain_stats(:).std]',N.vel/2,2),2);

Gain_min = Gain_med - Gain_std;
% Gain_med = [0.5 0.45 0.4 0.27 0.17];

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
ylim([-0.2 3])

end