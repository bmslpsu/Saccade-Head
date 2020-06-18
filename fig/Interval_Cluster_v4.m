function [] = Interval_Cluster_v4()
%% Interval_Cluster:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'GRAND','U','N','SACCADE')

clms = N.vel/2;
CC = repmat(hsv(clms),2,1);

Vel = U.vel{1};
Speed = Vel(1:N.vel/2);

%% Sort Saccades
clearvars -except clms CC Speed Vel COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

Ts = SACCADE.head_saccade{1}.Ts; % sampling time [s]
Fs = 1/Ts; % sampling frequency [Hz]

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

% Filter
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

% Combine directions
IntSort_Speed = [];
for jj = 1:N.vel/2
    for kk = 1:n_main_fields
        field_names = string(fieldnames(IntSort.(main_fields(kk))));
        for f = 1:length(field_names)
            cw  = IntSort.(main_fields(kk))(jj).(field_names(f));
            ccw = IntSort.(main_fields(kk))(jj + N.vel/2).(field_names(f));
            
            if strcmp(field_names(f),'time')
                flip = 1;
            else
                flip = -1;
            end
            
            dim = [size(cw,1) ; size(ccw,1)];
            sz_diff = diff(dim);
            if sz_diff == 0
                comb = [cw, flip*ccw];
            elseif sz_diff > 0
                comb = [padmatrix(cw,sz_diff,nan,2), flip*ccw];
            elseif sz_diff < 0
                comb = [cw, flip*padmatrix(ccw,-sz_diff,nan,2)];
            end
            
            IntSort_Speed.(main_fields(kk))(jj).(field_names(f)) = comb;
        end
        %disp(sz_diff)
        IntSort_Speed.(main_fields(kk))(jj) = struct_center(IntSort_Speed.(main_fields(kk))(jj), 0, false, 1);
    end
end

% Stats
IntSort_Speed_Stats = [];
int_times_speed = [];
for jj = 1:N.vel/2
    for kk = 1:n_main_fields
        field_names = string(fieldnames(IntSort.(main_fields(kk))));
        for f = 1:length(field_names)
            IntSort_Speed_Stats.(main_fields(kk))(jj).(field_names(f)) = ...
                basic_stats(IntSort_Speed.(main_fields(kk))(jj).(field_names(f)), 2);
        end
    end
    int_times_speed(jj).all = [int_times(jj).all , int_times(jj + N.vel/2).all];
    int_times_speed_stats(jj).all = basic_stats(int_times_speed(jj).all,2);
end

%% Find saturation and peak times
figure (100); clf
ax = subplot(1,1,1) ; hold on
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8)

move_win = round(0.1*Fs);
IntProp = [];
for jj = 1:N.vel/2
	time_all = IntSort_Speed.norm_interval(jj).time;
    vel_all = IntSort_Speed.norm_interval(jj).velocity;
  	move_vel_all = movmean(vel_all, move_win, 'omitnan');
    for ii = 1:size(time_all,2)
        time = time_all(:,ii);
        vel = vel_all(:,ii);
        last_val = find(isnan(vel), 1, 'first') - 1;
        
        [~,locs] = findpeaks(vel,'MinPeakHeight', 0.5, 'MinPeakProminence', 1);
        
        if sum(~isnan(vel))*Ts > 0.5
            move_vel = move_vel_all(:,ii);
        else
            move_vel = vel;
        end
        
        if ~isempty(locs)
            move_vel(1:locs(1)) = nan;
        end
        
        sat_point = find(move_vel < 3, 1, 'first');
        
        if isempty(sat_point)
            if ~isempty(locs)
                temp_vel = vel;
                temp_vel(1:locs(1)) = nan;
                [~,sat_point] = min(temp_vel);
            else
                [~,sat_point] = min(vel);
            end
            
        end
        
        if sat_point > last_val
            sat_point = last_val;
        end
        
        time_cut = time(sat_point);
        
        cut_vel = vel;
        cut_vel(sat_point+1:end) = nan;
        cut_time = time;
        cut_time(sat_point+1:end) = nan;
        
        [max_vel,max_idx] = max(cut_vel);
        peak_time = time(max_idx);
    
        IntProp(jj).Cut_Vel(:,ii)   = cut_vel;
        IntProp(jj).Cut_Time(:,ii)	= cut_time;
        IntProp(jj).Time_Cut(ii)    = time_cut;
        IntProp(jj).Peak_Vel(ii)    = max_vel;
        IntProp(jj).Peak_Time(ii)   = peak_time;
        IntProp(jj).Mean_Vel(ii)    = nanmean(cut_vel);
      	IntProp(jj).Mean_Gain(ii) 	= nanmean(cut_vel) / Speed(jj);
        IntProp(jj).Peak_Gain(ii) 	= max_vel / Speed(jj);
        
        if time_cut < 0.05
%             cla
%             plot(time, vel, 'Color', 'k', 'LineWidth', 1)
%             plot(time, move_vel, 'Color', [0.5 0.5 0.5], 'LineWidth', 1)
%             plot(cut_time, cut_vel, 'Color', CC(jj,:), 'LineWidth', 2)
%             plot(time_cut, vel(sat_point), '.m', 'MarkerSize', 15)
%             plot(peak_time, max_vel, '.c', 'MarkerSize', 15)
% 
%             pause()
        end
    end
end

% Stats
Mean_Gain_Stats = cellfun(@(x) basic_stats(x,2), {IntProp.Mean_Gain});
Peak_Gain_Stats = cellfun(@(x) basic_stats(x,2), {IntProp.Peak_Gain});
Temp_Freq = Speed / U.wave;

%% Save
savedir = 'C:\Users\BC\Box\Research\Manuscripts\Head Saccade\Data';
fname = ['Int_Gain_Stats_wave=' num2str(U.wave) '.mat'];
% save(fullfile(savedir,fname), 'IntProp', 'Mean_Gain_Stats', 'Peak_Gain_Stats', 'Temp_Freq', 'U')

%% Interval Position Norm
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = 1.0*[2 2 clms*(4/3) 1*(3/2)];
FIG.Name = 'Interval Time Visualization';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel/2,1);
clear h
for jj = 1:N.vel/2
    ax(jj) = subplot(1,clms,jj) ; hold on
    title([num2str(Vel(jj)) ' (°/s)']);
        time = IntSort_Speed.normstart_interval(jj).time;
        pos = IntSort_Speed.norm_interval(jj).position;
        med_time = IntSort_Speed_Stats.norm_interval(jj).time.mean;
        med_pos = IntSort_Speed_Stats.norm_interval(jj).position.mean;
        std_pos = IntSort_Speed_Stats.norm_interval(jj).position.std;
        
        h.trial = plot(time, pos, 'Color', [0.5 0.5 0.5 , 0.2], 'LineWidth', 0.5);
        
    	span_std = 1: Fs*(int_times_speed_stats(jj).all.mean + 1*int_times_speed_stats(jj).all.std);
        span_med = 1: Fs*(int_times_speed_stats(jj).all.mean);
        
    	[h.patch(jj,1),h.med(jj,1)] = PlotPatch(med_pos(span_std), std_pos(span_std), med_time(span_std), 1, 1, ...
                                        CC(jj,:), 0.5*CC(jj,:), 0.3, 1);
        h.med(jj,1).Color = 'none';
        
        plot(med_time(span_med), med_pos(span_med), 'Color', CC(jj,:), 'LineWidth', 1)
end
linkaxes(ax,'xy')
set(ax,'LineWidth',1,'YLim',30*[-1 1])
set(ax,'XLim',[0 2])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (s)')
YLabelHC = get(ax(1,:), 'YLabel');
set([YLabelHC], 'String', 'Head Position (°)')

%% Interval Velocity Norm
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = 1.0*[2 2 clms*(4/3) 1*(3/2)];
FIG.Name = 'Interval Time Visualization';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel/2,1);
clear h
for jj = 1:N.vel/2
    ax(jj) = subplot(1,clms,jj) ; hold on
    title([num2str(Vel(jj)) ' (°/s)']);
        time = IntSort_Speed.normstart_interval(jj).time;
        vel = IntSort_Speed.norm_interval(jj).velocity;
        med_time = IntSort_Speed_Stats.norm_interval(jj).time.mean;
        med_vel = IntSort_Speed_Stats.norm_interval(jj).velocity.mean;
        std_vel = IntSort_Speed_Stats.norm_interval(jj).velocity.std;
        
        h.trial = plot(time, vel, 'Color', [0.5 0.5 0.5 , 0.2], 'LineWidth', 0.5);
        
    	span_std = 1: Fs*(int_times_speed_stats(jj).all.mean + 1*int_times_speed_stats(jj).all.std);
        span_med = 1: Fs*(int_times_speed_stats(jj).all.mean);
        
    	[h.patch(jj,1),h.med(jj,1)] = PlotPatch(med_vel(span_std), std_vel(span_std), med_time(span_std), 1, 1, ...
                                        CC(jj,:), 0.5*CC(jj,:), 0.3, 1);
        h.med(jj,1).Color = 'none';
        
        plot(med_time(span_med), med_vel(span_med), 'Color', CC(jj,:), 'LineWidth', 1)
        
        plot([0 2], [0 0], '-k', 'LineWidth', 0.5)
end
linkaxes(ax,'xy')
set(ax,'LineWidth',1,'FontSize',8,'XLim',[0 2])
set(ax, 'YLim', [-50 110])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (s)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Head Velocity (°/s)')

%% Interval Times %%
fig = figure (3);
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 9 2])
ax = gobjects(N.vel/2,1);
edges = 0:0.1:7;
h = gobjects(N.vel/2,1);
for jj = 1:N.vel/2
    ax(jj) = subplot(1,N.vel/2,jj); hold on
        h(jj) = histogram(int_times_speed(jj).all, edges, 'FaceColor', CC(jj,:));
end
set(h, 'Normalization', 'Probability', 'EdgeColor', 'none', 'FaceAlpha', 0.9)
set(ax,'LineWidth', 1, 'FontSize', 8, 'YLim', [-0.01 0.4])
linkaxes(ax, 'xy')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Probability')
XLabelHC = get(ax(1), 'XLabel');
set([XLabelHC], 'String', 'Interval Time (s)')

%% Interval Stats
groups = (1:N.vel/2)';
Labels = Speed;
G = cellfun(@(x,y) y*ones(size(x)), {IntProp(:).Time_Cut}', num2cell(groups), ...
    'UniformOutput', false);
G = cat(2,G{:})';

box_data = ["Time_Cut","Peak_Time","Peak_Vel","Mean_Vel"];
n_plot = length(box_data);
ylim_list = [2 1 150 70];

FIG = figure (6) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 2*n_plot 1.5];
FIG.Name = 'Interval Times';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(n_plot,1);
clear h
for kk = 1:n_plot
    ax(kk) = subplot(1,n_plot,kk) ; hold on
        data = cat(2,IntProp(:).(box_data(kk)))';
        bx = boxplot(data,G,'Labels',{Labels},'Width',0.5,'Symbol','','Whisker',2);
        box off
        xlabel('Stimulus Speed (°/s)')
        ylabel(box_data(kk), 'interpreter', 'none')

        h = get(bx(5,:),{'XData','YData'});
        for jj = 1:size(h,1)
           patch(h{jj,1},h{jj,2},CC(jj,:));
        end

        set(findobj(ax(kk),'tag','Median'), 'Color', 'w','LineWidth',1.5);
        set(findobj(ax(kk),'tag','Box'), 'Color', 'none');
        set(findobj(ax(kk),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
        set(findobj(ax(kk),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
        ax(kk).Children = ax(kk).Children([end 1:end-1]);
        
        ax(kk).YLim = [-0.1*ylim_list(kk) ylim_list(kk)];
end
set(ax, 'LineWidth', 1, 'FontSize', 8)

%% Interval Integrated Velocity Error
FIG = figure (7) ; clf
FIG.Units = 'inches';
FIG.Position = 1.0*[2 2 clms*(4/3) 1*(3/2)];
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel/2,1);
clear h
for jj = 1:N.vel/2
    ax(jj) = subplot(1,clms,jj) ; hold on
    title([num2str(Vel(jj)) ' (°/s)']);
        time = IntSort_Speed.normstart_interval(jj).time;
        pos = IntSort_Speed.error(jj).velocity;
        med_time = IntSort_Speed_Stats.norm_interval(jj).time.mean;
        med_pos = IntSort_Speed_Stats.error(jj).velocity.mean;  
        std_pos = IntSort_Speed_Stats.error(jj).velocity.std;
        
        h.trial = plot(time, pos, 'Color', [0.5 0.5 0.5 , 0.2], 'LineWidth', 0.5);
        
    	span_std = 1: Fs*(int_times_speed_stats(jj).all.mean + 1*int_times_speed_stats(jj).all.std);
        span_med = 1: Fs*(int_times_speed_stats(jj).all.mean);
        
    	[h.patch(jj,1),h.med(jj,1)] = PlotPatch(med_pos(span_std), std_pos(span_std), med_time(span_std), 1, 1, ...
                                        CC(jj,:), 0.5*CC(jj,:), 0.3, 1);
        h.med(jj,1).Color = 'none';
        
        plot(med_time(span_med), med_pos(span_med), 'Color', CC(jj,:), 'LineWidth', 1)
end
linkaxes(ax,'xy')
set(ax,'LineWidth',1,'YLim',100*[-1 1])
set(ax,'XLim',[0 2])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (s)')
YLabelHC = get(ax(1,:), 'YLabel');
set([YLabelHC], 'String', 'Head Position (°)')

end