function [] = Process_Interval_Ramp()
%% Interval_Gain:
root = 'E:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'U','N','SACCADE','Stim')

%% Get filtered intervals
clearvars -except SACCADE U N Stim root

n_speed = N.vel/2;
Vel = U.vel{1};
Speed = Vel(1:N.vel/2);
CC = repmat(hsv(n_speed),2,1);
TempFreq = Speed / U.wave;

Fs = round(SACCADE.head_saccade{1}.Fs);
Fc = 7;
[b,a] = butter(3, Fc/(Fs/2), 'low');

clear Int
Int.pos = cell(N.fly,n_speed);
Int.vel = cell(N.fly,n_speed);
int_win = [];
for n = 1:N.file
    head_saccade = SACCADE.head_saccade{n};
    vel = SACCADE.vel(n);
    fly = SACCADE.fly(n);
    
    pos_filt = filtfilt(b, a, head_saccade.position);
    if vel > n_speed
        vel = vel - n_speed;
        pos_filt = -pos_filt;
    end
    vel_filt = central_diff(pos_filt, head_saccade.Ts); % velocity
    vel_filt = filtfilt(b, a, vel_filt);
    
    [~,int_pos] = getSaccade(head_saccade, pos_filt, [], int_win);
    [~,int_vel] = getSaccade(head_saccade, vel_filt, [], int_win);
    
    Int.pos{fly,vel} = [Int.pos{fly,vel}; int_pos(2:end)];
    Int.vel{fly,vel} = [Int.vel{fly,vel}; int_vel(2:end)];
end
Int.pos = cellfun(@(x) nancat(2, x{:}), Int.pos, 'UniformOutput', false);
Int.vel = cellfun(@(x) nancat(2, x{:}), Int.vel, 'UniformOutput', false);

Int.index_end = cellfun(@(x) sum(~isnan(x),1), Int.pos, 'UniformOutput', false);
Int.time_end = cellfun(@(x) x ./ Fs, Int.index_end, 'UniformOutput', false);

% Remove starts & ends of intervals where saccades have started/ended
cut_win = round(0.035 * Fs); % window for mean velocity at end of interval
start_thresh = -20; % cut below this velocity at end of interval
end_thresh = -20; % cut below this velocity at start of interva1
n_rep = 5; % how many times to repeat the process
for w = 1:n_rep
    for v = 1:n_speed
        for f = 1:N.fly
            %pos_int_all = Int.vel{f,v};
            vel_int_all = Int.vel{f,v};
            endI = sum(~isnan(vel_int_all),1);
            for k = 1:size(endI,2)
                %pos = pos_int_all(:,k);
                vel = vel_int_all(:,k);
                span_start = round(1:cut_win);
                span_end = round(endI(k)-cut_win:endI(k));

                % If we can get the window size
                if all(span_end >=1) && all(span_end <= length(vel)) && (span_start(end) <= length(vel))
                    startV = mean(vel(span_start));
                    endV = mean(vel(span_end));
                    
                    if endV < end_thresh
                        Int.pos{f,v}(span_end,k) = nan;
                        Int.vel{f,v}(span_end,k) = nan;
                    end
                    
                    if startV < start_thresh
                        shift_pos = Int.pos{f,v}(:,k);
                        shift_vel = Int.vel{f,v}(:,k);
                        sI = find(shift_vel > start_thresh, 1 , 'first');
                        eI = find(~isnan(shift_vel), 1 , 'last');
                        int_range = round(sI:eI);
                        new_range = round(1:length(int_range));
                        new_pos = nan(size(shift_pos));
                        new_vel = nan(size(shift_vel));
                        new_pos(new_range) = shift_pos(int_range);
                        new_vel(new_range) = shift_vel(int_range);
                        
                        Int.pos{f,v}(:,k) = new_pos;
                        Int.vel{f,v}(:,k) = new_vel;
                    end
                end
            end
        end
    end
end

%% Compute Error
Int.pos_error = cell(N.fly,n_speed);
Int.vel_error = cell(N.fly,n_speed);
Int.int_pos_error = cell(N.fly,n_speed);
Int.int_vel_error = cell(N.fly,n_speed);
Int.pos_error_end = cell(N.fly,n_speed);
Int.vel_error_end = cell(N.fly,n_speed);
Int.int_pos_error_end = cell(N.fly,n_speed);
Int.int_vel_error_end = cell(N.fly,n_speed);
for v = 1:n_speed
    vel = U.vel{1}(v);
    for f = 1:N.fly
        Int.vel_error{f,v} = vel - Int.vel{f,v};
        Int.time{f,v} = nan(size(Int.pos{f,v}));
        for t = 1:size(Int.pos{f,v},2)
            endI = sum(~isnan(Int.pos{f,v}(:,t)));
            Int.time{f,v}(1:endI,t) = (1/Fs) * (0:endI-1);
            Int.pos_error{f,v}(:,t) = (vel .* Int.time{f,v}(:,t)) - (Int.pos{f,v}(:,t) - Int.pos{f,v}(1,t));
            Int.pos_error_end{f,v}(:,t) = Int.pos_error{f,v}(endI,t);
            Int.vel_error_end{f,v}(:,t) = Int.vel_error{f,v}(endI,t);
            
            Int.int_pos_error{f,v}(:,t) = cumtrapz(1/Fs, Int.pos_error{f,v}(:,t), 1);
            Int.int_vel_error{f,v}(:,t) = cumtrapz(1/Fs, Int.vel_error{f,v}(:,t), 1);
            
            Int.int_pos_error_end{f,v}(:,t) = Int.int_pos_error{f,v}(endI,t);
            Int.int_vel_error_end{f,v}(:,t) = Int.int_vel_error{f,v}(endI,t);
        end
    end
end

%% Compute all stats
clear Int_all
clear Int_stats
Int_stats.fly_all = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    Int, 'UniformOutput', false);

fnames = string(fieldnames(Int_stats.fly_all));
n_field = length(fnames);
for f = 1:n_field
    for v = 1:n_speed
        Int_all.(fnames(f)){v} = nancat(2, Int.(fnames(f)){:,v});
        Int_stats.fly_mean.(fnames(f)){v} = nancat(2, Int_stats.fly_all.(fnames(f))(:,v).mean);
        Int_stats.fly_median.(fnames(f)){v} = nancat(2, Int_stats.fly_all.(fnames(f))(:,v).median);
        Int_stats.fly_std.(fnames(f)){v} = nancat(2, Int_stats.fly_all.(fnames(f))(:,v).std);
    end
end

Int_stats.all = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    Int_all, 'UniformOutput', false);

Int_stats.fly = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    Int_stats.fly_mean, 'UniformOutput', false);

%% Intervals all
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 8 3])
movegui(fig, 'center')
n_std = 1;
ax = gobjects(2,n_speed); clear h
for v = 1:n_speed
    span_std = 1:round(Int_stats.all.index_end(v).median + n_std*Int_stats.all.index_end(v).std);
    span_med = 1:round(Int_stats.all.index_end(v).median + 0*Int_stats.all.index_end(v).std);
    
    ax(1,v) = subplot(2,n_speed,v); cla ; hold on
        plot(Int_all.time{v}, Int_all.pos{v}, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
        %plot(Int_stats.fly_mean.time{v}(span_std,:), Int_stats.fly_mean.pos{v}(span_std,:), ...
                            %'Color', [0.5 0.5 0.5 0.7], 'LineWidth', 0.5)

        [h.pos_patch(v),h.pos_mean(v)] = PlotPatch(Int_stats.all.pos(v).mean(span_std), ...
                                                    Int_stats.all.pos(v).std(span_std), ...
                                                    Int_stats.all.time(v).mean(span_std), 1, 1, ...
                                                    'k', CC(v,:), 0.3, 1);
        plot(Int_stats.all.time(v).mean(span_med), Int_stats.all.pos(v).mean(span_med,:), ...
            'Color', CC(v,:), 'LineWidth', 1)
        
   ax(2,v) = subplot(2,n_speed,v + n_speed); cla ; hold on
        plot(Int_all.time{v}, Int_all.vel{v}, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
        %plot(Int_stats.fly_mean.time{v}(span_std,:), Int_stats.fly_mean.vel{v}(span_std,:), ...
                            %'Color', [0.5 0.5 0.5 0.7], 'LineWidth', 0.5)

        [h.vel_patch(v),h.vel_mean(v)] = PlotPatch(Int_stats.all.vel(v).mean(span_std), ...
                                                    Int_stats.all.vel(v).std(span_std), ...
                                                    Int_stats.all.time(v).mean(span_std), 1, 1, ...
                                                    'k', CC(v,:), 0.3, 1);
        plot(Int_stats.all.time(v).mean(span_med), Int_stats.all.vel(v).mean(span_med,:), ...
            'Color', CC(v,:), 'LineWidth', 1)
        plot([0 10], [0 0], 'k', 'LineWidth', 0.5)
      	plot([0 10], Speed(v)*[1 1], '--k', 'LineWidth', 0.5)
end
set(ax, 'LineWidth', 1)
set(ax, 'XLim', [-0.2 2])
linkaxes(ax,'x')
linkaxes(ax(1,:),'y')
linkaxes(ax(2,:),'y')
set(ax(1,:), 'YLim', 25*[-1 1])
set(ax(1,:), 'YTick', -20:10:20)
set(ax(2,:), 'YLim', [-200 250])
set(ax(2,:), 'YTick', -200:50:200)
set(ax(1:2,2:n_speed), 'YTickLabel', [])
set(ax(1:2,2:n_speed), 'YColor', 'none')
set(ax(1,:), 'XTickLabels', [])
set(ax(1,:), 'XColor', 'none')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Head Position (°)')
XLabelHC = get(ax(2,1), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Head Velocity (°/s)')

%% Intervals by fly
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 8 3])
movegui(fig, 'center')
n_std = 0;
ax = gobjects(2,n_speed); clear h
for v = 1:n_speed
    span_std = 1:round(Int_stats.all.index_end(v).median + n_std*Int_stats.all.index_end(v).std);
    span_med = 1:round(Int_stats.all.index_end(v).median + 0*Int_stats.all.index_end(v).std);
    
    ax(1,v) = subplot(2,n_speed,v); cla ; hold on
        %plot(Int_all.time{v}, Int_all.pos{v}, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
        plot(Int_stats.fly_mean.time{v}(span_std,:), Int_stats.fly_mean.pos{v}(span_std,:), ...
                            'Color', [0.5 0.5 0.5 0.7], 'LineWidth', 0.5)

        [h.pos_patch(v),h.pos_mean(v)] = PlotPatch(Int_stats.fly.pos(v).mean(span_std), ...
                                                    Int_stats.fly.pos(v).std(span_std), ...
                                                    Int_stats.fly.time(v).mean(span_std), 1, 1, ...
                                                    'k', CC(v,:), 0.3, 1);
%         plot(Int_stats.all.time(v).mean(span_med), Int_stats.all.pos(v).mean(span_med,:), ...
%             'Color', CC(v,:), 'LineWidth', 1)
        
   ax(2,v) = subplot(2,n_speed,v + n_speed); cla ; hold on
        %plot(Int_all.time{v}, Int_all.vel{v}, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
        plot(Int_stats.fly_mean.time{v}(span_std,:), Int_stats.fly_mean.vel{v}(span_std,:), ...
                            'Color', [0.5 0.5 0.5 0.7], 'LineWidth', 0.5)

        [h.pos_patch(v),h.pos_mean(v)] = PlotPatch(Int_stats.fly.vel(v).mean(span_std), ...
                                                    Int_stats.fly.vel(v).std(span_std), ...
                                                    Int_stats.fly.time(v).mean(span_std), 1, 1, ...
                                                    'k', CC(v,:), 0.3, 1);
%         plot(Int_stats.all.time(v).mean(span_med), Int_stats.all.vel(v).mean(span_med,:), ...
%             'Color', CC(v,:), 'LineWidth', 1)
        
        plot([0 10], [0 0], 'k', 'LineWidth', 0.5)
      	plot([0 10], Speed(v)*[1 1], '--k', 'LineWidth', 0.5)
end
set(ax, 'LineWidth', 1)
set(ax, 'XLim', [-0.2 1])
linkaxes(ax,'x')
linkaxes(ax(1,:),'y')
linkaxes(ax(2,:),'y')
set(ax(1,:), 'YLim', 15*[-1 1])
set(ax(1,:), 'YTick', -20:10:20)
set(ax(2,:), 'YLim', [-20 150])
set(ax(2,:), 'YTick', -200:50:200)
set(ax(1:2,2:n_speed), 'YTickLabel', [])
set(ax(1:2,2:n_speed), 'YColor', 'none')
set(ax(1,:), 'XTickLabels', [])
set(ax(1,:), 'XColor', 'none')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Head Position (°)')
XLabelHC = get(ax(2,1), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Head Velocity (°/s)')

%% Intervals Error
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 8 3])
movegui(fig, 'center')
n_std = 1;
ax = gobjects(2,n_speed); clear h
for v = 1:n_speed
    span_std = 1:round(Int_stats.all.index_end(v).median + n_std*Int_stats.all.index_end(v).std);
    span_med = 1:round(Int_stats.all.index_end(v).median + 0*Int_stats.all.index_end(v).std);
    
    ax(1,v) = subplot(2,n_speed,v); cla ; hold on
        %plot(Int_all.time{v}, Int_all.pos_error{v}, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
        [h.pos_patch(v),h.pos_mean(v)] = PlotPatch(Int_stats.all.pos_error(v).mean(span_std), ...
                                                    Int_stats.all.pos_error(v).std(span_std), ...
                                                    Int_stats.all.time(v).mean(span_std), 1, 1, ...
                                                    'k', CC(v,:), 0.3, 1);
        plot(Int_stats.all.time(v).mean(span_med), Int_stats.all.pos_error(v).mean(span_med,:), ...
            'Color', CC(v,:), 'LineWidth', 1)
        
   ax(2,v) = subplot(2,n_speed,v + n_speed); cla ; hold on
        %plot(Int_all.time{v}, Int_all.vel_error{v}, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
        [h.vel_patch(v),h.vel_mean(v)] = PlotPatch(Int_stats.all.vel_error(v).mean(span_std), ...
                                                    Int_stats.all.vel_error(v).std(span_std), ...
                                                    Int_stats.all.time(v).mean(span_std), 1, 1, ...
                                                    'k', CC(v,:), 0.3, 1);
        plot(Int_stats.all.time(v).mean(span_med), Int_stats.all.vel_error(v).mean(span_med,:), ...
            'Color', CC(v,:), 'LineWidth', 1)
        plot([0 10], [0 0], 'k', 'LineWidth', 0.5)
      	plot([0 10], Speed(v)*[1 1], '--k', 'LineWidth', 0.5)
end
set(ax, 'LineWidth', 1)
set(ax, 'XLim', [-0.2 2])
linkaxes(ax,'x')
linkaxes(ax(1,:),'y')
linkaxes(ax(2,:),'y')
set(ax(1,:), 'YLim', [-1 60])
set(ax(1,:), 'YTick', -20:10:20)
set(ax(2,:), 'YLim', [-200 250])
set(ax(2,:), 'YTick', -200:50:200)
set(ax(1:2,2:n_speed), 'YTickLabel', [])
set(ax(1:2,2:n_speed), 'YColor', 'none')
set(ax(1,:), 'XTickLabels', [])
set(ax(1,:), 'XColor', 'none')

YLabelHC = get(ax(1,1), 'YLabel');
set([YLabelHC], 'String', 'Head Position (°)')
XLabelHC = get(ax(2,1), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')
YLabelHC = get(ax(2,1), 'YLabel');
set([YLabelHC], 'String', 'Head Velocity (°/s)')

%% Interval Times
fig = figure (2); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 6 2])
clear ax h
ax = gobjects(1,1);
h = gobjects(n_speed,1);
edges = (-0.2:0.2:5)+0.001;
pp = 1;
for v = 1:n_speed
    ax(pp) = subplot(1,n_speed,v); hold on
        h(pp) = histogram(Int_all.time_end{v}, edges, 'FaceColor', CC(v,:));
    pp = pp + 1;
end
set(h, 'Normalization', 'Probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7)
linkaxes(ax, 'xy')

set(ax, 'LineWidth', 1)
axis tight
% set(ax,'YLim', [-0.01 0.6])
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Probability')
XLabelHC = get(ax(1), 'XLabel');
set([XLabelHC], 'String', 'Interval Time (s)')
set(ax(2:end),'XTickLabels',[],'YTickLabels',[])

%% Find saturation and peak times
showplot = false;
clear ax
Fc = 2;
[b,a] = butter(3, Fc/(Fs/2), 'low');
Int.fly_mean.Prop = [];
sat_thresh = 10;
for v = 1:n_speed
    for f = 1:N.fly
        % Get signal and filter
        time_nan = Int_stats.fly_mean.time{v}(:,f);
        pos_nan = Int_stats.fly_mean.pos{v}(:,f);
        vel_nan = Int_stats.fly_mean.vel{v}(:,f);
        
        endI = sum(~isnan(vel_nan));
        pos = pos_nan(1:endI);
        vel = vel_nan(1:endI);
        time = time_nan(1:endI);
        pos_filt = filtfilt(b, a, pos);
        vel_filt = filtfilt(b, a, vel);
        
        % Detect peaks
        [peaks,locs] = findpeaks(vel, 'MinPeakHeight', 10, ...
                                        'MinPeakProminence', 5, ...
                                        'SortStr', 'none');
        
       	% Find peak velocity and time in raw signal
        time_loc = time(locs);
        peak_vel = peaks(1);
        peak_time = time_loc(1);
        
        % Find peak velocity and time in filtered signal
        [peaks_filt,locs_filt] = findpeaks(vel_filt, 'MinPeakHeight', 10, ...
                                'MinPeakProminence', 5, ...
                                'SortStr', 'none');
      	time_loc_filt = time(locs_filt);
        peak_vel_filt = peaks_filt(1);
        peak_time_filt = time_loc_filt(1);
        
        % Get data after and before filtered peak
        idx_after = round(peak_time_filt*Fs):endI;
        time_after = time(idx_after);
        pos_after = pos_filt(idx_after);
        vel_filt_after = vel_filt(idx_after);
        
        % Find saturation time
        sat_idx = round(find(vel_filt_after < sat_thresh, 1, 'first'));
        if isempty(sat_idx)
            warning('true saturation not found')
            [~,sat_idx] = min(vel_filt_after);
        end
        sat_time = time_after(sat_idx);
        sat_pos = pos_after(sat_idx);
        sat_vel = vel_filt_after(sat_idx);
        
        % Find saturation time and initial time
        start_idx = round(find(vel > -5, 1, 'first'));
        start_time = time(start_idx);
        start_pos = pos(start_idx);
        start_vel = vel(start_idx);
        
        % Pull out stabilizing reigon and find mean velocity
        idx_stable = round(start_idx:sat_time*Fs);
        time_stable = nan(size(time_nan));
        vel_stable = nan(size(time_nan));
        time_stable(idx_stable) = time(idx_stable);
        vel_stable(idx_stable) = vel(idx_stable);
        mean_vel = nanmean(vel_stable);
        
        if range(pos_filt)==0
            disp('here')
        end
        
        if showplot
            figure (100)
            ax(1) = subplot(1,2,1) ; hold on ; cla
                plot(time, 0*time, 'Color', [0.5 0.5 0.5])
                plot(time, sat_thresh*ones(size(time)), 'g')

                plot(time, vel, 'k', 'LineWidth', 1)
                plot(time, vel_filt, 'c', 'LineWidth', 1)
                plot(time_after, vel_filt_after, 'b', 'LineWidth', 1)
                plot(time_stable, vel_stable, '--', 'Color', [0.7 0.1 0.2], 'LineWidth', 1)

                plot(peak_time, peak_vel, '.r', 'MarkerSize', 10)
                plot(peak_time_filt, peak_vel_filt, '.y', 'MarkerSize', 10)

                plot(start_time, start_vel, '.m', 'MarkerSize', 10)
                plot(sat_time, sat_vel, '.m', 'MarkerSize', 10)
                ylim([-50 100])

            ax(2) = subplot(1,2,2) ; hold on ; cla
                plot(time, pos, 'Color', 'k', 'LineWidth', 1)
                plot(time, pos_filt, 'Color', 'b', 'LineWidth', 1)
                plot(time_after, pos_after, '--', 'Color', 'r', 'LineWidth', 1)
                plot(start_time, start_pos, '.m', 'MarkerSize', 10)
                plot(sat_time, sat_pos, '.m', 'MarkerSize', 15)

                ylim(25*[-1 1])

            set(ax, 'LineWidth', 1, 'XLim', [0 3])

            pause
        end
      	
      	% Store properties
        Int_stats.fly_mean.Prop(v).peak_vel(f)          = peak_vel;
        Int_stats.fly_mean.Prop(v).peak_time(f)         = peak_time;
        Int_stats.fly_mean.Prop(v).sat_pos(f)           = sat_pos;
        Int_stats.fly_mean.Prop(v).sat_vel(f)           = sat_vel;
        Int_stats.fly_mean.Prop(v).sat_time(f)       	= sat_time;
        Int_stats.fly_mean.Prop(v).start_vel(f)         = start_vel;
        Int_stats.fly_mean.Prop(v).start_time(f)        = start_time;
        Int_stats.fly_mean.Prop(v).vel_stable(:,f)      = vel_stable;
        Int_stats.fly_mean.Prop(v).time_stable(:,f) 	= time_stable;
        Int_stats.fly_mean.Prop(v).pos_range(:,f)       = range(pos_filt);
        Int_stats.fly_mean.Prop(v).pos_amp(:,f)         = pos_filt(end) - start_pos;
        Int_stats.fly_mean.Prop(v).mean_vel(f)          = mean_vel;
        Int_stats.fly_mean.Prop(v).peak_gain(f)     	= peak_vel / Speed(v);
        Int_stats.fly_mean.Prop(v).mean_gain(f)     	= mean_vel / Speed(v);   
    end
end

fnames = string(fieldnames(Int_stats.fly_mean.Prop));
for v = 1:n_speed
    for f = 1:length(fnames)
        Int_stats.fly_mean.Prop_stats(v).(fnames(f)) = basic_stats(Int_stats.fly_mean.Prop(v).(fnames(f)), 2);
    end
end

%% Interval time at saturation
n_std = 1;
time_cut = nan(n_speed,1);
percent_sat = nan(n_speed,1);
for v = 1:n_speed
    int_times = Int_all.time_end{v};
    time_cut(v) = median(int_times) + n_std*std(int_times);
    percent_sat(v) = sum(int_times > time_cut(v)) / length(int_times);
end

%% Save Intervals
fname = ['Ramp_intervals_wave=' num2str(U.wave)];
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'Int', 'Int_all', 'Int_stats', 'U', 'N', 'TempFreq', 'Fc', ...
                            'cut_win', 'start_thresh', 'end_thresh', 'sat_thresh', 'Stim');

end