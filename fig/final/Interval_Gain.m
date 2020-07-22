function [] = Interval_Gain()
%% Interval_Gain:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'U','N','SACCADE')

%% Get filtered intervals
clearvars -except SACCADE U N

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
for n = 1:N.file
    head_saccade = SACCADE.head_saccade{n};
    vel = SACCADE.vel(n);
    fly = SACCADE.fly(n);
    
    pos_filt = filtfilt(b, a, head_saccade.position);
    if vel > n_speed
        vel = vel - n_speed;
        pos_filt = -pos_filt;
    end
    vel_filt = diff(pos_filt) *  head_saccade.Fs;
    vel_filt = [vel_filt(1) ; vel_filt];
    vel_filt = filtfilt(b, a, vel_filt);
    
    [~,int_pos] = getSaccade(head_saccade, pos_filt);
    [~,int_vel] = getSaccade(head_saccade, vel_filt);
    
    Int.pos{fly,vel} = [Int.pos{fly,vel}; int_pos(2:end)];
    Int.vel{fly,vel} = [Int.vel{fly,vel}; int_vel(2:end)];    
end
Int.pos = cellfun(@(x) nancat(2, x{:}), Int.pos,'UniformOutput', false);
Int.vel = cellfun(@(x) nancat(2, x{:}), Int.vel,'UniformOutput', false);
Int.time_end = cellfun(@(x) sum(~isnan(x),1)' ./ Fs, Int.pos,'UniformOutput', false);

%% Remove ends of intervals where saccades have started
win = round(0.035 * Fs); % window for mean velocity at end of interval
vthresh = -20; % cut below this velocity
n_rep = 2; %how many times to repeat the process
for w = 1:n_rep
    for v = 1:n_speed
        for f = 1:N.fly
            vel_int_all = Int.vel{f,v};
            endI = sum(~isnan(vel_int_all),1);
            for k = 1:size(endI,2)
                vel = vel_int_all(:,k);
                span = round(endI(k)-win:endI(k));

                % If we can get the window size
                if all(span >=1) && all(span <= length(vel)) 
                    endV = mean(vel(span));
                    if endV < vthresh
                        Int.pos{f,v}(span,k) = nan;
                        Int.vel{f,v}(span,k) = nan;
                    end
                end        
            end
        end
    end
end

%% Compute stats
Int.pos_stats = cellfun(@(x) basic_stats(x,2), Int.pos);
Int.vel_stats = cellfun(@(x) basic_stats(x,2), Int.vel);
Int.time_end_stats = cellfun(@(x) basic_stats(x,1), Int.time_end);

Int.fly_mean.pos = arrayfun(@(x) x.mean, Int.pos_stats,'UniformOutput',false);
Int.fly_mean.vel = arrayfun(@(x) x.mean, Int.vel_stats,'UniformOutput',false);
Int.fly_mean.time_end = arrayfun(@(x) x.mean, Int.time_end_stats,'UniformOutput',false);

% For all intervals (not by fly)
Int.time = cell(N.fly,n_speed);
Int.all.pos = cell(n_speed,1);
Int.all.vel = cell(n_speed,1);
Int.all.time_end = cell(n_speed,1);
Int.fly_mean.all.pos = cell(n_speed,1);
Int.fly_mean.all.vel = cell(n_speed,1);
Int.fly_mean.all.time_end = cell(n_speed,1);
Int.fly_mean.all.time = cell(n_speed,1);
for v = 1:n_speed
    Int.all.pos{v} = nancat(2, Int.pos{:,v});
    Int.all.vel{v} = nancat(2, Int.vel{:,v});
    Int.all.time_end{v} = nancat(2, Int.time_end{:,v});
    Int.all.time{v} = (0:size(Int.all.pos{v},1)-1)' ./ Fs;
    
    Int.fly_mean.all.pos{v} = nancat(2, Int.fly_mean.pos{:,v});
    Int.fly_mean.all.vel{v} = nancat(2, Int.fly_mean.vel{:,v});
    Int.fly_mean.all.time_end{v} = nancat(2, Int.fly_mean.time_end{:,v})';
    Int.fly_mean.all.time{v} = (0:size(Int.fly_mean.all.pos{v},1)-1)' ./ Fs;
    
    for f = 1:N.fly
        Int.time{f,v} = (0:size(Int.pos{f,v},1)-1)' ./ Fs;
    end
end
Int.all.pos_stats = cellfun(@(x) basic_stats(x,2), Int.all.pos);
Int.all.vel_stats = cellfun(@(x) basic_stats(x,2), Int.all.vel);
Int.all.time_end_stats = cellfun(@(x) basic_stats(x,1), Int.all.time_end);

Int.fly_mean.all.pos_stats = cellfun(@(x) basic_stats(x,2), Int.fly_mean.all.pos);
Int.fly_mean.all.vel_stats = cellfun(@(x) basic_stats(x,2), Int.fly_mean.all.vel);
Int.fly_mean.all.time_end_stats = cellfun(@(x) basic_stats(x,1), Int.fly_mean.all.time_end);

%% Intervals all
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 8 3])
movegui(fig, 'center')
n_std = 1;
ax = gobjects(2,n_speed); clear h
for v = 1:n_speed
    span_std = 1:Fs*(Int.all.time_end_stats(v).median + n_std*Int.all.time_end_stats(v).std);
    span_med = 1:Fs*(Int.all.time_end_stats(v).median + 0*Int.all.time_end_stats(v).std);
    
    ax(1,v) = subplot(2,n_speed,v); cla ; hold on
        plot(Int.all.time{v}, Int.all.pos{v}, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
        %plot(Int.fly_mean.all.time{v}(span_std), Int.fly_mean.all.pos{v}(span_std,:), ...
                            %'Color', [0.5 0.5 0.5 0.7], 'LineWidth', 0.5)

        [h.pos_patch(v),h.pos_mean(v)] = PlotPatch(Int.all.pos_stats(v).mean(span_std), ...
                                                    Int.all.pos_stats(v).std(span_std), ...
                                                    Int.time{v}(span_std), 1, 1, 'k', CC(v,:), 0.3, 1);
        %delete(h.pos_mean(v))
        plot(Int.all.time{v}(span_med), Int.all.pos_stats(v).mean(span_med,:), ...
            'Color', CC(v,:), 'LineWidth', 1.5)
        
   ax(2,v) = subplot(2,n_speed,v + n_speed); cla ; hold on
        plot(Int.all.time{v}, Int.all.vel{v}, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
        %plot(Int.fly_mean.all.time{v}(span_std), Int.fly_mean.all.vel{v}(span_std,:), ...
                            %'Color', [0.5 0.5 0.5 0.7], 'LineWidth', 0.5)

        [h.vel_patch(v),h.vel_mean(v)] = PlotPatch(Int.all.vel_stats(v).mean(span_std), ...
                                                    Int.all.vel_stats(v).std(span_std), ...
                                                    Int.time{v}(span_std), 1, 1, 'k', CC(v,:), 0.3, 1);
        %delete(h.vel_mean(v))
        plot(Int.all.time{v}(span_med), Int.all.vel_stats(v).mean(span_med,:), ...
            'Color', CC(v,:), 'LineWidth', 1.5)
        plot([0 10], [0 0], 'k', 'LineWidth', 0.5)
      	plot([0 10], Speed(v)*[1 1], 'k', 'LineWidth', 0.5)
end
set(ax, 'LineWidth', 1)
set(ax, 'XLim', [-0.2 2])
linkaxes(ax,'x')
linkaxes(ax(1,:),'y')
linkaxes(ax(2,:),'y')
set(ax(1,:), 'YLim', 25*[-1 1])
set(ax(1,:), 'YTick', -20:10:20)
set(ax(2,:), 'YLim', 350*[-1 1])
set(ax(2,:), 'YTick', -200:50:200)

XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (s)')
YLabelHC = get(ax(1,:), 'YLabel');
set([YLabelHC{:}], 'String', 'Head Position (°)')
YLabelHC = get(ax(2,:), 'YLabel');
set([YLabelHC{:}], 'String', 'Head Velocity (°/s)')

%% Intervals by fly
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 8 3])
movegui(fig, 'center')
n_std = 1;
ax = gobjects(2,n_speed); clear h
for v = 1:n_speed
    span_std = 1:Fs*(Int.fly_mean.all.time_end_stats(v).mean + n_std*Int.fly_mean.all.time_end_stats(v).std);
    span_med = 1:Fs*(Int.fly_mean.all.time_end_stats(v).mean + 0*Int.fly_mean.all.time_end_stats(v).std);
    
    ax(1,v) = subplot(2,n_speed,v); cla ; hold on
        %plot(Int.all.time{v}, Int.all.pos{v}, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
        plot(Int.fly_mean.all.time{v}(span_med), Int.fly_mean.all.pos{v}(span_med,:), ...
                            'Color', [0.5 0.5 0.5 0.7], 'LineWidth', 0.5)

        [h.pos_patch(v),h.pos_mean(v)] = PlotPatch(Int.fly_mean.all.pos_stats(v).mean(span_std), ...
                                                    Int.fly_mean.all.pos_stats(v).std(span_std), ...
                                                    Int.fly_mean.all.time{v}(span_std), ...
                                                    1, 1, 'k', CC(v,:), 0.3, 1);
        plot(Int.all.time{v}(span_med), Int.all.pos_stats(v).mean(span_med,:), ...
            'Color', CC(v,:), 'LineWidth', 1.5)
        
   ax(2,v) = subplot(2,n_speed,v + n_speed); cla ; hold on
        %plot(Int.all.time{v}, Int.all.vel{v}, 'Color', [0.5 0.5 0.5 0.2], 'LineWidth', 0.25)
        plot(Int.fly_mean.all.time{v}(span_med), Int.fly_mean.all.vel{v}(span_med,:), ...
                            'Color', [0.5 0.5 0.5 0.7], 'LineWidth', 0.5)

        [h.vel_patch(v),h.vel_mean(v)] = PlotPatch(Int.fly_mean.all.vel_stats(v).mean(span_std), ...
                                                    Int.fly_mean.all.vel_stats(v).std(span_std), ...
                                                    Int.fly_mean.all.time{v}(span_std), ...
                                                    1, 1, 'k', CC(v,:), 0.3, 1);
        plot(Int.all.time{v}(span_med), Int.all.vel_stats(v).mean(span_med,:), ...
            'Color', CC(v,:), 'LineWidth', 1.5)
        plot([0 10], [0 0], 'k', 'LineWidth', 0.5)
end
set(ax, 'LineWidth', 1)
set(ax, 'XLim', [-0.2 2])
linkaxes(ax,'x')
linkaxes(ax(1,:),'y')
linkaxes(ax(2,:),'y')
set(ax(1,:), 'YLim', 20*[-1 1])
set(ax(1,:), 'YTick', -20:10:20)
set(ax(2,:), 'YLim', 150*[-1 1])
set(ax(2,:), 'YTick', -150:50:150)

XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (s)')
YLabelHC = get(ax(1,:), 'YLabel');
set([YLabelHC{:}], 'String', 'Head Position (°)')
YLabelHC = get(ax(2,:), 'YLabel');
set([YLabelHC{:}], 'String', 'Head Velocity (°/s)')

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
        h(pp) = histogram(Int.all.time_end{v}, edges, 'FaceColor', CC(v,:));
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
Fc = 2;
[b,a] = butter(3, Fc/(Fs/2), 'low');
Int.fly_mean.Prop = [];
sat_thresh = 10;
for v = 1:n_speed
    for f = 1:N.fly
        % Get signal and filter
        time_nan = Int.fly_mean.all.time{v};
        vel_nan = Int.fly_mean.all.vel{v}(:,f);
        
        endI = sum(~isnan(vel_nan));
        vel = vel_nan(1:endI);
        time = time_nan(1:endI);
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
        vel_filt_after = vel_filt(idx_after);
        
        % Find saturation time
        sat_idx = round(find(vel_filt_after < sat_thresh, 1, 'first'));
        if isempty(sat_idx)
            warning('true saturation not found')
            [~,sat_idx] = min(vel_filt_after);
        end
        sat_time = time_after(sat_idx);
        sat_vel = vel_filt_after(sat_idx);
        
        % Find saturation time and initial time
        start_idx = round(find(vel > 0, 1, 'first'));
        start_time = time(start_idx);
        start_vel = vel(start_idx);
        
        % Pull out stabilizing reigon and find mean velocity
        idx_stable = round(start_idx:sat_time*Fs);
        time_stable = nan(size(time_nan));
        vel_stable = nan(size(time_nan));
        time_stable(idx_stable) = time(idx_stable);
        vel_stable(idx_stable) = vel(idx_stable);
        mean_vel = nanmean(vel_stable);
        
        if showplot
        figure (100)
        ax = subplot(1,1,1) ; hold on ; cla ; set(ax, 'LineWidth', 1)
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
            
            
            xlim([0 3])
            ylim([-50 100])
            
            pause
        end
      	
      	% Store properties
        Int.fly_mean.Prop(v).peak_vel(f)        = peak_vel;
        Int.fly_mean.Prop(v).peak_time(f)       = peak_time;
        Int.fly_mean.Prop(v).sat_vel(f)         = sat_vel;
        Int.fly_mean.Prop(v).sat_time(f)        = sat_time;
        Int.fly_mean.Prop(v).start_vel(f)       = start_vel;
        Int.fly_mean.Prop(v).start_time(f)      = start_time;
        Int.fly_mean.Prop(v).vel_stable(:,f)  	= vel_stable;
        Int.fly_mean.Prop(v).time_stable(:,f) 	= time_stable;
        Int.fly_mean.Prop(v).mean_vel(f)        = mean_vel;
        Int.fly_mean.Prop(v).peak_gain(f)     	= peak_vel / Speed(v);
        Int.fly_mean.Prop(v).mean_gain(f)     	= mean_vel / Speed(v);   
    end
end

fnames = string(fieldnames(Int.fly_mean.Prop));
for v = 1:n_speed
    for f = 1:length(fnames)
        Int.fly_mean.Prop_stats(v).(fnames(f)) = basic_stats(Int.fly_mean.Prop(v).(fnames(f)), 2);
    end
end

%% Interval Stats
stat_names = ["sat_time","mean_vel","peak_vel","mean_gain","peak_gain","peak_time"];
n_plot = length(stat_names);
ylim_list = [1.5 150 150 1 2 0.15];

fig = figure (4) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2*n_plot 2])
clear ax h
ax = gobjects(n_plot,1);
clear h
for kk = 1:n_plot
    ax(kk) = subplot(1,n_plot,kk) ; hold on
        data = cat(1,Int.fly_mean.Prop.(stat_names(kk)))';
        bx = boxplot(data, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
        
        ylabel(stat_names(kk), 'interpreter', 'none')

        h = get(bx(5,:),{'XData','YData'});
        for v = 1:size(h,1)
           patch(h{v,1},h{v,2},CC(v,:));
        end

        set(findobj(ax(kk),'tag','Median'), 'Color', 'w','LineWidth',1.5);
        set(findobj(ax(kk),'tag','Box'), 'Color', 'none');
        set(findobj(ax(kk),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
        set(findobj(ax(kk),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
        ax(kk).Children = ax(kk).Children([end 1:end-1]); 
        ax(kk).YLim(1) = 0;
        ax(kk).YLim(2) = ylim_list(kk);
end
set(ax, 'LineWidth', 1, 'Box', 'off')

%% Save Intervals
fname = ['Ramp_Intervals_wave=' num2str(U.wave)];
savedir = 'C:\Users\BC\Box\Research\Manuscripts\Head Saccade\Data';
save(fullfile(savedir, [fname '.mat']), 'Int', 'U', 'N', 'TempFreq');

end