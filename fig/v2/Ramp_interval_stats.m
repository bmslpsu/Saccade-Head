function [] = Ramp_interval_stats()
%% Ramp_interval_stats:
root = 'E:\DATA\Rigid_Data\Saccade\processed';

[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','on');
FILE = cellstr(FILE);

n_wave = length(FILE);
wave_data = cell(n_wave,1);
for n = 1:n_wave
    wave_data{n} = load(fullfile(PATH,FILE{n}));
end
n_speed = 5;
Speed = wave_data{1}.U.vel{1}(1:n_speed);
vel_color = hsv(n_speed);
wave_color = parula(n_wave);

%% Mean Gain
fig = figure(1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 2])
clear h s
wave = nan(n_wave,1);
ax = subplot(1,1,1); hold on
for n = 1:n_wave
    wave(n) = wave_data{n}.U.wave;
    x = wave_data{n}.TempFreq;
    %x = Speed;
    data = [wave_data{n}.Int_stats.fly_mean.Prop_stats.mean_gain];
    gain_mean = [data.mean];
    gain_std = [data.std];
    [h.patch(n),h.mean(n)] = PlotPatch(gain_mean, gain_std, x, 1, 1,  ...
                            wave_color(n,:), 0.5*wave_color(n,:), 0.3, 2);
	s(n)= scatter(x, gain_mean, 200, vel_color, '.');
                        
end
uistack(h.mean, 'top')
uistack(s, 'top')
leg = legend([h.mean], string(wave), 'Box', 'off');
leg.Title.String = 'Spatial Wavength (°)';

set(ax , 'LineWidth', 1, 'YLim', [0 1])
xlabel('Temporal Frequency (Hz)')
ylabel('Mean Gain (°s^{-1} / °s^{-1})')

%% Peak Gain
fig = figure(2); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 2])
clear h s
wave = nan(n_wave,1);
ax = subplot(1,1,1); hold on
for n = 1:n_wave
    wave(n) = wave_data{n}.U.wave;
    x = wave_data{n}.TempFreq;
    %x = Speed;
    data = [wave_data{n}.Int_stats.fly_mean.Prop_stats.peak_gain];
    gain_mean = [data.mean];
    gain_std = [data.std];
    [h.patch(n),h.mean(n)] = PlotPatch(gain_mean, gain_std, x, 1, 1,  ...
                            wave_color(n,:), 0.5*wave_color(n,:), 0.3, 2);
	s(n)= scatter(x, gain_mean, 200, vel_color, '.');
                        
end
uistack(h.mean, 'top')
uistack(s, 'top')
leg = legend([h.mean], string(wave), 'Box', 'off');
leg.Title.String = 'Spatial Wavength (°)';

set(ax , 'LineWidth', 1, 'YLim', [0 2.2])
xlabel('Temporal Frequency (Hz)')
ylabel('Peak Gain (°s^{-1} / °s^{-1})')

%% Mean Velocity
fig = figure(3); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 2])
clear h s
wave = nan(n_wave,1);
ax = subplot(1,1,1); hold on
for n = 1:n_wave
    wave(n) = wave_data{n}.U.wave;
    x = wave_data{n}.TempFreq;
    %x = Speed;
    data = [wave_data{n}.Int_stats.fly_mean.Prop_stats.mean_vel];
    gain_mean = [data.mean];
    gain_std = [data.std];
    [h.patch(n),h.mean(n)] = PlotPatch(gain_mean, gain_std, x, 1, 1,  ...
                            wave_color(n,:), 0.5*wave_color(n,:), 0.3, 2);
	s(n)= scatter(x, gain_mean, 200, vel_color, '.');
                        
end
uistack(h.mean, 'top')
uistack(s, 'top')
leg = legend([h.mean], string(wave), 'Box', 'off');
leg.Title.String = 'Spatial Wavength (°)';

set(ax , 'LineWidth', 1, 'YLim', [0 60])
ax.YTick = 0:30:60;
xlabel('Temporal Frequency (Hz)')
ylabel('Mean Velocity (°s^{-1})')

%% Peak Velocity
fig = figure(4); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 2])
clear h s
wave = nan(n_wave,1);
ax = subplot(1,1,1); hold on
for n = 1:n_wave
    wave(n) = wave_data{n}.U.wave;
    x = wave_data{n}.TempFreq;
    %x = Speed;
    data = [wave_data{n}.Int_stats.fly_mean.Prop_stats.peak_vel];
    gain_mean = [data.mean];
    gain_std = [data.std];
    [h.patch(n),h.mean(n)] = PlotPatch(gain_mean, gain_std, x, 1, 1,  ...
                            wave_color(n,:), 0.5*wave_color(n,:), 0.3, 2);
	s(n)= scatter(x, gain_mean, 200, vel_color, '.');
                        
end
uistack(h.mean, 'top')
uistack(s, 'top')
leg = legend([h.mean], string(wave), 'Box', 'off');
leg.Title.String = 'Spatial Wavength (°)';

set(ax , 'LineWidth', 1, 'YLim', [0 150])
ax.YTick = 0:30:150;
xlabel('Temporal Frequency (Hz)')
ylabel('Peak Velocity (°s^{-1})')

%% Interval Stats
stat_names = ["peak_time","peak_vel","sat_time","mean_vel","mean_gain","peak_gain","pos_amp"];
n_plot = length(stat_names);

clear int_stats
for kk = 1:n_plot
    int_stats.(stat_names(kk)) = [];
    for n = 1:n_wave
        data = cat(1, wave_data{n}.Int_stats.fly_mean.Prop.(stat_names(kk)))';
        int_stats.(stat_names(kk)) = [int_stats.(stat_names(kk)) ; data];
    end
end

ylim_list = [0 0.1;0 150;0 1.5;0 150;0 1;0 2.5;0 30];

fig = figure (5) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2*n_plot 1.5])
clear ax h
ax = gobjects(n_plot,1);
clear h
for kk = 1:n_plot
    ax(kk) = subplot(1,n_plot,kk) ; hold on
        bx = boxplot(int_stats.(stat_names(kk)), 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
        
        ylabel(stat_names(kk), 'interpreter', 'none')

        h = get(bx(5,:),{'XData','YData'});
        for v = 1:size(h,1)
           patch(h{v,1},h{v,2},vel_color(v,:), 'EdgeColor', 'none');
        end

        set(findobj(ax(kk),'tag','Median'), 'Color', 'k','LineWidth',0.5);
        set(findobj(ax(kk),'tag','Box'), 'Color', 'none');
        set(findobj(ax(kk),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
        set(findobj(ax(kk),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
        ax(kk).Children = ax(kk).Children([end 1:end-1]); 
        ax(kk).YLim = ylim_list(kk,:);
        
        if any(kk == [2,4])
            plot(1:n_speed, Speed, '--k')
        end
end
set(ax, 'LineWidth', 1, 'Box', 'off')

%% Stats
% [p,tb,stats] = anovan(int_stats.vel_error_end, G, 'varnames', {'Speed'});
% [p,tb,stats] = anovan(int_stats.int_vel_error_end, G, 'varnames', {'Speed'});
% [c,m] = multcompare(stats);

% % Kruskalwallis
% data = int_stats.peak_time;
% data = reshap(data, )
[p,tb,stats] = kruskalwallis(int_stats.pos_amp);
[c,m] = multcompare(stats, 'Alpha', 0.001);

%% Error Stats
stat_names = ["pos_error_end","vel_error_end","int_pos_error_end","int_vel_error_end"];
n_plot = length(stat_names);

clear int_stats
for kk = 1:n_plot
    int_stats.(stat_names(kk)) = [];
    G = [];
    for n = 1:n_wave
        data = wave_data{n}.Int_all.(stat_names(kk));
        g = cellfun(@(x,y) y*ones(length(x),1), data, num2cell(1:length(data)), 'UniformOutput', false);
        data = cat(2, data{:})';
        g = cat(1, g{:});
        int_stats.(stat_names(kk)) = [int_stats.(stat_names(kk)) ; data];
        G = [G ; g];
    end
end

ylim_list = [-50 200;-50 300;-10 40;-50 200];

fig = figure (6) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2*n_plot 1.5])
clear ax h
ax = gobjects(n_plot,1);
clear h
for kk = 1:n_plot
    ax(kk) = subplot(1,n_plot,kk) ; hold on
        bx = boxplot(int_stats.(stat_names(kk)), G, 'Labels', {Speed}, ...
            'Width', 0.5, 'Symbol', '', 'Whisker', 2, 'OutlierSize', 0.5);
        
        ylabel(stat_names(kk), 'interpreter', 'none')

        h = get(bx(5,:),{'XData','YData'});
        for v = 1:size(h,1)
           patch(h{v,1},h{v,2},vel_color(v,:), 'EdgeColor', 'none');
        end

        set(findobj(ax(kk),'tag','Median'), 'Color', 'k', 'LineWidth', 0.5);
        set(findobj(ax(kk),'tag','Box'), 'Color', 'none');
        set(findobj(ax(kk),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
        set(findobj(ax(kk),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
        ax(kk).Children = ax(kk).Children([end 1:end-1]); 
        ax(kk).YLim = ylim_list(kk,:);
        
        if any(kk == [2])
            plot(1:n_speed, Speed, '--k')
        end
end
set(ax, 'LineWidth', 1, 'Box', 'off')

%% Stats
% [p,tb,stats] = anovan(int_stats.vel_error_end, G, 'varnames', {'Speed'});
% [p,tb,stats] = anovan(int_stats.int_vel_error_end, G, 'varnames', {'Speed'});
% [c,m] = multcompare(stats);

% % Kruskalwallis
[p,tb,stats] = kruskalwallis(int_stats.int_vel_error_end, G);
[c,m] = multcompare(stats, 'Alpha', 0.001);

%% Interval Times Histogram
fig = figure (6); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 6 2])
clear ax h
ax = gobjects(1,1);
h = gobjects(n_speed,1);
edges = (-0.2:0.1:8)+0.001;
pp = 1;
int_times = cell(n_speed,1);
for v = 1:n_speed
    for n = 1:n_wave
        data = wave_data{n}.Int_all.time_end{v}(:);
        data = data(~isnan(data));
        int_times{v} = [int_times{v} ; data];
    end
    
    ax(pp) = subplot(1,n_speed,v); hold on
        h(pp) = histogram(int_times{v}, edges, 'FaceColor', vel_color(v,:));
    
    pp = pp + 1;
end
set(h, 'Normalization', 'Probability', 'EdgeColor', 'none', 'FaceAlpha', 0.7)
linkaxes(ax, 'xy')

set(ax, 'LineWidth', 1)
axis tight
set(ax,'YLim', [-0.005 0.11])
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Probability')
XLabelHC = get(ax(1), 'XLabel');
set([XLabelHC], 'String', 'Interval Time (s)')
set(ax(2:end),'XTickLabels',[],'YTickLabels',[])

%% Interval Times Violin
fig = figure (7); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 7 4])
clear ax h bw
ax = gobjects(1,1);
h = gobjects(n_speed,1);
bw = nan(n_speed,1);
pp = 1;
int_times = cell(n_speed,1);
for v = 1:n_speed
    for n = 1:n_wave
        data = wave_data{n}.Int_all.time_end{v}(:);
        data = data(~isnan(data));
        int_times{v} = [int_times{v} ; data];
    end
    
    ax(pp) = subplot(1,n_speed,v); hold on
        [h(pp),L,~,~,bw(pp)] = violin(int_times{v}, 'facecolor', vel_color(v,:), ...
                                    'mc', 'k', 'medc', 'g');
        if v > 1
            delete(L)
        end
        %h(pp).FaceColor = vel_color(v,:);
    
    pp = pp + 1;
end
set(h, 'EdgeColor', 'k', 'FaceAlpha', 0.7)
linkaxes(ax, 'xy')

set(ax, 'LineWidth', 1, 'Box', 'off')
axis tight
set(ax,'YLim', [-0.5 3])
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Interval Times (s)')
XLabelHC = get(ax(1), 'XLabel');
set([XLabelHC], 'String', 'Kernel Distribution')
set(ax(2:end),'XTickLabels',[],'YTickLabels',[])

%% Percent Saturate
clc
time_med = cellfun(@(x) median(x), int_times);
time_std = cellfun(@(x) std(x), int_times);
n_mad = 3;
time_mad = cellfun(@(x) mad(x), int_times);
time_cut = time_med + n_mad *time_mad;

time_sat = time_med + 2*time_std;
time_rmv = cellfun(@(x) rmoutliers(x, 'median', 'ThresholdFactor', n_mad), int_times, 'UniformOutput', false);
sat_ratio = 1 - ( cellfun(@(x) length(x), time_rmv) ./ cellfun(@(x) length(x), int_times) );
disp(sat_ratio)
time_cut =  cellfun(@(x) max(x), time_rmv);
% n_std = 3;
% time_cut = nan(n_speed,1);
% percent_sat = nan(n_speed,1);
% for v = 1:n_speed
%     I = int_times{v};
%     time_cut(v) = median(I) + n_std*std(I);
%     percent_sat(v) = sum(I > time_cut(v)) / length(I);
% end

fig = figure (17); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 8])
clear ax h
edges = 0:0.05:3;
for v = 1:n_speed
    ax(v) = subplot(n_speed,1,v); cla ; hold on ; title(Speed(v))
        histogram(int_times{v}, edges, 'Normalization', 'probability', ...
            'FaceColor', 'r', 'FaceAlpha', 1, 'EdgeColor', 'none')
        histogram(time_rmv{v}, edges, 'Normalization', 'probability', ...
            'FaceColor', 'k', 'FaceAlpha', 1, 'EdgeColor', 'none')
end
legend('Outlier', 'All')

%% Amplitude
int_amp = cell(n_wave,1);
int_amp_fly = cell(n_wave,1);
for w = 1:n_wave
    for v = 1:wave_data{w}.N.vel/2
        for f = 1:wave_data{w}.N.fly
            data = wave_data{w}.Int.pos{f,v};
            endI = sum(~isnan(data),1);
            endI = endI(endI ~= 0);
            for c = 1:size(endI,2)
                int_amp{w}{f,v}(1,c) = data(endI(c),c) - data(1,c);
                %int_amp{w}{f,v}(1,c) = data(endI(c),c);
            end
        end
    end
    int_amp_fly{w} = cellfun(@(x) mean(x,2), int_amp{w}, 'UniformOutput', true);
end
int_amp_fly = cat(1, int_amp_fly{:});

fig = figure (8) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2*1 1.5])
clear ax h

kk = 1;
ax(kk) = subplot(1,1,kk) ; hold on
    bx = boxplot(int_amp_fly, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    ylabel('Amplitude (°)')

    h = get(bx(5,:),{'XData','YData'});
    for v = 1:size(h,1)
       patch(h{v,1},h{v,2},vel_color(v,:), 'EdgeColor', 'none');
    end

    set(findobj(ax(kk),'tag','Median'), 'Color', 'w','LineWidth', 1);
    set(findobj(ax(kk),'tag','Box'), 'Color', 'none');
    set(findobj(ax(kk),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(kk),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(kk).Children = ax(kk).Children([end 1:end-1]);
    ylim([0 25])

set(ax, 'LineWidth', 1, 'Box', 'off')

%% Stats
% [p,tb,stats] = anovan(int_stats.vel_error_end, G, 'varnames', {'Speed'});
% [p,tb,stats] = anovan(int_stats.int_vel_error_end, G, 'varnames', {'Speed'});
% [c,m] = multcompare(stats);

% % Kruskalwallis
[p,tb,stats] = kruskalwallis(int_amp_fly);
[c,m] = multcompare(stats, 'Alpha', 0.001);

end