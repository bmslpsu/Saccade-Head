function [] = Int_Tuning_Curve()
%% Int_Tuning_Curve:
root = 'C:\Users\BC\Box\Research\Manuscripts\Head Saccade\Data';

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
    data = [wave_data{n}.Int.fly_mean.Prop_stats.mean_gain];
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
    data = [wave_data{n}.Int.fly_mean.Prop_stats.peak_gain];
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

set(ax , 'LineWidth', 1, 'YLim', [0 2])
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
    data = [wave_data{n}.Int.fly_mean.Prop_stats.mean_vel];
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
    data = [wave_data{n}.Int.fly_mean.Prop_stats.peak_vel];
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
stat_names = ["peak_time","peak_vel","sat_time","mean_vel","mean_gain","peak_gain",];
n_plot = length(stat_names);

clear int_stats
for kk = 1:n_plot
    int_stats.(stat_names(kk)) = [];
    for n = 1:n_wave
        data = cat(1, wave_data{n}.Int.fly_mean.Prop.(stat_names(kk)))';
        int_stats.(stat_names(kk)) = [int_stats.(stat_names(kk)) ; data];
    end
end

ylim_list = [0 0.15;0 150;0 1.5;0 150;0 1;0 2.5];

fig = figure (5) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2*n_plot 1.5])
clear ax h
ax = gobjects(n_plot,1);
clear h
for kk = 1:n_plot
    ax(kk) = subplot(1,n_plot,kk) ; hold on
        bx = boxplot(int_stats.(stat_names(kk)), 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
        
        ylabel(stat_names(kk), 'interpreter', 'none')

        h = get(bx(5,:),{'XData','YData'});
        for v = 1:size(h,1)
           patch(h{v,1},h{v,2},vel_color(v,:));
        end

        set(findobj(ax(kk),'tag','Median'), 'Color', 'w','LineWidth',0.5);
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
        data = wave_data{n}.Int.all.time_end{v}(:);
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
        data = wave_data{n}.Int.all.time_end{v}(:);
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

end