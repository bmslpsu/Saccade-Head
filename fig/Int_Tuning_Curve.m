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
n_vel = 5;
Vel = wave_data{1}.U.vel{1}(1:n_vel);
vel_color = hsv(n_vel);
wave_color = parula(n_wave);

%% Mean Gain
fig = figure(1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 2])
clear h s
wave = nan(n_wave,1);
ax = subplot(1,1,1); hold on
for n = 1:n_wave
    wave(n) = wave_data{n}.U.wave;
    x =  wave_data{n}.Temp_Freq;
    gain_mean = [wave_data{n}.Mean_Gain_Stats.mean];
    gain_std = [wave_data{n}.Mean_Gain_Stats.std]; 
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
    x =  wave_data{n}.Temp_Freq;
    gain_mean = [wave_data{n}.Peak_Gain_Stats.mean];
    gain_std = [wave_data{n}.Peak_Gain_Stats.std]; 
    [h.patch(n),h.mean(n)] = PlotPatch(gain_mean, gain_std, x, 1, 1,  ...
                            wave_color(n,:), 0.5*wave_color(n,:), 0.3, 2);
	s(n)= scatter(x, gain_mean, 200, vel_color, '.');
                        
end
uistack(h.mean, 'top')
uistack(s, 'top')
leg = legend([h.mean], string(wave), 'Box', 'off');
leg.Title.String = 'Spatial Wavength (°)';

set(ax , 'LineWidth', 1, 'YLim', [0 1.6])
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
    x =  wave_data{n}.Temp_Freq;
    vel_mean = [wave_data{n}.Mean_Vel_Stats.mean];
    vel_std = [wave_data{n}.Mean_Vel_Stats.std]; 
    [h.patch(n),h.mean(n)] = PlotPatch(vel_mean, vel_std, x, 1, 1,  ...
                            wave_color(n,:), 0.5*wave_color(n,:), 0.3, 2);
	s(n)= scatter(x, vel_mean, 200, vel_color, '.');
                        
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
    x =  wave_data{n}.Temp_Freq;
    vel_mean = [wave_data{n}.Peak_Vel_Stats.mean];
    vel_std = [wave_data{n}.Peak_Vel_Stats.std]; 
    [h.patch(n),h.mean(n)] = PlotPatch(vel_mean, vel_std, x, 1, 1,  ...
                            wave_color(n,:), 0.5*wave_color(n,:), 0.3, 2);
	s(n)= scatter(x, vel_mean, 200, vel_color, '.');
                        
end
uistack(h.mean, 'top')
uistack(s, 'top')
leg = legend([h.mean], string(wave), 'Box', 'off');
leg.Title.String = 'Spatial Wavength (°)';

set(ax , 'LineWidth', 1, 'YLim', [0 90])
ax.YTick = 0:30:90;
xlabel('Temporal Frequency (Hz)')
ylabel('Peak Velocity (°s^{-1})')

end