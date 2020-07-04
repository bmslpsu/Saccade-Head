%% Example detetction of wing saccades from WBA
clear ; close all ; clc
% load ('Example_1') % note this data has clear periods of fixation and saccades
% fly 11 tiral 25
root = 'H:\EXPERIMENTS\RIGID\Experiment_Asymmetry_Control_Verification\HighContrast\30\Vid\wing_filt\tracked_head_wing';
[FILE,PATH] = uigetfile({'*.csv'},'Select data file', root, 'MultiSelect','off');
benifly_data = ImportBenifly(fullfile(PATH, FILE)); % load head & wing angles

%%
time = linspace(0,10,size(benifly_data,1)); % time vector
Fs = round(1 / mean(diff(time))); % sampling frequency [Hz]
wing_left = rad2deg(hampel(time, benifly_data.LWing)); % left wing
wing_right = rad2deg(hampel(time, benifly_data.RWing)); % right wing
dwba = wing_left - wing_right; % delta-wba

% Filter
Fc = 3; % low-pass
[b,a] = butter(3, Fc/(Fs/2), 'low');
wing_left_filt = filtfilt(b, a, wing_left); % left wing
wing_right_filt = filtfilt(b, a, wing_right); % right wing
filt_dwba = filtfilt(b, a, wing_left_filt - wing_right_filt);

Fc = 0.3; % high-pass
[b,a] = butter(3, Fc/(Fs/2), 'high');
filt_dwba = filtfilt(b, a, filt_dwba);

wba_vel = diff(filt_dwba) * Fs;
wba_vel = [wba_vel(1) ; wba_vel];

wba_vel_norm = max(filt_dwba) * wba_vel / max(wba_vel);

%% Find peaks (need to tune these properties too)
min_peak_dist_time = 0.5; % [s]
min_peak_width_time = 5; % [s]
[vel_pks,pks_idx] = findpeaks(wba_vel, ... 
                        'MinPeakProminence', 50, ...
                        'MinPeakHeight', 25, ...
                        'MinPeakDistance', min_peak_dist_time*Fs, ...
                        'MaxPeakWidth', min_peak_width_time*Fs);

filt_dwba_pks = filt_dwba(pks_idx);
dwba_pks = dwba(pks_idx);
time_pks = time(pks_idx);

[vel_vly,vly_idx] = findpeaks(-wba_vel, ...
                        'MinPeakProminence', 5, ...
                        'MinPeakDistance', min_peak_dist_time*Fs, ...
                        'MaxPeakWidth', min_peak_width_time*Fs);

vel_vly = -vel_vly;
filt_dwba_vly = filt_dwba(vly_idx);
dwba_vly = dwba(vly_idx);
time_vly = time(vly_idx);

% Plot signals, red markers are saccade peaks
fig = figure (1);
set(fig, 'Color', 'w')
clear ax
ax(1) = subplot(3,1,1) ; cla ; hold on ; title('Raw')
    ylabel('\DeltaWBA (V)')
    plot([0 time(end)], [0 0], '--', 'Color', [0.5 0.5 0.5])
    plot(time, dwba, 'k')
    plot(time_pks, dwba_pks, '.r', 'MarkerSize', 15)
    plot(time_vly, dwba_vly, '.b', 'MarkerSize', 15)
    
ax(2) = subplot(3,1,2) ; cla ; hold on ; title('Filtered')
    ylabel('\DeltaWBA (V)')
    plot([0 time(end)], [0 0], '--', 'Color', [0.5 0.5 0.5])
    plot(time, filt_dwba, 'k', 'LineWidth', 1)
    plot(time_pks, filt_dwba_pks, '.r', 'MarkerSize', 15)
  	plot(time_vly, filt_dwba_vly, '.b', 'MarkerSize', 15)
    
ax(3) = subplot(3,1,3) ; cla ; hold on ; title('Absolute Value')
    ylabel('\DeltaWBA (V/s)')
    plot([0 time(end)], [0 0], '--', 'Color', [0.5 0.5 0.5])
    plot(time, wba_vel, 'k')
    plot(time_pks, vel_pks, '.r', 'MarkerSize', 15)
    plot(time_vly, vel_vly, '.b', 'MarkerSize', 15)
    
    xlabel('Time (s)')

set(ax, 'LineWidth', 1)
linkaxes(ax,'x')

% Extract saccades & inter-saccade intervals, align to peak time
thresh = 0;
amp_cut = 0;
direction = 0;
sacd_length = nan;
showplot = true;
obj = saccade(filt_dwba, time, thresh, amp_cut, direction, pks_idx, sacd_length, showplot);
