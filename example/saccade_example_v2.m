%% Example detetction of wing saccades from WBA
clear ; close all ; clc
% load ('Example_1') % note this data has clear periods of fixation and saccades

% time = t_p; % time vector
% Fs = 1 / mean(diff(time)); % sampling frequency [Hz]
% wing_left = data(:,4); % left wing
% wing_right = data(:,5); % right wing
% dwba = wing_left - wing_right; % delta-wba

root = 'H:\EXPERIMENTS\RIGID\Experiment_Asymmetry_Control_Verification\HighContrast\30\Vid\vid_filt\tracked_head_wing';
[FILE,PATH] = uigetfile({'*.csv'},'Select data file', root, 'MultiSelect','off');
benifly_data = ImportBenifly(fullfile(PATH, FILE)); % load head & wing angles
%%
time = linspace(0,10,size(benifly_data,1)); % time vector
Fs = 1 / mean(diff(time)); % sampling frequency [Hz]
wing_left = benifly_data.LWing; % left wing
wing_right = benifly_data.RWing; % right wing
dwba = wing_left - wing_right; % delta-wba

% Filter (tune these based on data properties)
Fc = 10; % low-pass
[b,a] = butter(3, Fc/(Fs/2), 'low');
filt_dwba = filtfilt(b, a, dwba);

Fc = 0.5; % high-pass
[b,a] = butter(3, Fc/(Fs/2), 'high');
filt_dwba = filtfilt(b, a, filt_dwba);

filt_dwba = filt_dwba - max(filt_dwba);
filt_dwba = 10 * filt_dwba / max(abs(filt_dwba));

abs_filt_wba = abs(filt_dwba); % take absolute value of dWBA

% Find peaks (need to tune these properties too)
min_peak_dist_time = 0.3; % [s]
min_peak_width_time = 0.3; % [s]
[abs_pks,pks_idx] = findpeaks(abs_filt_wba, ...  % https://www.mathworks.com/help/signal/ref/findpeaks.html#namevaluepairarguments
                        'MinPeakProminence', 5, ...
                        'MinPeakHeight', 5, ...
                        'MinPeakDistance', min_peak_dist_time*Fs, ...
                        'MaxPeakWidth', min_peak_width_time*Fs);

filt_dwba_pks = filt_dwba(pks_idx);
dwba_pks = dwba(pks_idx);
time_pks = time(pks_idx);

% Plot signals, red markers are saccade peaks
fig = figure (1);
set(fig, 'Color', 'w')
clear ax
ax(1) = subplot(3,1,1) ; cla ; hold on ; title('Raw')
    ylabel('\DeltaWBA (V)')
    plot([0 time(end)], [0 0], '--', 'Color', [0.5 0.5 0.5])
    plot(time, dwba, 'k')
    plot(time_pks, dwba_pks, '.r', 'MarkerSize', 15)
ax(2) = subplot(3,1,2) ; cla ; hold on ; title('Filtered')
    ylabel('\DeltaWBA (V)')
    plot([0 time(end)], [0 0], '--', 'Color', [0.5 0.5 0.5])
    plot(time, filt_dwba, 'b', 'LineWidth', 1)
    plot(time_pks, filt_dwba_pks, '.r', 'MarkerSize', 15)
ax(3) = subplot(3,1,3) ; cla ; hold on ; title('Absolute Value')
    ylabel('\DeltaWBA (V)')
    plot(time, abs_filt_wba, 'g', 'LineWidth', 1)
    plot(time_pks, abs_pks, '.r', 'MarkerSize', 15)
    xlabel('Time (s)')

set(ax, 'LineWidth', 1)
linkaxes(ax,'x')

% Extract saccades & inter-saccade intervals, align to peak time
thresh = 0;
amp_cut = 0;
direction = -1;
sacd_length = 0.5;
showplot = true;
obj = saccade(filt_dwba, time, thresh, amp_cut, direction, pks_idx, sacd_length, showplot);
