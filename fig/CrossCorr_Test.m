
clear ; close all ; clc

% Make smooth signals
Fs = 1000;
T = 10;
t = (0:(1/Fs):T)';
n = length(t);
f = 2;
phi_1 = deg2rad(0);
phi_2 = deg2rad(45);
phase_diff = phi_1 - phi_2;
time_diff_smooth = -(1/f) * (phase_diff/(2*pi));
x = 1*sin(2*pi*f*t + phi_1);
y = 1*sin(2*pi*f*t + phi_2);

% Add "saccades"
sacd_length = 0.5;
sacd_start = [2];
sacd_end = sacd_start + sacd_length;
td = -30e-3;
td_idx = Fs * td;
offset = 4;

for k = 1:length(sacd_start)
    idx_start = round(Fs * sacd_start(k));
    idx_end = round(Fs * sacd_end(k));
    x(idx_start:idx_end) = x(idx_start:idx_end) + offset;
    y(td_idx + (idx_start:idx_end) ) = y(td_idx + (idx_start:idx_end) ) + offset;
end

% Compute cross-corr and set window for max value
[acor,lags] = xcorr(x, y, 'normalized'); % full cross-corr
timelags = lags' / Fs; % time lags
maxlag_time = 0.1; % window time to find max cross-corr [s]
maxlag = round(maxlag_time * Fs); % window in samples to find max cross-corr
cent = ceil(length(lags) / 2); % center index
lag_win = cent-maxlag:cent+maxlag; % max lag window

% Set all outside window to nan
timelags_win = nan(size(timelags));
acor_lag_win = nan(size(timelags));
timelags_win(lag_win) = timelags(lag_win);
acor_lag_win(lag_win) = acor(lag_win);

% Find max cross-corr & time difference
[maxcor,idx] = max(acor_lag_win);
time_diff_cross = timelags(idx);

% Plot
fig = figure (1);
set(fig, 'Color', 'w')
xcolor = [0 0 0.7];
ycolor = [0.7 0 0];
ax(1) = subplot(2,1,1) ; cla ; hold on
        plot(t, x, 'Color', xcolor, 'LineWidth', 1)
        plot(t, y, 'Color', ycolor, 'LineWidth', 1)
        xlabel('Time (s)')

ax(2) = subplot(2,1,2) ; cla ; hold on ; ylim([-1 1])
        title(['TD_{cc} = ' num2str(1000*time_diff_cross) , '      TD_{smooth} = ' num2str(1000*time_diff_smooth)])
        plot(timelags, acor, 'k', 'LineWidth', 1)
        plot([time_diff_cross time_diff_cross],[ax(2).YLim(1) maxcor], 'm', 'LineWidth', 2)
        xlim(4*maxlag_time*[-1 1])
        
        xlabel('Time Lag (s)')
    
    
    
    
    