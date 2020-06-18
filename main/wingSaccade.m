function [pks_idx,time_pks,dwba_pks] = wingSaccade(dwba,time,dir,showplot)
%% wingSaccade: find wing saccades in delta wing-beat-amplitude signal
%
%   INPUT:
%       dwba        : delta wing-beat-amplitude signal
%       time       	: time signal
%       dir         : direction of saccades to detect
%                        1 = CW
%                       -1 = CCW
%                        0 = All
%       showplot    : show debug plot (boolean)
%
%   OUTPUT:
%       peaks       : peak locations (indicies)
%

Fs = 1 / mean(diff(time)); % sampling frequency [Hz]
filt_dwba = dwba - median(dwba);

% Filter (tune these based on data properties)
Fc = 5; % low-pass
[b,a] = butter(3, Fc/(Fs/2), 'low');
filt_dwba = filtfilt(b, a, filt_dwba);

Fc = 0.5; % high-pass
[b,a] = butter(3, Fc/(Fs/2), 'high');
filt_dwba = filtfilt(b, a, filt_dwba);

% Find peaks (need to tune these properties too)
min_peak_dist_time = 0.5; % [s]
min_peak_width_time = 1; % [s]
MinPeakProminence = 0.1;
MinPeakHeight = 0.4;
MinPeakDistance = min_peak_dist_time * Fs;
MaxPeakWidth = min_peak_width_time * Fs;
if dir == 1
    % Normalize
    norm_dwba = filt_dwba / abs(max(filt_dwba));
    %norm_dwba = norm_dwba - max(norm_dwba);
    [norm_pks,pks_idx] = findpeaks(norm_dwba, ...
                            'MinPeakProminence', MinPeakProminence, ...
                            'MinPeakHeight', MinPeakHeight, ...
                            'MinPeakDistance', MinPeakDistance, ...
                            'MaxPeakWidth', MaxPeakWidth);
elseif dir == -1
	% Normalize
    norm_dwba = filt_dwba / abs(min(filt_dwba));
    norm_dwba = -norm_dwba;
    %norm_dwba = norm_dwba - max(norm_dwba);
    [norm_pks,pks_idx] = findpeaks(norm_dwba, ...
                            'MinPeakProminence', MinPeakProminence, ...
                            'MinPeakHeight', MinPeakHeight, ...
                            'MinPeakDistance', MinPeakDistance, ...
                            'MaxPeakWidth', MaxPeakWidth);
	norm_pks = -norm_pks;
elseif dir == 0
    [abs_pks_cw,pks_idx_cw] = findpeaks(abs_filt_wba, ...
                            'MinPeakProminence', MinPeakProminence, ...
                            'MinPeakHeight', MinPeakHeight, ...
                            'MinPeakDistance', MinPeakDistance, ...
                            'MaxPeakWidth', MaxPeakWidth);
    [abs_pks_ccw,pks_idx_ccw] = findpeaks(-abs_filt_wba, ...
                            'MinPeakProminence', MinPeakProminence, ...
                            'MinPeakHeight', MinPeakHeight, ...
                            'MinPeakDistance', MinPeakDistance, ...
                            'MaxPeakWidth', MaxPeakWidth);
                        
     pks_idx = sort([pks_idx_cw , pks_idx_ccw]);
     norm_pks = [abs_pks_cw , -abs_pks_ccw];
end

% Get peaks values in signals
time_pks = time(pks_idx);
dwba_pks = dwba(pks_idx);
norm_dwba_pks = norm_dwba(pks_idx);

% Plot signals, red markers are saccade peaks
if showplot
    fig = figure (5);
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
        plot(time, norm_dwba, 'b', 'LineWidth', 1)
        plot(time_pks, norm_dwba_pks, '.r', 'MarkerSize', 15)
    ax(3) = subplot(3,1,3) ; cla ; hold on ; title('Absolute Value')
        ylabel('\DeltaWBA (V)')
%         plot(time, abs_filt_wba, 'g', 'LineWidth', 1)
%         plot(time_pks, abs_pks, '.r', 'MarkerSize', 15)
        xlabel('Time (s)')

    set(ax, 'LineWidth', 1)
    linkaxes(ax,'x')
end

% % Extract saccades & inter-saccade intervals, align to peak time
% thresh = 0;
% amp_cut = 0;
% direction = -1;
% sacd_length = 0.5;
% showplot = true;
% obj = saccade(filt_dwba, time, thresh, amp_cut, direction, pks_idx, sacd_length, showplot);
end
