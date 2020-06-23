function [head2wing] = head_wing_saccade(head_saccade,wing_saccade,win,showplot)
%% head_wing_saccade: find wing saccades in delta wing-beat-amplitude signal
%
%   INPUT:
%       head_saccade	: head saccade object
%       wing_saccade  	: wing saccade object
%       win             : window size (+/-) for cross-correlation aroudn each head saccade [s]
%       showplot        : show debug plot (boolean)
%
%   OUTPUT:
%       head2wing       : structure containing head-wing relationship properties
%

% Get signal data
time = head_saccade.time; % time vector [s]
Fs = head_saccade.Fs; % sampling frequency [Hz]
span = round(win * Fs); % window in samples around each head saccade
head = head_saccade.position; % head signal
wing = wing_saccade.position; % wing signal
head_color = [0 0 0.7]; % head color
wing_color = [0.7 0 0]; % wing color

% Compute cross-corr and set window for max value
[acor,lags] = xcorr(head, wing, 'normalized'); % full cross-corr
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

% Find max cross-corr
[maxcor,idx] = max(acor_lag_win);
timediff = timelags(idx);

% Compute cross-corr around each head saccade
lag_size = 4*span + 1;
int_acor = nan(lag_size,head_saccade.count);
int_lags = nan(lag_size,head_saccade.count);
sync = head_saccade.SACD.PeakIdx; % sync to end of head saccade
for s = 1:head_saccade.count
    interval = (sync(s) - span):(sync(s) + span); % interval around head saccade
    
    if (max(interval) > length(head)) || (min(interval) < 1)
        % skip because interval is not complete
    	int_acor(:,s) = nan;
        int_lags(:,s) = nan;
    else
        hint = head(interval); % head interval
        wint = wing(interval); % wing interval
        [acor_int,lags_int] = xcorr(hint, wint, 'normalized'); % cross-corr within interval
        timelags_int = lags_int / Fs; % time lags within interval

        % Store cross-corr in matric columns
        int_acor(:,s) = acor_int;
        int_lags(:,s) = timelags_int;
    end
end

% Find mean cross-corr acorss all saccades and find max
if isempty(int_lags)
    int_lags = nan;
    int_acor_stats = nan;
    int_maxcor = nan;
    int_timediff = nan;
else
    int_lags = int_lags(:,1);
    int_acor_stats = basic_stats(int_acor, 2);
    [int_maxcor,idx] = max(int_acor_stats.mean);
    int_timediff = int_lags(idx);
end

% Assign output properties
head2wing.acor          = acor;
head2wing.timelags      = timelags;

head2wing.maxcor     	= maxcor;
head2wing.timediff      = timediff;

head2wing.int_acor      = int_acor;
head2wing.int_lags      = int_lags;

head2wing.int_maxcor    = int_maxcor;
head2wing.int_timediff  = int_timediff;

if showplot
    fig = figure (10);
    set(fig, 'Color', 'w')
    ax = gobjects(2,1);
        ax(1) = subplot(2,1,1) ; cla ; hold on
            yyaxis left ; hold on ; cla ; ylabel('Head (°)')
                plot(time, head, 'Color', head_color, 'LineWidth', 1)
                plot(time(sync), head(sync), '.b', ...
                    'LineWidth', 1, 'MarkerSize', 20)
                %ylim(max(abs(ax(1).YLim))*[-1 1])
            yyaxis right ; hold on ; cla ; ylabel('\DeltaWBA (°)')
                plot(time, wing, 'Color', wing_color, 'LineWidth', 1)
                plot(wing_saccade.SACD.PeakTime, wing_saccade.SACD.PeakPos, '.r', ...
                    'LineWidth', 1, 'MarkerSize', 20)
                %ylim(max(abs(ax(1).YLim))*[-1 1])
                xlabel('Time (s)')
            ax(1).YAxis(1).Color = head_color;
            ax(1).YAxis(2).Color = wing_color;

            for s = 1:head_saccade.count
                interval = (sync(s) - span):(sync(s) + span);
                if (max(interval) > length(head)) || (min(interval) < 1)
                    % skip
                else
                    tint = time(interval);
                    xx = [tint(1), tint(1), tint(end), tint(end)];
                    pp = [ax(1).YLim(1), ax(1).YLim(2), ax(1).YLim(2), ax(1).YLim(1)];
                    patch(xx, pp, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                end
            end

        ax(2) = subplot(2,1,2) ; cla ; hold on
            title(['TD_{all} = ' num2str(1000*timediff) , '        TD_{scd} = ' num2str(1000*int_timediff)])

            plot(1000*int_lags, int_acor, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1)
            plot(1000*int_lags, int_acor_stats.mean,  'g', 'LineWidth', 3)
            plot(1000*[int_timediff int_timediff],[ax(2).YLim(1) int_maxcor], 'g', 'LineWidth', 2)

            plot(1000*timelags, acor, 'k', 'LineWidth', 2)
            plot(1000*timelags_win, acor_lag_win, 'm', 'LineWidth', 3)
            plot(1000*[timediff timediff],[ax(2).YLim(1) maxcor], 'm', 'LineWidth', 2)

            xlim(1000*4*maxlag_time*[-1 1])
            xlabel('Time Lag (ms)')

        set(ax, 'LineWidth', 1.5)
end

end
