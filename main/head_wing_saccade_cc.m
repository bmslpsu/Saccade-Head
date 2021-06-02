function [head2wing] = head_wing_saccade_cc(head_saccade, wing_saccade, win, win_head, win_wing, align_wing, showplot)
%% head_wing_saccade_cc: compare head and wing saccades and compute cross-correlation
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
span = round(win * Fs); % window in samples around each head saccade for detection
span_head = round(win_head * Fs); % window in samples around each head saccade
span_wing = round(win_wing * Fs); % window in samples around each wing saccade
head = head_saccade.position; % head signal

wing = wing_saccade.position; % wing signal
head_color = [0 0 1]; % head color
wing_color = [1 0 0]; % wing color

% Compute cross-corr and set window for max value
[acor,lags] = xcorr(head, wing_saccade.extra.dwba, 'normalized'); % full cross-corr
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
sync = head_saccade.SACD.PeakIdx; % sync to peak of head saccade
tint = nan(2*span_head+1,head_saccade.count); % head intervals
hint = nan(2*span_head+1,head_saccade.count); % head intervals
wint = nan(2*span_head+1,head_saccade.count); % wing intervals
hint_vel = nan(2*span_head+1,head_saccade.count); % head intervals
wint_vel = nan(2*span_head+1,head_saccade.count); % wing intervals
for s = 1:head_saccade.count
    interval = (sync(s) - span):(sync(s) + span); % interval around head saccade
    interval_keep = (sync(s) - span_head):(sync(s) + span_head); % interval around head saccade to keep
    
    if (max(interval_keep) > length(head)) || (min(interval_keep) < 1)
        % skip because interval is not complete
    else
        temp_head = head(interval); % head interval
        temp_wing = wing_saccade.extra.dwba(interval); % wing interval
        
        hint(:,s) = head(interval_keep); % head interval
        wint(:,s) = wing_saccade.extra.dwba(interval_keep); % wing interval
        hint_vel(:,s) = head_saccade.velocity(interval_keep); % head interval
        wint_vel(:,s) = wing_saccade.extra.dwba_vel(interval_keep); % wing interval
        tint(:,s) = time(interval_keep) - time(sync(s)); % time interval
        
        [acor_int,lags_int] = xcorr(temp_head, temp_wing, 'normalized'); % cross-corr within interval
        timelags_int = lags_int / Fs; % time lags within interval

        % Store cross-corr in matric columns
        int_acor(:,s) = acor_int;
        int_lags(:,s) = timelags_int;
    end
end

% Find mean cross-corr across all saccades and find max
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

% Head/wing dynamics and time differenence
wpeak_all = wing_saccade.SACD.PeakIdx;
sync_head = head_saccade.SACD.PeakIdx;
sync_wing = wing_saccade.SACD.PeakIdx;
TD = nan(head_saccade.count,3);
Amp = nan(head_saccade.count,2);
Dur = nan(head_saccade.count,2);
PkVel = nan(head_saccade.count,2);

% Synchronized head/wing saccades
hsacd.time = nan(2*span_head+1,head_saccade.count);
hsacd.pos  = nan(2*span_head+1,head_saccade.count);
hsacd.vel  = nan(2*span_head+1,head_saccade.count);
wsacd.time = nan(2*span_wing+1,head_saccade.count);
wsacd.pos  = nan(2*span_wing+1,head_saccade.count);
wsacd.vel  = nan(2*span_wing+1,head_saccade.count);

% Desynchronized head/wing saccades
hsacd.pos_desync  = nan(2*span_head+1,head_saccade.count);
hsacd.vel_desync  = nan(2*span_head+1,head_saccade.count);
wsacd.pos_desync  = nan(2*span_wing+1,head_saccade.count);
wsacd.vel_desync  = nan(2*span_wing+1,head_saccade.count);

% Head synchronization rates
n_sync_head = 0;

% Find time difference based on start, peak, & end times when head & wing saccades occur within window
for s = 1:head_saccade.count
    hpeak = head_saccade.SACD.PeakIdx(s);
    wpeak = find( (wpeak_all <= (hpeak + span)) & (wpeak_all >= (hpeak - span)) );

    if length(wpeak) > 1 % if multiple wing saccades in window around head saccade ==> use closest one
        [~,close_peak] = min( abs(hpeak - wpeak));
        wpeak = wpeak(close_peak);
    end
    
    % Make sure head and wing saccade are the same direction
    if wing_saccade.SACD.Direction(wpeak) ~= head_saccade.SACD.Direction(s)
        wpeak = [];
    end
    
    % Pull out time difference and dynamics
    if length(wpeak) == 1 % wing saccade found in window around head saccade
        n_sync_head = n_sync_head + 1;
        TD(s,1) = head_saccade.SACD.StartTime(s) - wing_saccade.SACD.StartTime(wpeak);
        TD(s,2) = head_saccade.SACD.PeakTime(s) - wing_saccade.SACD.PeakTime(wpeak);
        TD(s,3) = head_saccade.SACD.EndTime(s) - wing_saccade.SACD.EndTime(wpeak);
        Amp(s,1) = head_saccade.SACD.Amplitude(s);
        Amp(s,2) = wing_saccade.SACD.Amplitude(wpeak);
      	Dur(s,1) = head_saccade.SACD.Duration(s);
        Dur(s,2) = wing_saccade.SACD.Duration(wpeak);
      	PkVel(s,1) = head_saccade.SACD.PeakVel(s);
        PkVel(s,2) = wing_saccade.SACD.PeakVel(wpeak);
        
        shift_idx = 0*round(TD(s,2) * Fs);
        
        % Get intervals around head & wing saccades relative to peak time
      	interval_head = (sync_head(s) - span_head):(sync_head(s) + span_head); % interval around head saccade
        
        if ~align_wing % wing saccades around head saccade, not aligned
            interval_wing = (sync_head(s) - span_wing):(sync_head(s) + span_wing); % interval around wing saccade
        else % wing saccades aligned to peak and shifted by peak time difference
            interval_wing = shift_idx + ((sync_wing(wpeak) - span_wing):(sync_wing(wpeak) + span_wing)); % interval around wing saccade
        end
        
        out_range = any( (interval_wing < 1) | (interval_wing > head_saccade.n) ) || ...
                    any( (interval_head < 1) | (interval_head > head_saccade.n) );
        
        if ~out_range
            % Sync to head peak time and wing offset
            hsacd.time(:,s) = head_saccade.time(interval_head) - head_saccade.time(sync_head(s));
            hsacd.pos(:,s)  = head_saccade.position(interval_head);
            hsacd.vel(:,s)  = head_saccade.velocity(interval_head);
            
            if ~align_wing % wing saccades around head saccade, not aligned
                wsacd.time(:,s) = wing_saccade.time(interval_wing) - head_saccade.time(sync_head(s));
            else % wing saccades aligned to peak and shifted by peak time difference
                wsacd.time(:,s) = wing_saccade.time(interval_wing) - wing_saccade.time(sync_wing(wpeak)) ...
                                                                - shift_idx/Fs;
            end

            %wsacd.pos(:,s)	= wing_saccade.position(interval_wing);
            %wsacd.vel(:,s)	= wing_saccade.velocity(interval_wing);
            wsacd.pos(:,s)	= wing_saccade.extra.dwba(interval_wing);
            wsacd.vel(:,s)  = wing_saccade.extra.dwba_vel(interval_wing);
           	%wsacd.pos(:,s)	= wing_saccade.position_filt_detect(interval_wing);
            %wsacd.vel(:,s)  = wing_saccade.velocity_filt_detect(interval_wing);
        end
	elseif isempty(wpeak) % no wing saccade found in window
        interval_all = (sync_head(s) - span_head):(sync_head(s) + span_head); % interval around head saccade
        out_range =  any( (interval_all < 1) | (interval_all > head_saccade.n) );
        if ~out_range
            hsacd.time(:,s)         = head_saccade.time(interval_all) - head_saccade.time(sync_head(s));
            hsacd.pos_desync(:,s)   = head_saccade.position(interval_all);
            hsacd.vel_desync(:,s)   = head_saccade.velocity(interval_all);
            %wsacd.pos_desync(:,s)  = wing_saccade.position(interval_all);
            %wsacd.vel_desync(:,s)  = wing_saccade.velocity(interval_all);
            wsacd.pos_desync(:,s)	= wing_saccade.extra.dwba(interval_all);
            wsacd.vel_desync(:,s)	= wing_saccade.extra.dwba_vel(interval_all);
        end
    end
end
hsacd.stats = structfun(@(x) basic_stats(x,2), hsacd, 'UniformOutput', false);
wsacd.stats = structfun(@(x) basic_stats(x,2), wsacd, 'UniformOutput', false);
sync_head_rate = n_sync_head / head_saccade.count;
sync_wing_rate = (wing_saccade.count - n_sync_head) / wing_saccade.count;

% Assign output properties
head2wing.hint          = hint;
head2wing.wint          = wint;
head2wing.hint_vel      = hint_vel;
head2wing.wint_vel      = wint_vel;
head2wing.tint          = tint;

head2wing.acor          = acor;
head2wing.timelags      = timelags;

head2wing.maxcor     	= maxcor;
head2wing.timediff      = timediff;

head2wing.int_acor      = int_acor;
head2wing.int_lags      = int_lags;

head2wing.int_maxcor    = int_maxcor;
head2wing.int_timediff  = int_timediff;

head2wing.TimeDiff  	= 1000*TD; % [ms]
head2wing.Amplitude   	= Amp;
head2wing.Duration     	= Dur;
head2wing.PeakVel     	= PkVel;
head2wing.hsacd         = hsacd;
head2wing.wsacd         = wsacd;

head2wing.sync_head_rate = sync_head_rate;
head2wing.sync_wing_rate = sync_wing_rate;

if showplot && all(~isnan(sync))
    fig = figure (10);
    set(fig, 'Color', 'w')
    ax = gobjects(6,1);
        ax(1) = subplot(4,2,1:2) ; cla ; hold on
            yyaxis left ; hold on ; cla ; ylabel('Head (°)')
                plot(time, head, 'Color', [head_color 0.5], 'LineWidth', 1)
                for ww = 1:head_saccade.count
                   plot(head_saccade.saccades{ww}.Time, head_saccade.saccades{ww}.Position,...
                        '-','LineWidth', 2, 'Color', head_color)
                end
%              	plot(time(sync), head(sync), '.b', ...
%                     'LineWidth', 1, 'MarkerSize', 20, 'MarkerEdgeColor', head_color, 'MarkerFaceColor', 'none')
                ylim(25*[-1 1])
                %ylim(max(abs(ax(1).YLim))*[-1 1])
            yyaxis right ; hold on ; cla ; ylabel('\DeltaWBA (°)')
                plot(time, wing_saccade.extra.dwba, '-', 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
                plot(time, wing, 'Color', [wing_color 0.5], 'LineWidth', 1)
                for ww = 1:wing_saccade.count
                   plot(wing_saccade.saccades{ww}.Time, wing_saccade.saccades{ww}.Position,...
                        '-','LineWidth', 2, 'Color', wing_color)
                end
%                 plot(wing_saccade.SACD.PeakTime, wing_saccade.SACD.PeakPos, '.r', ...
%                     'LineWidth', 1,  'MarkerSize', 20, 'MarkerEdgeColor', wing_color, 'MarkerFaceColor', 'none')
                ylim(25*[-1 1])    
                %ylim(max(abs(ax(1).YLim))*[-1 1])
                xlabel('Time (s)')
            ax(1).YAxis(1).Color = head_color;
            ax(1).YAxis(2).Color = wing_color;
            %ax(1).YLim = ceil(max(abs(ax(1).YLim)))*[-1 1];

            for s = 1:head_saccade.count
                interval = (sync(s) - span):(sync(s) + span);
                if (max(interval) > length(head)) || (min(interval) < 1)
                    % skip
                else
                    tint = time(interval);
                    xx = [tint(1), tint(1), tint(end), tint(end)];
                    pp = [ax(1).YLim(1), ax(1).YLim(2), ax(1).YLim(2), ax(1).YLim(1)];
                    if isnan(TD(s,1))
                        facecolor = [0.5 0.5 0.5];
                    else
                        facecolor = 'g';
                    end
                    patch(xx, pp, facecolor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                end
            end

        % Synchronized head/wing saccades
        ax(2) = subplot(4,2,3) ; cla ; hold on
            yyaxis left ; hold on ; cla ; ylabel('Head (°)')
                plot(hsacd.time, hsacd.pos, '-','Color', [0.7*head_color 0.5])
                plot(hsacd.stats.time.mean, hsacd.stats.pos.mean, ...
                     '-', 'Color', head_color, 'LineWidth', 2)
              	ylim(25*[-1 1])

         	yyaxis right ; hold on ; cla ; ylabel('\DeltaWBA (°)')
                plot(wsacd.time, wsacd.pos, '-','Color', [0.7*wing_color 0.5])
                plot(wsacd.stats.time.mean, wsacd.stats.pos.mean, ...
                     '-', 'Color', wing_color, 'LineWidth', 2)
             	ylim(25*[-1 1])
          	
          	ax(2).YAxis(1).Color = head_color;
            ax(2).YAxis(2).Color = wing_color;
            
        ax(3) = subplot(4,2,4) ; cla ; hold on
            yyaxis left ; hold on ; cla ; ylabel('Head (°/s)')
                plot(hsacd.time, hsacd.vel, '-','Color', [0.7*head_color 0.5])
                plot(hsacd.stats.time.mean, hsacd.stats.vel.mean, ...
                     '-', 'Color', head_color, 'LineWidth', 2)

         	yyaxis right ; hold on ; cla ; ylabel('\DeltaWBA (°/s)')
                plot(wsacd.time, wsacd.vel, '-','Color', [0.7*wing_color 0.5])
                plot(wsacd.stats.time.mean, wsacd.stats.vel.mean, ...
                     '-', 'Color', wing_color, 'LineWidth', 2)
          	
          	ax(3).YAxis(1).Color = head_color;
            ax(3).YAxis(2).Color = wing_color;
            
       	% Desynchronized head/wing saccades
        ax(4) = subplot(4,2,5) ; cla ; hold on
            yyaxis left ; hold on ; cla ; ylabel('Head (°)')
                plot(hsacd.time, hsacd.pos_desync, '-','Color', [0.7*head_color 0.5])
                plot(hsacd.stats.time.mean, hsacd.stats.pos_desync.mean, ...
                     '-', 'Color', head_color, 'LineWidth', 2)

         	yyaxis right ; hold on ; cla ; ylabel('\DeltaWBA (°)')
                plot(hsacd.time, wsacd.pos_desync, '-','Color', [0.7*wing_color 0.5])
                plot(hsacd.stats.time.mean, wsacd.stats.pos_desync.mean, ...
                     '-', 'Color', wing_color, 'LineWidth', 2)
          	
          	ax(4).YAxis(1).Color = head_color;
            ax(4).YAxis(2).Color = wing_color;
            
        ax(5) = subplot(4,2,6) ; cla ; hold on
            yyaxis left ; hold on ; cla ; ylabel('Head (°/s)')
                plot(hsacd.time, hsacd.vel_desync, '-','Color', [0.7*head_color 0.5])
                plot(hsacd.stats.time.mean, hsacd.stats.vel_desync.mean, ...
                     '-', 'Color', head_color, 'LineWidth', 2)

         	yyaxis right ; hold on ; cla ; ylabel('\DeltaWBA (°/s)')
                plot(hsacd.time, wsacd.vel_desync, '-','Color', [0.7*wing_color 0.5])
                plot(hsacd.stats.time.mean, wsacd.stats.vel_desync.mean, ...
                     '-', 'Color', wing_color, 'LineWidth', 2)
          	
          	ax(5).YAxis(1).Color = head_color;
            ax(5).YAxis(2).Color = wing_color;
        
        ax(6) = subplot(4,2,7:8) ; cla ; hold on
            title(['TD_{all} = ' num2str(1000*timediff) , '        TD_{scd} = ' num2str(1000*int_timediff)])

            plot(1000*int_lags, int_acor, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 1)
            plot(1000*int_lags, int_acor_stats.mean,  'g', 'LineWidth', 3)
            plot(1000*[int_timediff int_timediff],[ax(6).YLim(1) int_maxcor], 'g', 'LineWidth', 2)

            plot(1000*timelags, acor, 'k', 'LineWidth', 2)
            plot(1000*timelags_win, acor_lag_win, 'm', 'LineWidth', 3)
            plot(1000*[timediff timediff],[ax(6).YLim(1) maxcor], 'm', 'LineWidth', 2)

            xlim(1000*4*maxlag_time*[-1 1])
            xlabel('Time Lag (ms)')

        set(ax, 'LineWidth', 1.5)
        linkaxes(ax(2:5), 'x')
        
        linkaxes(ax([2,4]),'xy')
        linkaxes(ax([3,5]),'xy')
end

end
