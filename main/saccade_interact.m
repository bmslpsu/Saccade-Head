classdef saccade_interact
    % saccade_interact: class to investigate interacting saccades in seperate signals
    %   
 	
    properties (SetAccess = private, Hidden = false)
        % Input signal attributes
        name                % variable inputs
        in                  % input saccade object
        out                 % output saccade object
        win_signal          % window around saccade peaks in time units
        win_sync           	% window around saccade peaks for correlation in time units
        time                % full time vector for input & output
     	n                   % # of data points
        Ts                  % sampling time
        Fs                  % sampling frequency
        
        % Cross-correlation
        acor                % normalized cross-correlation
        timelags            % time lags
        zrcorr              % cross-correlation at 0 lag
        maxcor              % max cross-correlation in window
        timediff            % time difference from max cross-correlation       
        
        % Intervals around saccades
        interval            %
        interval_rel        %
        interval_rel_sync   %
        
        % Synchronous dynamics
        TimeDiff            %
        PeakVel             %
        Amplitude           %
        Duration            %
        
    end
    
  	properties (SetAccess = public, Hidden = false)
    	showplot        % show plots if true
        extra         	% anything else to store
    end
    
 	properties (SetAccess = private, Hidden = true)
    end
     
	properties (Transient) 
    end
    
    methods
        function obj = saccade_interact(in, out, win_signal, win_sync, showplot)
            % saccade_interact: Construct an instance of this class
            %
            %   INPUTS:
            %       in              : input saccade object
            %       out            	: output saccade object
            %       win_signal    	: window around saccade peaks in time units
            %       win_cc          : window around saccade peaks for correlation in time units
            %       showplot        : show plots if true
            %
            
            % Get input/output variable names
            obj.name.in = inputname(1);
            obj.name.out = inputname(2);
            obj.in = in;
            obj.out = out;
            
            % Assign properties
            obj.win_signal  = win_signal;
            obj.win_sync  	= win_sync;
            
            % Make sure signals are aligned
            assert(all(in.time == out.time), 'saccade object must be synchronized to have same times')
            obj.time        = in.time;
           	obj.n           = in.n;
            obj.Ts          = in.Ts;
            obj.Fs          = in.Fs;
            
            % Compute cross-correlation
            maxlag = obj.time(end);
            obj = crosscorr(obj, maxlag);
            
            % Extraxt raw intervals
            [obj.interval.in.pos, obj.interval.in.vel, ~, obj.interval.time_align] = ...
                        pull_interval(obj, obj.in, obj.win_signal);
            [obj.interval.out.pos, obj.interval.out.vel] = ...
                        pull_interval(obj, obj.out, obj.win_signal);

            % Extraxt intervals relative to each other
            [obj.interval_rel.out.pos, obj.interval_rel.out.vel, ~, obj.interval_rel.out.time_align] = ...
                        pull_interval_relative(obj, obj.in, obj.out, obj.win_signal);
            [obj.interval_rel.in.pos, obj.interval_rel.in.vel, ~, obj.interval_rel.in.time_align] = ...
                        pull_interval_relative(obj, obj.out, obj.in, obj.win_signal);
                    
            % Seperate synchrononous and desynchronous saccades
            dirflag = false;
            %[obj] = pull_sync(obj, obj.in, obj.out, obj.win_sync, obj.win_sync, dirflag);
            
        end
        
        function obj = crosscorr(obj, maxlag)
            % crosscorr: cross-correlation between input & output
            %
            %   INPUTS:
            %       maxlag	: max lag window in time
            %
                        
            % Compute cross-corr and set window for max value
            x = obj.in.position;
            y = obj.out.position;
            
            [obj.acor,lags] = xcorr(x, y, 'normalized'); % full cross-correlation
            obj.timelags = lags' / obj.Fs; % time lags
            maxlag = round(maxlag * obj.Fs); % window in samples to find max cross-corr
            cent = ceil(length(lags) / 2); % center index
            lag_win = cent - maxlag : cent + maxlag; % max lag window

            % Set all outside window to nan
            %timelags_win = nan(size(obj.timelags));
            %timelags_win(lag_win) = obj.timelags(lag_win);
            acor_lag_win = nan(size(obj.timelags));
            acor_lag_win(lag_win) = obj.acor(lag_win);
            
            % Find zero & max cross-corr
            obj.zrcorr = obj.acor(lags==0);
            [obj.maxcor,idx] = max(acor_lag_win);
            obj.timediff = obj.timelags(idx);
        end
        
        function [pos,vel,time,time_align] = pull_interval(obj, saccade_obj, win)
            % pull_interval: extract interval around saccade peaks and align to peak
            %
            %   INPUTS:
            %       saccade_obj  	: saccade object
            %       win            	: windows around saccade peaks (+/-)
            %
            
            span        = ceil(win * obj.Fs);
            int_sz      = 2*span + 1;
            pos         = nan(int_sz,saccade_obj.count);
            vel         = nan(int_sz,saccade_obj.count);
            time        = nan(int_sz,saccade_obj.count);
            time_align  = nan(int_sz,1);
            
            % Get intervals
            for s = 1:saccade_obj.count
                pkI = saccade_obj.SACD.PeakIdx(s);
                I = (pkI - span) : (pkI + span);
                
                in_range = (I >= 1) & (I <= obj.n);
                I = I(in_range);
                
                if all(in_range)
                    pos(in_range,s) 	= saccade_obj.position(I) - mean(saccade_obj.position(I));
                    %pos(in_range,s) 	= saccade_obj.position(I) - saccade_obj.SACD.PeakPos(s);
                    vel(in_range,s)   	= saccade_obj.velocity(I);
                    time(in_range,s)  	= saccade_obj.time(I);

                    time_align(in_range,1)	= saccade_obj.time(I) - saccade_obj.SACD.PeakTime(s);
                end
            end
        end
        
       function [pos_out,vel_out,time,time_align] = pull_interval_relative(obj, in, out, win)
            % pull_interval_relative: extract interval around saccade peaks and align to 
            %                         peak for input, align output realitve to input
            %
            %   INPUTS:
            %       in          : input saccade object
           	%       out       	: out saccade object
            %       win       	: windows around saccade peaks (+/-)
            %
            
            span        = ceil(win * obj.Fs);
            int_sz      = 2*span + 1;
            pos_out   	= nan(int_sz,in.count);
            vel_out   	= nan(int_sz,in.count);
            time        = nan(int_sz,in.count);
            time_align  = nan(int_sz,1);
                        
            [b,a] = butter(3, 10 / (obj.Fs/2), 'low');
            vel_filt = filtfilt(b, a, out.velocity);
            vel_out_filt = nan(int_sz,in.count);
            
            % Get output intervals aligned to input
            for s = 1:in.count
                pkI = in.SACD.PeakIdx(s);
                I = (pkI - span) : (pkI + span);
                
                in_range = (I >= 1) & (I <= obj.n);
                I = I(in_range);
                if all(in_range)
                    pos_out(in_range,s) 	= out.position(I) - mean(out.position(I));
                    vel_out(in_range,s)   	= out.velocity(I);
                    time(in_range,s)        = out.time(I);

                    time_align(in_range,s)	= out.time(I) - in.SACD.PeakTime(s);

                    vel_out_filt(in_range,s) = vel_filt(I);
                end
            end
            
%             % Find peaks to make sure there are no false detetcions
%             locs = nan(in.count,1);
%             peaks = nan(in.count,1);
%             for s = 1:in.count
%                 [~,loc] = findpeaks(out.direction*vel_out_filt(:,s), ...
%                                         'MinPeakProminence', 0.7*out.min_pkprom, ...
%                                         'MinPeakHeight', out.threshold, ...
%                                         'MinPeakDistance', round(out.min_pkdist*obj.Fs), ...
%                                         'MinPeakWidth', 0.7*round(out.min_pkwidth*obj.Fs), ...
%                                         'SortStr','descend', 'NPeaks', 1);
%                 if ~isempty(loc)
%                     locs(s)= loc;
%                     peaks(s) = vel_out_filt(locs(s),s);
%                 end
%             end
            
%             % Remove false detections
%             reverse_cond = isnan(peaks) | (out.direction*peaks < 100);
%             if any(reverse_cond)
%                 warning('Possible false detection')
%             end
%             pos_out(:,reverse_cond) = nan;
%             vel_out(:,reverse_cond) = nan;
%             time(:,reverse_cond) = nan;
%             time_align(:,reverse_cond) = nan;
       end
        
       function [obj] = pull_sync(obj, in, out, win_sync, win, dirflag)
            % pull_sync: extract interval around saccade peaks and align to 
            %                         peak for input, align output realitve to input
            %
            %   INPUTS:
            %       in          : input saccade object
           	%       out       	: out saccade object
            %       win       	: windows around saccade peaks (+/-)
            %
            
            %dirflag = false;
            
            span_sync  	= ceil(win_sync * obj.Fs);
            span        = ceil(win * obj.Fs);
            int_sz      = 2*span + 1;
            pos_out   	= nan(int_sz,in.count);
            vel_out   	= nan(int_sz,in.count);
            time        = nan(int_sz,in.count);
            time_align  = nan(int_sz,1);
            
            % Get output intervals aligned to input
            in_peak = in.SACD.PeakIdx;
            out_peak = out.SACD.PeakIdx;
            in_peak_sync = nan(in.count,1);
            out_peak_sync = nan(in.count,1);
            obj.TimeDiff = nan(in.count,3);
            obj.PeakVel = nan(in.count,3);
            obj.Amplitude = nan(in.count,3);
            obj.Duration = nan(in.count,3);
            for s = 1:in.count
                I_sync = (in_peak(s) - span_sync) : (in_peak(s) + span_sync); % search window            
                w = find( (out_peak <= I_sync(end)) & (out_peak  >= I_sync(1)) ); % get synchronous peak if it exists
                
                if length(w) > 1 % if multiple wing saccades in window around head saccade ==> use closest one
                    [~,close_peak] = min( abs(in_peak(s) - w));
                    w(s) = w(close_peak);
                end
                
                % Make sure head and wing saccade are the same direction
                if dirflag
                    if out.SACD.Direction(w) ~= in.SACD.Direction(s)
                        w = [];
                    end
                end
                
                if ~isempty(w)
                  	out_peak_sync(s) = out_peak(w);
                    in_peak_sync(s)  = in_peak(s);
                    
                    % Get dynamics
                  	obj.TimeDiff(s,1) = in.SACD.StartTime(s) - out.SACD.StartTime(w);
                    obj.TimeDiff(s,2) = in.SACD.PeakTime(s) - out.SACD.PeakTime(w);
                    obj.TimeDiff(s,3) = in.SACD.EndTime(s) - out.SACD.EndTime(w);
                    obj.PeakVel(s,1) = in.SACD.PeakVel(s);
                    obj.PeakVel(s,2) = out.SACD.PeakVel(w);
                    obj.Amplitude(s,1) = in.SACD.Amplitude(s);
                    obj.Amplitude(s,2) = out.SACD.Amplitude(w);
                    obj.Duration(s,1) = in.SACD.Duration(s);
                    obj.Duration(s,2) = out.SACD.Duration(w);
                    
                    % Get Intervals
                   	I = (in_peak(s) - span) : (in_peak(s) + span);
                    in_range = (I >= 1) & (I <= obj.n);
                    I = I(in_range);

                    pos_out(in_range,s) 	= out.position(I) - mean(out.position(I));
                    vel_out(in_range,s)   	= out.velocity(I);
                    time(in_range,s)        = out.time(I);

                    time_align(in_range,s)	= out.time(I) - in.SACD.PeakTime(s);
                else
                    % pass
                end

            end
            
            % Dynamics ratios
            obj.PeakVel(:,3) = obj.PeakVel(:,1) ./ obj.PeakVel(:,2);
            obj.Amplitude(:,3) = obj.Amplitude(:,1) ./ obj.Amplitude(:,2);
            obj.Duration(:,3) = obj.Duration(:,1) ./ obj.Duration(:,2);
       end
    	
    end
end

