classdef saccade
    % saccade: Summary of this class goes here
    %   Detailed explanation goes here
    
	properties (SetAccess = public, Hidden = false)

    end
       
    properties (SetAccess = private, Hidden = false)
     	time
        position
    	velocity
        ismanual
        SACD
        removed

    	threshold
      	nstd_thresh
    	count
        rate
        Ts
        Fs
        n           % # of data points
        
     	saccade_position
        saccade_velocity
        saccade_time
    end
    
 	properties (SetAccess = private, Hidden = false)
        allthresh
        abs_velocity
    	median
        std
        
      	peaks
        starts
        ends
        amplitudes
        durations
        directions
        
    end
    
	properties (Transient)
        
    end
    
    
    methods
        function obj = saccade(position,time,threshold,debug)
            % saccade: Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.position            = position(:);                      % position vector
            obj.time                = time(:);                          % time vector
         	obj.n                   = length(obj.time);                 % # of data points
            obj.Ts                  = mean(diff(obj.time));             % sampling time
            obj.Fs                  = 1/obj.Ts;                         % sampling frequency
            obj.velocity            = [0 ; diff(obj.position)/obj.Ts]; 	% velocity
            obj.abs_velocity      	= abs(obj.velocity);                % absolute velocity

            obj.median.abs_velocity	= median(obj.abs_velocity);	% median absolute velocity
            obj.std.abs_velocity 	= std(obj.abs_velocity);   	% STD absolute velocity
            
            if threshold <= 6 % calculate threshold based on # of standard deviations
                obj = calculateThreshold(obj,threshold);
            else % specify threshold manually
                obj.threshold = threshold; % detection threshold
            end
            
            obj = detectSaccades(obj);
            
            if debug
                plotSaccade(obj)
            end
            
        end
        
        function obj = calculateThreshold(obj,N)
            % calculateThreshold: compute detetcion threshold
            %   N   : # of stds for threshold
            obj.threshold = obj.median.abs_velocity + N*obj.std.abs_velocity;
        end
        
        function obj = detectSaccades(obj,peaks)
            % detetcSaccades: compute detection threshold
            %   peaks   : saccade peak locations [indicies]
            
            if nargin < 2 % if peaks are not given
                obj.ismanual = false;
                
                % Limit peak detetction to absolute velocity above threshold
                svel = obj.abs_velocity;
                svel(obj.abs_velocity<obj.threshold) = 0;

                % Find peaks in absolute of velocity
                [~, locs] = findpeaks(svel,'MINPEAKDISTANCE',40); % find local maxima
                window = obj.n*0.015;  % window length at start & end to ignore saccades [samples]
                I = (locs > window) & (locs < (obj.n - (window+1))); % saccades inside middle window
                obj.peaks.index = locs(I);
                
            elseif nargin == 2 % user specifies peaks
                obj.peaks.index = peaks;
                obj.threshold = nan;
                obj.ismanual = true;
            end
            
            obj.peaks.time      = obj.time(obj.peaks.index);
            obj.peaks.position  = obj.position(obj.peaks.index);
            obj.peaks.velocity  = obj.velocity(obj.peaks.index);
            obj.count           = length(obj.peaks.index);
            obj.rate            = obj.count / (obj.time(end) - obj.time(1));
            
            % Names for saccade statistics table
            tablenames = {'Duration'  ,  'Amplitude'  , 'Direction'  , ...
                          'StartIdx'  ,  'PeakIdx'    , 'EndIdx'     ,...
                       	  'StartTime' ,  'PeakTime'   , 'EndTime'    ,...
                          'StartPos'  ,  'PeakPos'    , 'EndPos'     , ...
                          'StartVel'  ,  'PeakVel'    , 'EndVel'};
           
            if ~isempty(obj.peaks.index) % if any saccades are detected
                obj.allthresh = repmat(obj.threshold,obj.count,1); % same threshold used for each saccades
                boundThresh = 1/4; % saccades start & end at 25% of peak velocity
                for ww = 1:obj.count % every saccade
                    % Find start & end of saccade
                    obj.starts.index(ww,1) = find(obj.abs_velocity(1:obj.peaks.index(ww)) <= ... % saccade start index
                        obj.peaks.index(ww)*boundThresh,1,'last');
                    
                    Eind = find(obj.abs_velocity <= abs(obj.peaks.velocity(ww))*boundThresh); % all values below 1/4 peak
                    Es = find(Eind > obj.peaks.index(ww),1,'first'); % first value under 25% of peak after the start index is the end index
                    if ~isempty(Es) % make sure saccade ends
                        obj.ends.index(ww,1) = Eind(Es); % saccade end index
                    end
                end
                
                obj.starts.time   	= obj.time      (obj.starts.index);
                obj.starts.position	= obj.position  (obj.starts.index);
                obj.starts.velocity	= obj.velocity  (obj.starts.index);
                obj.ends.time       = obj.time      (obj.ends.index);
                obj.ends.position   = obj.position  (obj.ends.index);
                obj.ends.velocity   = obj.velocity  (obj.ends.index);

                obj.durations     	= obj.starts.time - obj.ends.time;
                obj.amplitudes     	= obj.starts.position - obj.ends.position;
                obj.directions      = sign(obj.peaks.velocity);

                % Collect data
                STATS = [obj.durations        , obj.amplitudes      , obj.peaks.velocity    , ...
                         obj.starts.index     , obj.peaks.index     , obj.ends.index        , ...
                         obj.starts.time	  , obj.peaks.time      , obj.ends.time         , ...
                         obj.starts.position  , obj.peaks.position  , obj.ends.position     , ...
                         obj.starts.velocity  , obj.peaks.velocity  , obj.ends.velocity     ];
                
            else
                STATS = nan(size(varnames));
                obj.count = nan;
                obj.rate = nan;
            end
                
            % Saccade table
            obj.SACD = splitvars(table(STATS));
            obj.SACD.Properties.VariableNames = tablenames;
            
        end
        
      	function obj = removeSaccades(obj)
      	% detetcSaccades: compute detection threshold
     	%   peaks   : saccade peak locations [indicies]
            
            % obj.removed.time = time.data
            if ~isnan(obj.count) % if any saccades are detected
                
            end
            
        end
        
        function plotSaccade(obj)
            % plotSaccade: plots data & detected saccades
            %   
            
            FIG = figure ; clf
            FIG.Color = 'w';

            ax(1) = subplot(3,1,1) ; hold on
                ylabel('Position')
                h(1) = plot(obj.time,obj.position,'k');
                plot(obj.time, zeros(obj.n,1),'--','Color',[0.5 0.5 0.5])
                plot(obj.starts.time , obj.starts.position ,'*g')
              	plot(obj.peaks.time  , obj.peaks.position  ,'*b')
                plot(obj.ends.time   , obj.ends.position   ,'*r')
                ax(1).YLim = max(abs(ax(1).YLim))*[-1 1];

            ax(2) = subplot(3,1,2) ; hold on
                ylabel('Velocity')
                h(2) = plot(obj.time,obj.velocity,'k');
                plot(obj.time,  obj.threshold*ones(obj.n,1),'--m')
              	plot(obj.time, -obj.threshold*ones(obj.n,1),'--m')

                
                plot(obj.starts.time , obj.starts.velocity ,'*g')
              	plot(obj.peaks.time  , obj.peaks.velocity  ,'*b')
                plot(obj.ends.time   , obj.ends.velocity   ,'*r')
                ax(2).YLim = max(abs(ax(2).YLim))*[-1 1];

            ax(3) = subplot(3,1,3) ; hold on
                ylabel('Removed)')
                xlabel('Time')

            set(ax,'LineWidth',1)
            linkaxes(ax,'x')
            set(h,'LineWidth',0.5)
        end
    end
end

