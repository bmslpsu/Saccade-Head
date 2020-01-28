classdef saccade
    % saccade: class to detect and analyze saccade in a data set
    %   
    
	properties (SetAccess = public, Hidden = false)

    end
       
    properties (SetAccess = private, Hidden = false)
     	time                % raw time
        position            % raw position
    	velocity            % raw velocity
        index               % indicies
        ismanual            % boolean: if saccades were detetcted manually
        direction           % direction of saccades to detect (- = -1, all = 0, + = 1)
        count               % # of saccades deteteced
        rate                % saccades / time unit
        SACD                % table fo saccade statistics
        saccades            % saccade only index/time/position/velocity tables in cells
        intervals           % inter-saccade intervals
        shift               % data with interpolation between saccades
        normpeak_saccade    % normalized saccades, aligned to peak time
     	normstart_interval 	% normalized intervals, aligned to start time
     	normend_interval 	% normalized intervals, aligned to end time

        saccades_all        % all saccade data
     	removed_all         % data with saccades removed & nan's in their place
    	threshold           % saccade detcteion threshold (nan if manual)
      	nstd_thresh         % # of std's used for detection threshold
        boundThresh         % ratio of peak velocity to bound saccades on both sides
        Ts                  % sampling time
        Fs                  % sampling frequency
        n                   % # of data points
    end
    
 	properties (SetAccess = private, Hidden = true)
        abs_velocity
    	median
        std
      	peaks
        starts
        ends
        amplitudes
        durations
        directions
        cmap        % colormap
    end
    
	properties (Transient)
        
    end
    
    
    methods
        function obj = saccade(position,time,threshold,direction,debug)
            % saccade: Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.position            = position(:);                      % position vector
            obj.time                = time(:);                          % time vector
            obj.direction           = sign(direction);                	% saccade direction to detect (-1,0,1)
         	obj.n                   = length(obj.time);                 % # of data points
            obj.index            	= (1:obj.n)';                     	% index vector
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
                obj.nstd_thresh = nan;
            end
            
            obj = detectSaccades(obj);
            obj = removeSaccades(obj);
            
            if debug
                plotSaccade(obj)
            end
            
            obj = normSaccade(obj);
            
        end
        
        function obj = calculateThreshold(obj,N)
            % calculateThreshold: compute detetcion threshold
            %   N   : # of stds for threshold
            
            obj.nstd_thresh = N;
            obj.threshold = obj.median.abs_velocity + obj.nstd_thresh*obj.std.abs_velocity;
        end
        
        function obj = detectSaccades(obj,peaks)
            % detetcSaccades: detetc saccades in data
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
             	obj.ismanual = true;
                obj.peaks.index = peaks;
                obj.threshold = nan;
            end
            
            % Only use saccades in the direction specified
            if obj.direction==0
                % analyze all saccades
            else % remove saccades that don't match the input direction
                dir_include = sign(obj.velocity(obj.peaks.index)) == obj.direction;
                obj.peaks.index = obj.peaks.index(dir_include);
            end
            
            obj.peaks.time      = obj.time(obj.peaks.index);
            obj.peaks.position  = obj.position(obj.peaks.index);
            obj.peaks.velocity  = obj.velocity(obj.peaks.index);
            obj.count           = length(obj.peaks.index);
            obj.rate            = obj.count / (obj.time(end) - obj.time(1));
            obj.cmap            = hsv(obj.count);
            
            % Names for saccade statistics table
            tablenames = {'Duration'  ,  'Amplitude'  , 'Direction'  , ...
                          'StartIdx'  ,  'PeakIdx'    , 'EndIdx'     ,...
                       	  'StartTime' ,  'PeakTime'   , 'EndTime'    ,...
                          'StartPos'  ,  'PeakPos'    , 'EndPos'     , ...
                          'StartVel'  ,  'PeakVel'    , 'EndVel'};
            obj.boundThresh = 0.15; % saccades start & end at boundThresh*peak velocity
            if ~isempty(obj.peaks.index) % if any saccades are detected                
                for ww = 1:obj.count % every saccade
                    % Find start & end of saccade
                    obj.starts.index(ww,1) = find(obj.abs_velocity(1:obj.peaks.index(ww)) <= ... % saccade start index
                        abs(obj.peaks.velocity(ww))*obj.boundThresh,1,'last');
                    
                    Eind = find(obj.abs_velocity <= abs(obj.peaks.velocity(ww))*obj.boundThresh); % all values below 1/4 peak
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
      	% removeSaccades: extract saccade kinematics from data & grab
      	% intervals between saccades
     	% 
            if ~isempty(obj.count) % if any saccades are detected
                obj.saccades  = cell(obj.count,1); % pull out each saccade
                obj.intervals = cell(obj.count,1); % pull out each interval
                for ww = 1:obj.count % every saccade
                    % Create saccade table
                    n_saccade_idx = obj.SACD.EndIdx(ww) - obj.SACD.StartIdx(ww) + 1;
                    obj.saccades{ww} = splitvars(table(nan(n_saccade_idx,4)));
                    obj.saccades{ww}.Properties.VariableNames = {'Index','Time','Position','Velocity'};
                    
                    % Assign variables
                    obj.saccades{ww}.Index      = (obj.SACD.StartIdx(ww):obj.SACD.EndIdx(ww))';
                 	obj.saccades{ww}.Time   	= obj.time(obj.saccades{ww}.Index);
                    obj.saccades{ww}.Position 	= obj.position(obj.saccades{ww}.Index);
                    obj.saccades{ww}.Velocity 	= obj.velocity(obj.saccades{ww}.Index);
                    
                    % Create interval table
                    if ww~=1
                        n_interval_idx = obj.SACD.StartIdx(ww) - obj.SACD.EndIdx(ww-1) + 1;
                        obj.intervals{ww} = splitvars(table(nan(n_interval_idx,4)));
                        obj.intervals{ww}.Properties.VariableNames = {'Index','Time','Position','Velocity'};
                        
                        obj.intervals{ww}.Index     = (obj.SACD.EndIdx(ww-1):obj.SACD.StartIdx(ww))';
                     	obj.intervals{ww}.Time      = obj.time(obj.intervals{ww}.Index);
                        obj.intervals{ww}.Position  = obj.position(obj.intervals{ww}.Index);
                        obj.intervals{ww}.Velocity  = obj.velocity(obj.intervals{ww}.Index);
                    else % don't include 1st interval
                        n_interval_idx = obj.SACD.StartIdx(ww);
                        obj.intervals{ww} = splitvars(table(nan(n_interval_idx,4)));
                     	obj.intervals{ww}.Properties.VariableNames = {'Index','Time','Position','Velocity'};
                        
                        obj.intervals{ww}.Index     = (1:obj.SACD.StartIdx(ww))';
                        obj.intervals{ww}.Time      = nan*obj.time(obj.intervals{ww}.Index);
                        obj.intervals{ww}.Position  = nan*obj.position(obj.intervals{ww}.Index);
                        obj.intervals{ww}.Velocity  = nan*obj.velocity(obj.intervals{ww}.Index);
                    end
                end
                
                % Get only saccade data;
                obj.saccades_all = cat(1,obj.saccades{:});
                saccade_idx = ismember(obj.index,obj.saccades_all.Index)';
                
                % Remove saccades & replace with nan's
                obj.removed_all = table(obj.index,obj.time,obj.position,obj.velocity,...
                    'VariableNames',{'Index','Time','Position','Velocity'});
                
                obj.removed_all.Index(saccade_idx)      = nan;
                obj.removed_all.Time(saccade_idx)       = nan;
                obj.removed_all.Position(saccade_idx)   = nan;
                obj.removed_all.Velocity(saccade_idx)   = nan;
                
                % Setup table to store shifted data with removed saccades
                obj.shift = obj.removed_all;
                obj.shift.Index = obj.index;
                obj.shift.Time = obj.time;
                
                % Shift position between saccades so data after a saccade
                % starts where the last saccades ends
                for ww = 1:obj.count % every saccade
                    obj.shift.Position(obj.SACD.EndIdx(ww):end) = obj.shift.Position(obj.SACD.EndIdx(ww):end) ...
                        - ( obj.SACD.EndPos(ww) - obj.SACD.StartPos(ww) );
                end
                
                % Remove nan's for interpolation
                time_shift = obj.removed_all.Time(~isnan(obj.removed_all.Time));
                pos_shift  = obj.shift.Position(~isnan(obj.removed_all.Time));
                vel_shift  = obj.shift.Velocity(~isnan(obj.removed_all.Time));
                
                % Interpolate between saccades
                intrp_pos = table( interp1(time_shift, pos_shift, obj.shift.Time, 'pchip'),... % interpolate
                    'VariableNames',{'IntrpPosition'} );
                intrp_vel= table( interp1(time_shift, vel_shift, obj.shift.Time, 'pchip'),...  % interpolate
                    'VariableNames',{'IntrpVelocity'} );
                
                % intrp_vel = table( [0;diff(intrp_pos{:,1})./diff(obj.time)], ...
                %    'VariableNames',{'IntrpVelocity'} );
                
                % Add to shift table
                obj.shift = [obj.shift,intrp_pos,intrp_vel];

            end
        end
        
        function [obj] = normSaccade(obj)
        % plotSaccade: plots data & detected saccades
        %  
%         obj.normpeak_saccade  	= cell(obj.count,1); % saccades aligned to peak time
%         obj.normstart_interval	= cell(obj.count,1); % intervals aligned to start time
%         obj.normend_interval    = cell(obj.count,1); % intervals aligned to end time
        
%         obj.normpeak_saccade	= obj.saccades;  % saccades aligned to peak time
%         obj.normstart_interval	= obj.intervals; % intervals aligned to start time
%         obj.normend_interval	= obj.intervals; % intervals aligned to end time
        
        obj.normpeak_saccade.time = cellfun(@(x,y) x.Time - y, obj.saccades, num2cell(obj.SACD.PeakTime), ...
            'UniformOutput', false); % aligned to peak times of saccades
        
        [obj.normpeak_saccade.time ,~,~,~,align_shift,~] = nancat_center(obj.normpeak_saccade.time ,0,1);
        
        
        
        
        end
        
        function plotSaccade(obj)
            % plotSaccade: plots data & detected saccades
            %   
            
            FIG = figure ; clf
            FIG.Color = 'w';

            ax(1) = subplot(4,1,1) ; hold on
                ylabel('Position')
                h(1) = plot(obj.time,obj.position,'k');
                plot(obj.time, zeros(obj.n,1),'--','Color',[0.5 0.5 0.5])
                for ww = 1:obj.count
                   plot(obj.saccades{ww}.Time,obj.saccades{ww}.Position,...
                       'LineWidth', 1, 'Color', obj.cmap(ww,:))
                   plot(obj.intervals{ww}.Time,obj.intervals{ww}.Position,...
                       'LineWidth', 1, 'Color', 0.5*obj.cmap(ww,:))
                end
              	plot(obj.starts.time , obj.starts.position , '*g')
              	plot(obj.peaks.time  , obj.peaks.position  , '*b')
                plot(obj.ends.time   , obj.ends.position   , '*r')
                
                ax(1).YLim = max(abs(ax(1).YLim))*[-1 1];

            ax(2) = subplot(4,1,2) ; hold on
                ylabel('Velocity')
                h(2) = plot(obj.time,obj.velocity,'k');                
                for ww = 1:obj.count
                   plot(obj.saccades{ww}.Time,obj.saccades{ww}.Velocity,...
                       'LineWidth', 1, 'Color', obj.cmap(ww,:))
                   plot(obj.intervals{ww}.Time,obj.intervals{ww}.Velocity,...
                       'LineWidth', 1, 'Color', 0.5*obj.cmap(ww,:))
                end
                plot(obj.starts.time , obj.starts.velocity  , '*g')
              	plot(obj.peaks.time  , obj.peaks.velocity   , '*b')
                plot(obj.ends.time   , obj.ends.velocity    , '*r')
                plot(obj.time,  obj.threshold*ones(obj.n,1) , '--m')
              	plot(obj.time, -obj.threshold*ones(obj.n,1) , '--m')
                
                ax(2).YLim = max(abs(ax(2).YLim))*[-1 1];

            ax(3) = subplot(4,1,3) ; hold on
                ylabel('Removed Position')                
                plot(obj.shift.Time, obj.shift.IntrpPosition,      	 'r', 'LineWidth', 1);
                plot(obj.shift.Time, obj.shift.Position,             'k', 'LineWidth', 1);
                plot(obj.removed_all.Time, obj.removed_all.Position, 'c', 'LineWidth', 1);
                
            ax(4) = subplot(4,1,4) ; hold on
                ylabel('Removed Velocity')
                xlabel('Time')
                plot(obj.shift.Time, obj.shift.IntrpVelocity,   'r', 'LineWidth', 0.5);
                plot(obj.shift.Time, obj.shift.Velocity,        'k', 'LineWidth', 0.5);
                
            set(ax,'LineWidth',1)
            linkaxes(ax,'x')
            set(h,'LineWidth',0.5)
            % align_Ylabels(FIG)
        end
    end
end

