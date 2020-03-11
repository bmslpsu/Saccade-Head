classdef saccade
    % saccade: class to detect and analyze saccade in a data set
    %   
 	
    properties (SetAccess = private, Hidden = false)
        n                   % # of data points
        Ts                  % sampling time
        Fs                  % sampling frequency
     	index               % indicies
        time                % raw time
        position            % raw position
        velocity            % raw velocity
        abs_velocity        % absolute value of velocity
        ismanual            % boolean: if saccades were detetected manually
        direction           % direction of saccades to detect (- = -1, all = 0, + = 1)
     	nstd_thresh         % # of std's used for detection threshold
        threshold           % saccade detcteion threshold (nan if manual)
        boundThresh         % ratio of peak velocity to bound saccades on both sides
        count               % # of saccades deteteced
        rate                % saccades / time unit
        amp_cutoff = 2;    	% cut off for saccade amplitudes (< amp_cutoff are removed)
        SACD                % table of saccade statistics
        
        saccades            % saccade only index/time/position/velocity tables in cells
        intervals           % inter-saccade intervals
        shift               % data with interpolation between saccades
        normpeak_saccade    % normalized saccades , aligned to peak time
        norm_interval       % normalized intervals, aligned to start time
        normstart_interval	% normalized intervals, aligned to start time and initial position
        normend_interval    % normalized intervals, aligned to end time and end position

        saccades_all        % all saccade data
        removed_all         % data with saccades removed & nan's in their place

        stimlus_position    % stimulus position
        stimlus_velocity    % stimulus veloicty
        stimulus            % stimulus inter-saccade intervals
        normstart_stimulus 	% normalized stimulus intervals, aligned to start time
        normend_stimulus  	% normalized stimulus intervals, aligned to end time
        error               % error (stimulus - position) for inter-saccade intervals
        int_error        	% integrated error (stimulus - position) for inter-saccade intervals
    end
    
  	properties (SetAccess = public, Hidden = false)
        cmap                % colormap
    end
    
 	properties (SetAccess = private, Hidden = true)
    	median
        std
      	peaks
        starts
        ends
        amplitudes
        durations
        directions
        tablenames
        
        align_saccade_peak
        align_interval_start
        align_interval_end
    end
     
	properties (Transient)
        
    end
    
    methods
        function obj = saccade(position, time, threshold, direction, pks, debug)
            % saccade: Construct an instance of this class
            %   Saccade detector
            
            if nargin==0
                obj = setDeafults(obj);
                return
            end
            
            if nargin<6
                debug = true; % default show plots
                if nargin<5
                   pks = []; % default no peaks given
                end
            end
            
            obj.position     	= position(:);                      % position vector
            obj.time          	= time(:);                          % time vector
            obj.direction     	= sign(direction);                	% saccade direction to detect (-1,0,1)
         	obj.n             	= length(obj.time);                 % # of data points
            obj.index       	= (1:obj.n)';                     	% index vector
            obj.Ts           	= mean(diff(obj.time));             % sampling time
            obj.Fs           	= 1/obj.Ts;                         % sampling frequency
            obj.velocity      	= diff(obj.position)/obj.Ts;        % velocity
            obj.velocity      	= [obj.velocity(1) ; obj.velocity];	% velocity
            obj.abs_velocity	= abs(obj.velocity);                % absolute velocity

            obj.median.abs_velocity	= median(obj.abs_velocity);	% median absolute velocity
            obj.std.abs_velocity 	= std(obj.abs_velocity);   	% STD absolute velocity
            
            % Calculate saccade detcetion threshold
            if threshold <= 6 % calculate threshold based on # of standard deviations
                obj = calculateThreshold(obj,threshold);
            else % specify threshold manually
                obj.threshold = threshold; % detection threshold
                obj.nstd_thresh = nan;
            end
            
            % Detect, isolate, & normalize sacades & inter-saccade intervals
            obj = detectSaccades(obj,pks);
            obj = removeSaccades(obj);
            obj = normSaccade(obj);
            
            if debug
                plotSaccade(obj)
                if obj.count~=0
                    plotInterval(obj)
                else
                    warning('No saccades detetected')
                end
            end
            
        end
        
    	function obj = setDeafults(obj)
            % setDeafults: sets values to default
            % 
            
            empty_feilds = ["time","position","velocity"];
            
            for f = 1:length(empty_feilds)
                obj.normpeak_saccade.(empty_feilds(f))      = [];
                obj.norm_interval.(empty_feilds(f))         = [];
                obj.normstart_interval.(empty_feilds(f))   	= [];
                obj.normend_interval.(empty_feilds(f))      = [];
                obj.normstart_stimulus.(empty_feilds(f))  	= [];
                obj.normend_stimulus.(empty_feilds(f))   	= [];
                obj.error.(empty_feilds(f))                 = [];
                obj.int_error.(empty_feilds(f))             = [];
            end
        end
        
        function obj = calculateThreshold(obj,N)
            % calculateThreshold: compute detetcion threshold
            %   N   : # of stds for threshold
            
            obj.nstd_thresh = N;
            obj.threshold = obj.median.abs_velocity + obj.nstd_thresh*obj.std.abs_velocity;
            
            if obj.threshold < 350
                obj.threshold = 350;
            end
        end
        
        function obj = detectSaccades(obj,pks)
            % detetcSaccades: detetc saccades in data
            %   pks   : saccade peak locations [indicies]
            
            if (nargin < 2) || isempty(pks) % if peaks are not given
                obj.ismanual = false;
                
                % Limit peak detetction to absolute velocity above threshold
                svel = obj.abs_velocity;
                svel(obj.abs_velocity<obj.threshold) = 0;

                % Find peaks in absolute of velocity
                [~, locs] = findpeaks(svel,'MINPEAKDISTANCE',15); % find local maxima
                window = obj.n*0.015; % window length at start & end to ignore saccades [samples]
                I = (locs > window) & (locs < (obj.n - (window+1))); % saccades inside middle window
                obj.peaks.index = locs(I);
                
            elseif nargin == 2 % user specifies peaks
             	obj.ismanual = true;
                obj.peaks.index = pks;
                obj.threshold = nan;
            end
            
            % Only use saccades in the direction specified
            if obj.direction==0
                % analyze all saccades
            else % remove saccades that don't match the input direction
                if obj.direction~=0
                    dir_include = sign(obj.velocity(obj.peaks.index)) == obj.direction;
                    obj.peaks.index = obj.peaks.index(dir_include);
                end
            end
            
            obj.peaks.time      = obj.time(obj.peaks.index);
            obj.peaks.position  = obj.position(obj.peaks.index);
            obj.peaks.velocity  = obj.velocity(obj.peaks.index);
            obj.count           = length(obj.peaks.index);
            
            % Names for saccade statistics table
            obj.tablenames = {'Duration'  ,  'Amplitude'  , 'Direction'  , ...
                          'StartIdx'  ,  'PeakIdx'    , 'EndIdx'     ,...
                       	  'StartTime' ,  'PeakTime'   , 'EndTime'    ,...
                          'StartPos'  ,  'PeakPos'    , 'EndPos'     , ...
                          'StartVel'  ,  'PeakVel'    , 'EndVel'};
                      
            obj.boundThresh = 0.25; % saccades start & end at boundThresh*peak velocity
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
                obj.amplitudes     	= obj.ends.position - obj.starts.position;
                
                % Remove saccades below position amplitude threshold
              	rmv_amp = abs(obj.amplitudes) < obj.amp_cutoff;
                if any(rmv_amp)
                    obj.amplitudes = obj.amplitudes(~rmv_amp,:);

                    obj.peaks.index  	= obj.peaks.index(~rmv_amp);
                    obj.peaks.time      = obj.peaks.time(~rmv_amp);
                    obj.peaks.position  = obj.peaks.position(~rmv_amp);
                    obj.peaks.velocity  = obj.peaks.velocity(~rmv_amp);
                    
                    obj.starts.index  	= obj.starts.index(~rmv_amp);
                    obj.starts.time   	= obj.starts.time(~rmv_amp);
                    obj.starts.position = obj.starts.position(~rmv_amp);
                    obj.starts.velocity = obj.starts.velocity(~rmv_amp);

                    obj.ends.index  	= obj.ends.index(~rmv_amp);
                    obj.ends.time       = obj.ends.time(~rmv_amp);
                    obj.ends.position   = obj.ends.position(~rmv_amp);
                    obj.ends.velocity   = obj.ends.velocity(~rmv_amp);

                    obj.count           = length(obj.peaks.index);
                end
                
                obj.durations     	= obj.ends.time - obj.starts.time;
                obj.directions      = sign(obj.peaks.velocity);
             	obj.rate            = obj.count / (obj.time(end) - obj.time(1));
                obj.cmap            = hsv(obj.count);
                
                if obj.count == 0
                    STATS = nan(size(obj.tablenames));
                else
                    % Collect data
                    STATS = [obj.durations        , obj.amplitudes      , obj.directions    , ...
                             obj.starts.index     , obj.peaks.index     , obj.ends.index        , ...
                             obj.starts.time	  , obj.peaks.time      , obj.ends.time         , ...
                             obj.starts.position  , obj.peaks.position  , obj.ends.position     , ...
                             obj.starts.velocity  , obj.peaks.velocity  , obj.ends.velocity     ,];
                end
                     
                % Saccade table
                obj.SACD = splitvars(table(STATS));
                obj.SACD.Properties.VariableNames = obj.tablenames;
                
            else
                STATS = nan(size(obj.tablenames));
                obj.count = 0;
                obj.rate = 0;
                
             	% Saccade table
                obj.SACD = splitvars(table(STATS));
                obj.SACD.Properties.VariableNames = obj.tablenames;
            end
            
        end
        
      	function obj = removeSaccades(obj)
      	% removeSaccades: extract saccade kinematics from data & grab
      	% intervals between saccades
     	% 
            if obj.count~=0 % if any saccades are detected
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
                    else % store 1st interval as nan's because we don't know the start of this interval
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
        % normSaccade: align saccades to peak time & align intervals to
        % start and end time.
        %  
            if obj.count~=0 % if there are any saccades
                % Normalize saccade times to saccade peak times
                obj.normpeak_saccade.time = cellfun(@(x,y) x.Time - y, obj.saccades, ...
                    num2cell(obj.SACD.PeakTime),'UniformOutput', false);

                [obj.normpeak_saccade.time ,~,~,~,obj.align_saccade_peak,~] = ...
                    nancat_center(obj.normpeak_saccade.time, 0, 1);

                % Align saccade positions to saccade peak times
                obj.normpeak_saccade.position = cellfun(@(x,y) padmat(x.Position,y,nan,1), obj.saccades, ...
                    obj.align_saccade_peak, 'UniformOutput', false);
                obj.normpeak_saccade.position = cat(2,obj.normpeak_saccade.position{:});

                % Align saccade velocities to saccade peak times
                obj.normpeak_saccade.velocity = cellfun(@(x,y) padmat(x.Velocity,y,nan,1), obj.saccades, ...
                    obj.align_saccade_peak, 'UniformOutput', false);
                obj.normpeak_saccade.velocity = cat(2,obj.normpeak_saccade.velocity{:});

                % Compute saccade stats
                obj.normpeak_saccade.time_stats = basic_stats(obj.normpeak_saccade.time,2);
                obj.normpeak_saccade.position_stats = basic_stats(obj.normpeak_saccade.position,2);
                obj.normpeak_saccade.velocity_stats = basic_stats(obj.normpeak_saccade.velocity,2);

                %-------------------------------------------

                % Normalize interval times to intervals start times
                obj.normstart_interval.time = cellfun(@(x) x.Time - x.Time(1), obj.intervals, ...
                    'UniformOutput', false);

                [obj.normstart_interval.time ,~,~,~,obj.align_interval_start,~] = ...
                    nancat_center(obj.normstart_interval.time, 0 ,1);
               
                % Align interval start times without normalizing positions
                obj.norm_interval.time = obj.normstart_interval.time;
                
                obj.norm_interval.position = cellfun(@(x,y) padmat(x.Position, y,nan,1), ...
                    obj.intervals, obj.align_interval_start, 'UniformOutput', false);
                obj.norm_interval.position = cat(2,obj.norm_interval.position{:});

                % Align interval positions to interval start times
                obj.normstart_interval.position = cellfun(@(x,y) padmat(x.Position - x.Position(1),y,nan,1), ...
                    obj.intervals, obj.align_interval_start, 'UniformOutput', false);
                obj.normstart_interval.position = cat(2,obj.normstart_interval.position{:});

                % Align interval velocities to interval start times
                obj.normstart_interval.velocity = cellfun(@(x,y) padmat(x.Velocity,y,nan,1), obj.intervals, ...
                    obj.align_interval_start, 'UniformOutput', false);
                obj.normstart_interval.velocity = cat(2,obj.normstart_interval.velocity{:});
                
                obj.norm_interval.velocity = obj.normstart_interval.velocity;

                % Compute interval stats
                obj.normstart_interval.time_stats = basic_stats(obj.normstart_interval.time,2);
                obj.normstart_interval.position_stats = basic_stats(obj.normstart_interval.position,2);
                obj.normstart_interval.velocity_stats = basic_stats(obj.normstart_interval.velocity,2);
                
                obj.norm_interval.position_stats = basic_stats(obj.norm_interval.position,2);

                % Get interval times
                obj.normstart_interval.endidx = sum(~isnan(obj.normstart_interval.time));
                obj.normstart_interval.endidx(1) = obj.normstart_interval.endidx(1) + 1;
                IntTime = table(max(obj.normstart_interval.time(obj.normstart_interval.endidx,:))',...
                    'VariableNames',{'IntTime'});
                obj.SACD = [obj.SACD , IntTime];

                %-------------------------------------------

                % Normalize interval times to intervals end times
                obj.normend_interval.time = cellfun(@(x) x.Time - x.Time(end), obj.intervals, ...
                    'UniformOutput', false);

                [obj.normend_interval.time ,~,~,~,obj.align_interval_end,~] = ...
                    nancat_center(obj.normend_interval.time, 0 ,1);

                % Align interval positions to interval end times
                obj.normend_interval.position = cellfun(@(x,y) padmat(x.Position,y,nan,1), obj.intervals, ...
                    obj.align_interval_end, 'UniformOutput', false);
                obj.normend_interval.position = cat(2,obj.normend_interval.position{:});

                % Align interval velocities to interval end times
                obj.normend_interval.velocity = cellfun(@(x,y) padmat(x.Velocity,y,nan,1), obj.intervals, ...
                    obj.align_interval_end, 'UniformOutput', false);
                obj.normend_interval.velocity = cat(2,obj.normend_interval.velocity{:});

                % Compute interval stats
                obj.normend_interval.time_stats = basic_stats(obj.normstart_interval.time,2);
                obj.normend_interval.position_stats = basic_stats(obj.normstart_interval.position,2);
                obj.normend_interval.velocity_stats = basic_stats(obj.normstart_interval.velocity,2);
            else
             	IntTime = table(nan,'VariableNames',{'IntTime'});
                obj.SACD = [obj.SACD , IntTime];
                obj = setDeafults(obj);
            end
        end
        
        function [obj] = stimSaccade(obj,stim,debug)
        % stimSaccade: computes saccade/interval relationship to visual
        % motion stimulus (stimulus position, stimulus velocity, error,integrated error)
        %  
        
         	% Get stimulus data
            obj.stimlus_position = stim(:);
            obj.stimlus_velocity = diff(obj.stimlus_position)./diff(obj.time);
            obj.stimlus_velocity = [obj.stimlus_velocity(1) ; obj.stimlus_velocity];
            obj.stimulus = cell(obj.count,1); % pull out each saccade
        
            if obj.count~=0 % if there are any saccades
                if nargin < 3
                    debug = false; % default
                end

                % Extract stimulus intervals
                if ~isempty(stim) % if any saccades are detected and a stimulus if given
                    for ww = 1:obj.count % every saccade                   
                        % Create stimulus table
                        if ww~=1
                            n_stim_idx = obj.SACD.StartIdx(ww) - obj.SACD.EndIdx(ww-1) + 1;
                            obj.stimulus{ww} = splitvars(table(nan(n_stim_idx,4)));
                            obj.stimulus{ww}.Properties.VariableNames = {'Index','Time','Position','Velocity'};

                            obj.stimulus{ww}.Index     = (obj.SACD.EndIdx(ww-1):obj.SACD.StartIdx(ww))';
                            obj.stimulus{ww}.Time      = obj.time(obj.stimulus{ww}.Index);
                            obj.stimulus{ww}.Position  = obj.stimlus_position(obj.stimulus{ww}.Index);
                            obj.stimulus{ww}.Velocity  = obj.stimlus_velocity(obj.stimulus{ww}.Index);
                        else % store 1st interval as nan's because we don't know the start of this interval
                            n_stim_idx = obj.SACD.StartIdx(ww);
                            obj.stimulus{ww} = splitvars(table(nan(n_stim_idx,4)));
                            obj.stimulus{ww}.Properties.VariableNames = {'Index','Time','Position','Velocity'};

                            obj.stimulus{ww}.Index     = (1:obj.SACD.StartIdx(ww))';
                            obj.stimulus{ww}.Time      = nan*obj.time(obj.stimulus{ww}.Index);
                            obj.stimulus{ww}.Position  = nan*obj.stimlus_position(obj.stimulus{ww}.Index);
                            obj.stimulus{ww}.Velocity  = nan*obj.stimlus_velocity(obj.stimulus{ww}.Index);
                        end
                    end
                    % Assign times
                    obj.normstart_stimulus.time = obj.normstart_interval.time;
                    obj.normend_stimulus.time   = obj.normend_interval.time;
                    %-------------------------------------------
                    
                    % Align stimulus interval positions to interval start times
                    obj.normstart_stimulus.position = cellfun(@(x,y) padmat(x.Position - x.Position(1),y,nan,1), ...
                        obj.stimulus, obj.align_interval_start, 'UniformOutput', false);
                    obj.normstart_stimulus.position = cat(2,obj.normstart_stimulus.position{:});

                    % Align stimulus interval velocities to interval start times
                    obj.normstart_stimulus.velocity = cellfun(@(x,y) padmat(x.Velocity,y,nan,1), ...
                        obj.stimulus, obj.align_interval_start, 'UniformOutput', false);
                    obj.normstart_stimulus.velocity = cat(2,obj.normstart_stimulus.velocity{:});
                    
                    %-------------------------------------------
                    
                    % Align stimulus interval positions to interval end times
                    obj.normend_stimulus.position = cellfun(@(x,y) padmat(x.Position,y,nan,1), ...
                        obj.stimulus, obj.align_interval_end, 'UniformOutput', false);
                    obj.normend_stimulus.position = cat(2,obj.normend_stimulus.position{:});

                    % Align stimulus interval velocities to interval end times
                    obj.normend_stimulus.velocity = cellfun(@(x,y) padmat(x.Velocity,y,nan,1), ...
                        obj.stimulus, obj.align_interval_end, 'UniformOutput', false);
                    obj.normend_stimulus.velocity = cat(2,obj.normend_stimulus.velocity{:});
                    
                    %-------------------------------------------

                    % Get time vectors for integration & compute error & integrated error within intervals
                    [~,max_time] = max(obj.normstart_interval.endidx);
                    tvector = obj.normstart_interval.time(:,max_time);

                    % Compute error within intervals
                    obj.error.time  	= obj.normstart_interval.time;
                    obj.error.position  = obj.normstart_stimulus.position - obj.normstart_interval.position;
                    obj.error.velocity  = obj.normstart_stimulus.velocity - obj.normstart_interval.velocity;

                    obj.error.position_end = arrayfun(@(x,y)  obj.error.position(y,x), ...
                                1:size(obj.error.position,2), obj.normstart_interval.endidx, ...
                                'UniformOutput',true);
                    obj.error.velocity_end = arrayfun(@(x,y)  obj.error.velocity(y,x), ...
                                1:size(obj.error.velocity,2), obj.normstart_interval.endidx, ...
                                'UniformOutput',true);

                    % Compute integrated error within intervals
                    obj.int_error.time      = obj.normstart_interval.time;
                    obj.int_error.position  = cumtrapz(tvector, obj.error.position, 1);
                    obj.int_error.velocity  = cumtrapz(tvector, obj.error.velocity, 1);

                    obj.int_error.position_end = arrayfun(@(x,y)  obj.int_error.position(y,x), ...
                                1:size(obj.int_error.position,2), obj.normstart_interval.endidx, ...
                                'UniformOutput',true);
                    obj.int_error.velocity_end = arrayfun(@(x,y)  obj.int_error.velocity(y,x), ...
                                1:size(obj.int_error.velocity,2), obj.normstart_interval.endidx, ...
                                'UniformOutput',true);
                    obj.int_error.position_end(1) = nan; % integrated error is nan for 1st saccade
                    obj.int_error.velocity_end(1) = nan; % integrated error is nan for 1st saccade

                    % Make table of finals erros & combine with SACD table
                    ERR = table(obj.error.position_end',     obj.error.velocity_end', ...
                                obj.int_error.position_end', obj.int_error.velocity_end', ...
                                'VariableNames', {'ErrorPos','ErrorVel','IntErrorPos','IntErrorVel'});
                    obj.SACD = [obj.SACD , ERR];

                    % Compute interval stats
                    obj.error.position_stats = basic_stats(obj.error.position,2);
                    obj.error.velocity_stats = basic_stats(obj.error.velocity,2);

                    obj.int_error.position_stats = basic_stats(obj.int_error.position,2);
                    obj.int_error.velocity_stats = basic_stats(obj.int_error.velocity,2);
                end

                % Plot if specified
                if debug
                    if ~obj.count==0
                        plotStimulus(obj) 
                    else
                        warning('No saccades detetected')
                    end
                end
                
            else
                % Saccade table
                obj.SACD = splitvars(table([obj.SACD , splitvars(table(nan(1,4)))]));
                obj.SACD.Properties.VariableNames = [obj.tablenames, ...
                    {'IntTime','ErrorPos','ErrorVel','IntErrorPos','IntErrorVel'}];
            end
        end
        
        function plotSaccade(obj)
            % plotSaccade: plots data & detected saccades
            %   
            
            FIG = figure ; clf
            FIG.Color = 'w';
            FIG.Name = 'Saccades';
            set(FIG,'WindowStyle','docked')

            ax(1) = subplot(4,1,1) ; hold on
                ylabel('Position')
                h(1) = plot(obj.time,obj.position,'k');
                plot(obj.time, zeros(obj.n,1),'--','Color',[0.5 0.5 0.5])
                if obj.count~=0
                    for ww = 1:obj.count
                       plot(obj.saccades{ww}.Time,obj.saccades{ww}.Position,...
                           'LineWidth', 1, 'Color', obj.cmap(ww,:))
                       plot(obj.intervals{ww}.Time,obj.intervals{ww}.Position,...
                           'LineWidth', 1, 'Color', 0.7*obj.cmap(ww,:))
                    end
                    plot(obj.starts.time , obj.starts.position , '*g')
                    plot(obj.peaks.time  , obj.peaks.position  , '*b')
                    plot(obj.ends.time   , obj.ends.position   , '*r')
                end
                
                ax(1).YLim = max(abs(ax(1).YLim))*[-1 1];

            ax(2) = subplot(4,1,2) ; hold on
                ylabel('Velocity')
                h(2) = plot(obj.time,obj.velocity,'k'); 
                plot(obj.time, -obj.threshold*ones(obj.n,1) , '--m')
                plot(obj.time,  obj.threshold*ones(obj.n,1) , '--m')
                if obj.count~=0
                    for ww = 1:obj.count
                       plot(obj.saccades{ww}.Time,obj.saccades{ww}.Velocity,...
                           'LineWidth', 1, 'Color', obj.cmap(ww,:))
                       plot(obj.intervals{ww}.Time,obj.intervals{ww}.Velocity,...
                           'LineWidth', 1, 'Color', 0.7*obj.cmap(ww,:))
                    end
                    plot(obj.starts.time , obj.starts.velocity  , '*g')
                    plot(obj.peaks.time  , obj.peaks.velocity   , '*b')
                    plot(obj.ends.time   , obj.ends.velocity    , '*r')
                end
                ax(2).YLim = max(abs(ax(2).YLim))*[-1 1];

            ax(3) = subplot(4,1,3) ; hold on
                if obj.count~=0
                    ylabel('Removed Position')                
                    plot(obj.shift.Time, obj.shift.IntrpPosition,      	 'r', 'LineWidth', 1);
                    plot(obj.shift.Time, obj.shift.Position,             'k', 'LineWidth', 1);
                    plot(obj.removed_all.Time, obj.removed_all.Position, 'c', 'LineWidth', 1);
                end
            ax(4) = subplot(4,1,4) ; hold on
                if obj.count~=0
                    ylabel('Removed Velocity')
                    xlabel('Time')
                    plot(obj.shift.Time, obj.shift.IntrpVelocity,   'r', 'LineWidth', 0.5);
                    plot(obj.shift.Time, obj.shift.Velocity,        'k', 'LineWidth', 0.5);
                end
            set(ax,'LineWidth',1,'FontWeight','bold')
            linkaxes(ax,'x')
            set(h,'LineWidth',0.5)
            align_Ylabels(FIG)
        end
        
        function plotInterval(obj)
            % plotInterval: plots extracted normalized saccades & intervals
            %   
            
            FIG = figure ; clf
            FIG.Color = 'w';
            FIG.Name = 'Normalized Saccades & Intervals';
            set(FIG,'WindowStyle','docked')
            ax(1) = subplot(2,2,1) ; hold on ; title('Saccades')
                ylabel('Position')
                h.sacdpos = plot(obj.normpeak_saccade.time,obj.normpeak_saccade.position);
                [hstd(1),~] = PlotPatch(obj.normpeak_saccade.position_stats.median, ...
                                obj.normpeak_saccade.position_stats.std, ...
                                obj.normpeak_saccade.time_stats.median,...
                                1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(1),'bottom')
                                
          	ax(2) = subplot(2,2,3) ; hold on
                ylabel('Velocity')
                h.sacdvel = plot(obj.normpeak_saccade.time,obj.normpeak_saccade.velocity);
               	[hstd(2),~] = PlotPatch(obj.normpeak_saccade.velocity_stats.median, ...
                                obj.normpeak_saccade.velocity_stats.std, ...
                                obj.normpeak_saccade.time_stats.median,...
                                1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(2),'bottom')
                                
          	ax(3) = subplot(2,2,2) ; hold on ; title('Intervals')
                h.intpos = plot(obj.normstart_interval.time,obj.norm_interval.position);
                [hstd(3),~] = PlotPatch(obj.norm_interval.position_stats.median, ...
                                obj.norm_interval.position_stats.std, ...
                                obj.normstart_interval.time_stats.median,...
                                1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(3),'bottom')
                            
          	ax(4) = subplot(2,2,4) ; hold on
                h.intvel = plot(obj.normstart_interval.time,obj.normstart_interval.velocity);
                [hstd(4),~] = PlotPatch(obj.normstart_interval.velocity_stats.median, ...
                                obj.normstart_interval.velocity_stats.std, ...
                                obj.normstart_interval.time_stats.median,...
                                1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(4),'bottom')
                            
            set(h.sacdpos, {'color'}, num2cell(obj.cmap,2))
            set(h.sacdvel, {'color'}, num2cell(obj.cmap,2))
            set(h.intpos,  {'color'}, num2cell(obj.cmap,2))
            set(h.intvel,  {'color'}, num2cell(obj.cmap,2))
            
          	set([h.sacdpos h.sacdvel h.intpos h.intvel],'LineWidth',1)
            set(ax,'LineWidth',1,'FontWeight','bold')
            linkaxes(ax(1:2),'x')
            linkaxes(ax(3:4),'x')
            align_Ylabels(FIG)
        end
        
        function plotStimulus(obj)
            % plotStimulus: plots extracted stimulus, error, & integrated error intervals
            %   
            
            FIG = figure ; clf
            FIG.Color = 'w';
            FIG.Name = 'Stimulus & Error';
            set(FIG,'WindowStyle','docked')

            ax(1) = subplot(2,3,1) ; hold on ; title('Stimulus')
                ylabel('Position')
                h.stimpos = plot(obj.normstart_interval.time, obj.normstart_stimulus.position);
                
          	ax(2) = subplot(2,3,4) ; hold on
                ylabel('Velocity')
                h.stimvel = plot(obj.normstart_interval.time, obj.normstart_stimulus.velocity);
                ax(2).YLim = 1.1*abs(max(max(obj.normstart_stimulus.velocity)))*[-1 1];
                
          	ax(3) = subplot(2,3,2) ; hold on ; title('Error')
                h.errpos = plot(obj.normstart_interval.time,obj.error.position);
               	[hstd(1),~] = PlotPatch(obj.error.position_stats.median, ...
                                        obj.error.position_stats.std, ...
                                        obj.normstart_interval.time_stats.median,...
                                        1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(1),'bottom')
                            
          	ax(4) = subplot(2,3,5) ; hold on
                h.errvel = plot(obj.normstart_interval.time,obj.error.velocity);
               	[hstd(2),~] = PlotPatch(obj.error.velocity_stats.median, ...
                                        obj.error.velocity_stats.std, ...
                                        obj.normstart_interval.time_stats.median,...
                                        1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(2),'bottom')
                                    
          	ax(5) = subplot(2,3,3) ; hold on ; title('Integrated Error')
                h.interrpos = plot(obj.normstart_interval.time,obj.int_error.position);
               	[hstd(3),~] = PlotPatch(obj.int_error.position_stats.median, ...
                                        obj.int_error.position_stats.std, ...
                                        obj.normstart_interval.time_stats.median,...
                                        1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(3),'bottom')
                                    
          	ax(6) = subplot(2,3,6) ; hold on
                h.interrvel = plot(obj.normstart_interval.time,obj.int_error.velocity);
               	[hstd(4),~] = PlotPatch(obj.int_error.velocity_stats.median, ...
                                        obj.int_error.velocity_stats.std, ...
                                        obj.normstart_interval.time_stats.median,...
                                        1,1,'k',[0.4 0.4 0.6],0.3,2); uistack(hstd(4),'bottom')
                                    
            set(h.stimpos,      {'color'}, num2cell(obj.cmap,2))
            set(h.stimvel,      {'color'}, num2cell(obj.cmap,2))
            set(h.errpos,       {'color'}, num2cell(obj.cmap,2))
            set(h.errvel,       {'color'}, num2cell(obj.cmap,2))
            set(h.interrpos,    {'color'}, num2cell(obj.cmap,2))
            set(h.interrvel,    {'color'}, num2cell(obj.cmap,2))
            
          	set([h.stimpos h.stimvel h.errpos h.errvel h.interrpos h.interrvel],'LineWidth',1)
            set(ax,'LineWidth',1,'FontWeight','bold')
            linkaxes(ax,'x')
            linkaxes(ax([1,3,6]),'y')
            align_Ylabels(FIG)
        end
        
    end
end

