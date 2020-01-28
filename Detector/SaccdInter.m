function [Saccade,Interval,Stimulus,Error,IntError,matchFlag] = saccdinter(xx,tt,SACD,match,STIM,debug)
%% SaccdInter: extracts saccade & inter-saccade intervals
%   INPUTS:
%       xx          :   position
%       tt          :   time vector
%       SACD      	:   saccade table
%       match      	:   only analyze saccades of a specified direction
%                         0 = all-directions
%                      	  1 = co-directional
%                        -1 = anti-directional
%                         2 = positive
%                        -2 = negative
%       Stim      	:   stimulus data (same size as "xx")
%       debug       :   showplot boolean
%
%   OUTPUTS:
%       Saccade     :   saccade intervals
%       Interval   	:   inter-saccade intervals
%       Stimulus   	:   stimulus intervals
%       Error       :   error intervals
%       IntError  	:   integrated error intervals
%       matchFlag  	:   true if there are no saccades with this condition
%

% Check if there are any saccades in the data
if isnan(SACD.Duration(1))
    SACD(1,:) = [];
    warning('No saccades detected')
end

if nargin>3
  	matchFlag = false;
    % Remove saccades in specified direction
    if match~=0
        if (match==1) || (match==-1)
            [dir] = find(SACD.Match~=match);
            if ~isempty(dir)
                SACD(dir,:) = [];
            end

            if isempty(SACD)
                warning('No saccades fit the match condition')
                matchFlag = true;
            end
        elseif (match==2) || (match==-2)
          	[dir] = find(SACD.Direction~=sign(match));
            if ~isempty(dir)
                SACD(dir,:) = [];
            end

            if isempty(SACD)
                warning('No saccades fit the match condition')
                matchFlag = true;
            end
        else
            error('match input incorrect')
        end
    end
end

xx = xx(:);
tt = tt(:);
Ts  = mean(diff(tt));           % sampling time
vel(:,1) = [0 ; diff(xx)/Ts];  	% derivative of data

if ~isempty(STIM) % we know the stimulus
    nn = 4;
    vel(:,2) = [0 ; diff(STIM)/Ts];	% derivative of stimulus
else
    nn = 2;
end

n = size(SACD,1);               % # of saccades
Inter = cell(n,nn);          	% store saccde & inter-sacade intervals
Time = cell(n,2);            	% store saccde & inter-sacade time intervals
PeakTime = cell(n,1);           % peak velocity time for each saccade
StartTime = cell(n,1);          % start times for each inter-saccade interval
EndTime = cell(n,1);            % start times for each inter-saccade interval
for kk = 1:n % all saccades
    % Saccade Intervals
    SI = SACD.StartIdx(kk);             % saccade interval start index
    EI = SACD.EndIdx(kk);               % saccade interval end index
    Inter{kk,1}(:,1) = xx(SI:EI);       % saccade interval position
 	Inter{kk,1}(:,2) = vel(SI:EI);      % saccade interval velocity
    Time{kk,1} = tt(SI:EI);             % saccade interval times
    PeakTime{kk} = SACD.PeakTime(kk);   % peak velocity time
    StartTime{kk} = SACD.StartTime(kk);	% start times for each inter-saccade interval (based on 
    EndTime{kk} = SACD.EndTime(kk);     % start times for each inter-saccade interval
    if ~isempty(STIM) % we know the stimulus
        Inter{kk,3}(:,1) = STIM(SI:EI);
        Inter{kk,3}(:,2) = vel(SI:EI,2);
    end
    
    % Inter-Saccde Intervals
    if kk==1 % 1st saccade does not have valid inter-saccade interval
    	Time{kk,2}(:,1) = tt(1:SI);
      	Inter{kk,2}(:,1:2) = nan(length(Time{kk,2}(:,1)),2);
        if ~isempty(STIM) % we know the stimulus
            Inter{kk,4}(:,1:2) = nan(length(Time{kk,2}(:,1)),2);
        end
    else
        SI = SACD.EndIdx(kk-1);             % inter-saccade interval start index
        EI = SACD.StartIdx(kk);             % inter-saccade interval end index
        Inter{kk,2}(:,1) = xx(SI:EI);       % inter-saccade interval position
    	Inter{kk,2}(:,2) = vel(SI:EI,1);    % inter-saccade interval velocity
    	Time{kk,2}(:,1) = tt(SI:EI);        % inter-saccade interval time
        if ~isempty(STIM) % we know the stimulus
            Inter{kk,4}(:,1) = STIM(SI:EI);
            Inter{kk,4}(:,2) = vel(SI:EI,2);
        end        
    end
end

TimeNorm = cellfun(@(x) x - x(1), Time, 'UniformOutput', false); % normalized times
TimePeak = cellfun(@(x,y) x - y, Time(:,1), PeakTime, 'UniformOutput', false); % aligned to peak times of saccades
% TimeEnd = cellfun(@(x,y) x - y, Time(:,2), StartTime, ... % aligned to end time of inter-saccde intervals (saccde start)
%     'UniformOutput', false);

[Saccade.Time,~,~,~,dR,~] = nancat_center(TimePeak,0,1);

maxInter = cellfun(@(x) length(x), Inter(:,2), 'UniformOutput', true);
[Interval.Time,~,~,~,~,~] = nancat_center(TimeNorm(:,2),[],1,1);

Saccade.Pos    = nan(size(Saccade.Time));
Saccade.Vel    = nan(size(Saccade.Time));
Interval.Pos   = nan(size(Interval.Time));
Interval.Vel   = nan(size(Interval.Time));

Stimulus.Saccade.Pos 	= nan(size(Saccade.Time));
Stimulus.Saccade.Vel 	= nan(size(Saccade.Time));
Stimulus.Interval.Pos 	= nan(size(Interval.Time));
Stimulus.Interval.Vel 	= nan(size(Interval.Time));

for kk = 1:n
    Saccade.Pos(:,kk)      = cat_pad(Inter{kk,1}(:,1), dR{kk},nan);
	Saccade.Vel(:,kk)      = cat_pad(Inter{kk,1}(:,2), dR{kk},nan);
    
    interOff = [0,max(maxInter)] - [0,maxInter(kk)];
    
    Interval.Pos(:,kk)     = cat_pad(Inter{kk,2}(:,1), interOff, nan);
	Interval.Vel(:,kk)     = cat_pad(Inter{kk,2}(:,2), interOff, nan);
    
    Interval.Pos(:,kk) = Interval.Pos(:,kk) - Interval.Pos(1,kk);
    
    if ~isempty(STIM) % we know the stimulus
        Stimulus.Saccade.Pos(:,kk)      = cat_pad(Inter{kk,3}(:,1), dR{kk},nan);
        Stimulus.Saccade.Vel(:,kk)      = cat_pad(Inter{kk,3}(:,2), dR{kk},nan); 
        Stimulus.Interval.Pos(:,kk)     = cat_pad(Inter{kk,4}(:,1), interOff, nan);
        Stimulus.Interval.Vel(:,kk)   	= cat_pad(Inter{kk,4}(:,2), interOff, nan);
        
        Stim_Sacd_Pos_Norm = Stimulus.Saccade.Pos(~isnan(Stimulus.Saccade.Pos(:,kk)),kk);
        Stim_Sacd_Pos_Norm = Stim_Sacd_Pos_Norm(1);
        Stim_Inter_Pos_Norm = Stimulus.Interval.Pos(~isnan(Stimulus.Interval.Pos(:,kk)),kk);
        if all(isnan(Stim_Inter_Pos_Norm))
            Stim_Inter_Pos_Norm = nan;
        else
            Stim_Inter_Pos_Norm = Stim_Inter_Pos_Norm(1);
        end
        
        Stimulus.Saccade.Pos(:,kk)     = Stimulus.Saccade.Pos(:,kk) - Stim_Sacd_Pos_Norm;
        Stimulus.Interval.Pos(:,kk)    = Stimulus.Interval.Pos(:,kk) - Stim_Inter_Pos_Norm;
    end
    
end

Error.Saccade.Pos = Stimulus.Saccade.Pos - Saccade.Pos;
Error.Saccade.Vel = Stimulus.Saccade.Vel - Saccade.Vel;
Error.Interval.Pos = Stimulus.Interval.Pos - Interval.Pos;
Error.Interval.Vel = Stimulus.Interval.Vel - Interval.Vel;

if n==0
    IntError.Saccade.Pos(:,kk) = nan(1,1);
    IntError.Saccade.Vel(:,kk) = nan(1,1);

    IntError.Interval.Pos = nan(1,1);
    IntError.Interval.Vel = nan(1,1);
else
    for kk = 1:n
        err_saccd_time = Saccade.Time(:,kk);
        saccd_nanIdx = ~isnan(err_saccd_time);
        err_saccd_time = err_saccd_time(saccd_nanIdx);

        err_saccd_pos = Error.Saccade.Pos(:,kk);   
        err_saccd_pos = err_saccd_pos(saccd_nanIdx);

        err_saccd_vel = Error.Saccade.Vel(:,kk);   
        err_saccd_vel = err_saccd_vel(saccd_nanIdx);

        err_inter_time = Interval.Time(:,kk);
        inter_nanIdx = ~isnan(err_inter_time);
        err_inter_time = err_inter_time(inter_nanIdx);

        err_inter_pos = Error.Interval.Pos(:,kk);   
        err_inter_pos = err_inter_pos(inter_nanIdx);

        err_inter_vel= Error.Interval.Vel(:,kk);   
        err_inter_vel = err_inter_vel(inter_nanIdx);

        IntError.Saccade.Pos(:,kk) = nan(length(saccd_nanIdx),1);
        IntError.Saccade.Pos(saccd_nanIdx,kk) = cumtrapz(err_saccd_pos,err_saccd_time);
        IntError.Saccade.Vel(:,kk) = nan(length(saccd_nanIdx),1);
        IntError.Saccade.Vel(saccd_nanIdx,kk) = cumtrapz(err_saccd_vel,err_saccd_time);

        IntError.Interval.Pos(:,kk) = nan(length(inter_nanIdx),1);
        IntError.Interval.Pos(inter_nanIdx,kk) = cumtrapz(err_inter_pos,err_inter_time);
        IntError.Interval.Vel(:,kk)  = nan(length(inter_nanIdx),1);
        IntError.Interval.Vel(inter_nanIdx,kk) = cumtrapz(err_inter_vel,err_inter_time);
    end
end

Sacd_Pos_Mean = nanmean(Saccade.Pos,2);
Sacd_Vel_Mean = nanmean(Saccade.Vel,2);
Sacd_Time_Mean = nanmean(Saccade.Time,2);
Inter_Pos_Mean = nanmean(Interval.Pos,2);
Inter_Vel_Mean = nanmean(Interval.Vel,2);
Inter_Time_Mean = nanmean(Interval.Time,2);

if matchFlag
    return
end

if debug
    FIG = figure ; clf % data with saccade & inter-saccde intervals highlighted
    FIG.Color = 'w';
    FIG.Units = 'inches';
    FIG.Position = [2 2 6 4];
        clear ax
        ax(1) = subplot(2,1,1) ; hold on
        ax(1).XLabel.String = 'Time';
        ax(1).YLabel.String = 'Data';
        
        plot(tt,xx,'--k','LineWidth',0.5)
        for kk = 1:n
            h.sacd  = plot(Time{kk,1},Inter{kk,1}(:,1),'-','LineWidth',2);
            h.inter = plot(Time{kk,2},Inter{kk,2}(:,1),'-k','LineWidth',0.5);
        end
%         legend([h.sacd h.inter],'Saccade','Interval')

        ax(2) = subplot(2,1,2) ; hold on
        ax(2).XLabel.String = 'Time';
        ax(2).YLabel.String = 'Diff Data';
        plot(tt,vel(:,1),'--k','LineWidth',0.5)
       
        for kk = 1:n
            h.sacd  = plot(Time{kk,1},Inter{kk,1}(:,2),'-','LineWidth',2);
            h.inter = plot(Time{kk,2},Inter{kk,2}(:,2),'-k','LineWidth',0.5);
        end
        
        if ~isempty(STIM)
            plot(tt,vel(:,2),'--','Color','c','LineWidth',1)
        end
        
        cc = 'k';
        set(ax,'Color','w','YColor',cc,'XColor',cc)
    
	FIG = figure ; clf % normalized saccade & inter-saccde intervals
    FIG.Color = 'w';
    FIG.Units = 'inches';
    FIG.Position = [2+6.5 2 6 4];
        clear ax
        ax(1) = subplot(2,2,1) ; hold on
        ax(1).YLabel.String = 'Saccade';
     	h.sacd  = plot(Sacd_Time_Mean,Saccade.Pos,'-','LineWidth',1);
        plot(Sacd_Time_Mean,Sacd_Pos_Mean,'-k','LineWidth',2);
        plot([0 0],ax(1).YLim,'--k','LineWidth',1.5);
                
        ax(2) = subplot(2,2,2) ; hold on
        ax(2).YLabel.String = 'Interval-Sacadde Interval';
       	h.inter = plot(Interval.Time,Interval.Pos,'-','LineWidth',1);
        plot(Inter_Time_Mean,Inter_Pos_Mean,'-k','LineWidth',2);
        
        if ~isempty(STIM)
            plot(Interval.Time,Stimulus.Interval.Pos,'--r','LineWidth',2)
        end
        
        ax(3) = subplot(2,2,3) ; hold on
        ax(3).XLabel.String = 'Normalized Time';
        ax(3).YLabel.String = 'Saccade Diff';
       	h.sacd  = plot(Saccade.Time,Saccade.Vel,'-','LineWidth',1);
        plot(Sacd_Time_Mean,Sacd_Vel_Mean,'k','LineWidth',2);
        plot([0 0],ax(3).YLim,'--k','LineWidth',1.5);

        ax(4) = subplot(2,2,4) ; hold on
        ax(4).XLabel.String = 'Normalized Time';
        ax(4).YLabel.String = 'Interval-Sacadde Interval Diff';
       	h.inter = plot(Interval.Time,Interval.Vel,'-','LineWidth',1);
       	plot(Inter_Time_Mean,Inter_Vel_Mean,'-k','LineWidth',2);
        
        if ~isempty(STIM)
            plot(Interval.Time,Stimulus.Interval.Vel,'--r','LineWidth',2)
        end

        cc = 'k';
        set(ax,'Color','w','YColor',cc,'XColor',cc)
end

end