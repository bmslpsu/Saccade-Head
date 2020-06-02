function [SACD,thresh,count,rate,SaccdRmv_X,vel,svel] = sacd_detect_manual(xx,tt,debug)
%% sacd_detect_manual.: user defines saccades
%   INPUTS:
%       xx          :   position
%       tt          :   time vector
%       debug       :   showplot boolean
%   OUTPUTS:
%       SACD        :   saccade data table
%       count       :   # of saccades
%       rate        :   saccades/s
%       SaccdRmv_X  :   positon data without saccades
%

n = length(xx); % # of data points
IDX = (1:n)'; % index vector
if nargin<3
	debug = false;
    if nargin<2
        tt = IDX; % use index if no time data
    end
end

if isempty(tt)
    tt = IDX;
end

time.data       = tt;                               % time vector
Ts              = mean(diff(time.data));        	% sampling time
pos.data        = xx;                               % position
vel.data        = [0 ; diff(pos.data)/Ts];          % velocity
avel.data       = abs(vel.data);                    % absolute velocity
vel.mean      	= mean(avel.data);                  % mean absolute velocity
vel.std       	= std(vel.data);                    % std absolute velocity
svel            = avel.data;                        % preallocate velocity data above threshold

% Variable names
varnames = {'Duration','Amplitude','Direction','StartIdx','PeakIdx','EndIdx','StartTime','PeakTime','EndTime','StartPos',...
                'PeakPos','EndPos','StartVel','PeakVel','EndVel','Threshold'};

% Maually ind velocity spikes
% [avel.pks, loc.pks] = findpeaks(svel,'MINPEAKDISTANCE',40); % find local maxima
fig(1) = figure; clf
set(fig,'Color','w','Units','inches','Position',[2,2,11,9])
movegui(fig(1),'center')
CC = [0 0 0];
ax(1) = subplot(2,1,1);
    plot(IDX(4:end),pos.data(4:end),'Color',CC)
    ylabel('Position')
ax(2) = subplot(2,1,2);
    plot(IDX(4:end),vel.data(4:end),'Color',CC)
    ylabel('Velocity')
    xlabel('Index')
    ylim(max(abs(ax(2).YLim))*[-1 1])
set(ax,'XLim',[min(IDX),max(IDX)])
linkaxes(ax,'x')

[loc.pks,pts] = getpts(ax(2));
loc.pks = round(loc.pks);
loc.pks = sort(loc.pks);
loc.pks = unique(loc.pks);
count = length(loc.pks);
tol = 0.2;
tol = round(tol/Ts);
for ww = 1:count % every saccade
    tpk = loc.pks(ww);
    tlow = tpk - tol;
    thigh = tpk + tol;
    if tlow<1
        tlow = 1;
    end
    if thigh<1
        tlow = 1;
    end
    trange = tlow:thigh;
    idxrange = zeros(size(IDX));
    idxrange(trange) = 1;
    velrange = vel.data;
    if pts(ww)<0
        velrange = -velrange;
    end
    velrange(~idxrange) = 0;
    [avel.pks(ww),loc.pks(ww)] = max(velrange);
end
window      = n*0.015;  % window length (in samples)
I           = find(loc.pks > (window) & loc.pks < ((length(vel.data) - (window+1)))); % ignore saccades at start & end
time.pks    = time.data(loc.pks);
pos.pks     = pos.data(loc.pks);
vel.pks     = vel.data(loc.pks);

ax(1) = subplot(2,1,1); hold on
    plot(loc.pks,pos.pks,'.b','MarkerSize',15)
ax(2) = subplot(2,1,2); hold on
    plot(loc.pks,vel.pks,'.b','MarkerSize',15)

rate = count/(time.data(end)-time.data(1));
boundThresh = 1/4;
if ~isempty(I)
    thresh      = min(avel.pks);
    for ww = 1:count % every saccade
        % Find start & end velocity
        loc.start(ww,1) = find(avel.data(1:loc.pks(ww)) <= avel.pks(ww)/4,1,'last'); % saccade start index               
        Eind = find(avel.data <= avel.pks(ww)*boundThresh); % all values below 1/4 peak
        Es = find(Eind > loc.pks(ww,1),1,'first'); % first value after the start index is the end index
        if ~isempty(Es) % make sure data did not start above threshold
            loc.end(ww,1) = Eind(Es); % saccade end index
        else
            loc.end(ww,1) = 1;
        end
    end
    
    % Stats
    avel.start  	= avel.data(loc.start);
    avel.end        = avel.data(loc.end);
    vel.start       = vel.data(loc.start);
    vel.end         = vel.data(loc.end);
    pos.start       = pos.data(loc.start);
    pos.end         = pos.data(loc.end);
    time.start   	= time.data(loc.start);
    time.end      	= time.data(loc.end);
    pos.amp         = pos.end - pos.start;
    time.dur        = time.end - time.start;
    dir             = sign(vel.pks);
    
    ax(1) = subplot(2,1,1); hold on
        plot(loc.start,pos.start,'.g','MarkerSize',15)
        plot(loc.end,pos.end,'.r','MarkerSize',15)
    ax(2) = subplot(2,1,2); hold on
        plot(loc.start,vel.start,'.g','MarkerSize',15)
    	plot(loc.end,vel.end,'.r','MarkerSize',15)

    
    THRESH = repmat(thresh,count,1); % assign threshold used to all saccades
    
    % All data
    DATA = [time.dur , pos.amp , dir, loc.start , loc.pks , loc.end , time.start , time.pks , time.end , pos.start , pos.pks, ...
                pos.end , vel.start , vel.pks , vel.end , THRESH];
else % no saccades
    DATA    = nan(1,length(varnames));
    count   = nan;
    rate    = nan;
    thresh  = nan;
end

% Saccade table
SACD = splitvars(table(DATA));
SACD.Properties.VariableNames = varnames;

% Remove Saccades
time.rmv = time.data;
pos.rmv = pos.data;
vel.rmv = vel.data;
if ~isempty(I)
    loc.cell_idx    = cell(count,1);
    loc.cell_all    = cell(count,1);
    time.cell_rmv   = cell(count,1);
 	pos.cell_rmv    = cell(count,1);
 	vel.cell_rmv    = cell(count,1);
    for ww = 1:count % every saccade
        loc.cell_idx{ww}    = (loc.start(ww):loc.end(ww))';
        loc.cell_all{ww}    = IDX( ismember( IDX, loc.cell_idx{ww} ) );
        time.cell_rmv{ww}   = time.data(loc.cell_all{ww});
       	pos.cell_rmv{ww}    = pos.data(loc.cell_all{ww});
       	vel.cell_rmv{ww}    = vel.data(loc.cell_all{ww});
    end
    loc.idx = cat(1,loc.cell_idx{:});
    loc.all = ismember(IDX, loc.idx );
    
    time.rmv(loc.idx) = nan;
    pos.rmv(loc.idx)  = nan;
    vel.rmv(loc.idx)  = nan;
    
    time.nan = time.rmv;
    pos.rmv_n = pos.rmv;
    
    for ww = 1:count
        pos.rmv_n(loc.end(ww):end) = pos.rmv_n(loc.end(ww):end) - (pos.end(ww) - pos.start(ww));
    end
  	pos.nan         = pos.rmv_n;
    pos.rmv_n       = pos.rmv_n(~isnan(pos.rmv_n));
    time.rmv        = time.rmv(~isnan(time.rmv));
    pos.rmv_intrp   = interp1(time.rmv,pos.rmv_n,time.data,'linear');

else
    pos.rmv_intrp = pos.data;
end

SaccdRmv_X = pos.rmv_intrp;

% Debug plots
if debug
    FIG = figure ; clf
    FIG.Color = 'w';
    
    ax1 = subplot(3,1,1) ; hold on
        ylabel('Velocity')
        plot(time.data(4:end),vel.data(4:end),'k')
        plot(time.data, zeros(n,1)*thresh,'--','Color',[0.5 0.5 0.5])
        plot(time.data, ones(n,1)*thresh,'--m')
        plot(time.data,-ones(n,1)*thresh,'--m')
        for kk = 1:count
            plot(time.cell_rmv{kk},vel.cell_rmv{kk},'-c','MarkerSize',5)
        end
        ax1.YLim = max(abs(ax1.YLim))*[-1 1];

        plot(time.pks,vel.pks,'*b')
        plot(time.start,vel.start,'*g')
        plot(time.end,vel.end,'*r')
    
    ax2 = subplot(3,1,2) ; hold on
        ylabel('Position')
        plot(time.data,pos.data,'k')
        plot(time.data, zeros(n,1)*thresh,'--','Color',[0.5 0.5 0.5])
        for kk = 1:count
            plot(time.cell_rmv{kk},pos.cell_rmv{kk},'-c','MarkerSize',5)
        end
        plot(time.pks,pos.pks,'*b')
        plot(time.start,pos.start,'*g')
        plot(time.end,pos.end,'*r')
        ax2.YLim = max(abs(ax2.YLim))*[-1 1];
        
    ax3 = subplot(3,1,3) ; hold on
        plot(time.data,pos.rmv_intrp,'c','LineWidth',1)
        plot(time.nan ,pos.nan,'k','LineWidth',1.5)
        ylabel('Position (Removed Saccades)')
        xlabel('Time')
end


end