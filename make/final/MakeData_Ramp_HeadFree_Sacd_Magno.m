function [] = MakeData_Ramp_HeadFree_Sacd_Magno()
%% MakeData_Ramp_HeadFree_Sacd_HeadWing:
%   INPUTS:
%       -
%
%   OUTPUTS:
%       -
%

% Data location
rootdir = 'H:\EXPERIMENTS\MAGNO\Experiment_Ramp';

% Output file name
filename = 'Ramp_HeadFree_SACCD_Magno';

% Setup Directories 
PATH.daq  = rootdir; % DAQ data location
PATH.body = fullfile(PATH.daq,'tracked_body');
PATH.reg  = fullfile(PATH.daq,'registered');
PATH.head = fullfile(PATH.reg,'tracked_head');
% PATH.wing = fullfile(PATH.reg,'wing_filt', 'tracked_head_wing');

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.head,'*.mat',false,'fly','trial','vel');
% [D,I,N,U,T,~,~,basename] = GetFileData(PATH.wing,'*.csv',false,'fly','trial','vel','wave');

%% Get Data %%
clc
disp('Loading...')
Fs = 100; % sampling frequency [s]
tintrp = (0:(1/Fs):15)'; % time vector for interpolation

% HEAD saccade detection parameters
head.showplot = false;
head.Fc_detect = [10 nan];
head.Fc_ss = [30 nan];
head.amp_cut = 4;
head.thresh = 70;
head.sacd_length = nan;
head.pks = [];
head.min_pkdist = 0.05;
head.min_pkwidth = 0.02;
head.min_pkprom = 50;
head.min_pkthresh = 0;
head.boundThresh = 0.40;

% BODY saccade detection parameters
body.showplot = false;
body.Fc_detect = [10 nan];
body.Fc_ss = [30 nan];
body.amp_cut = 7;
body.thresh = 50;
body.sacd_length = nan;
body.pks = [];
body.min_pkdist = 0.05;
body.min_pkwidth = 0.02;
body.min_pkprom = 20;
body.min_pkthresh = 0;
body.boundThresh = 0.25;

Vel = 3.75*U.vel{1}; % velocities
Stim = (Vel*tintrp')'; % stimuli
SACCADE = [I , splitvars(table(num2cell(zeros(N.file,3))))]; % store saccade objects
SACCADE.Properties.VariableNames(4:6) = {'head_saccade','body_saccade','head2body'};
SACCADE_STATS.head = []; % store saccade stats
SACCADE_STATS.body = []; % store saccade stats
HEAD_DATA = cell(N.fly,N.vel);
for kk = 1:N.file
    disp(kk)
    % Load HEAD, BODY, & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p');
    load(fullfile(PATH.head, [basename{kk} '.mat']),'hAngles');
    load(fullfile(PATH.body, [basename{kk} '.mat']),'bAngles');
 
  	% Sync frames and get pattern data
	[trig,~] = sync_pattern_trigger(t_p, data(:,2), 15, data(:,1), true, [], false, false);
    if length(hAngles) ~= length(trig.time_sync)
        trig.time_sync = [trig.time_sync ; trig.time_sync(end) + (1/mean(diff(trig.time_sync)))];
    end
    
    % Get data
    head.pos = interp1(trig.time_sync, hAngles, tintrp, 'pchip');
    body.pos = interp1(trig.time_sync, bAngles, tintrp, 'pchip');
       
    % Extract saccades
    direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    %direction = 0;
    head_saccade = saccade_all(head.pos, tintrp, head.thresh, head.Fc_detect, head.Fc_ss, head.amp_cut, ...
                                direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, head.showplot);
    head_saccade = stimSaccade(head_saccade, Stim(:,I.vel(kk)), false);
    
    % Extract body saccades
    body_saccade = saccade_all(body.pos, tintrp, body.thresh, body.Fc_detect, body.Fc_ss, ...
                            body.amp_cut, direction, body.pks, body.sacd_length, ...
                            body.min_pkdist, body.min_pkwidth, body.min_pkprom, ...
                            body.min_pkthresh, body.boundThresh, body.showplot);
 	body_saccade = stimSaccade(body_saccade, Stim(:,I.vel(kk)), false);
    
    % Store saccade objects
    SACCADE.head_saccade(kk) = {head_saccade};
    SACCADE.body_saccade(kk)  = {body_saccade};

    % Create STAT tables
    if head_saccade.count == 0
        rep = 1;
    else
        rep = head_saccade.count;
        HEAD_DATA{I.fly(kk),I.vel(kk)}(end+1,1) = head_saccade;
    end
    VTable = table(D.vel(kk),'VariableNames',{'Vel'});
    ITable = [I(kk,:),VTable];
    ITable = repmat(ITable,rep,1);
    SACCADE_STATS.head = [SACCADE_STATS.head ; [ITable , head_saccade.SACD]];
    
    if head.showplot
        figure (1)
        pause
        close all
    end
            
    if body_saccade.count == 0
        rep = 1;
    else
        rep = body_saccade.count;
    end
    VTable = table(D.vel(kk),'VariableNames',{'Vel'});
    ITable = [I(kk,:),VTable];
    ITable = repmat(ITable,rep,1);
    SACCADE_STATS.body = [SACCADE_STATS.body  ; [ITable , body_saccade.SACD]];

    if body.showplot
        pause
        close all
    end
    
    % Head ===> Body
    if head_saccade.count > 0
        head2body = saccade_interact(head_saccade, body_saccade, 0.5, 0.1, true);
        SACCADE.head2body(kk) = {head2body};
    end
end
disp('Done')

%%
files = 1;
keepI = cellfun(@(x) isstruct(x) | isobject(x), SACCADE.head2body);
Saccade = SACCADE(keepI,:);
Saccade = Saccade(:,:);

velI = Saccade.vel;
cw = velI <= 5;
clear Head Body Time

Time.align = cellfun(@(x) x.interval.time_align, Saccade.head2body, 'UniformOutput', false);
Time.align = nanmean(cat(2,Time.align{:}),2);

Head.pos = cellfun(@(x) x.interval_rel.in.pos , Saccade.head2body, 'UniformOutput', false);
Head.vel = cellfun(@(x) x.interval_rel.in.vel , Saccade.head2body, 'UniformOutput', false);
Head.pos(cw) = cellfun(@(x) -x , Head.pos(cw), 'UniformOutput', false);
Head.vel(cw) = cellfun(@(x) -x , Head.vel(cw), 'UniformOutput', false);
Head = structfun(@(x) cat(2,x{:}), Head, 'UniformOutput', false);
Head.stats = structfun(@(x) basic_stats(x,2), Head, 'UniformOutput', false);

Body.pos = cellfun(@(x) x.interval.out.pos , Saccade.head2body, 'UniformOutput', false);
Body.vel = cellfun(@(x) x.interval.out.vel , Saccade.head2body, 'UniformOutput', false);
Body.pos(cw) = cellfun(@(x) -x , Body.pos(cw), 'UniformOutput', false);
Body.vel(cw) = cellfun(@(x) -x , Body.vel(cw), 'UniformOutput', false);
Body = structfun(@(x) cat(2,x{:}), Body, 'UniformOutput', false);
Body.stats = structfun(@(x) basic_stats(x,2), Body, 'UniformOutput', false);

fig = figure (102); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 5])
movegui(fig, 'center')
clear ax
hcolor = [0 0 1 ];
bcolor = [0.1 1 0.5];
alpha = 0.3;
ax(1) = subplot(2,1,1); cla ; hold on ; ylabel('(°)')
    plot(Time.align, Body.pos, 'Color', [0.7*bcolor alpha], 'LineWidth', 0.25)
    plot(Time.align, Head.pos, 'Color', [0.7*hcolor alpha], 'LineWidth', 0.25)
    [~,h.mean(1)] = PlotPatch(Body.stats.pos.mean, Body.stats.pos.std, ...
        Time.align, 1, 1, bcolor, 0.7*bcolor, 0.3, 1);
    [~,h.mean(2)] = PlotPatch(Head.stats.pos.mean, Head.stats.pos.std, ...
        Time.align, 1, 1, hcolor, 0.7*hcolor, 0.3, 1);
    
ax(2) = subplot(2,1,2); cla ; hold on ; ylabel('(°/s)')
    plot(Time.align, Body.vel, 'Color', [0.7*bcolor alpha], 'LineWidth', 0.25)
    plot(Time.align, Head.vel, 'Color', [0.7*hcolor alpha], 'LineWidth', 0.25)
    [~,~] = PlotPatch(Body.stats.vel.mean, Body.stats.vel.std, ...
        Time.align, 1, 1, bcolor, 0.7*bcolor, 0.3, 1);
    [~,~] = PlotPatch(Head.stats.vel.mean, Head.stats.vel.std, ...
        Time.align, 1, 1, hcolor, 0.7*hcolor, 0.3, 1);
    
    xlabel('Time (s)')
    legend(h.mean, 'Body', 'Head', 'Box', 'off')
    
set(ax, 'LineWidth', 1, 'XLim', 0.2*[-1 1], 'XTick', -0.2:0.05:0.2, 'Box', 'on')
set(ax(1), 'YLim', 40*[-1 1])
set(ax(2), 'YLim', [-600 1600])
linkaxes(ax,'x')

%%
TD = cellfun(@(x) x.TimeDiff, Saccade.head2body, 'UniformOutput', false);
TD = cat(1,TD{:});


%% Extract & group saccades & intervals by speed & by fly
fields = {'normpeak_saccade','norm_interval','normstart_interval','normend_interval',...
    'normstart_stimulus','error','int_error'};
nfield = length(fields);
norm_fields = {'time','position','velocity'};
norm_fields_stats = string(norm_fields) + "_stats";

center = 0; % normalization center for saccades & inter-saccade intervals
dim = 1; % dimension to center
dim_stats = 2; % dimension for statistics

FLY = []; % all trials per speed per fly
GRAND = []; % all trials per speed
for kk = 1:nfield % for each field in saccade structure
    if kk < 2 % even padding around center only for saccades (not intervals)
        even = true;
    else
        even = false;
    end
    
    FLY.(fields{kk}) = cellfun(@(x) struct_center([x.(fields{kk})], center, even, dim, norm_fields), ...
                            HEAD_DATA, 'UniformOutput', true);
    for jj = 1:N.vel
        GRAND.(fields{kk})(jj) = struct_center(FLY.(fields{kk})(:,jj), center, even , dim, norm_fields);
        
        % Calculate stats by fly
        for ff = 1:length(norm_fields)            
            for ww = 1:N.fly
                FLY.(fields{kk})(ww,jj).(norm_fields_stats(ff)) = basic_stats( ...
                    FLY.(fields{kk})(ww,jj).(norm_fields{ff}), dim_stats);
            end
        end
    end
end

% Calculate stats for all
for kk = 1:nfield % for each field in saccade structure
    for jj = 1:N.vel
        for ff = 1:length(norm_fields)
            GRAND.(fields{kk})(jj).(norm_fields_stats(ff)) = basic_stats( ...
                GRAND.(fields{kk})(jj).(norm_fields{ff}), dim_stats);
        end
    end
end

%% SAVE %%
% disp('Saving...')
% save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
%       'PATH','SACCADE_STATS','SACCADE','FLY','GRAND',...
%       'Stim','D','I','U','N','T','-v7.3')
% disp('SAVING DONE')
end