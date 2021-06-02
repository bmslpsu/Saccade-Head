function [] = Make_Ramp_HeadFree_head_wing_glue_iron()
%% Make_Ramp_HeadFree_head_wing_leg:
%   INPUTS:
%       wave    :   spatial wavelength of data
%
%   OUTPUTS:
%       -
%
warning('off', 'signal:findpeaks:largeMinPeakHeight')

% Data location
% rootdir = 'E:\EXPERIMENTS\RIGID\Experiment_ramp_glue';
rootdir = 'E:\EXPERIMENTS\RIGID\Experiment_ramp_glue_head';

% Output file name
filename = 'Ramp_HeadFree_head_glue';

% Setup Directories
PATH.daq  = rootdir;
PATH.vid  = fullfile(PATH.daq,'');
PATH.head = fullfile(PATH.vid,'tracked_head');
PATH.wing = fullfile(PATH.vid,'filt', 'tracked_head_wing');

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.head,'*.mat',false,'fly','trial','vel','wave');
% [D,I,N,U,T,~,~,basename] = GetFileData(PATH.wing,'*.csv',false,'fly','trial','vel','wave');

%% Get Data
clc
close all
disp('Loading...')
Fs = 200; % sampling frequency [s]
tintrp = (0:(1/Fs):(10 - 1/Fs))'; % time vector for interpolation

% HEAD saccade detection parameters
head.showplot = true;
head.Fc_detect = [10 nan];
head.Fc_ss = [nan nan];
head.amp_cut = 4;
head.dur_cut = inf;
head.thresh = [70, 2, 1, 0];
head.true_thresh = 250;
head.sacd_length = nan;
head.pks = [];
head.min_pkdist = 0.2;
head.min_pkwidth = 0.03;
head.min_pkprom = 75;
head.min_pkthresh = 0;
head.boundThresh = [0.2 60];
head_carry.Fc = 40;
[head_carry.b, head_carry.a] = butter(3, head_carry.Fc / (Fs/2) ,'low');

% WING saccade detection parameters
wing.showplot = false;
wing.Fc_detect = [5 nan];
wing.Fc_ss = [5 nan];
wing.amp_cut = 4;
wing.dur_cut = inf;
wing.thresh = [50, 2, 0, 0];
wing.true_thresh = 70;
wing.sacd_length = nan;
wing.pks = [];
wing.min_pkdist = 0.3;
wing.min_pkwidth = 0.03;
wing.min_pkprom = 20;
wing.min_pkthresh = 0;
wing.boundThresh = 0.35;
wing.Fc = 8;
[wing.b, wing.a] = butter(3, wing.Fc / (Fs/2) ,'low');
wing_carry.Fc = 40;
[wing_carry.b, wing_carry.a] = butter(3, wing_carry.Fc / (Fs/2) ,'low');

Vel = U.vel{1}; % velocities
Stim = (Vel*tintrp')'; % stimuli
SACCADE = [D , splitvars(table(num2cell(zeros(N.file,10))))]; % store saccade objects
SACCADE.Properties.VariableNames(5:14) = {'head_saccade','head_scd_pos','head_scd_vel', 'head_scd_accel','scd_time', ...
                                          'wing_saccade', 'wing_scd_pos', 'wing_scd_vel','head2wing','leg'};
HEAD_DATA = cell(N.fly,N.vel);
HEAD_SACCADE_STATS = [];
WING_SACCADE_STATS = [];
for kk = 1:N.file
    disp(kk)
    basename{kk}
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p');
    %load(fullfile(PATH.head, [basename{kk} '.mat']),'head_data');
    load(fullfile(PATH.head, [basename{kk} '.mat']),'hAngles');
 
  	% Sync frames and get pattern data
	[trig,pat] = sync_pattern_trigger(t_p, data(:,2), 10, data(:,1), true, 1, false, false);

    % Get head data
    head.pos = interp1(trig.time_sync, hAngles, tintrp, 'pchip');
    
	% Head with filter
 	head.pos_filt = filtfilt(head_carry.b, head_carry.a, head.pos);
       
    % Extract head saccades
    direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    head_saccade = saccade_all(head.pos_filt, tintrp, head.thresh, head.true_thresh, head.Fc_detect, ...
                                head.Fc_ss, head.amp_cut, head.dur_cut, direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, head.showplot);
    head_saccade = stimSaccade(head_saccade, Stim(:,I.vel(kk)), false); % with approximate pattern position
    SACCADE.head_saccade(kk) = {head_saccade}; % store data in cell

    if any(head_saccade.peaks.velocity > 1200)
       disp('here') 
    end
    
    if head_saccade.count == 0
        rep = 1;
    else
        HEAD_DATA{I.fly(kk),I.vel(kk)}(end+1,1) = head_saccade;
        rep = head_saccade.count;
        
        [scds,~,~,~] = getSaccade(head_saccade, head_saccade.position, 0.5, 0.5, true);
        SACCADE.head_scd_pos(kk) = {cat(2,scds{:})};
        [scds,~,~,~] = getSaccade(head_saccade, head_saccade.velocity, 0.5, 0.5, true);
        SACCADE.head_scd_vel(kk) = {cat(2,scds{:})};
        [scds,~,scd_time,~] = getSaccade(head_saccade, head_saccade.acceleration, 0.5, 0.5, true);
        SACCADE.head_scd_accel(kk) = {cat(2,scds{:})};
        SACCADE.scd_time(kk) = scd_time(1);
    end
    VTable = table(D.vel(kk),'VariableNames',{'Vel'});
    ITable = [I(kk,:),VTable];
    ITable = repmat(ITable,rep,1);
    HEAD_SACCADE_STATS = [HEAD_SACCADE_STATS ; [ITable , head_saccade.SACD]];
    
    if head.showplot
        figure (1)
        pause
        close all
    end
    
	% Load WING data if avialable
    wfile = fullfile(PATH.wing, [basename{kk} '.csv']);
    if exist(wfile,'file') == 2
        benifly = ImportBenifly(wfile); % load wing angles
        [b_h,a_h] = butter(3, 0.1 / (Fs/2), 'high');
        
        % Get wing data
        hampel_dx = round(0.05*Fs);
        wing.left       = rad2deg(benifly.LWing);
        wing.right      = rad2deg(benifly.RWing);
     	wing.left       = hampel(trig.time_sync, wing.left, hampel_dx);
        wing.right      = hampel(trig.time_sync, wing.right, hampel_dx);
        wing.left       = interp1(trig.time, wing.left, tintrp, 'pchip');
        wing.right      = interp1(trig.time, wing.right, tintrp, 'pchip');

        wing.dwba       = wing.left - wing.right;
        %wing.dwba       = filtfilt(b_h, a_h, wing.dwba);
        wing.dwba       = filtfilt(wing_carry.b, wing_carry.a, wing.dwba);
        wing.dwba_vel   = central_diff(wing.dwba, 1/Fs);
        wing.left_filt	= filtfilt(wing.b, wing.a, wing.left);
        wing.right_filt	= filtfilt(wing.b, wing.a, wing.right);
        wing.dwba_filt	= wing.left_filt - wing.right_filt;
        
%         cla ; hold on
%             plot(tintrp, head.pos_filt, 'k', 'LineWidth', 1)
%             %plot(tintrp, wing.dwba, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5)
%             plot(tintrp, wing.dwba_filt - mean(wing.dwba_filt), 'r', 'LineWidth', 1)
%             ylim(20*[-1 1])
%             disp(basename{kk})
%             pause
        
      	% Extract wing saccades
        wing_saccade = saccade_all(wing.dwba_filt, tintrp, wing.thresh, wing.true_thresh, wing.Fc_detect, ...
                                wing.Fc_ss, wing.amp_cut, wing.dur_cut, direction, wing.pks, wing.sacd_length, ...
                                wing.min_pkdist, wing.min_pkwidth, wing.min_pkprom, ...
                                wing.min_pkthresh, wing.boundThresh, wing.showplot);
     	wing_saccade.extra.dwba = wing.dwba; % carry unfiltered dwba singal
        wing_saccade.extra.dwba_vel = wing.dwba_vel; % carry unfiltered dwba velocity singal
        wing_saccade.extra.lwing = wing.left; % carry unfiltered left wba singal
        wing_saccade.extra.rwing = wing.right; % carry unfiltered right wba singal
        
%         if wing.showplot
%             figure (1) ; subplot(4,1,1)
%             plot(tintrp, head.pos_filt, 'm', 'LineWidth', 0.25)
%             plot(tintrp, wing.dwba, 'b', 'LineWidth', 0.25)
%         end
        
        if wing_saccade.count == 0
            rep = 1;
        else
            rep = wing_saccade.count;
        end
        VTable = table(D.vel(kk),'VariableNames',{'Vel'});
        ITable = [I(kk,:),VTable];
        ITable = repmat(ITable,rep,1);
        WING_SACCADE_STATS = [WING_SACCADE_STATS ; [ITable , wing_saccade.SACD]];
        
        SACCADE.wing_saccade(kk) = {wing_saccade};
        
        if wing.showplot
            pause
            close all
        end

        if head_saccade.count > 0
            [scds,~,~,~] = getSaccade(head_saccade, wing.dwba, 0.5, 0.5, true);
            SACCADE.wing_scd_pos(kk) = {cat(2,scds{:})};
            [scds,~,~,~] = getSaccade(head_saccade, wing.dwba_vel, 0.5, 0.5, true);
            SACCADE.wing_scd_vel(kk) = {cat(2,scds{:})};
        end
    end
end

% Fill in empty saccade trials
empty_idx = cellfun(@(x) isempty(x), HEAD_DATA);
HEAD_DATA(empty_idx) = {saccade_all(0*head.pos, tintrp, head.thresh, head.true_thresh, head.Fc_detect, ...
                                head.Fc_ss, head.amp_cut, head.dur_cut, direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, false)};

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

%% Distributions
head_pos = cellfun(@(x) x.position, SACCADE.head_saccade, 'UniformOutput', false);
DATA = [SACCADE , table(head_pos, 'VariableNames', {'head_pos'})];
[stats] = table_fly_stats(DATA, 3, 15, false);

stat_name = 'head_pos';
n_speed = N.vel/2;
vel = U.vel{1};

hbins = -21.5:1:21.5;
% hbins = -21.25:0.5:22.15;

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 N.fly*2 3])
movegui(fig, 'center')
ax = gobjects(n_speed,N.fly);
h = gobjects(n_speed,N.fly,2);
pp = 1;
P = nan(n_speed,N.fly);
for v = 1:n_speed
    for n = 1:N.fly
        ax(v,n) = subplot(n_speed,N.fly,pp); cla ; hold on
            %hdata = [stats.(stat_name){n,v} , stats.(stat_name){n,v+n_speed}];
        
        h(v,n,1) = histogram(stats.(stat_name){n,v}, hbins, 'Normalization', 'probability', ...
            'EdgeColor', 'none', 'FaceColor', 'b');
        h(v,n,2) = histogram(stats.(stat_name){n,v+n_speed}, hbins, 'Normalization', 'probability', ...
            'EdgeColor', 'none', 'FaceColor', 'r');
        xline(0, '--', 'Color', [0.5 0.5 0.5])
        
        pp = pp + 1;
    end
end
linkaxes(ax, 'xy')
set(ax, 'Color', 'none', 'LineWidth', 1, 'XColor', 'k')
% set(ax, 'YLim', [-0.1 2])
set(ax(:,2:end), 'YColor', 'none')

% set(b(1,:), 'BoxFaceColor', 'r')
% set(b(2,:), 'BoxFaceColor', 'b')

YLabelHC = get(ax(:,1), 'YLabel');
set([YLabelHC{:}], 'String', 'Probability', 'Interpreter', 'none')

%% Rate by fly
scd_rate = cellfun(@(x) x.rate, SACCADE.head_saccade);
DATA = [SACCADE , table(num2cell(scd_rate), 'VariableNames', {'rate'})];
[stats] = table_fly_stats(DATA, 3, 15, false);

stat_name = 'rate';
n_speed = N.vel/2;
vel = U.vel{1};

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 N.fly*2 3])
movegui(fig, 'center')
ax = gobjects(n_speed,N.fly);
b = gobjects(n_speed,N.fly);
pp = 1;
P = nan(n_speed,N.fly);
for v = 1:n_speed
    for n = 1:N.fly
        ax(v,n) = subplot(n_speed,N.fly,pp); cla ; hold on
        x = stats.(stat_name){n,v};
        y = stats.(stat_name){n,v+n_speed};
        z = [x , y];
        G = [vel(v)*ones(size(stats.(stat_name){n,v})) , ...
            vel(v+n_speed)*ones(size(stats.(stat_name){n,v+n_speed}))];
        [~,p] = ttest2(x,y);
        P(v,n) = p;
        b(v,n) = boxchart(categorical(G), z, 'LineWidth', 1, 'MarkerStyle', '.', ...
            'JitterOutliers', 'on', 'Notch', 'off');
        
        pp = pp + 1;
    end
end
linkaxes(ax, 'y')
set(ax, 'Color', 'none', 'LineWidth', 1, 'XColor', 'k')
set(ax, 'YLim', [-0.1 1.5])
set(ax(:,2:end), 'YColor', 'none')

set(b(1,:), 'BoxFaceColor', 'r')
set(b(2,:), 'BoxFaceColor', 'b')

YLabelHC = get(ax(:,1), 'YLabel');
set([YLabelHC{:}], 'String', stat_name)

%% Dynamics by fly
[stats] = table_fly_stats(HEAD_SACCADE_STATS, 5, [6 7 8 18 20 22 24 25], false);

stat_name = 'StartPos';
n_speed = N.vel/2;
vel = U.vel{1};

fig = figure (3) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 N.fly*2 3])
movegui(fig, 'center')
ax = gobjects(n_speed,N.fly);
b = gobjects(n_speed,N.fly);
pp = 1;
for v = 1:n_speed
    for n = 1:N.fly
        ax(v,n) = subplot(n_speed,N.fly,pp); cla ; hold on
        y = [stats.(stat_name){n,v} , stats.(stat_name){n,v+n_speed}];
        dir = -[stats.Direction{n,v} , stats.Direction{n,v+n_speed}];
        
        z = y .* dir;
        G = [vel(v)*ones(size(stats.(stat_name){n,v})) , ...
            vel(v+n_speed)*ones(size(stats.(stat_name){n,v+n_speed}))];
        
        b(v,n) = boxchart(categorical(G), z, 'LineWidth', 1, 'MarkerStyle', '.', ...
            'JitterOutliers', 'on', 'Notch', 'off');
        
        pp = pp + 1;
    end
end
linkaxes(ax, 'y')
set(ax, 'Color', 'none', 'LineWidth', 1, 'XColor', 'k')
% set(ax, 'YLim', [-0.1 2])
set(ax(:,2:end), 'YColor', 'none')

set(b(1,:), 'BoxFaceColor', 'r')
set(b(2,:), 'BoxFaceColor', 'b')

YLabelHC = get(ax(:,1), 'YLabel');
set([YLabelHC{:}], 'String', stat_name)

%% SAVE %%
disp('Saving...')
save(['E:\DATA\Rigid_Data\Saccade\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','HEAD_DATA','SACCADE','HEAD_SACCADE_STATS','WING_SACCADE_STATS','FLY','GRAND',...
      'Stim','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end