function [] = MakeData_Static_HeadFree_Sacd_HeadWing()
%% MakeData_Static_HeadFree_Sacd_HeadWing:
%   INPUTS:
%       -
%
%   OUTPUTS:
%       -
%

direction = 0; % get saccades in all directions

switch direction
    case 0
        dirlabel = 'ALL';
    case 1
        dirlabel = 'CW';
    case -1
        dirlabel = 'CCW';
    otherwise
        error('direction must be 0,1, or -1')
end

% Data location
rootdir = 'H:\EXPERIMENTS\RIGID\Experiment_Static_Wave';

% Output file name
filename = ['Static_HeadFree_SACCD_' dirlabel '_HeadWingALL_Wave'];

% Setup Directories 
PATH.daq  = rootdir;
PATH.vid  = fullfile(PATH.daq,'');
PATH.head = fullfile(PATH.vid,'tracked_head');
PATH.wing = fullfile(PATH.vid,'wing_filt', 'tracked_head_wing');

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.head,'*.mat',false,'fly','trial','wave','vel');
% [D,I,N,U,T,~,~,basename] = GetFileData(PATH.wing,'*.csv',false,'fly','trial','wave','vel');

%% Get Data %%
clc
disp('Loading...')
Fs = 200; % sampling frequency [s]
tintrp = (0:(1/Fs):(10 - 1/Fs))'; % time vector for interpolation

% HEAD saccade detection parameters
head.showplot = false;
head.Fc_detect = [10 nan];
head.Fc_ss = [40 nan];
head.amp_cut = 4;
head.thresh = 80;
head.true_thresh = 250;
head.sacd_length = nan;
head.pks = [];
head.min_pkdist = 0.1;
head.min_pkwidth = 0.03;
head.min_pkprom = 75;
head.min_pkthresh = 0;
head.boundThresh = [0.25 50];
head_carry.Fc = 60;
[head_carry.b, head_carry.a] = butter(3, head_carry.Fc / (Fs/2) ,'low');

% WING saccade detection parameters
wing.showplot = false;
wing.Fc_detect = [5 nan];
wing.Fc_ss = [5 nan];
wing.amp_cut = 7;
wing.thresh = 40;
wing.true_thresh = [];
wing.sacd_length = nan;
wing.pks = [];
wing.min_pkdist = 0.5;
wing.min_pkwidth = 0.03;
wing.min_pkprom = 20;
wing.min_pkthresh = 0;
wing.boundThresh = 0.35;
wing.Fc = 8;
[wing.b, wing.a] = butter(3, wing.Fc / (Fs/2) ,'low');
wing_carry.Fc = 40;
[wing_carry.b, wing_carry.a] = butter(3, wing_carry.Fc / (Fs/2) ,'low');

SACCADE = [I , splitvars(table(num2cell(zeros(N.file,4))))]; % store saccade objects
SACCADE.Properties.VariableNames(5:8) = {'head_saccade','wing_saccade','head2wing','head2wing_align_wing'};
HEAD_DATA = cell(N.fly,N.wave);
HEAD_SACCADE_STATS = []; % store saccade stats
WING_SACCADE_STATS = []; % store saccade stats
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    load(fullfile(PATH.head, [basename{kk} '.mat']),'hAngles'); % load head angles % time arrays
 
  	% Sync frames and get pattern data
	[trig,~] = sync_pattern_trigger(t_p, data(:,2), 10, data(:,1), true, 1, false, false);
    
    % Get head data
    head.pos = interp1(trig.time_sync, hAngles, tintrp, 'pchip');
    
	% Head with filter
 	head.pos_filt = filtfilt(head_carry.b, head_carry.a, head.pos);
       
    % Extract head saccades
    head_saccade = saccade_all(head.pos_filt, tintrp, head.thresh, head.true_thresh, head.Fc_detect, ...
                                head.Fc_ss, head.amp_cut, direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, head.showplot);
    SACCADE{kk,5} = {head_saccade}; % store data in cell
    
    if head_saccade.count == 0
        rep = 1;
    else
        HEAD_DATA{I.fly(kk),I.wave(kk)}(end+1,1) = head_saccade;
        rep = head_saccade.count;
    end
    VTable = table(D.wave(kk),'VariableNames',{'Wave'});
    ITable = [I(kk,:),VTable];
    ITable = repmat(ITable,rep,1);
    HEAD_SACCADE_STATS = [HEAD_SACCADE_STATS ; [ITable , head_saccade.SACD]];
    
    if head.showplot
        figure (1)
        pause
        close all
    end
    
    if any(head_saccade.SACD.Duration > 0.09)
        plotSaccade(head_saccade)
        plotInterval(head_saccade)
        pause
        close all
    end
    
	% Load WING data if we have it
    wfile = fullfile(PATH.wing, [basename{kk} '.csv']);
    if exist(wfile,'file') == 2
        benifly = ImportBenifly(wfile); % load wing angles
        
        % Get wing data
        wing.left       = hampel(trig.time_sync, benifly.LWing);
        wing.right      = hampel(trig.time_sync, benifly.RWing);
        wing.left       = rad2deg(interp1(trig.time, wing.left, tintrp, 'pchip'));
        wing.right      = rad2deg(interp1(trig.time, wing.right, tintrp, 'pchip'));
        wing.dwba       = wing.left - wing.right;
        wing.dwba       = filtfilt(wing_carry.b, wing_carry.a, wing.dwba);
        wing.dwba_vel   = diff(wing.dwba) * Fs;
        wing.dwba_vel   = [wing.dwba_vel(1) ; wing.dwba_vel];
        wing.left_filt	= filtfilt(wing.b, wing.a, wing.left);
        wing.right_filt	= filtfilt(wing.b, wing.a, wing.right);
        wing.dwba_filt	= wing.left_filt - wing.right_filt;
        
      	% Extract wing saccades
        wing_saccade = saccade_all(wing.dwba_filt, tintrp, wing.thresh, wing.true_thresh, wing.Fc_detect, ...
                                wing.Fc_ss, wing.amp_cut, direction, wing.pks, wing.sacd_length, ...
                                wing.min_pkdist, wing.min_pkwidth, wing.min_pkprom, ...
                                wing.min_pkthresh, wing.boundThresh, wing.showplot);
     	wing_saccade.extra.dwba = wing.dwba; % carry unfiltered dwba singal
        wing_saccade.extra.dwba_vel = wing.dwba_vel; % carry unfiltered dwba velocity singal
        
        if wing_saccade.count == 0
            rep = 1;
        else
            rep = wing_saccade.count;
        end
        VTable = table(D.wave(kk),'VariableNames',{'Vel'});
        ITable = [I(kk,:),VTable];
        ITable = repmat(ITable,rep,1);
        WING_SACCADE_STATS = [WING_SACCADE_STATS ; [ITable , wing_saccade.SACD]];
        
        SACCADE{kk,6} = {wing_saccade};
        
        if wing.showplot
            pause
            close all
        end
        
        if head_saccade.count > 0
            head2wing = head_wing_saccade_cc(head_saccade, wing_saccade, 0.15, 0.5, 0.5, false, false);
            %head2wing_align_wing = head_wing_saccade_cc(head_saccade, wing_saccade, 0.15, 0.5, 0.5, true, false);
            %pause
        	SACCADE{kk,7} = {head2wing};
            %SACCADE{kk,8} = {head2wing_align_wing};
        end
    end
end

% Fill in empty saccade trials
empty_idx = cellfun(@(x) isempty(x), HEAD_DATA);
HEAD_DATA(empty_idx) = {saccade_all(0*head.pos, tintrp, head.thresh, head.true_thresh, head.Fc_detect, ...
                                head.Fc_ss, head.amp_cut, direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, false)};

%% Extract & group saccades & intervals by speed & by fly
fields = {'normpeak_saccade','norm_interval','normstart_interval','normend_interval'};
nfield = length(fields);
norm_fields = {'time','position','velocity'};
norm_fields_stats = string(norm_fields) + "_stats";

center = 0; % normalization center for saccades & inter-saccade intervals
dim = 1; % dimension to center
dim_stats = 2; % dimension for statisticsHEAD_SACCADE_STATS

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
    for jj = 1:N.wave
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
    for jj = 1:N.wave
        for ff = 1:length(norm_fields)
            GRAND.(fields{kk})(jj).(norm_fields_stats(ff)) = basic_stats( ...
                GRAND.(fields{kk})(jj).(norm_fields{ff}), dim_stats);
        end
    end
end
                        
%% SAVE %%
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','HEAD_DATA','SACCADE','HEAD_SACCADE_STATS','WING_SACCADE_STATS','FLY','GRAND',...
      'D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end