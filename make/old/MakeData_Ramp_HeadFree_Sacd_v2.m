function [] = MakeData_Ramp_HeadFree_Sacd_v2(wave)
%% MakeData_Ramp_HeadFree_Sacd:
%   INPUTS:
%       wave    :   spatial wavelength of data
%
%   OUTPUTS:
%       -
%

wave = 30;

% Data location
rootdir = ['H:\EXPERIMENTS\RIGID\Experiment_Asymmetry_Control_Verification\HighContrast\' num2str(wave)];

% Output file name
filename = ['Ramp_HeadFree_SACCD_Anti_Wave=' num2str(wave)];

% Setup Directories 
PATH.daq  = rootdir; % DAQ data location
PATH.vid  = fullfile(PATH.daq,'Vid'); % video data location
PATH.head = fullfile(PATH.vid,'tracked_head'); % tracked kinematic data location
PATH.wing = fullfile(PATH.vid,'wing_filt', 'tracked_head_wing'); % tracked kinematic data location

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.wing,'*.csv',false,'fly','trial','vel','wave');

%% Get Data %%
clc
disp('Loading...')
Fs = 200; % sampling frequency [s]
tintrp = (0:(1/Fs):(10 - 1/Fs))'; % time vector for interpolation

% HEAD saccade detection parameters
head.showplot = false;
head.Fc_detect = [8 0.5];
head.Fc_ss = [nan nan];
head.amp_cut = 4;
head.thresh = 100;
head.sacd_length = nan;
head.pks = [];
head.min_pkdist = 0.05;
head.min_pkwidth = 0.01;
head.min_pkprom = 20;
head.min_pkthresh = 0;
head.boundThresh = 0.25;

% WING saccade detection parameters
wing.showplot = false;
wing.Fc_detect = [4 0.5];
wing.Fc_ss = [4 nan];
wing.amp_cut = 7;
wing.thresh = 25;
wing.sacd_length = nan;
wing.pks = [];
wing.min_pkdist = 0.5;
wing.min_pkwidth = 0.1;
wing.min_pkprom = 20;
wing.min_pkthresh = 0;
wing.boundThresh = 0.35;
wing.Fc = 8;
[wing.b, wing.a] = butter(3, wing.Fc / (Fs/2) ,'low');

Vel = U.vel{1}; % velocities
Stim = (Vel*tintrp')'; % stimuli
SACCADE = [I , splitvars(table(num2cell(zeros(N.file,3))))]; % store saccade objects
SACCADE.Properties.VariableNames(5:7) = {'head_saccade','wing_saccade','head2wing'};
HEAD_DATA = cell(N.fly,N.vel);
HEAD_SACCADE_STATS = []; % store saccade stats
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    load(fullfile(PATH.head, [basename{kk} '.mat']),'hAngles'); % load head angles % time arrays
 
  	% Sync frames and get pattern data
	[trig,~] = sync_pattern_trigger(t_p, data(:,2), 10, data(:,1), true, 1, true, false);
    
    % Get head data
    head.pos = interp1(trig.time_sync, hAngles, tintrp, 'pchip');
    
	% Head with same filter as wings
 	head.pos_filt = filtfilt(wing.b, wing.a, head.pos);
       
    % Extract head saccades
    direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    head_saccade = saccade_all(head.pos, tintrp, head.thresh, head.Fc_detect, head.Fc_ss, head.amp_cut, ...
                                direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, head.showplot);
    head_saccade = stimSaccade(head_saccade, Stim(:,I.vel(kk)), false); % with approximate pattern position
    SACCADE{kk,5} = {head_saccade}; % store data in cell
    
    if head_saccade.count == 0
        rep = 1;
    else
        HEAD_DATA{I.fly(kk),I.vel(kk)}(end+1,1) = head_saccade;
        rep = head_saccade.count;
    end
    VTable = table(D.vel(kk),'VariableNames',{'Vel'});
    ITable = [I(kk,:),VTable];
    ITable = repmat(ITable,rep,1);
    HEAD_SACCADE_STATS = [HEAD_SACCADE_STATS ; [ITable , head_saccade.SACD]];
    
    if head.showplot
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
        
        wing.left_filt  = filtfilt(wing.b, wing.a, wing.left);
        wing.right_filt = filtfilt(wing.b, wing.a, wing.right);
        wing.dwba_filt 	= wing.left_filt - wing.right_filt;
        
      	% Extract wing saccades
        wing_saccade = saccade_all(wing.dwba_filt, tintrp, wing.thresh, wing.Fc_detect, wing.Fc_ss, ...
                                wing.amp_cut, direction, wing.pks, wing.sacd_length, ...
                                wing.min_pkdist, wing.min_pkwidth, wing.min_pkprom, ...
                                wing.min_pkthresh, wing.boundThresh, wing.showplot);
     	wing_saccade.extra.dwba = wing.dwba; % carry unfiltered dwba singal
        SACCADE{kk,6} = {wing_saccade};
        
        if wing.showplot
            pause
            close all
        end
        
        if head_saccade.count > 0
            head2wing = head_wing_saccade_cc(head_saccade, wing_saccade, 0.15, 0.03, 0.5, false);
        	SACCADE{kk,7} = {head2wing};
        end 
    end
end

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

%%
keepI = cellfun(@(x) isstruct(x) | isobject(x), SACCADE.head2wing);
SACCADE = SACCADE(keepI,:);
clear head wing
head.fly = [];
wing.fly = [];
head.grand = [];
wing.grand = [];
for v = 1:N.vel
    velI = v == SACCADE.vel;
    for a = 1:N.fly
        flyI = a == SACCADE.fly;
        sI = velI & flyI;
        hsacd = cellfun(@(x) x.hsacd, SACCADE.head2wing(sI,:), 'UniformOutput', true);
        wsacd = cellfun(@(x) x.wsacd, SACCADE.head2wing(sI,:), 'UniformOutput', true);
        get_fields = string(fieldnames(hsacd));
        get_fields = get_fields(1:end-1);
        stats_fields = get_fields + "_stats";
        for f = 1:length(get_fields)
           head.fly(a,v).(get_fields(f)) = cat(2, hsacd.(get_fields(f)));
           wing.fly(a,v).(get_fields(f)) = cat(2, wsacd.(get_fields(f)));
           head.fly(a,v).(stats_fields(f)) = basic_stats(head.fly(a,v).(get_fields(f)),2);
           wing.fly(a,v).(stats_fields(f)) = basic_stats(wing.fly(a,v).(get_fields(f)),2);
        end
    end
end

for v = 1:N.vel
    for f = 1:length(stats_fields)
        temp = cat(2,head.fly(:,v).(stats_fields(f)));
        head.grand(v).(get_fields(f)) = cat(2,temp.mean);
        head.grand(v).(stats_fields(f)) = basic_stats(head.grand(v).(get_fields(f)),2);
        
        temp = cat(2,wing.fly(:,v).(stats_fields(f)));
        wing.grand(v).(get_fields(f)) = cat(2,temp.mean);
        syncI = round(1000*nanmean(wing.grand(v).('time'),2)) == round(1000*-0.2);
      	for k = 1:size(wing.grand(v).(get_fields(f)),2)
            wing.grand(v).(get_fields(f))(:,k) = wing.grand(v).(get_fields(f))(:,k) - ...
                wing.grand(v).(get_fields(f))(syncI,k);
        end
        wing.grand(v).(stats_fields(f)) = basic_stats(wing.grand(v).(get_fields(f)),2);
    end
end

%%

TD = cellfun(@(x) x.TimeDiff, SACCADE.head2wing, 'UniformOutput', false);
TD = cat(1,TD{:});
TD = TD(~isnan(TD(:,1)),:);

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = 1*[2 2 1.5 2];
movegui(FIG,'center')
ww = 1;
ax(ww) = subplot(1,1,ww); axis tight
CC = parula(size(TD,2));
bx = boxplot(TD, 'Labels', {'Start','Peak','End'}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
ylabel('Time Difference (ms)')

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1},h{kk,2},CC(kk,:));
end

set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(ww).Children = ax(ww).Children([end 1:end-1]);
ax(ww).YLim(1) = -300;
ax(ww).YLim(2) = 300;

set(ax,'LineWidth', 1, 'Box', 'on')
set(ax(ww), 'YTick', -300:100:300)

%% SAVE %%
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','HEAD_DATA','SACCADE','HEAD_SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end