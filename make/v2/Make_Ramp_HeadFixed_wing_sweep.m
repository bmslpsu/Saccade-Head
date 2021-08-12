function [] = Make_Ramp_HeadFixed_wing_sweep(wing,tag)
%% Make_Ramp_HeadFixed_wing_sweep:
%   INPUTS:
%       wing    : wing saccade detetcion parameters
%       tag     : last part of file name
%
%   OUTPUTS:
%       -
%
warning('off', 'signal:findpeaks:largeMinPeakHeight')

% Data location
rootdir = 'E:\EXPERIMENTS\RIGID\Experiment_Ramp_30_HeadFixed';

% Output file name
filename = 'Ramp_HeadFixed_wing';

% Setup Directories 
PATH.daq  = rootdir;
PATH.vid  = fullfile(PATH.daq,''); %
PATH.wing = fullfile(PATH.vid,'wing_filt', 'tracked_wing_new');

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.wing,'*.csv',false,'fly','trial','vel');

%% Get Data
clc
disp('Loading...')
Fs = 100; % sampling frequency [s]
tintrp = (0:(1/Fs):(10 - 1/Fs))'; % time vector for interpolation

wing_carry.Fc = 40;
[wing_carry.b, wing_carry.a] = butter(3, wing_carry.Fc / (Fs/2) ,'low');

direction = -1;
Vel = U.vel{1};
Stim = (Vel*tintrp')';
SACCADE = [I , splitvars(table(num2cell(zeros(N.file,1))))]; % store saccade objects
SACCADE.Properties.VariableNames(4) = {'wing_saccade'};
WING_SACCADE_STATS = [];
for kk = 1:N.file
    disp(kk)
    % Load WING & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    wfile = fullfile(PATH.wing, [basename{kk} '.csv']);
 	benifly = ImportBenifly(wfile); % load wing angles
    
  	% Sync frames and get pattern data
	[trig,~] = sync_pattern_trigger(t_p, data(:,2), 10, data(:,1), true, [], false, false);
    
  	% Get wing data
    n_detrend = 2;
    wing.left       = rad2deg(benifly.LWing);
    wing.right      = rad2deg(benifly.RWing);
    wing.left       = hampel(trig.time_sync, wing.left);
    wing.right      = hampel(trig.time_sync, wing.right);
    wing.left       = detrend(interp1(trig.time, wing.left, tintrp, 'pchip'), n_detrend);
    wing.right      = detrend(interp1(trig.time, wing.right, tintrp, 'pchip'), n_detrend);
    wing.dwba       = wing.left - wing.right;
    wing.dwba       = filtfilt(wing_carry.b, wing_carry.a, wing.dwba);
    wing.dwba_vel   = diff(wing.dwba) * Fs;
    wing.dwba_vel   = [wing.dwba_vel(1) ; wing.dwba_vel];
    wing.left_filt	= filtfilt(wing.b, wing.a, wing.left);
    wing.right_filt	= filtfilt(wing.b, wing.a, wing.right);
    wing.dwba_filt	= wing.left_filt - wing.right_filt;
    
    
    wing.dwba_filt = detrend(wing.dwba_filt,3);

    % Extract wing saccades
    wing_saccade = saccade_all(wing.dwba_filt, tintrp, wing.thresh, wing.true_thresh, wing.Fc_detect, ...
                            wing.Fc_ss, wing.amp_cut, wing.dur_cut, direction, wing.pks, wing.sacd_length, ...
                            wing.min_pkdist, wing.min_pkwidth, wing.min_pkprom, ...
                            wing.min_pkthresh, wing.boundThresh, wing.showplot);
    wing_saccade.extra.dwba = wing.dwba; % carry unfiltered dwba singal
    wing_saccade.extra.dwba_vel = wing.dwba_vel; % carry unfiltered dwba velocity singal

    if wing_saccade.count == 0
        rep = 1;
    else
        rep = wing_saccade.count;
    end
    VTable = table(D.vel(kk),'VariableNames',{'Vel'});
    ITable = [I(kk,:),VTable];
    ITable = repmat(ITable,rep,1);
    WING_SACCADE_STATS = [WING_SACCADE_STATS ; [ITable , wing_saccade.SACD]];

    SACCADE{kk,4} = {wing_saccade};

    if wing.showplot
        figure (1)
        pause
        close all
    end
end

%% SAVE
disp('Saving...')
save(fullfile(root, 'dataset', [filename '_' tag '.mat']),...
      'PATH','SACCADE','WING_SACCADE_STATS',...
      'Stim','D','I','U','N','T','wing','-v7.3')
disp('SAVING DONE')
end