function [] = MakeData_Ramp_HeadFixed_Sacd_v2()
%% MakeData_Ramp_HeadFixed_Sacd_v2:
%   INPUTS:
%       Fc      :   head data low-pass filter cutoff frequency [Hz]
%
%   OUTPUTS:
%       -
%

% Data location
rootdir = 'H:\EXPERIMENTS\RIGID\Experiment_Ramp_30_HeadFixed';

% Output file name
filename = 'Ramp_HeadFixed_SACCD_Anti';

% Setup Directories
PATH.daq  = rootdir; % DAQ data location
PATH.vid  = fullfile(PATH.daq,''); % video data location
PATH.wing = fullfile(PATH.vid,'wing_filt', 'tracked_wing_new'); % tracked kinematic data location

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.wing,'*.csv',false,'fly','trial','vel');

%% Get Data %%
clc
disp('Loading...')
Fs = 100; % sampling frequency [s]
tintrp = (0:(1/Fs):(10 - 1/Fs))'; % time vector for interpolation

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

Vel = U.vel{1}; % velocities
Stim = (Vel*tintrp')'; % stimuli
SACCADE = [I , splitvars(table(num2cell(zeros(N.file,1))))]; % store saccade objects
SACCADE.Properties.VariableNames(4) = {'wing_saccade'};
WING_DATA = cell(N.fly,N.vel);
WING_SACCADE_STATS = []; % store saccade stats
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
 
  	% Sync frames and get pattern data
	[trig,~] = sync_pattern_trigger(t_p, data(:,2), 10, data(:,1), true, [], false, false);
    
	% Load WING data if we have it
	direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    wfile = fullfile(PATH.wing, [basename{kk} '.csv']);
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
    
    if wing.showplot
        pause
        close all
    end
    
%     figure (100) ; cla ; hold on
%     plot(wing.dwba, 'k', 'LineWidth', 0.5)
%     plot(wing.dwba_filt, 'r', 'LineWidth', 1)
%     ylim(50*[-1 1])
%     xlim([1 1000])
%     pause

    % Store data in cells
    SACCADE{kk,4} = {wing_saccade};
    
    if wing_saccade.count==0
        rep = 1;
    else
        WING_DATA{I.fly(kk),I.vel(kk)}(end+1,1) = wing_saccade;
        rep = wing_saccade.count;
    end
    VTable = table(D.vel(kk),'VariableNames',{'Vel'});
    ITable = [I(kk,:),VTable];
    ITable = repmat(ITable,rep,1);
    WING_SACCADE_STATS = [WING_SACCADE_STATS ; [ITable , wing_saccade.SACD]];
end


%% SAVE %%
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','WING_DATA','SACCADE','WING_SACCADE_STATS','Stim','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end