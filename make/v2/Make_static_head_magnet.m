function [] = Make_static_head_magnet()
%% Make_static_head_magnet:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%
warning('off', 'signal:findpeaks:largeMinPeakHeight')

direction = 0; % get saccades in all directions

% Data location
rootdir = 'E:\EXPERIMENTS\RIGID\Experiment_head_magnet';

% Output file name
filename = 'Magnet_static';

% Setup Directories 
PATH.daq  = rootdir; % DAQ data location
% PATH.head = fullfile(PATH.daq,'tracked_head_edge');
PATH.filt = fullfile(PATH.daq,'filt');
PATH.wing = fullfile(PATH.filt,'tracked_wing');

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.wing,'*.csv',false,'fly','trial');

%% Get Data
clc
close all
disp('Loading...')
Fs = 100; % sampling frequency [s]
tintrp = (0:(1/Fs):(20 - 1/Fs))'; % time vector for interpolation

% HEAD saccade detection parameters
head.showplot = false;
head.Fc_detect = [10 nan];
head.Fc_ss = [nan nan];
head.amp_cut = 4;
head.dur_cut = 0.1;
head.thresh = [60 , 1, 1, 0];
head.true_thresh = 250;
head.sacd_length = nan;
head.pks = [];
head.min_pkdist = 0.5;
head.min_pkwidth = 0.03;
head.min_pkprom = 20;
head.min_pkthresh = 0;
head.boundThresh = [0.2 20];
head_carry.Fc = 40;
[head_carry.b, head_carry.a] = butter(3, head_carry.Fc / (Fs/2) ,'low');

% WING saccade detection parameters
wing.showplot = false;
wing.Fc_detect = [5 nan];
wing.Fc_ss = [5 nan];
wing.amp_cut = 7;
wing.dur_cut = inf;
wing.thresh = [10 , 1, 1.2, 0];
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
wing_carry.Fc = 25;
[wing_carry.b, wing_carry.a] = butter(3, wing_carry.Fc / (Fs/2) ,'low');

SACCADE = [I , splitvars(table(num2cell(zeros(N.file,4))))]; % store saccade objects
SACCADE.Properties.VariableNames(3:end) = {'head_saccade','wing_saccade','wing_on','wing_off'};
for kk = 1:N.file
    disp(kk)
    %basename{kk}
    % Load data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p');
    %load(fullfile(PATH.head, [basename{kk} '.mat']),'hAngles');
    benifly = ImportBenifly(fullfile(PATH.wing, [basename{kk} '.csv']));

  	% Sync frames and get pattern data
	[trig,pat] = sync_pattern_trigger(t_p, data(:,2), 20, data(:,1), true, 1, false, false);
    
    % Get magnet signal
    magnet = round(data(:,7));
    %magnet = magnet(pat.sync);
    magnet = interp1(pat.time_sync, magnet, tintrp, 'nearest');

    % Get head data
    %head.pos = interp1(trig.time_sync, hAngles, tintrp, 'pchip');
    
	% Head with filter
 	%head.pos_filt = filtfilt(head_carry.b, head_carry.a, head.pos);
    
    % Get wing data
    hampel_dx = round(0.05*Fs);
    wing.left       = rad2deg(benifly.LWing);
    wing.right      = rad2deg(benifly.RWing);
    wing.left       = hampel(trig.time_sync, wing.left, hampel_dx);
    wing.right      = hampel(trig.time_sync, wing.right, hampel_dx);
    wing.left       = interp1(trig.time, wing.left, tintrp, 'pchip');
    wing.right      = interp1(trig.time, wing.right, tintrp, 'pchip');

    wing.dwba       = wing.left - wing.right;
    wing.dwba       = filtfilt(wing_carry.b, wing_carry.a, wing.dwba);
    wing.dwba_vel   = central_diff(wing.dwba, 1/Fs);
    wing.left_filt	= filtfilt(wing.b, wing.a, wing.left);
    wing.right_filt	= filtfilt(wing.b, wing.a, wing.right);
    wing.dwba_filt	= wing.left_filt - wing.right_filt;

    % Extract wing saccades
    wing_saccade = saccade_all(wing.dwba_filt, tintrp, wing.thresh, wing.true_thresh, wing.Fc_detect, ...
                            wing.Fc_ss, wing.amp_cut, wing.dur_cut, direction, wing.pks, wing.sacd_length, ...
                            wing.min_pkdist, wing.min_pkwidth, wing.min_pkprom, ...
                            wing.min_pkthresh, wing.boundThresh, wing.showplot);
    wing_saccade.extra.dwba = wing.dwba; % carry unfiltered dwba singal
    wing_saccade.extra.dwba_vel = wing.dwba_vel; % carry unfiltered dwba velocity singal
    wing_saccade.extra.lwing = wing.left; % carry unfiltered left wba singal
    wing_saccade.extra.rwing = wing.right; % carry unfiltered right wba singal
    
%     figure (1) ; cla ; hold on
%     plot(tintrp, wing.dwba, 'k')
%     plot(tintrp, wing.dwba_filt, 'r')
%     plot(tintrp, 10*magnet, 'm')
%     pause
    
    win = logical(magnet);
    test = wing.dwba_filt(win);
    
    [W,~] = pull_windows(wing.dwba, win, 1.99*Fs, Fs, false);
    SACCADE.wing_on{kk} = W{1};
    SACCADE.wing_off{kk} = W{2};

    if wing.showplot
        plot(tintrp, wing.dwba, 'k', 'LineWidth', 0.25)
        plot(tintrp, wing.dwba_filt, 'r')
        plot(tintrp, magnet, 'm')
        pause
        close all
    end
end

%% SAVE %%
disp('Saving...')
save(['E:\DATA\Rigid_Data\Saccade\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','SACCADE','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end