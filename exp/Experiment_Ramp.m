function [] = Experiment_Ramp(Fn)
%% Experiment_Ramp: runs an experiment using LED arena and fly panel controller
%   Written for Panel Controller v3 and NiDAQ seesion mode
%   Set the stimulus to a constant speed
%
%   INPUTS:
%       Fn      :  	fly number
%

%% Set directories & experimental parameters
daqreset
imaqreset
%root = 'C:\BC\Rigid_data\Experiment_Ramp_Head_Fixed';
root = 'C:\BC\Rigid_data\Experiment_Ramp_forRoll';

% EXPERIMENTAL PARAMETERS
speed = 30;                 % stimulus speed [°/s]
ypos = 5;                   % 30° spatial wavelength
n_tracktime = 10 + 1;       % length(func)/fps; seconds for each EXPERIMENT
n_resttime = 1;             % seconds for each REST
n_pause = 0.2;              % seconds for each pause between panel commands
n_rep = 15;                	% number of cycles through spatial frequencies for each fly
FPS = 200;                  % camera frame rate
nFrame = FPS*n_tracktime;   % # of frames to log
Gain = 550;               	% camera gain
Fs = 5000;                  % DAQ sampling rate [Hz]
AI = 1:6;                	% Analog input channels
AO = 0;                     % Analog output channels

%% Set up data acquisition & camera
[s,~] = NI_USB_6212(Fs,AI,AO);

% Camera Trigger Signal
t = 0:1/s.Rate:n_tracktime;
TRIG = ((1/2)*(square(2*pi*FPS*t,50) - 1)');
TRIG(TRIG==-1) = 4;

[vid,~] = Basler_acA640_120gm(FPS,Gain,nFrame);

%% Set variable to control stimlus speed
vel = [speed ; -speed];         % [°/s] CW & CCW speeds
n_vel = length(vel);            % # of velocities
vel_all = repmat(vel,n_rep,1);	% reshuffle randomly
gain_all = vel_all/3.75;        % stimulus arena gain [panel/s]
n_trial = n_rep * n_vel;        % total # trials

%% EXPERIMENT LOOP
disp('Start Experiment:')
disp('--------------------------------------------------')
for kk = 1:n_trial
    disp('Trial')
    disp(num2str(kk));  % print counter to command line
    preview(vid);       % open video preview window
	start(vid)          % start video buffer
    
    % CLOSED LOOP BAR TRACKING %
    Arena_CL(1,'X',-15)
    pause(n_resttime)
    Panel_com('stop')
    
    pause(1) % pause between closed-loop & experiment 
    
    % EXPERIMENT SETUP %
    disp('Play Stimulus: ')
    disp(['Vel = ' num2str(vel_all(kk))])
    
	Panel_com('set_pattern_id', 2);                         % set pattern
        pause(n_pause)
    Panel_com('set_position',[randi([1,96]), ypos])         % set starting position (xpos,ypos)
        pause(n_pause)           
	Panel_com('set_funcX_freq', 50)                         % update rate for x-channel
        pause(n_pause)
    Panel_com('set_funcY_freq', 50)                         % update rate for y-channel
        pause(n_pause)
    Panel_com('set_mode', [0,0])                            % 0=open,1=closed,2=fgen,3=vmode,4=pmode
        pause(n_pause)
	Panel_com('send_gain_bias',[gain_all(kk) 0 0 0]); ...   % [xgain,xoffset,ygain,yoffset]

    % START EXPERIMENT & DATA COLLECTION %
    queueOutputData(s,TRIG) % set trigger AO signal
    T = timer('StartDelay',0.7,'TimerFcn',@(src,evt) Panel_com('start'));
    start(T)
    tic
        [data, t_p ] = s.startForeground; % data collection
        stop(vid) % stop video buffer
        Panel_com('stop')
   	toc
 	[vidData, t_v] = getdata(vid, vid.FramesAcquired); % get video data
    Fs = 1/mean(diff(t_v)); % check FPS of video
  	disp(['Fs = ' num2str(Fs)])
    
    % CLOSED LOOP BAR TRACKING %
    Arena_CL(1,'X',-15)
    
    % SAVE DATA %
    disp('Saving...')
    disp('----------------------------------------------------------------------')
    
    fname = ['Fly_' num2str(Fn) '_Trial_' num2str(kk) '_Vel_' num2str(vel_all(kk)) '.mat'];
    
    save(fullfile(root,fname),'-v7.3','data','t_p','vidData','t_v');
end

delete(vid)
disp('Done');
daqreset
imaqreset
PControl
end