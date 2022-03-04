function [] = Experiment_head_ramp(Fn)
%% Experiment_head_electro_magnet_control: runs a experiment using the LED arena and fly panel
% Fn is the fly number
daqreset
imaqreset
% Fn = 0;
%% Set directories & experimental parameters
root = 'C:\BC\Rigid_data\Experiment_ramp_glue_head';

%% EXPERIMENTAL PARAMETERS
n_tracktime = 10 + 1;     	% length(func)/fps; seconds for each EXPERIMENT
n_pause = 0.2;              % seconds for each pause between panel commands
n_rep = 10;                 % # of repetitions
patID = 3;                  % pattern ID
yPos = 5;                   % spatial frequency
FPS = 200;                  % camera frame rate
nFrame = FPS*n_tracktime;   % # of frames to log
Fs = 5000;                  % DAQ sampling rate [Hz]
AI = 1:6;                	% Analog input channels
AO = 0;                  	% Analog output channels
Gain = 950;                 % Camera gain

speed = [30 -30 60 -60];
arena_gain = speed ./ 3.75;

%% DAQ Setup
[s,~] = NI_USB_6212(Fs,AI,AO);

%% Camera Setup
% [vid,src] = Basler_acA640_750um(nFrame);
[vid,~] = Basler_acA640_120gm(FPS,Gain,nFrame);

%% Camera Trigger Signal & magnet timing
off = 0.1;
t = (0:1/s.Rate:n_tracktime + off)';
TRIG = ((1/2)*(square(2*pi*FPS*t,5) - 1));
TRIG(TRIG==-1) = 4;
end_off = round(Fs*off);
TRIG(end-end_off:end) = 0;

%% Randomized speeds
val = arena_gain;
n_func = length(val); % # of functions
func = (1:n_func)'; % position function indicies

% Create sequence of randomly shuffled functions
func_all = nan(n_func*n_rep,1);
pp = 0;
for kk = 1:n_rep
    func_rand = func(randperm(n_func),:);    % reshuffle randomly
    func_all(pp+1:pp+n_func,1) = func_rand;  % add rep
    pp = kk*n_func;
end
val_all = val(func(func_all));
n_trial = n_rep * n_func;

vel_all = val_all * 3.75;

%% EXPERIMENT LOOP
disp('Start Experiment:')
for ii = 1:n_trial
    disp(['Trial: ' num2str(ii) '     Vel: ' num2str(vel_all(ii))])
    preview(vid) % open video preview window
    
    Panel_com('stop')
    
    % Set AO trigger to 0
 	queueOutputData(s,zeros(5000,1)) % set trigger AO signal to 0
    [~,~] = s.startForeground; % data collection
    
    pause(1) % pause between buffer & experiment
    
    % EXPERIMENT SETUP
    disp('Play Stimulus:')
    Panel_com('set_pattern_id', patID);	% set pattern
    pause(n_pause)
    Panel_com('set_position', [1, yPos]); % set starting position (xpos,ypos)
    pause(n_pause)
	Panel_com('set_funcX_freq', 50); % update rate for x-channel
    pause(n_pause)
    Panel_com('set_funcY_freq', 50); % update rate for y-channel
    pause(n_pause)
    Panel_com('send_gain_bias', [val_all(ii) 0 0 0])
    pause(n_pause)
    Panel_com('set_mode', [0,0]); % 0=open,1=closed,2=fgen,3=vmode,4=pmode
	    
    % START EXPERIMENT & DATA COLLECTION
    start(vid) % start video buffer
    queueOutputData(s, TRIG) % set AO signals
    T = timer('StartDelay',0.5,'TimerFcn',@(src,evt) Panel_com('start'));
    start(T)
    tic
        [data, t_p ] = s.startForeground; % data collection
        stop(vid) % stop video buffer
        Panel_com('stop') % stop stimulus
        [vidData, t_v] = getdata(vid, vid.FramesAcquired); % get video data
    toc
    
    Fs = 1/mean(diff(t_v)); % check FPS of video
  	disp(['Fs = ' num2str(Fs)])
    
    % CL BUFFER
    Arena_CL(2,'X',-10)
    
    % SAVE DATA
    disp('Saving...')
    disp('-----------------------------------------------------------------')
    fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) ...
        '_vel_' num2str(vel_all(ii)) '_wave_' num2str(30) '.mat'];
    save(fullfile(root,fname),'-v7.3','data','t_p','vidData','t_v');
    Panel_com('stop')
end

delete(vid)
disp('Done');
daqreset
imaqreset
end