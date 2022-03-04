function [] = Experiment_head_electro_magnet_control(Fn)
%% Experiment_head_electro_magnet_control: runs a experiment using the LED arena and fly panel
% Fn is the fly number
daqreset
imaqreset
% Fn = 0;
%% Set directories & experimental parameters
root = 'C:\BC\Rigid_data\Experiment_head_electro_magnet_background';

%% EXPERIMENTAL PARAMETERS
n_tracktime = 20 + 1;     	% length(func)/fps; seconds for each EXPERIMENT
n_pause = 0.2;              % seconds for each pause between panel commands
n_rep = 20;                 % # of repetitions
patID = 3;                  % pattern ID
yPos = 5;                   % spatial frequency
FPS = 100;                  % camera frame rate
nFrame = FPS*n_tracktime;   % # of frames to log
Fs = 5000;                  % DAQ sampling rate [Hz]
AI = 1:7;                	% Analog input channels
AO = 0:1;                  	% Analog output channels
Gain = 800;                 % Camera gain

%% Set up data acquisition on MCC (session mode)
% DAQ Setup
[s,~] = NI_USB_6212_head_magnet(Fs,AI,AO);

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

on_time = 2; % magnet ON for this amount of time [s]
off_time = 2; % magnet OFF for this amount of time [s]
seq_time = on_time + off_time; % time for one cycle [s]
n_seq = (n_tracktime + off) ./ seq_time; % number of sequences
pad_end = round(seq_time*s.Rate*( n_seq - floor(n_seq)));

on_samp = on_time * s.Rate; % ON time in samples
off_samp = off_time * s.Rate; % OFF time in samples

seq = [0*ones(off_samp,1) ; 5*ones(on_samp,1)]; % one on-off sequence
seq_all = repmat(seq, [floor(n_seq) 1]); % sequences for length of experiment
seq_all = [seq_all ; seq_all(end)*ones(pad_end,1)]; % pad end in case sequences don't match with experiment time
seq_all = [seq_all ; seq_all(end)];% add one last point to match TRIGGER signal
seq_all(end-50:end) = 0; % stay off at end

AO_all = [TRIG, seq_all]; % all alanlog output signals

cntrl = repmat([1 0]', [n_rep 1]);

%% EXPERIMENT LOOP
disp('Start Experiment:')
for ii = 11:20
    disp(['Trial: ' num2str(ii) '     Control: ' num2str(cntrl(ii))])
    preview(vid) % open video preview window
    
    Panel_com('stop')
    
    % Set AO trigger to 0
 	queueOutputData(s,zeros(5000,2)) % set trigger AO signal to 0
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
    Panel_com('send_gain_bias', [0 0 0 0])
    pause(n_pause)
    Panel_com('set_mode', [0,0]); % 0=open,1=closed,2=fgen,3=vmode,4=pmode
	    
    % START EXPERIMENT & DATA COLLECTION
    start(vid) % start video buffer
    AO_temp = AO_all;
    AO_temp(:,2) = cntrl(ii)*AO_temp(:,2);
    queueOutputData(s, AO_temp) % set AO signals
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
    fname = ['fly_' num2str(Fn) '_trial_' num2str(ii) '_control_' num2str(cntrl(ii)) '.mat'];
    save(fullfile(root,fname),'-v7.3','data','t_p','vidData','t_v', 'off_time', 'on_time');
    Panel_com('stop')
end

delete(vid)
disp('Done');
daqreset
imaqreset
end