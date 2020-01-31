function [] = MakeData_Ramp_HeadFree_Sacd(wave,match,Fc)
%% MakeData_Ramp_HeadFree_Sacd:
%   INPUTS:
%       wave    :   spatial wavelength of data
%       match   :   what saccades direction to 
%       Fc      :   head data low-pass filter cutoff frequency [Hz]
%
%   OUTPUTS:
%       -
%
wave = 30;
match = 0;
Fc = 40;

% Data location
rootdir = ['H:\EXPERIMENTS\RIGID\Experiment_Asymmetry_Control_Verification\HighContrast\' num2str(wave)];

% What saccades we will get
if      match==0
    clss = 'All';
elseif  match==1
    clss = 'CO';
elseif  match==-1
    clss = 'Anti';
elseif  match==2
    clss = 'Positive';
elseif  match==-2
    clss = 'Negative';
else
    error('Invalid match condition')
end

% Output file name
filename = ['NewRamp_HeadFree_SACCD_' clss '_filt=' num2str(Fc) '_Wave=' num2str(wave)];

% Setup Directories 
PATH.daq = rootdir; % DAQ data location
PATH.vid = fullfile(PATH.daq,'\Vid'); % video data location
PATH.ang = fullfile(PATH.daq,'\Vid\tracked'); % tracked kinematic data location

% Select files
[D,I,N,U,T,FILES,~,basename] = GetFileData(PATH.ang,'*.csv',false,'fly','trial','vel','wave');

%% Get Data %%
disp('Loading...')
showplot = false;
tintp = (0:(1/200):(10 - 1/200))';
Vel = U.vel{1};
Stim = (Vel*tintp')';
SACCADE = [I , table(num2cell(zeros(N.file,1)))]; % store saccade objects
SACCADE.Properties.VariableNames{5} = 'saccade';
ALL_DATA = cell(N.vel/2,2);
for jj = 1:N.vel
    ALL_DATA{jj} = cell(N.fly,1);
end
SACCADE_STATS = []; % store saccade stats
for kk = 1:N.file
    disp(kk)   
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    % load(fullfile(PATH.vid, [basename{kk} '.mat']),'t_v'); % load head angles % time arrays
    benifly = ImportBenifly(fullfile(PATH.ang, FILES{kk}));
    
    % Sync video with trigger & pattern
    trig.raw_time   = t_p; % DAQ raw times for trigger
    trig.pos        = round(data(:,1)); % trigger values
    trig.diff       = diff(trig.pos); % trigger derivative (rising edge triggers frame)
    [~,trig.locs]   = findpeaks(trig.diff); % where each frame starts
    trig.time       = [0;trig.raw_time(trig.locs+1)]; % where each frame starts
   	
    % Get head data
    benifly.Head(1) = benifly.Head(2);
    Head = process_signal(trig.time, rad2deg(benifly.Head), Fc, tintp, [4 8 16 32 64]);
    
    % Get Saccade Stats
    peaks = [];
    direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    head_saccade = saccade(Head.X(:,1), Head.Time, 3.5, direction, peaks,  showplot);
    head_saccade = stimSaccade(head_saccade, Stim(:,I.vel(kk)), showplot);
    
    if isnan(head_saccade.count)
        rep = 1;
        SACCADE{kk,5} = {[]};
        % ALL_DATA{I.vel(kk)}{I.fly(kk)}(end+1,1) = {[]};
    else
        SACCADE{kk,5} = {head_saccade};
        ALL_DATA{I.vel(kk)}{I.fly(kk)}(end+1,1) = {head_saccade};
        rep = head_saccade.count;
    end
    ITable = repmat(I(kk,:),rep,1);
    SACCADE_STATS = [SACCADE_STATS ; [ITable , head_saccade.SACD]];
end

%%
% Extract & normalize saccades
clear NORM
NORM.saccade_all  = cell(N.vel,1);
NORM.saccade_fly_all = cell(N.vel,1);
for jj = 1:N.vel
    for kk = 1:N.fly
        % Get all saccades for all trials for each fly & calculate fly
        % statistics
        time = cellfun(@(x) x.normpeak_saccade.time, ALL_DATA{jj}{kk}, 'UniformOutput', false);
        [time,~,~,~,dR] = nancat_center(time, 0, 1, [], true);
        
    	position = cellfun(@(x,y) padmat(x.normpeak_saccade.position,y,nan,1), ALL_DATA{jj}{kk}, ...
            dR, 'UniformOutput', false);
        position = cat(2,position{:});
        
    	velocity = cellfun(@(x,y) padmat(x.normpeak_saccade.velocity,y,nan,1), ALL_DATA{jj}{kk}, ...
            dR, 'UniformOutput', false);
        velocity = cat(2,velocity{:});
        
        NORM.saccade_all{jj}(kk).time = time;
        NORM.saccade_all{jj}(kk).position = position;
        NORM.saccade_all{jj}(kk).velocity = velocity;
        
        NORM.saccade_fly_all{jj}(kk).time = basic_stats(time,2);
        NORM.saccade_fly_all{jj}(kk).position = basic_stats(position,2);
        NORM.saccade_fly_all{jj}(kk).velocity = basic_stats(velocity,2);
    end
    % Get median fly saccades & calculate stats for each speed
    time = [NORM.saccade_fly_all{jj}.time];
    position = [NORM.saccade_fly_all{jj}.position];
    velocity = [NORM.saccade_fly_all{jj}.velocity];
    
	[time,~,~,~,dR] = nancat_center({time.median}',0, 1, [], true);
    position = cellfun(@(x,y)  padmat(x,y,nan,1), {position.median}', dR, 'UniformOutput', false);
    velocity = cellfun(@(x,y)  padmat(x,y,nan,1), {velocity.median}', dR, 'UniformOutput', false);
    position = cat(2,position{:});
    velocity = cat(2,velocity{:});
    
    NORM.saccade_fly_median(jj).time = time;
    NORM.saccade_fly_median(jj).position = position;
    NORM.saccade_fly_median(jj).velocity = velocity;
    
    NORM.saccade_vel(jj).time = basic_stats(time,2);
    NORM.saccade_vel(jj).position = basic_stats(position,2);
	NORM.saccade_vel(jj).velocity = basic_stats(velocity,2);
end
% Get median saccade by speed & calculate stats for each direction (CW,CCW)
time = [NORM.saccade_vel.time];
position = [NORM.saccade_vel.position];
velocity = [NORM.saccade_vel.velocity];

[time,~,~,~,dR] = nancat_center({time.median}',0, 1, [], true);

position = cellfun(@(x,y)  padmat(x,y,nan,1), {position.median}', dR, 'UniformOutput', false);
velocity = cellfun(@(x,y)  padmat(x,y,nan,1), {velocity.median}', dR, 'UniformOutput', false);
position = cat(2,position{:});
velocity = cat(2,velocity{:});

NORM.saccade_vel_median.time = time;
NORM.saccade_vel_median.position = position;
NORM.saccade_vel_median.velocity = velocity;

NORM.saccade_ALL(1).time = basic_stats(time(1:N.vel/2),2);
NORM.saccade_ALL(1).position = basic_stats(position(1:N.vel/2),2);
NORM.saccade_ALL(1).velocity = basic_stats(velocity(1:N.vel/2),2);

NORM.saccade_ALL(2).time = basic_stats(time(N.vel/2+1:end),2);
NORM.saccade_ALL(2).position = basic_stats(position(N.vel/2+1:end),2);
NORM.saccade_ALL(2).velocity = basic_stats(velocity(N.vel/2+1:end),2);

% Extract & normalize inter-saccade intervals
NORM.interval_all  = cell(N.vel,1);
NORM.interval_fly_all = cell(N.vel,1);
for jj = 1:N.vel
    for kk = 1:N.fly
        % Get all saccades for all trials for each fly & calculate fly
        % statistics
        time = cellfun(@(x) x.normstart_interval.time, ALL_DATA{jj}{kk}, 'UniformOutput', false);
        [time,~,~,~,dR] = nancat_center(time, 0, 1, [], false);
        
    	position = cellfun(@(x,y) padmat(x.normstart_interval.position,y,nan,1), ALL_DATA{jj}{kk}, ...
            dR, 'UniformOutput', false);
        position = cat(2,position{:});
        
    	velocity = cellfun(@(x,y) padmat(x.normstart_interval.velocity,y,nan,1), ALL_DATA{jj}{kk}, ...
            dR, 'UniformOutput', false);
        velocity = cat(2,velocity{:});
        
        NORM.interval_all{jj}(kk).time = time;
        NORM.interval_all{jj}(kk).position = position;
        NORM.interval_all{jj}(kk).velocity = velocity;
        
        NORM.interval_fly_all{jj}(kk).time = basic_stats(time,2);
        NORM.interval_fly_all{jj}(kk).position = basic_stats(position,2);
        NORM.interval_fly_all{jj}(kk).velocity = basic_stats(velocity,2);
    end
    % Get median fly saccades & calculate stats for each speed
    time = [NORM.interval_fly_all{jj}.time];
    position = [NORM.interval_fly_all{jj}.position];
    velocity = [NORM.interval_fly_all{jj}.velocity];
    
	[time,~,~,~,dR] = nancat_center({time.median}',0, 1, [], false);
    position = cellfun(@(x,y)  padmat(x,y,nan,1), {position.median}', dR, 'UniformOutput', false);
    velocity = cellfun(@(x,y)  padmat(x,y,nan,1), {velocity.median}', dR, 'UniformOutput', false);
    position = cat(2,position{:});
    velocity = cat(2,velocity{:});
    
    NORM.interval_fly_median(jj).time = time;
    NORM.interval_fly_median(jj).position = position;
    NORM.interval_fly_median(jj).velocity = velocity;
    
    NORM.interval_vel(jj).time = basic_stats(time,2);
    NORM.interval_vel(jj).position = basic_stats(position,2);
	NORM.interval_vel(jj).velocity = basic_stats(velocity,2);
end
% Get median saccade by speed & calculate stats for each direction (CW,CCW)
time = [NORM.interval_vel.time];
position = [NORM.interval_vel.position];
velocity = [NORM.interval_vel.velocity];

[time,~,~,~,dR] = nancat_center({time.median}',0, 1, [], false);

position = cellfun(@(x,y)  padmat(x,y,nan,1), {position.median}', dR, 'UniformOutput', false);
velocity = cellfun(@(x,y)  padmat(x,y,nan,1), {velocity.median}', dR, 'UniformOutput', false);
position = cat(2,position{:});
velocity = cat(2,velocity{:});

NORM.interval_vel_median.time = time;
NORM.interval_vel_median.position = position;
NORM.interval_vel_median.velocity = velocity;

NORM.interval_ALL(1).time = basic_stats(time(1:N.vel/2),2);
NORM.interval_ALL(1).position = basic_stats(position(1:N.vel/2),2);
NORM.interval_ALL(1).velocity = basic_stats(velocity(1:N.vel/2),2);

NORM.interval_ALL(2).time = basic_stats(time(N.vel/2+1:end),2);
NORM.interval_ALL(2).position = basic_stats(position(N.vel/2+1:end),2);
NORM.interval_ALL(2).velocity = basic_stats(velocity(N.vel/2+1:end),2);


%% SAVE %%
disp('Saving...')
% save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
%     'SACD','SACCADE','INTERVAL','Stim','TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')

save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','SACD','SACCADE','INTERVAL','Stim','D','I','U','N','T','-v7.3')

disp('SAVING DONE')
end