function [] = MakeData_Ramp_HeadFree_Sacd(wave,Fc)
%% MakeData_Ramp_HeadFree_Sacd:
%   INPUTS:
%       wave    :   spatial wavelength of data
%       Fc      :   head data low-pass filter cutoff frequency [Hz]
%
%   OUTPUTS:
%       -
%

wave = 30;
Fc = 40;

% Data location
rootdir = ['H:\EXPERIMENTS\RIGID\Experiment_Asymmetry_Control_Verification\HighContrast\' num2str(wave)];

% Output file name
filename = ['NewRamp_HeadFree_SACCD_Anti_filt=' num2str(Fc) '_Wave=' num2str(wave)];

% Setup Directories 
PATH.daq = rootdir; % DAQ data location
PATH.vid = fullfile(PATH.daq,'\Vid'); % video data location
% PATH.ang = fullfile(PATH.daq,'\Vid\tracked'); % tracked kinematic data location
PATH.ang = fullfile(PATH.daq,'\Vid\tracked_head'); % tracked kinematic data location

% Select files
[D,I,N,U,T,FILES,~,basename] = GetFileData(PATH.ang,'*.mat',false,'fly','trial','vel','wave');

%% Get Data %%
disp('Loading...')
showplot = false;
tintp = (0:(1/200):(10 - 1/200))';
Vel = U.vel{1};
Stim = (Vel*tintp')';
SACCADE = [I , table(num2cell(zeros(N.file,1)))]; % store saccade objects
SACCADE.Properties.VariableNames{5} = 'saccade'; 
ALL_DATA = cell(N.fly,N.vel);
COUNT = cell(N.fly,N.vel);
SACCADE_STATS = []; % store saccade stats
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    load(fullfile(PATH.ang, [basename{kk} '.mat']),'hAngles'); % load head angles % time arrays
    % benifly = ImportBenifly(fullfile(PATH.ang, FILES{kk})); % load head angles
    
    % Sync video with trigger & pattern
    trig.raw_time   = t_p; % DAQ raw times for trigger
    trig.pos        = round(data(:,1)); % trigger values
    trig.diff       = diff(trig.pos); % trigger derivative (rising edge triggers frame)
    [~,trig.locs]   = findpeaks(trig.diff); % where each frame starts
    trig.time       = [0 ; trig.raw_time(trig.locs+1)]; % where each frame starts
   	
    % Get head data
    % benifly.Head(1) = benifly.Head(2);
    Head = process_signal(trig.time, hAngles, Fc, tintp, [4 8 16 32 64]);
    
    % Get Saccade Stats
    peaks = [];
    direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    head_saccade = saccade(Head.X(:,1), Head.Time, 350, direction, peaks, showplot);
    head_saccade = stimSaccade(head_saccade, Stim(:,I.vel(kk)), false);
    % figure (1) ; title(num2str(D.vel(kk)))
    
    COUNT{I.fly(kk),I.vel(kk)}(end+1,1) = head_saccade.count;
    
%     if head_saccade.count < 3
%      	temp = saccade(Head.X(:,1), Head.Time, 350, direction, peaks, true);
%         pause
%         close all
%     end
    
    if head_saccade.count==0
        rep = 1;
        SACCADE{kk,5} = {head_saccade};
        % ALL_DATA{I.fly(kk),I.vel(kk)}(end+1,1) = saccade();
    else
        SACCADE{kk,5} = {head_saccade};
        ALL_DATA{I.fly(kk),I.vel(kk)}(end+1,1) = head_saccade;
        rep = head_saccade.count;
    end
    VTable = table(D.vel(kk),'VariableNames',{'Vel'});
    ITable = [I(kk,:),VTable];
    ITable = repmat(ITable,rep,1);
    SACCADE_STATS = [SACCADE_STATS ; [ITable , head_saccade.SACD]];
    
    if showplot
        pause
        close all
    end
end

%% Extract & group saccades & intervals by speed & by fly
fields = {'normpeak_saccade','norm_interval','normstart_interval','normend_interval',...
    'normstart_stimulus','error','int_error'};
nfield = length(fields);
norm_fields = {'time','position','velocity'};

center = 0; % normalization cenetr for saccades & inter-saccade intervals
dim = 1; % dimension to center

FLY = []; % all trials per speed per fly
GRAND = []; % all trials per speed
for kk = 1:nfield % for each field in the saccade structure
    if kk < 2 % even padding around center only for saccades (not intervals)
        even = true;
    else
        even = false;
    end
    
    FLY.(fields{kk}) = cellfun(@(x) struct_center([x.(fields{kk})], center, even, dim, norm_fields), ...
                            ALL_DATA, 'UniformOutput', true);
    for jj = 1:N.vel
        GRAND.(fields{kk})(:,jj) = struct_center(FLY.(fields{kk})(:,jj),center, even , dim, norm_fields);
    end
end

%% SAVE %%
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end