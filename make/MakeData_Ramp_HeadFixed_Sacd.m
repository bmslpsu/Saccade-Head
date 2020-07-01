function [] = MakeData_Ramp_HeadFixed_Sacd(Fc)
%% MakeData_Ramp_HeadFree_Sacd:
%   INPUTS:
%       Fc      :   head data low-pass filter cutoff frequency [Hz]
%
%   OUTPUTS:
%       -
%

Fc = 20;

% Data location
rootdir = 'H:\EXPERIMENTS\RIGID\Experiment_Ramp_30_HeadFixed';

% Output file name
filename = ['Ramp_HeadFixed_SACCD_Anti_filt=' num2str(Fc)];

% Setup Directories 
PATH.daq  = rootdir; % DAQ data location
PATH.wing = fullfile(PATH.daq,'wing_filt','tracked_wing'); % tracked kinematic data location

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.wing,'*.csv',false,'fly','trial','vel');

%% Get Data %%
disp('Loading...')
showplot = true;
Fs = 100;

pat_total_time = 10;
tintrp = (0:(1/Fs):(pat_total_time - 1/Fs))';
reg = true; % use interp times
start_idx = [];
add1 = false; % add 1st frame
[b_high,a_high] = butter(3,0.5/(Fs/2), 'high');
[wing.b,wing.a] = butter(2,Fc/(Fs/2),'low');
Vel = U.vel{1};
Stim = (Vel*tintrp')';
SACCADE = [I , splitvars(table(num2cell(zeros(N.file,2))))]; % store saccade objects
SACCADE.Properties.VariableNames(4:5) = {'head_saccade','dWBA'}; 
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    benifly = ImportBenifly(fullfile(PATH.wing, [basename{kk} '.csv'])); % load head & wing angles
    
    % Sync video with trigger & pattern  
    [trig,~] = sync_pattern_trigger(t_p, data(:,2), pat_total_time, ... 
                                data(:,1), reg, start_idx, add1, false);
    
  	% Get wing data
    wing.left  = rad2deg(interp1(trig.time_sync, benifly.LWing, tintrp, 'pchip'));
    wing.right = rad2deg(interp1(trig.time_sync, benifly.RWing, tintrp, 'pchip'));
    wing.left  = hampel(tintrp, wing.left);
    wing.right = hampel(tintrp, wing.right);
    wing.left  = filtfilt(wing.b, wing.a, wing.left);
    wing.right = filtfilt(wing.b, wing.a, wing.right);
    wing.wba = wing.left - wing.right;
    wing.wba = filtfilt(b_high, a_high, wing.wba);
    
    fig = figure (1) ; cla ; hold on ; title(num2str(D.vel(kk)))
    xlabel('Time (s)')
    ylabel('WBA')
    set(fig, 'Color', 'w')
%     plot(tintrp, wing.left, 'Color', [0.7 0 0], 'LineWidth', 1)
%     plot(tintrp, wing.right, 'Color', [0 0 0.7], 'LineWidth', 1)
    plot(tintrp, wing.wba, 'Color', [0 0 0], 'LineWidth', 1)
    ylim(20*[-1 1])
    pause
    
    % Get Saccade Stats
%     SACCADE{kk,6} = {wing.wba};
%     VTable = table(D.vel(kk),'VariableNames',{'Vel'});
%     ITable = [I(kk,:),VTable];
%     ITable = repmat(ITable,rep,1);
%     SACCADE_STATS = [SACCADE_STATS ; [ITable , head_saccade.SACD]];
    
%     if showplot
%         pause
%         close all
%     end
end

% Fill in empty saccade trials
% empty_idx = cellfun(@(x) isempty(x), ALL_DATA);
% ALL_DATA(empty_idx) = {saccade(nan*head.pos, nan*head.pos, 300, 0, 0, [], [], false)};

%% Extract & group saccades & intervals by speed & by fly
fields = {'normpeak_saccade','norm_interval','normstart_interval','normend_interval',...
    'normstart_stimulus','error','int_error'};
nfield = length(fields);
norm_fields = {'time','position','velocity'};

center = 0; % normalization center for saccades & inter-saccade intervals
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
        GRAND.(fields{kk})(:,jj) = struct_center(FLY.(fields{kk})(:,jj), center, even , dim, norm_fields);
    end
end

%% SAVE %%
disp('Saving...')
save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','ALL_DATA','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end