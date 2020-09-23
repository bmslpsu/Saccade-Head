function [] = MakeData_Sine_HeadFree_Sacd(amp)
%% MakeData_Ramp_HeadFree_Sacd:
%   INPUTS:
%       amp         :   amplitude of sineusoid stimulus
%       direction   :   only get saccades in this direction
%
%   OUTPUTS:
%       -
%
warning('off', 'signal:findpeaks:largeMinPeakHeight')

amp = 11.25;
direction = 0; % get saccades in all directions

% Data location
rootdir = ['H:\EXPERIMENTS\RIGID\Experiment_Sinusoid\' num2str(amp)];

% Output file name
filename = ['Sine_HeadFree_SACCD_Amp=' num2str(amp)];

% Setup Directories 
PATH.daq = rootdir;
PATH.vid = fullfile(PATH.daq,'\Vid');
PATH.head = fullfile(PATH.daq,'\Vid\tracked_head');

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.head,'*.mat',false);

%% Get Data %%
clc
close all
disp('Loading...')
Fs = 200; % sampling frequency [s]
tintrp = (0:(1/Fs):(10 - 1/Fs))'; % time vector for interpolation

% HEAD saccade true parameters
head.showplot = false;
head.Fc_detect = [10 nan];
head.Fc_ss = [nan nan];
head.amp_cut = 4;
head.dur_cut = 0.1;
head.thresh = [50 , 1, 1, 70];
head.true_thresh = 300;
head.sacd_length = nan;
head.pks = [];
head.min_pkdist = 0.1;
head.min_pkwidth = 0.03;
head.min_pkprom = 75;
head.min_pkthresh = 0;
head.boundThresh = [0.25 50];
head_carry.Fc = 40;
[head_carry.b, head_carry.a] = butter(3, head_carry.Fc / (Fs/2) ,'low');

AmpT = table(amp*ones(size(I,1),1), 'VariableNames', {'amp'});
SACCADE = [I , AmpT, splitvars(table(num2cell(zeros(N.file,1))))]; % store saccade objects
SACCADE.Properties.VariableNames(5) = {'head_saccade'};
HEAD_DATA = cell(N.fly,N.freq);
HEAD_SACCADE_STATS = []; % store saccade stats
for kk = 1:N.file
    disp(kk)
    disp(basename{kk})
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    load(fullfile(PATH.head, [basename{kk} '.mat']),'hAngles'); % load head angles % time arrays
 
  	% Sync frames and get pattern data
	[trig,~] = sync_pattern_trigger(t_p, data(:,2), 10, data(:,1), true, 1, true, false);
    
  	% Get pattern data
    pat.time = t_p;
    pat.Fs = 1 / mean(diff(pat.time));
    pat.Fc = 2.5*D.freq(kk);
    [pat.b,pat.a] = butter(3, pat.Fc/(pat.Fs/2),'low');
    pat.pos = 3.75*round((96/10)*data(:,2));
    pat.pos = filtfilt(pat.b, pat.a, pat.pos);
    pat.pos = interp1(pat.time, pat.pos, tintrp, 'pchip');
   	pat.pos = pat.pos - mean(pat.pos);
    
    % Get head data
    head.pos = interp1(trig.time_sync, hAngles, tintrp, 'pchip');
    
	% Head with filter
 	head.pos_filt = filtfilt(head_carry.b, head_carry.a, head.pos);
       
    % Extract head saccades
    head_saccade = saccade_all(head.pos_filt, tintrp, head.thresh, head.true_thresh, head.Fc_detect, ...
                                head.Fc_ss, head.amp_cut, head.dur_cut, direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, head.showplot);
    head_saccade = stimSaccade(head_saccade, pat.pos, false); % with approximate pattern position
    SACCADE{kk,5} = {head_saccade};

%     cla ; hold on
%     plot(tintrp, pat.pos, 'k')
%     plot(tintrp, head.pos_filt, 'b')
%     pause

%     if head_saccade.count > 5
%         plotSaccade(head_saccade)
%         plotInterval(head_saccade)
%         pause
%         close all
%     end
        
    if head_saccade.count == 0
        rep = 1;
    else
        HEAD_DATA{I.fly(kk),I.freq(kk)}(end+1,1) = head_saccade;
        rep = head_saccade.count;
    end
    VTable = table(D.freq(kk),'VariableNames',{'Freq'});
    ITable = [I(kk,:),VTable, table(amp,'VariableNames',{'amp'})];
    ITable = repmat(ITable,rep,1);
    HEAD_SACCADE_STATS = [HEAD_SACCADE_STATS ; [ITable , head_saccade.SACD]];
    
    if head.showplot
        figure (1)
        pause
        close all
    end
end

% Fill in empty saccade trials
empty_idx = cellfun(@(x) isempty(x), HEAD_DATA);
HEAD_DATA(empty_idx) = {saccade_all(0*head.pos, tintrp, head.thresh, head.true_thresh, head.Fc_detect, ...
                                head.Fc_ss, head.amp_cut, head.dur_cut, direction, head.pks, head.sacd_length, ...
                                head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, false)};

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
    for jj = 1:N.freq
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
    for jj = 1:N.freq
        for ff = 1:length(norm_fields)
            GRAND.(fields{kk})(jj).(norm_fields_stats(ff)) = basic_stats( ...
                GRAND.(fields{kk})(jj).(norm_fields{ff}), dim_stats);
        end
    end
end

%% SAVE %%
disp('Saving...')
save(['H:\DATA\Rigid_Data\Saccade\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','SACCADE','HEAD_SACCADE_STATS','FLY','GRAND',...
      'D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end