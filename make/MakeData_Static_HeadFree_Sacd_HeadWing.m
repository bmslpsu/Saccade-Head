function [] = MakeData_Static_HeadFree_Sacd_HeadWing(Fc)
%% MakeData_StaticHeadFree_Sacd_HeadWing:
%   INPUTS:
%       wave    :   spatial wavelength of data
%       Fc      :   head data low-pass filter cutoff frequency [Hz]
%
%   OUTPUTS:
%       -
%

Fc = 20;

% Data location
rootdir = 'H:\EXPERIMENTS\RIGID\Experiment_Static_Wave';

% Setup Directories 
PATH.daq  = rootdir; % DAQ data location
PATH.vid  = fullfile(PATH.daq,'wing_filt'); % video data location
PATH.head = fullfile(PATH.daq,'tracked_head'); % tracked kinematic data location
PATH.wing = fullfile(PATH.vid,'tracked_head_wing'); % tracked kinematic data location

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.wing,'*.csv',false,'fly','trial','wave','vel');

%% Get Data %%
disp('Loading...')
showplot = false;
amp_cut = 7;
Fs = 200;
pattern_total_time = 10; % [s]
reg = true; % use interp times
start_idx = Fs*1; % 1 second in
add1 = false; % add 1st frame

tintrp = (0:(1/Fs):(10 - 1/Fs))';
[wing.b,wing.a] = butter(2,Fc/(Fs/2),'low');
[head.b,head.a] = butter(2,40/(Fs/2),'low');
SACCADE = [I , splitvars(table(num2cell(zeros(N.file,3))))]; % store saccade objects
SACCADE.Properties.VariableNames(5:7) = {'head_saccade','wing_saccade','head2wing'}; 
Head_Wing = cell(N.fly,N.vel);
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    load(fullfile(PATH.head, [basename{kk} '.mat']),'hAngles'); % load head angles % time arrays
    benifly = ImportBenifly(fullfile(PATH.wing, [basename{kk} '.csv'])); % load head & wing angles
    
    % Sync video with trigger & pattern  
    [trig,~] = sync_pattern_trigger(t_p, data(:,2), pattern_total_time, ... 
                                data(:,1), reg, start_idx, add1, false);
    
    % Get head data
    head.pos = hAngles;
    head.pos = interp1(trig.time_sync, head.pos, tintrp, 'pchip');
    head.pos = filtfilt(head.b, head.a, head.pos);
    
  	% Get wing data
    wing.left  = benifly.LWing;
    wing.right = benifly.RWing;
  	wing.left  = hampel(trig.time_sync, wing.left);
    wing.right = hampel(trig.time_sync, wing.right);
    wing.left  = interp1(trig.time_sync, wing.left, tintrp, 'pchip');
    wing.right = interp1(trig.time_sync, wing.right, tintrp, 'pchip');
    wing.left  = filtfilt(wing.b, wing.a, wing.left);
    wing.right = filtfilt(wing.b, wing.a, wing.right);
    wing.dwba  = rad2deg(wing.left - wing.right);
    [b_high,a_high] = butter(3,0.5/(Fs/2), 'high');
    wing.wba = filtfilt(b_high, a_high, wing.dwba);
    
    % Get head saccades
    peaks = [];
    head_saccade = saccade(head.pos, tintrp, 250, amp_cut, 0, peaks, nan, showplot);
    
    % Get wing saccades
    [dwba_pks,~,~] = wingSaccade(wing.dwba, tintrp, 1, false);
    wing_saccade = saccade(wing.dwba, tintrp, 0, 0, 1, dwba_pks, 0.5, false);
    
    if head_saccade.count==0
        % skip
     	SACCADE{kk,5} = {[]};
        SACCADE{kk,6} = {[]};
        SACCADE{kk,7} = {[]};
    else
        % head2wing
        win = 0.1; % [s]
        head2wing = head_wing_saccade(head_saccade, wing_saccade, win, true);
        pause
        Head_Wing{I.fly(kk),I.vel(kk)}(end+1,1) = head2wing;
        SACCADE{kk,5} = {head_saccade};
        SACCADE{kk,6} = {wing_saccade};
        SACCADE{kk,7} = {head2wing};
    end

    if showplot
        pause
        close all
    end
end

% Fill in empty saccade trials
empty_idx = cellfun(@(x) isempty(x), Head_Wing);
Head_Wing(empty_idx) = {saccade(nan*head.pos, nan*head.pos, 300, 0, 0, [], [], false)};

empty_idx = cellfun(@(x) isempty(x), SACCADE.head2wing);
SACCADE = SACCADE(~empty_idx,:);

%% Wing saccade interval
Wing_all = cellfun(@(x) x.wint, SACCADE.head2wing, 'UniformOutput', false);
Wing_all = cat(2,Wing_all{:});

%% Cross-Correlation
Tlag_all = cellfun(@(x) x.timelags, SACCADE.head2wing, 'UniformOutput', false);
Tlag_all = mean(cat(2,Tlag_all{:}),2);
Tlag_sacd_all = cellfun(@(x) x.int_lags, SACCADE.head2wing, 'UniformOutput', false);
Tlag_sacd_all = mean(cat(2,Tlag_sacd_all{:}),2);

% By trial
Acor_all = cellfun(@(x) x.acor, SACCADE.head2wing, 'UniformOutput', false);
Acor_all = cat(2,Acor_all{:});
Acor_stats = basic_stats(Acor_all,2);
[Acor_max, Acor_maxI] = max(Acor_stats.mean);
Acor_td = Tlag_all(Acor_maxI);

Acor_sacd_all = cellfun(@(x) x.int_acor, SACCADE.head2wing, 'UniformOutput', false);
Acor_sacd_all = cat(2,Acor_sacd_all{:});
Acor_sacd_stats = basic_stats(Acor_sacd_all,2);
[Acor_sacd_max, Acor_sacd_maxI] = max(Acor_sacd_stats.mean);
Acor_sacd_td = Tlag_sacd_all(Acor_sacd_maxI);

% By fly
Acor_fly = cellfun(@(x) mean(cat(2,x(:).acor),2), Head_Wing, 'UniformOutput', false);
Acor_fly = cat(2,Acor_fly{:});
Acor_fly_stats = basic_stats(Acor_fly,2);
[Acor_fly_max, Acor_fly_maxI] = max(Acor_fly_stats.mean);
Acor_fly_td = Tlag_all(Acor_fly_maxI);

Acor_sacd_fly = cellfun(@(x) mean(cat(2,x(:).int_acor),2), Head_Wing, 'UniformOutput', false);
Acor_sacd_fly = cat(2,Acor_sacd_fly{:});
Acor_sacd_fly_stats = basic_stats(Acor_sacd_fly,2);
[Acor_sacd_fly_max, Acor_sacd_fly_maxI] = max(Acor_sacd_fly_stats.mean);
Acor_sacd_fly_td = Tlag_sacd_all(Acor_sacd_fly_maxI);

trial_color = [0 0 0];
sacd_color = [0 0 0.7];

%% Cross-Correlation for all trials and saccades (not by fly)
fig = figure (1);
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 4])
movegui(fig, 'center')
    ax(1) = subplot(2,1,1) ; cla ; hold on ; ylabel('Cross Correlation')
            plot(Tlag_all, Acor_all, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            plot(Tlag_all, Acor_stats.mean, 'Color', trial_color, 'LineWidth', 2)
            plot([Acor_td Acor_td], [0 Acor_max], 'r', 'LineWidth', 1)
    ax(2) = subplot(2,1,2) ; cla ; hold on ; ylabel('Cross Correlation')
            plot(Tlag_sacd_all, Acor_sacd_all, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            plot(Tlag_sacd_all, Acor_sacd_stats.mean, 'Color', sacd_color, 'LineWidth', 2)
            plot([Acor_sacd_td Acor_sacd_td], [0 Acor_sacd_max], 'r', 'LineWidth', 1)
            xlabel('Time Lag (s)')

 	set(ax, 'LineWidth', 1, 'YLim', [-1 1])
    set(ax(1), 'XLim', 5*[-1 1]) 
    set(ax(2), 'XLim', 0.5*[-1 1])
    
%% Cross-Correlation for all trials and saccades (by fly)
fig = figure (2);
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 1.5 3])
movegui(fig, 'center')
    ax(1) = subplot(2,1,1) ; cla ; hold on ; title(['TD = ' num2str(1000*Acor_fly_td) 'ms'])
            ylabel('Cross Correlation')
            plot(Tlag_all, Acor_fly, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            %plot(Tlag_all, Acor_fly_stats.mean, 'Color', trial_color, 'LineWidth', 2)
            [~] = PlotPatch(Acor_fly_stats.mean, Acor_fly_stats.std, Tlag_all, 1, 1, ...
                               trial_color, 0.5*trial_color, 0.3, 1);
            plot([Acor_fly_td Acor_fly_td], [0 Acor_fly_max], 'r', 'LineWidth', 1)
    ax(2) = subplot(2,1,2) ; cla ; hold on ; title(['TD = ' num2str(1000*Acor_sacd_fly_td) 'ms'])
            ylabel('Cross Correlation')
            plot(Tlag_sacd_all, Acor_sacd_fly, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
            %plot(Tlag_sacd_all, Acor_sacd_fly_stats.mean, 'Color', sacd_color, 'LineWidth', 2)
            [~] = PlotPatch(Acor_sacd_fly_stats.mean, Acor_sacd_fly_stats.std, Tlag_sacd_all, 1, 1, ...
                               sacd_color, 0.5*sacd_color, 0.3, 1);
            plot([Acor_sacd_fly_td Acor_sacd_fly_td], [0 Acor_sacd_fly_max], 'r', 'LineWidth', 1)
            xlabel('Time Lag (s)')

 	set(ax, 'LineWidth', 1, 'YLim', [-1 1])
    set(ax(1), 'XLim', 2*[-1 1]) 
    set(ax(2), 'XLim', 0.2*[-1 1])
    
    %set(ax(1), 'XLim', 0.1*[-1 1], 'YLim', [0.45 0.75]) 
    %set(ax(1), 'XLim', 2*[-1 1], 'YLim', [-0.2 0.8]) 

end