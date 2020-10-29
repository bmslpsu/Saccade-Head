function [] = Make_Passive_head(wave)
%% Make_Passive_head:
%   INPUTS:
%       wave    :   spatial wavelength of data
%
%   OUTPUTS:
%       -
%
warning('off', 'signal:findpeaks:largeMinPeakHeight')

% Data location
rootdir = 'H:\EXPERIMENTS\RIGID\Passive Head Displacement\Vid';

% Output file name
filename = 'Passive_head';

% Setup Directories
PATH.vid  = rootdir; % video data location
PATH.head = fullfile(PATH.vid,'tracked_head'); % tracked kinematic data location

% Select files
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.head,'*.mat',false);

%% Get Data %%
clc
close all
disp('Loading...')
% Fs = 200; % sampling frequency [s]
% tintrp = (0:(1/Fs):(10 - 1/Fs))'; % time vector for interpolation

% HEAD saccade detection parameters
head.showplot = false;
head.Fc_detect = [10 nan];
head.Fc_ss = [nan nan];
head.amp_cut = 4;
head.dur_cut = inf;
head.thresh = [60 , 0, 0, 0];
head.true_thresh = 100;
head.direction = 0;
head.sacd_length = nan;
head.pks = [];
head.min_pkdist = 0.2;
head.min_pkwidth = 0.03;
head.min_pkprom = 75;
head.min_pkthresh = 0;
head.boundThresh = [0 50];
head.Fc = 60;

t_win = 0.2;
% end_win = 1;
DATA = [I , splitvars(table(num2cell(zeros(N.file,8))))];
DATA.Properties.VariableNames(3:end) = {'head_saccade','time', 'shift_time', 'pos', 'vel', 'accel', ...
                                            'fit', 'gof'};
DATA = [DATA , splitvars(table(zeros(N.file,8)))];
DATA.Properties.VariableNames(11:end) = {'r2', 'tau', 'rise_time', 'peak_time', 'amp', 'peak_vel', ...
                                            'init_pos', 'init_vel'};
for kk = 1:N.file
    disp(kk)
    %disp(basename{kk})
    % Load HEAD & DAQ data
    load(fullfile(PATH.head, [basename{kk} '.mat']),'hAngles','t_v'); % load head angles % time arrays
 
    % Get head data
    head.time = t_v;
    head.Fs = round(1 / mean(diff(head.time)));
    [head.b, head.a] = butter(3, head.Fc / (head.Fs/2) ,'low');
    head.pos = hAngles;
    
 	% Head with filter
 	head.pos_filt = filtfilt(head.b, head.a, head.pos);
    
    % Extract head pulse
    head_saccade = saccade_all(head.pos_filt, head.time, head.thresh, head.true_thresh, head.Fc_detect, ...
                                head.Fc_ss, head.amp_cut, head.dur_cut, head.direction, head.pks, ...
                                head.sacd_length, head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                                head.min_pkthresh, head.boundThresh, head.showplot);
    
  	assert(head_saccade.count == 1, 'No step dectected')
    scdI = head_saccade.saccades{1}.Index;
    winI = scdI(1) : scdI(end) + ( round(t_win*head.Fs) - length(scdI) );
    tt = (1/head.Fs)*(0:length(winI)-1)';
    peak_time = (1/head.Fs)*( head_saccade.SACD.PeakIdx - scdI(1) );
    %end_winI = scdI(end) : scdI(end) + round(end_win*head.Fs);
    %head_temp = head.pos_filt - median(head.pos_filt(end_winI));
    head_temp = head.pos_filt(winI);
    head_temp = head_temp - head_temp(end);
    
    DATA.time{kk} = tt;
    DATA.shift_time{kk} = tt - peak_time;
    DATA.pos{kk} = head_temp;
    DATA.vel{kk} = head_saccade.velocity(winI);
    DATA.accel{kk} = head_saccade.acceleration(winI);
    DATA.head_saccade{kk} = head_saccade;
    
    [fitresult, gof] = exp_fit(tt, DATA.pos{kk}, false);
    tau = 1000 * (1 / fitresult.b);
    
  	[fitresult, gof] = test_fit(tt, DATA.pos{kk}, false);

    %disp([1000/fitresult.b])
    %disp([1000/fitresult.d])
    %pause
  	%clc
    
    DATA.fit{kk} = fitresult;
    DATA.gof{kk} = gof;
    DATA.r2(kk) = gof.rsquare;
    DATA.tau(kk) = tau;
  	DATA.rise_time(kk) = 1000*head_saccade.SACD.RiseTime;
    DATA.peak_time(kk) = peak_time;
    DATA.peak_vel(kk) = head_saccade.SACD.PeakVel;
    DATA.amp(kk) = DATA.pos{kk}(end) - DATA.pos{kk}(1);
    DATA.init_pos(kk) = DATA.pos{kk}(1);
    DATA.init_vel(kk) = DATA.vel{kk}(1);
    
%     hold on ; cla
%     plot(DATA.shift_time{kk}, DATA.vel{kk}, 'k')
%     plot(0, DATA.peak_vel(kk), 'ro')
%     pause
    
    if head.showplot
        figure (1)
        pause
        close all
    end
end

%% Analytical Response
% Fit coefficents
% A = mean(cellfun(@(x) x.a, DATA.fit));
% B = mean(cellfun(@(x) x.b, DATA.fit));
% C = mean(cellfun(@(x) x.b, DATA.fit));
% 
% % Expressions
% clear F dF d2F
% syms t
% F(t) = A*exp(-B*t) + C;
% dF(t) = diff(F,t,1);
% d2F(t) = diff(F,t,2);
% 
% % Output
% F = double(F(tt));
% F = F - F(end);
% dF = double(dF(tt));
% d2F = double(d2F(tt));

%% Time constant
% fig = figure (1) ; clf
% set(fig, 'Color', 'w','Units', 'inches', 'Position', 1*[2 2 2 2])
% clear ax
% ax(1) = subplot(1,1,1); hold on
% ylim([0 50])
% ylabel('Time Constant (ms)')
% % bx = boxplot(cat(2,DATA.tau{:}), I.Fly, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
% bx = boxplot(DATA.tau, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
% 
% h = get(bx(5,:),{'XData','YData'});
% for kk = 1:size(h,1)
%    patch(h{kk,1}, h{kk,2}, 'k',  'EdgeColor', 'none');
% end
% ww = 1;
% set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth', 1);
% set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
% set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
% set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
% ax(ww).Children = ax(ww).Children([end 1:end-1]);
% ylim([0 60])
% set(ax , 'LineWidth', 1, 'Box', 'off')

%% Intitial Condition Response
% pos_all = cat(2, DATA.pos{:});
% pos_stats = basic_stats(pos_all,2);
% vel_all = cat(2, DATA.vel{:});
% vel_stats = basic_stats(vel_all,2);
% accel_all = cat(2, DATA.accel{:});
% accel_stats = basic_stats(accel_all,2);
% 
% fig = figure (2) ; clf
% set(fig, 'Color', 'w','Units', 'inches', 'Position', 1*[3 1 2 4])
% clear ax
% ax(1) = subplot(3,1,1); hold on
%     ylim([-40 10])
%     plot(tt, pos_all, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
%     [~,~] = PlotPatch(pos_stats.mean, pos_stats.std, tt, 0, 1, 'k', 'r', 0.2, 1.5);
%     plot(tt, F, '--', 'Color', 'r', 'LineWidth', 1)
%     ylabel('Position (°)')
% 
% ax(2) = subplot(3,1,2); hold on
%     ylim([-250 2500])
%     plot(tt, vel_all, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
%     [~,~] = PlotPatch(vel_stats.mean, vel_stats.std, tt, 0, 1, 'k', 'r', 0.2, 1.5);
%     %plot(tt, -dF, 'Color', 'r', 'LineWidth', 1)
%     ylabel('Velocity (°/s)')
%     
% ax(3) = subplot(3,1,3); hold on
%     ylim(2.6e5*[-1 1])
%     plot(tt, accel_all, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.5)
%     [~,~] = PlotPatch(accel_stats.mean, accel_stats.std, tt, 0, 1, 'k', 'r', 0.2, 1.5);
%     %plot(tt, -d2F, 'Color', 'r', 'LineWidth', 1)
%     ylabel('Acceleration (°/s^{2})')
%     xlabel('Time (s)')
% 
% set(ax , 'LineWidth', 1, 'XLim', [-0.005 0.06])
% linkaxes(ax, 'x')

%% Fit Mean
% xx = -pos_stats.mean;
% [fitresult, gof] = exp_fit(tt, xx, true);
% tau = 1000 * (1 / fitresult.b);

%% SAVE %%
disp('Saving...')
save(['H:\DATA\Rigid_Data\Saccade\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','DATA','D','I','U','N','T','-v7.3')
disp('SAVING DONE')
end