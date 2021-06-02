function [] = Make_static_head_electromagnet()
%% Make_static_head_electromagnet:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%
warning('off', 'signal:findpeaks:largeMinPeakHeight')

% Data location
rootdir = 'E:\EXPERIMENTS\RIGID\Experiment_head_electro_magnet_background';

% Output file name
filename = 'Magnet_static';

% Setup Directories 
PATH.raw  = rootdir;
PATH.head = fullfile(PATH.raw,'tracked_bodypart_tip');
PATH.filt = fullfile(PATH.raw,'filt');
PATH.wing = fullfile(PATH.filt,'tracked_wing');

% Select files
% [D,I,N,U,T,~,~,basename] = GetFileData(PATH.wing,'*.csv',false,'fly','trial');
[D,I,N,U,T,~,~,basename] = GetFileData(PATH.head,'*.mat',false,'fly','trial');

%% Get Data
clc
close all
disp('Loading...')
Fs = 100; % sampling frequency [s]
tintrp = (0:(1/Fs):(20 - 1/Fs))'; % time vector for interpolation

head.Fc = 20;
[head.b, head.a] = butter(3, head.Fc / (Fs/2) ,'low');

wing.Fc = 20;
[wing.b, wing.a] = butter(3, wing.Fc / (Fs/2) ,'low');

wing_carry.Fc = 25;
[wing_carry.b, wing_carry.a] = butter(3, wing_carry.Fc / (Fs/2) ,'low');

DATA = [I , splitvars(table(num2cell(zeros(N.file,8))))]; % store saccade objects
DATA.Properties.VariableNames(4:end) = ...
    {'head_on','head_off','wing_on','wing_off','on_time','off_time','magnet_on','magnet_off'};
showplot = false;
for kk = 1:N.file
    disp(kk)
    %basename{kk}
    % Load data
	load(fullfile(PATH.raw, [basename{kk} '.mat']),'data','t_p','on_time','off_time');
    load(fullfile(PATH.head, [basename{kk} '.mat']),'bodypart');
    benifly = ImportBenifly(fullfile(PATH.wing, [basename{kk} '.csv']));

  	% Sync frames and get pattern data
	[trig,pat] = sync_pattern_trigger(t_p, data(:,2), 20, data(:,1), true, 1, false, false);
    vid_time = trig.time_sync + 0.01;
    
    % Get head data
    head.pos = interp1(vid_time, bodypart.angle, tintrp, 'pchip');
    
	% Head with filter
 	head.pos_filt = filtfilt(head.b, head.a, head.pos);
    
    % Get wing data
    hampel_dx       = round(0.05*Fs);
    T               = 2.5;
    wing.left       = rad2deg(benifly.LWing);
    wing.right      = rad2deg(benifly.RWing);
    wing.left       = hampel(vid_time, wing.left, hampel_dx, T);
    wing.right      = hampel(vid_time, wing.right, hampel_dx, T);
    wing.left       = interp1(vid_time, wing.left, tintrp, 'pchip');
    wing.right      = interp1(vid_time, wing.right, tintrp, 'pchip');
    
    wing.dwba       = wing.left - wing.right;
    wing.dwba       = filtfilt(wing_carry.b, wing_carry.a, wing.dwba);
    wing.dwba_vel   = central_diff(wing.dwba, 1/Fs);
    
    wing.left_filt	= filtfilt(wing.b, wing.a, wing.left);
    wing.right_filt	= filtfilt(wing.b, wing.a, wing.right);
    wing.dwba_filt	= wing.left_filt - wing.right_filt;

    % Get magnet signal
    switch D.control(kk)
        case 1 % magnet on
            magnet = round(data(:,7));
            magnet = interp1(pat.time_sync, magnet, tintrp, 'nearest');
        case 0 % magnet off
            dutyC = 100*(on_time / (off_time + on_time));
            f = 1 ./ (off_time + on_time);
            magnet = (1 - square(2*pi*f*tintrp, dutyC))/2;
        otherwise
            error('control must be 1 or 0')
    end
    
    % Get head & wing windows
    win_length = [0.1 1.99].*Fs;
    norm = [true true];
    [W,tt] = pull_windows(wing.dwba_filt, magnet, win_length, Fs, norm, false);
    normI = find(tt{1}(:,1) == 0.05):find(tt{1}(:,1) == 0.2);
    DATA.wing_on{kk} = W{1};
    DATA.wing_off{kk} = W{2} + mean(DATA.wing_on{kk}(normI,:),'all');
    DATA.on_time{kk} = mean(tt{1},2);
    DATA.off_time{kk} = mean(tt{2},2);
    
    [W,~] = pull_windows(head.pos_filt, magnet, win_length, Fs, norm, false);
    DATA.head_on{kk} = W{1};
    DATA.head_off{kk} = W{2} + mean(DATA.head_on{kk}(normI,:),'all');
        
    [W,~] = pull_windows(magnet, magnet, win_length, Fs, false, false);
    DATA.magnet_on{kk} = mean(W{1},2);
    DATA.magnet_off{kk} = mean(W{2},2);
    
    if showplot
        fig = figure (1); clf
        set(fig, 'Color', 'w', 'Units', 'inches')
        switch D.control(kk)
            case 1
                suptitle('Magnet: ON')
            case 0
                suptitle('Magnet: OFF')
        end
        ax = gobjects(2,2);
        clear h_on h_off
        subplot(2,2,1)
            yyaxis left ; hold on ; ylabel('Head (°)')
            ax(1,1) = gca;
                h_on(1,:) = plot(DATA.on_time{kk}, DATA.head_on{kk}, '-');
                yline(0, '--', 'Color', [0.5 0.5 0.5])
            yyaxis right
                plot(DATA.on_time{kk}, DATA.magnet_on{kk}, '-k', 'LineWidth', 1)
            
        subplot(2,2,3)
            yyaxis left ; hold on ; ylabel('\DeltaWBA (°)') ; xlabel('Time (s)')
            ax(2,1) = gca;
                h_on(2,:) = plot(DATA.on_time{kk}, DATA.wing_on{kk}, '-');
                yline(0, '--', 'Color', [0.5 0.5 0.5])
            yyaxis right
                plot(DATA.on_time{kk}, DATA.magnet_on{kk}, '-k', 'LineWidth', 1)

        subplot(2,2,2)
            yyaxis left ; hold on ; ylabel('Head (°)')
            ax(1,2) = gca;
                h_off(1,:) = plot(DATA.off_time{kk}, DATA.head_off{kk}, '-');
                yline(0, '--', 'Color', [0.5 0.5 0.5])
            yyaxis right
                plot(DATA.off_time{kk}, DATA.magnet_off{kk}, '-k', 'LineWidth', 1)

        subplot(2,2,4)
            yyaxis left ; hold on ; ylabel('\DeltaWBA (°)') ; xlabel('Time (s)')
            ax(2,2) = gca;
                h_off(2,:) = plot(DATA.off_time{kk}, DATA.wing_off{kk}, '-');
                yline(0, '--', 'Color', [0.5 0.5 0.5])
            yyaxis right
                ylim([-0.1 1.1])
                plot(DATA.off_time{kk}, DATA.magnet_off{kk}, '-k', 'LineWidth', 1)
                
        set(ax, 'Color', 'none', 'LineWidth', 1, 'Box', 'off', 'XLim', [-0.1 1])
        linkaxes(ax, 'x')
        
        for a = [1 2]
            yy = [ax(a).YAxis(1).Limits ; ax(a+2).YAxis(1).Limits];
            ymax = [min(yy(:,1)) , max(yy(:,2))];
            set([ax(a).YAxis(1) ax(a+2).YAxis(1)], 'Limits', ymax)
        end
        
        for a = 1:4
           	yy = ax(a).YAxis(1).Limits;
            rr = yy(1) / yy(2);
            ax(a).YAxis(2).Limits = [rr 1];
            set(ax(a).YAxis(1), 'Color', 'k')
            set(ax(a).YAxis(2), 'Color', 'none') 
        end
        
        cc = distinguishable_colors(size(h_on,2), [0.5 0.5 0.5]);
        set(h_on(1,:), {'color'}, num2cell(cc,2));
        set(h_on(2,:), {'color'}, num2cell(cc,2));
        
        cc = distinguishable_colors(size(h_off,2), [0.5 0.5 0.5]);
        set(h_off(1,:), {'color'}, num2cell(cc,2));
        set(h_off(2,:), {'color'}, num2cell(cc,2));
        
        pause
    end
end

%% All
stats = table_fly_stats(DATA, 3, 4:11, false);

flyI = 3;
cntrl = 2;
on_time = stats.val_stats.on_time(1).mean;
off_time = stats.val_stats.off_time(1,1).mean;

fig = figure (100); clf
set(fig, 'Color', 'w', 'Units', 'inches')
ax = gobjects(2,2);
clear h_on h_off h_on_med h_off_med
yaxis_cc = 'b';
subplot(2,2,1)
    yyaxis left ; hold on ; ylabel('Head (°)')
    ax(1,1) = gca;
        h_on(:,1) = plot(on_time, stats.comb.val_all.head_on{cntrl}, '-', 'Color', [0.5 0.5 0.5]);
        h_on_med(:,1) = plot(on_time, stats.val_stats.head_on(cntrl).mean, 'k-');
        yline(0, '--', 'Color', yaxis_cc)
    yyaxis right
        plot(on_time, DATA.magnet_on{kk}, '-r', 'LineWidth', 1)

subplot(2,2,3)
    yyaxis left ; hold on ; ylabel('\DeltaWBA (°)') ; xlabel('Time (s)')
    ax(2,1) = gca;
        h_on(:,2) = plot(on_time, stats.comb.val_all.wing_on{1,cntrl}, '-', 'Color', [0.5 0.5 0.5]);
        h_on_med(:,2) = plot(on_time, stats.val_stats.wing_on(cntrl).mean, 'k-');
        yline(0, '--', 'Color', yaxis_cc)
    yyaxis right
        plot(on_time, DATA.magnet_on{kk}, '-r', 'LineWidth', 1)

subplot(2,2,2)
    yyaxis left ; hold on
    ax(1,2) = gca;
        h_off(:,1) = plot(on_time, stats.comb.val_all.head_off{1,cntrl}, '-', 'Color', [0.5 0.5 0.5]);
        h_off_med(:,1) = plot(on_time, stats.val_stats.head_off(cntrl).mean, 'k-');
        yline(0, '--', 'Color', yaxis_cc)
    yyaxis right
        plot(off_time, DATA.magnet_off{kk}, '-r', 'LineWidth', 1)

subplot(2,2,4)
    yyaxis left ; hold on ; xlabel('Time (s)')
    ax(2,2) = gca; set(ax(2,2), 'YColor', 'k')
        h_off(:,2) = plot(on_time, stats.comb.val_all.wing_off{1,cntrl}, '-', 'Color', [0.5 0.5 0.5]);
        h_off_med(:,2) = plot(on_time, stats.val_stats.wing_off(cntrl).mean, 'k-');
        yline(0, '--', 'Color', yaxis_cc)
    yyaxis right
        ylim([-0.1 1.1])
        plot(off_time, DATA.magnet_off{kk}, '-r', 'LineWidth', 1)

set(ax, 'Color', 'none', 'LineWidth', 1, 'Box', 'off', 'XLim', [-0.11 0.20001])
linkaxes(ax, 'x')

% set(ax(1,1).YAxis(1), 'Limits', [-20 60])
set(ax(2,1).YAxis(1), 'Limits', [-15 20])
set(ax(2,2).YAxis(1), 'Limits', [-15 20])

for a = [1 2]
    yy = [ax(a).YAxis(1).Limits ; ax(a+2).YAxis(1).Limits];
    ymax = [min(yy(:,1)) , max(yy(:,2))];
    set([ax(a).YAxis(1) ax(a+2).YAxis(1)], 'Limits', ymax)
end

for a = 1:4
    yy = ax(a).YAxis(1).Limits;
    rr = yy(1) / yy(2);
    ax(a).YAxis(2).Limits = [rr 1];
    set(ax(a).YAxis(1), 'Color', 'k')
    set(ax(a).YAxis(2), 'Color', 'none') 
end

set([h_on_med h_off_med], 'LineWidth', 1, 'Color', 'k')
set([h_on ; h_off], 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.25)

set([ax(1,2).YAxis(1) ax(2,2).YAxis(1)], 'Color', 'none')
set(ax(1,:), 'XColor', 'none')

%% Fly
stats = table_fly_stats(DATA, 3, 4:11, false);

cntrl = 2;
on_time = stats.val_stats.on_time(1).mean;
off_time = stats.val_stats.off_time(1,1).mean;
cc_fly = distinguishable_colors(N.fly, [0.5 0.5 0.5]);

fig = figure (100); clf
set(fig, 'Color', 'w', 'Units', 'inches')
fig.Position(3:4) = [4 4];
ax = gobjects(2,2);
clear h_on h_off h_on_med h_off_med
yaxis_cc = 'k';
ax(1,1) = subplot(2,2,1) ; hold on ; ylabel('Head (°)')
    for n = 1:N.fly
        h_on{n}(:,1) = plot(on_time, stats.head_on{n,cntrl}, '-', ...
            'Color', [0.8*cc_fly(n,:) 0.2]);
        h_on_med{n}(:,1) = plot(on_time, stats.fly_stats.head_on(n,cntrl).mean, '-', ...
            'Color', cc_fly(n,:));
    end
    cellfun(@(x) uistack(x,'top'), h_on_med)
    yline(0, '--', 'Color', yaxis_cc)

ax(2,1) = subplot(2,2,3); hold on ; ylabel('\DeltaWBA (°)') ; xlabel('Time (s)')
    for n = 1:N.fly
        h_on{n}(:,2) = plot(on_time, stats.wing_on{n,cntrl}, '-', ...
            'Color', [0.8*cc_fly(n,:) 0.2]);
        h_on_med{n}(:,2) = plot(on_time, stats.fly_stats.wing_on(n,cntrl).mean, '-', ...
            'Color', cc_fly(n,:));
    end
    cellfun(@(x) uistack(x(:,2),'top'), h_on_med)
    yline(0, '--', 'Color', yaxis_cc)

ax(1,2) = subplot(2,2,2); hold on
    for n = 1:N.fly
        h_off{n}(:,1) = plot(off_time, stats.head_off{n,cntrl}, '-', ...
            'Color', [0.8*cc_fly(n,:) 0.2]);
        h_off_med{n}(:,1) = plot(off_time, stats.fly_stats.head_off(n,cntrl).mean, '-', ...
            'Color', cc_fly(n,:));
    end
    cellfun(@(x) uistack(x,'top'), h_off_med)
    yline(0, '--', 'Color', yaxis_cc)

ax(2,2) = subplot(2,2,4) ; hold on ; xlabel('Time (s)')
    for n = 1:N.fly
        h_off{n}(:,2) = plot(off_time, stats.wing_off{n,cntrl}, '-', ...
            'Color', [0.8*cc_fly(n,:) 0.2]);
        h_off_med{n}(:,2) = plot(off_time, stats.fly_stats.wing_off(n,cntrl).mean, '-', ...
            'Color', cc_fly(n,:));
    end
    cellfun(@(x) uistack(x(:,2),'top'), h_off_med)
    yline(0, '--', 'Color', yaxis_cc)

set(ax, 'Color', 'none', 'LineWidth', 0.7, 'Box', 'off', 'XLim', [-0.11 0.30001])
linkaxes(ax, 'x')

% set(ax(2,1).YAxis(1), 'Limits', [-15 20])
% set(ax(2,2).YAxis(1), 'Limits', [-15 20])
set(ax(1,:), 'YLim', [-40 30])
set(ax(1,:), 'YTick', -30:10:30)
set(ax(2,:), 'YLim', [-20 20])

for a = [1 2]
    yy = [ax(a).YAxis(1).Limits ; ax(a+2).YAxis(1).Limits];
    ymax = [min(yy(:,1)) , max(yy(:,2))];
    set([ax(a).YAxis(1) ax(a+2).YAxis(1)], 'Limits', ymax)
end

cellfun(@(x) set(x, 'LineWidth', 2), h_on_med)
cellfun(@(x) set(x, 'LineWidth', 2), h_off_med)
cellfun(@(x) set(x, 'LineWidth', 0.25), h_on)
cellfun(@(x) set(x, 'LineWidth', 0.25), h_off)

set([ax(1,2).YAxis(1) ax(2,2).YAxis(1)], 'Color', 'none')
set(ax(1,:), 'XColor', 'none')

xx = [0 0 0.2 0.2];
for n = 1:4
   subplot(2,2,n)
   yy = [-50 50 50 -50];
   h_patch = patch(xx, yy, 'c', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
   %uistack(h_patch, 'bottom')
end

%% SAVE %%
% disp('Saving...')
% save(['E:\DATA\Rigid_Data\Saccade\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
%       'PATH','SACCADE','D','I','U','N','T','-v7.3')
% disp('SAVING DONE')
end