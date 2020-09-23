function [] = Ramp_example_trial()
%% Ramp_example_trial:
root = 'H:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat'},'Select head angle trials', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'SACCADE','Stim','D','I','U','N')
n_speed = N.vel/2;

%% Active Tracking
active_30 = 100;
active_60 = 150;
active_90 = 265;
methods_30 = 36;
methods_90 = 60;

idx = active_60;
trial = SACCADE.head_saccade{idx};

FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 4 3];
FIG.Name = 'Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h
ax(1) = subplot(4,1,1) ; hold on ; title(['Stimulus: ' num2str(D.vel(idx)) ' (°/s)'])
    ylabel('Position (°)')
    h(1) = plot(trial.time, trial.position, 'k', 'LineWidth', 0.25);
    if trial.count~=0
        for ww = 1:trial.count
           plot(trial.saccades{ww}.Time,trial.saccades{ww}.Position,...
               'LineWidth', 0.5, 'Color', 'c')
        end
    end
    ax(1).YLim = 20*[-1 1];

ax(2) = subplot(4,1,2) ; hold on
    ylabel('Velocity (°/s)')
    h(2) = plot(trial.time, trial.velocity, 'k', 'LineWidth', 0.25);
    if trial.count~=0
        for ww = 1:trial.count
           plot(trial.saccades{ww}.Time, trial.saccades{ww}.Velocity,...
               'LineWidth', 0.5, 'Color', 'c')
        end
    end
    ax(2).YLim = 800*[-1 1];

ax(3) = subplot(4,1,3) ; hold on
    if trial.count~=0
        ylabel({'Removed', 'Position (°)'})                
        plot(trial.shift.Time, trial.shift.IntrpPosition,      	 'k', 'LineWidth', 1);
        plot(trial.shift.Time, trial.shift.Position,             'k', 'LineWidth', 1);
        plot(trial.removed_all.Time, trial.removed_all.Position, 'c', 'LineWidth', 1);
    end
    plot(trial.time, trial.stimlus_position, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    
ax(4) = subplot(4,1,4) ; hold on
    if trial.count~=0
        ylabel({'Removed', 'Velocity (°/s)'})
        xlabel('Time (s)')
        [b,a] = butter(2,2/(100/2),'low');
        filtvel = filtfilt(b, a, trial.shift.IntrpVelocity);
        plot(trial.shift.Time, trial.shift.IntrpVelocity,   'r', 'LineWidth', 0.5);
        plot(trial.shift.Time, trial.shift.Velocity,        'k', 'LineWidth', 0.5);
        plot(trial.shift.Time, filtvel,                     'c', 'LineWidth', 1);
    end
    ylim(400*[-1 1])
    
set(ax, 'LineWidth', 1,'XLim', [-0.2 10], 'XTick', 0:10)
linkaxes(ax,'x')
set(h,'LineWidth',0.25)
align_Ylabels(FIG)

%% Head Methods
idx = 229; % 150 deg/s
trial = SACCADE.head_saccade{idx};
                     
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 7 7])
movegui(fig, 'center')
clear ax h hh
mrk_sz = 5;
Lw = 0.5;
int_dull = 0.7;

ax(1) = subplot(6,4,1:2) ; hold on ; title('Position (°)')
    h.main(1) = plot(trial.time, trial.position, 'Color', 'k');

ax(2) = subplot(6,4,3:4) ; hold on ; title('Velocity (°)')
    h.main(2) = plot(trial.time, trial.velocity, 'Color', 'k');
    
ax(3) = subplot(6,4,5:6) ; hold on
    h.main(3) = plot(trial.time, trial.position_filt_detect, 'Color', 'k');
    
ax(4) = subplot(6,4,7:8) ; hold on
    h.main(4) = plot(trial.time, trial.velocity_filt_detect, 'Color', 'k');
 	plot([0 10],  [trial.threshold(1) trial.threshold(1)], 'Color', 'm');
    plot([0 10],  [trial.threshold(2) trial.threshold(2)], 'Color', 'm');
	h.main(5) = plot(trial.SACD.PeakTime, trial.peaks.velocity_detect, '.', 'MarkerFaceColor', 'none', ...
        'MarkerSize', mrk_sz, 'Color', 'm');
    
ax(5) = subplot(6,4,9:10) ; hold on
    h.main(6) = plot(trial.time, trial.position, 'Color', 'k');
    h.sacd(1) = plot(trial.SACD.StartTime, trial.SACD.StartPos, '.', 'MarkerFaceColor', 'none', ...
        'MarkerSize', mrk_sz, 'Color', 'g');
    h.sacd(2) = plot(trial.SACD.PeakTime, trial.SACD.PeakPos, '.', 'MarkerFaceColor', 'none', ...
        'MarkerSize', mrk_sz, 'Color', 'b');
    h.sacd(3) = plot(trial.SACD.EndTime, trial.SACD.EndPos, '.', 'MarkerFaceColor', 'none', ...
        'MarkerSize', mrk_sz, 'Color', 'r');

ax(6) = subplot(6,4,11:12) ; hold on
    h.main(7) = plot(trial.time, trial.velocity, 'Color', 'k');
    h.sacd(4) = plot(trial.SACD.StartTime, trial.SACD.StartVel, '.', 'MarkerFaceColor', 'none', ...
        'MarkerSize', mrk_sz, 'Color', 'g');
    h.sacd(5) = plot(trial.SACD.PeakTime, trial.SACD.PeakVel, '.', 'MarkerFaceColor', 'none', ...
        'MarkerSize', mrk_sz, 'Color', 'b');
    h.sacd(6) = plot(trial.SACD.EndTime, trial.SACD.EndVel, '.', 'MarkerFaceColor', 'none', ...
        'MarkerSize', mrk_sz, 'Color', 'r');
    plot([0 10], [trial.true_thresh trial.true_thresh], 'Color', 'b');
    plot([0 10], -[trial.true_thresh trial.true_thresh], 'Color', 'b');
    
ax(7) = subplot(6,4,13:14) ; hold on
    h.main(8) = plot(trial.time, trial.position, 'Color', 'k');
    for ww = 1:trial.count
       h.sacd(7) = plot(trial.saccades{ww}.Time, trial.saccades{ww}.Position,...
           'LineWidth', Lw, 'Color', trial.cmap(ww,:));
       h.sacd(8) = plot(trial.intervals{ww}.Time, trial.intervals{ww}.Position,...
           'LineWidth', Lw, 'Color', int_dull*trial.cmap(ww,:));
    end

ax(8) = subplot(6,4,15:16) ; hold on
    h.main(9) = plot(trial.time, trial.velocity, 'Color', 'k');
    for ww = 1:trial.count
       hh.sacd(9) = plot(trial.saccades{ww}.Time, trial.saccades{ww}.Velocity,...
           'LineWidth', Lw, 'Color', trial.cmap(ww,:));
       hh.sacd(10) = plot(trial.intervals{ww}.Time, trial.intervals{ww}.Velocity,...
           'LineWidth', Lw, 'Color', int_dull*trial.cmap(ww,:));
    end
    
ax(9) = subplot(6,4,17) ; hold on
    hh.sacdpos = plot(trial.normpeak_saccade.time, trial.normpeak_saccade.position);

ax(10) = subplot(6,4,18) ; hold on
    hh.intpos = plot(trial.normstart_interval.time, trial.norm_interval.position);
    
ax(11) = subplot(6,4,19) ; hold on
    hh.sacdvel = plot(trial.normpeak_saccade.time, trial.normpeak_saccade.velocity);

ax(12) = subplot(6,4,20) ; hold on
    [b, a] = butter(3, 10 / (trial.Fs/2), 'low');
    for ww = 1:trial.count
        endI = sum(~isnan(trial.norm_interval.velocity(:,ww)));
        if endI > 0
            filt_int_vel = filtfilt(b, a, trial.norm_interval.velocity(1:endI,ww));
        else
            filt_int_vel = trial.norm_interval.velocity(:,ww);
            endI = size(trial.norm_interval.velocity, 1);
        end
        tt = trial.Ts*(0:endI-1)';
        hh.intvel(ww) = plot(tt, filt_int_vel);
    end
    
ax(13) = subplot(6,4,21:22) ; hold on ; ylim([-300 20]) ; xlim([-0.5 10])
    plot(trial.shift.Time, trial.shift.IntrpPosition, 'r', 'LineWidth', 0.5);
    plot(trial.shift.Time, trial.shift.Position, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    plot(trial.removed_all.Time, trial.removed_all.Position, 'k', 'LineWidth', 0.5);
    %plot(trial.time, stim_pos, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    for ww = 1:trial.count
       hh.sacd(11) = plot(trial.intervals{ww}.Time, trial.intervals{ww}.Position,...
           'LineWidth', 1, 'Color', int_dull*trial.cmap(ww,:));
    end
    
ax(14) = subplot(6,4,23:24) ; hold on ; ylim(300*[-1 1]) ; yticks(300*[-1 0 1]) ; xlim([-0.5 10])
%ylim([-300 20])
    [b, a] = butter(3, 1 / (trial.Fs/2), 'low');
    vel_intrp = diff(trial.shift.IntrpPosition) * trial.Fs; vel_intrp = [vel_intrp(1) ; vel_intrp];
    filt_vel_intrp = filtfilt(b, a, vel_intrp);
 	plot(trial.shift.Time, trial.shift.IntrpVelocity, 'k', 'LineWidth', 0.5);
    plot(trial.shift.Time, filt_vel_intrp, 'r', 'LineWidth', 0.5);
%     plot(trial.shift.Time, trial.shift.Position, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
%     plot(trial.removed_all.Time, trial.removed_all.Position, 'g', 'LineWidth', 0.5);   

set(hh.sacdpos, {'color'}, num2cell(trial.cmap,2), 'LineWidth', Lw)
set(hh.sacdvel, {'color'}, num2cell(trial.cmap,2), 'LineWidth', Lw)
set(hh.intpos,  {'color'}, num2cell(int_dull*trial.cmap,2), 'LineWidth', Lw)
set(hh.intvel,  {'color'}, num2cell(int_dull*trial.cmap,2), 'LineWidth', Lw)
set(h.main, 'LineWidth', Lw)
set(h.sacd, 'LineWidth', Lw)

set(ax, 'LineWidth', 0.75)
set(ax(1:8), 'XLim', [-0.5 10]) 
set(ax([1:2:9,10]), 'YLim', 20*[-1 1], 'YTick', -20:10:20)
set(ax([2:2:8,11]), 'YLim', [-400 800])
set(ax([1:2:5]), 'XColor', 'none')
set(ax([2:2:6]), 'XColor', 'none')
set(ax([3:2:7]), 'YTickLabels', [])
set(ax([4:2:8]), 'YTickLabels', [])
set(ax([9,11]), 'XLim', 0.04*[-1 1])
set(ax([10,12]), 'XLim', [0 0.7])

%% Wing Methods (static)
mrk_sz = 5;
Lw = 0.5;
% int_dull = 0.7;

idx = 2;
head_trial = SACCADE.head_saccade{idx};
wing_trial = SACCADE.wing_saccade{idx};

[b, a] = butter(3, 10 / (head_trial.Fs/2), 'low');
head_pos = filtfilt(b, a, head_trial.position);
wing_pos = filtfilt(b, a, wing_trial.extra.dwba);

fig = figure (4) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 7 4])
movegui(fig, 'center')
clear ax h hh

ax(1) = subplot(1,1,1) ; hold on ; ylim(40*[-1 1]) ; xlim([-0.2 10])
    %plot(head_trial.time, head_trial.position, 'b', 'LineWidth', Lw)
    plot(head_trial.time, head_pos, 'b', 'LineWidth', Lw)
    %plot(wing_trial.time, wing_trial.extra.dwba, 'r', 'LineWidth', Lw)
    plot(wing_trial.time, wing_pos, 'r', 'LineWidth', Lw)

    for ww = 1:wing_trial.count
      	plot(wing_trial.saccades{ww}.Time, wing_trial.saccades{ww}.Position,...
           'LineWidth', 1, 'Color', 'k');
    end

end