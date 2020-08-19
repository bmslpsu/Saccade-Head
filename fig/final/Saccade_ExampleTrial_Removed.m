function [] = Saccade_ExampleTrial_Removed()
%% Saccade_ExampleTrials_Removed:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat'},'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE','Stim','D','I','U','N')

n_speed = N.vel/2;

%% Active Tracking vs Landing Behavior %%
% clearvars -except clms CC Vel PATH COUNT SACCADE

FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 7];
FIG.Name = 'Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h
sacd_color = [0.5 0.5 0.5];

active_switch = 309;
land = 24;
active_30 = 100;
active_60 = 150;
active_90 = 265;

idx = 116; % 6-9, 2-58
trial = SACCADE.head_saccade{idx};
trial_wing = SACCADE.wing_saccade{idx};
% [wing_sacd, wing_ints] = getSaccade(trial_wing, trial_wing.extra.dwba);

ax(1) = subplot(4,1,1) ; hold on ; title(['Stimulus: ' num2str(D.vel(idx)) ' (°/s)'])
    ylabel('Position (°)')
    h(1) = plot(trial.time,trial.position,'k');
    plot(trial.time, zeros(trial.n,1),'--','Color',[0.5 0.5 0.5])
    if trial.count~=0
        for ww = 1:trial.count
           plot(trial.saccades{ww}.Time,trial.saccades{ww}.Position,...
               'LineWidth', 1, 'Color', 'c')
        end
    end
    ax(1).YLim = 20*[-1 1];

ax(2) = subplot(4,1,2) ; hold on
    ylabel('Velocity (°/s)')
    h(2) = plot(trial.time,trial.velocity,'k');

ax(3) = subplot(4,1,3) ; hold on
%     ylim([-600 20])
    if trial.count~=0
        ylabel({'Removed', 'Position (°)'})
        h(1) = plot(trial.time, trial.position,'b');
        for ww = 1:trial.count
           plot(trial.saccades{ww}.Time, trial.saccades{ww}.Position,...
               'LineWidth', 1, 'Color', sacd_color)
        end
        plot(trial.shift.Time, trial.shift.IntrpPosition, 'Color', sacd_color, 'LineWidth', 1);
        plot(trial.shift.Time, trial.shift.Position, 'Color', 'b', 'LineWidth', 1);
        %plot(trial.removed_all.Time, trial.removed_all.Position, 'c', 'LineWidth', 1);
    end
    
    if trial_wing.count~=0
        %h(1) = plot(trial_wing.time, trial_wing.extra.dwba,'r');
        h(1) = plot(trial_wing.time, trial_wing.position,'r');
        for ww = 1:trial_wing.count
           plot(trial_wing.saccades{ww}.Time, trial_wing.saccades{ww}.Position,...
               'LineWidth', 1, 'Color', sacd_color)
        end
        plot(trial_wing.shift.Time, trial_wing.shift.IntrpPosition, 'Color', sacd_color, 'LineWidth', 1);
        plot(trial_wing.shift.Time, trial_wing.shift.Position, 'Color', 'r', 'LineWidth', 1);
        %plot(trial.removed_all.Time, trial.removed_all.Position, 'c', 'LineWidth', 1);
    end
    plot(trial.time, trial.stimlus_position, '--', 'Color', 'g', 'LineWidth', 1);
    
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
    % ax(4).YLim = ax(2).YLim;
    ylim(400*[-1 1])
    
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8)
linkaxes(ax,'x')
set(h,'LineWidth',0.5)
align_Ylabels(FIG)

end