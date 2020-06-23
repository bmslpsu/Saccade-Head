function [] = SS_Saccade_ExampleTrials()
%% SS_Saccade_ExampleTrials:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','D','I','U','N')

clearvars -except clms CC Vel PATH COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

%% Example Trial
clearvars -except clms CC Vel PATH COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 5];
FIG.Name = 'Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h

% 15
% idx_1 = 10;
% idx_2 = 5;
% idx_35 = 9;
% idx_65 = 3;
% idx_12 = 4;

% 3.75
idx_05 = 2;
idx_1 = 5;
idx_2 = 4;
idx_35 = 70;
idx_65 = 6;
idx_12 = 1;

% 18.75
% idx_1 = 10;
% idx_2 = 102;
% idx_35 = 9;
% idx_65 = 3;
% idx_12 = 4;

idx = idx_2;
trial = SACCADE.saccade{idx};

ax(1) = subplot(4,1,1) ; hold on ; title(['Stimulus: ' num2str(D.freq(idx)) ' (Hz)'])
    ylabel('Position (°)')
    stim = trial.stimlus_position;
    stim = stim - mean(stim);
    h(1) = plot(trial.time, stim, 'Color', [0.0 0.9 0.0]);
    plot(trial.time, trial.position - mean(trial.position),'k');
%     plot(trial.time, zeros(trial.n,1),'--','Color',[0.5 0.5 0.5])
%     if trial.count~=0
%         for ww = 1:trial.count
%            plot(trial.saccades{ww}.Time,trial.saccades{ww}.Position,...
%                'LineWidth', 1, 'Color', trial.cmap(ww,:))
%            plot(trial.intervals{ww}.Time,trial.intervals{ww}.Position,...
%                'LineWidth', 1, 'Color', 0.7*trial.cmap(ww,:))
%         end
%         plot(trial.starts.time , trial.starts.position , '*g')
%         plot(trial.peaks.time  , trial.peaks.position  , '*b')
%         plot(trial.ends.time   , trial.ends.position   , '*r')
%     end

    ax(1).YLim = 20*[-1 1];

ax(2) = subplot(4,1,2) ; hold on
    ylabel('Velocity (°/s)')
    h(2) = plot(trial.time,trial.velocity,'k'); 
    plot(trial.time, -trial.threshold*ones(trial.n,1) , '--m')
    plot(trial.time,  trial.threshold*ones(trial.n,1) , '--m')
    if trial.count~=0
        for ww = 1:trial.count
           plot(trial.saccades{ww}.Time,trial.saccades{ww}.Velocity,...
               'LineWidth', 1, 'Color', trial.cmap(ww,:))
           plot(trial.intervals{ww}.Time,trial.intervals{ww}.Velocity,...
               'LineWidth', 1, 'Color', 0.7*trial.cmap(ww,:))
        end
        plot(trial.starts.time , trial.starts.velocity  , '*g')
        plot(trial.peaks.time  , trial.peaks.velocity   , '*b')
        plot(trial.ends.time   , trial.ends.velocity    , '*r')
    end
    ax(2).YLim = 700*[-1 1];

ax(3) = subplot(4,1,3) ; hold on
    if trial.count~=0
        ylabel({'Removed', 'Position (°)'})                
        plot(trial.shift.Time, trial.shift.IntrpPosition,      	 'r', 'LineWidth', 1);
        plot(trial.shift.Time, trial.shift.Position,             'k', 'LineWidth', 1);
        plot(trial.removed_all.Time, trial.removed_all.Position, 'c', 'LineWidth', 1);
    end
    
ax(4) = subplot(4,1,4) ; hold on
    if trial.count~=0
        ylabel({'Removed', 'Velocity (°/s)'})
        xlabel('Time')
        plot(trial.shift.Time, trial.shift.IntrpVelocity,   'r', 'LineWidth', 0.5);
        plot(trial.shift.Time, trial.shift.Velocity,        'k', 'LineWidth', 0.5);
    end
    ax(4).YLim = ax(2).YLim;
    
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8)
linkaxes(ax,'x')
set(h,'LineWidth',0.5)
align_Ylabels(FIG)

end