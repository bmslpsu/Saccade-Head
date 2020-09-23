function [] = Sine_example_trial()
%% Sine_example_trial:
root = 'H:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'SACCADE','D','I','U','N')

%% Example Trial
% Amp = 3.75
% idx_05 = 2;
idx_1 = 5;
% idx_2 = 4;
% idx_35 = 70;
% idx_65 = 6;
% idx_12 = 1;

% Amp = 15
% idx_1 = 10;
% idx_2 = 5;
% idx_35 = 9;
% idx_65 = 3;
% idx_12 = 4;

% Amp = 18.75
% idx_1 = 46;
% idx_2 = 102;
% idx_35 = 9;
% idx_65 = 3;
% idx_12 = 4;

idx = idx_1;
trial = SACCADE.head_saccade{idx};

OR = mean(abs(trial.position - mean(trial.position)));

fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 6 2.5])
movegui(fig,'center')
clear ax h
ax(1) = subplot(2,1,1) ; hold on ; title(['Stimulus: ' num2str(D.freq(idx)) ' (Hz)'])
    ylabel('Position (°)')
    stim = 0.92*trial.stimlus_position;
    stim = stim - mean(stim);
    plot(trial.time, stim, 'Color', [0.0 0.9 0.0])
    plot(trial.time, trial.position - 0*mean(trial.position),'k', 'LineWidth', 0.5)
    if trial.count~=0
        for ww = 1:trial.count
           plot(trial.saccades{ww}.Time,trial.saccades{ww}.Position,...
               'LineWidth', 1, 'Color', 'c')
           %plot(trial.intervals{ww}.Time,trial.intervals{ww}.Position,...
               %'LineWidth', 1, 'Color', 0.7*trial.cmap(ww,:))
        end
        %plot(trial.starts.time , trial.starts.position , '*g')
        %plot(trial.peaks.time  , trial.peaks.position  , '*b')
        %plot(trial.ends.time   , trial.ends.position   , '*r')
    end
    plot([0 10],  OR*[1 1], '--r')
    plot([0 10], -OR*[1 1], '--r')
    
    ax(1).YLim = 20*[-1 1];

ax(2) = subplot(2,1,2) ; hold on
    ylabel('Velocity (°/s)')
    plot(trial.time,trial.velocity,'k', 'LineWidth', 0.5)
    %plot(trial.time, trial.threshold(1)*ones(trial.n,1) , '--m')
    %plot(trial.time, trial.threshold(2)*ones(trial.n,1) , '--m')
    if trial.count~=0
        for ww = 1:trial.count
           plot(trial.saccades{ww}.Time,trial.saccades{ww}.Velocity,...
               'LineWidth', 1, 'Color', 'c')
           %plot(trial.intervals{ww}.Time,trial.intervals{ww}.Velocity,...
               %'LineWidth', 1, 'Color', 0.7*trial.cmap(ww,:))
        end
        %plot(trial.starts.time , trial.starts.velocity  , '*g')
        %plot(trial.peaks.time  , trial.peaks.velocity   , '*b')
        %plot(trial.ends.time   , trial.ends.velocity    , '*r')
    end
    ax(2).YLim = 800*[-1 1];

set(ax, 'LineWidth', 1)
linkaxes(ax,'x')

end