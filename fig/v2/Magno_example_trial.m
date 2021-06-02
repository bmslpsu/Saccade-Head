function [] = Magno_example_trial()
%% Magno_example_trial:
%   INPUTS:
%       -
%   OUTPUTS:
%       -
%

root = 'E:\EXPERIMENTS\MAGNO\Experiment_Ramp';
[FILE,PATH] = uigetfile({'*.mat'}, 'Select file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'SACCADE','D','I','U','N')

%%

fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 6 2.5])
movegui(fig,'center')
clear ax h
ax(1) = subplot(2,1,1) ; hold on ; title(['Stimulus: ' num2str(D.freq(idx)) ' (Hz)'])
    ylabel('Position (°)')













end