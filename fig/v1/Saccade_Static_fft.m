function [] = Saccade_Static_fft()
%% Saccade_Static_fft:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE','D','I','U','N')

% clearvars -except SACCADE  D I U N

%% FFT of all trials
clc
n_trial = size(SACCADE,1);
ALL = cell(N.wave,1);
for n = 1:n_trial
    t = SACCADE.head_saccade{n}.time;
    x = SACCADE.head_saccade{n}.position;
    [Fv, Mag , ~ , ~] = FFT(t, x);
    ALL{I.wave(n)}(:,end+1) = Mag;
end

%% 
fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 7 2])
ax = gobjects(N.wave,1);
pp = 1;
[b, a] = butter(3, 0.1, 'low');
% ALL_filt = cellfun(@(x) filtfilt(b, a, x), ALL, 'UniformOutput', false);
for w = [2 3 4 6 1 5]
    ax(w) = subplot(1,N.wave,pp); cla ; hold on ; title(num2str(U.wave{1}(w)))
        mag_mean = mean(ALL{w},2);
        mag_mean_filt = filtfilt(b, a, mag_mean);
        [pks,locs] = findpeaks(mag_mean_filt, 'MinPeakDistance', 5, 'MinPeakHeight', 0.05, ...
                                        'MinPeakProminence', 0.05, 'MinPeakWidth', 20, 'SortStr', 'descend');
%         freq_peak = Fv(locs(1));
%         mag_peak = pks(1);
        
        plot(Fv, ALL{w}, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.25)
        %plot(Fv, mag_mean, 'k', 'LineWidth', 1)
        plot(Fv, mag_mean_filt, 'b', 'LineWidth', 1)
        %plot([freq_peak freq_peak], [0 mag_peak], 'r', 'LineWidth', 1)
        
    pp = pp + 1;
end

set(ax , 'LineWidth', 1, 'XLim', [-5 40], 'YLim', [-0.05 1])
set(ax(2:end), 'YTickLabels', [], 'XTickLabels', [])
linkaxes(ax)

end