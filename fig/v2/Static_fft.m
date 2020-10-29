function [] = Static_fft()
%% Static_fft:
root = 'H:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE','D','I','U','N')

%% FFT of all wavelengths by fly
warning('off', 'signal:findpeaks:largeMinPeakHeight')
clearvars -except SACCADE D I U N root
clc

keepI = cellfun(@(x) isobject(x) || isstruct(x), SACCADE.head_saccade);
Saccade = SACCADE(keepI,:);

n_trial = size(Saccade,1);
freq_wave_fly = cell(N.fly,N.wave);
for n = 1:n_trial
    t = Saccade.head_saccade{n}.time;
    x = Saccade.head_saccade{n}.position;
    %x = Saccade.wing_saccade{n}.position;
    [Fv, Mag , ~ , ~] = FFT(t, x);
    freq_wave_fly{Saccade.fly(n),Saccade.wave(n)}(:,end+1) = Mag;
end
emptyI = cellfun(@(x) isempty(x), freq_wave_fly, 'UniformOutput', true);
freq_wave_fly(emptyI) = {nan(1001,1)};

temp = cellfun(@(x) nanmean(x,2), freq_wave_fly, 'UniformOutput', false);
freq_wave = cell(1,N.wave);
freq_wave_fly_mean = cell(1,N.wave);
for w = 1:N.wave
    freq_wave{w} = cat(2, freq_wave_fly{:,w});
    freq_wave_fly_mean{w} = cat(2, temp{:,w});
end
freq_wave_fly_mean_gran = cellfun(@(x) nanmean(x,2), freq_wave_fly_mean, 'UniformOutput', false);
freq_wave_fly_mean_gran = cat(2, freq_wave_fly_mean_gran{:});

fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 7 2])
ax = gobjects(N.wave,1);
pp = 1;
[b, a] = butter(3, 0.1, 'low');
freq_peak = nan(1,N.wave);
mag_peak = nan(1,N.wave);
peak_range = 51:251;
for w = [2 3 4 6 1 5]
    ax(pp) = subplot(1,N.wave,pp); cla ; hold on ; title(num2str(U.wave{1}(w)))
        mag_mean = freq_wave_fly_mean_gran(:,w);
        %mag_mean = nanmean(freq_wave{w},2);
        mag_mean_filt = filtfilt(b, a, mag_mean);
        [pks,locs,~,~] = findpeaks(mag_mean_filt(peak_range), 'MinPeakDistance', 5, 'MinPeakHeight', 0.05, ...
                                        'MinPeakProminence', 0.01, 'MinPeakWidth', 20, 'SortStr', 'descend');
        
        plot(Fv, freq_wave{w}, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.25)
        %plot(Fv, freq_wave_fly_mean{w}, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.25)
        plot(Fv, mag_mean, 'k', 'LineWidth', 0.75)
        %plot(Fv, mag_mean_filt, 'b', 'LineWidth', 1)
        if ~isempty(locs)
            freq_peak(w) = Fv(locs(1) + peak_range(1));
            mag_peak(w) = pks(1);
            plot([freq_peak(w)  freq_peak(w) ], [0 mag_peak(w)], 'r', 'LineWidth', 1)
        end
        xline(11,'--r')
        xline(18,'--r')
        
    pp = pp + 1;
end

set(ax , 'LineWidth', 1, 'XLim', [-5 40], 'YLim', [-0.05 1])
set(ax(2:end), 'YTickLabels', [], 'XTickLabels', [])
linkaxes(ax)

%% Magnitude in time domain
clc
Fs = round(SACCADE.head_saccade{1}.Fs);
Fc = [11 18];
[b, a] = butter(3, Fc/(Fs/2), 'bandpass');

n_trial = size(Saccade,1);
mag_all = nan(n_trial,1);
pos_all = cell(1,N.wave);
for n = 1:n_trial
%     x = Saccade.head_saccade{n}.shift;
%     if isempty(x)
%         x = Saccade.head_saccade{n}.position;
%     else
%         x = x.IntrpPosition;
%     end
    x = Saccade.head_saccade{n}.position;
    x = filtfilt(b, a, x);
    mag_all(n) = 3*std(x);
    pos_all{I.wave(n)}(:,end+1) = x;
end
% mag_wave = cellfun(@(x) 3*std( filtfilt(b, a, x.position), 1)', ...
%     SACCADE.head_saccade, 'UniformOutput', true);

% Split into uniform & non-uniform wavelengths
wave_group_split = Saccade.wave;
wave_group_split(any(Saccade.wave==[1 5],2)) = 1;
wave_group_split(all(Saccade.wave~=[1 5],2)) = 2;

fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 2 2])
clear ax
ax(1) = subplot(1,1,1); hold on
set(ax , 'LineWidth', 1, 'XLim', [-5 40], 'YLim', [-0.05 0.5])

% bx = boxplot(Mag_all, Saccade.wave, 'Labels', U.wave, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
bx = boxplot(mag_all, wave_group_split, 'Labels', {'Uniform', 'Grating'}, ...
    'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
xlabel('Wavelength (°)')
ylabel('Amplitude (°)')

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1}, h{kk,2}, 'k', 'EdgeColor', 'none');
end

ww = 1;
set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth', 1);
set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(ww).Children = ax(ww).Children([end 1:end-1]);
ylim([-0.1 4.2])

%% Anova
% [p,tbl,stats] = anova1(mag_wave,wave_group_split);
[p,tb,stats] = anovan(mag_all, {wave_group_split Saccade.fly}, 'model','interaction', ...
    'varnames', {'Wave','Fly'});
[c,m] = multcompare(stats);

%% For all wavlengths
[p,tb,stats] = anovan(mag_all, {Saccade.wave}, 'model','interaction', ...
    'varnames', {'Wave'});
[c,m] = multcompare(stats, 'Alpha', 0.05);

%% Kruskalwallis
[p,tb,stats] = kruskalwallis(mag_all, Saccade.wave);
[c,m] = multcompare(stats);

%% Save magnitude data
fname = 'Static_freq_mag';
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'Fs', 'Fc', 'wave_group_split', 'mag_wave', 'I', 'U', 'N');

end