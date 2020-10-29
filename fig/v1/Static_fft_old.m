function [] = Static_fft()
%% Static_fft:
root = 'H:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE','D','I','U','N')

%% FFT of all wavelengths
clc
Fs = round(SACCADE.head_saccade{1}.Fs);
[b, a] = butter(10, 2/(Fs/2), 'high');
n_trial = size(SACCADE,1);
ALL = cell(1,N.wave);
for n = 1:n_trial
    t = SACCADE.head_saccade{n}.time;
    x = SACCADE.head_saccade{n}.shift;
    if isempty(x)
        x = SACCADE.head_saccade{n}.position;
    else
        x = x.IntrpPosition;
    end
    %x = filtfilt(b, a, x);
    [Fv, Mag , ~ , ~] = FFT(t, x);
    ALL{I.wave(n)}(:,end+1) = Mag;
end

fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 7 2])
ax = gobjects(N.wave,1);
pp = 1;
[b, a] = butter(3, 0.1, 'low');
% ALL_filt = cellfun(@(x) filtfilt(b, a, x), ALL, 'UniformOutput', false);
freq_peak = nan(1,N.wave);
mag_peak = nan(1,N.wave);
peak_range = 51:251;
for w = [2 3 4 6 1 5]
    ax(pp) = subplot(1,N.wave,pp); cla ; hold on ; title(num2str(U.wave{1}(w)))
        mag_mean = mean(ALL{w},2);
        mag_mean_filt = filtfilt(b, a, mag_mean);
        [pks,locs,~,~] = findpeaks(mag_mean_filt(peak_range), 'MinPeakDistance', 5, 'MinPeakHeight', 0.05, ...
                                        'MinPeakProminence', 0.01, 'MinPeakWidth', 20, 'SortStr', 'descend');
        
        plot(Fv, ALL{w}, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.25)
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

%% Magnitude of head movements
peak_width = 111:181;
Fv_peak = Fv(peak_width);
MAG = cellfun(@(x) mean(sum(x(peak_width,:),1)), ALL, 'UniformOutput', true);
MAG = nanmean(MAG,1);
MAG = mean(MAG([2:4,6])) - mean(MAG([1,5]));

fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 2 2])
clear ax
ax(1) = subplot(1,1,1); hold on ; title([num2str(MAG) '°'])
set(ax , 'LineWidth', 1, 'XLim', [-5 40], 'YLim', [-0.05 1])
clear mag_mean
for w = [2 3 4 6 1 5]
    %mag_mean = GRAND(:,w);
    mag_mean(:,w) = mean(ALL{w},2);
    %plot(Fv, mag_mean, 'k', 'LineWidth', 0.75)
end
mag_peak_mean = mean(mag_mean(:,[2:4,6]), 2);
mag_noise_mean = mean(mag_mean(:,[1,5]), 2);
mag_peak_only = mag_peak_mean(peak_width,:);
mag_noise_only = mag_noise_mean(peak_width,:);

plot(Fv, mag_noise_mean, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.75)
plot(Fv, mag_peak_mean, 'k', 'LineWidth', 0.75)

% plot(Fv_peak, mag_peak_only, 'r', 'LineWidth', 0.75)
% plot(Fv_peak, mag_noise_only, 'r', 'LineWidth', 0.75)
xx = [Fv_peak ; flipud(Fv_peak)];
yy = [mag_noise_only ; flipud(mag_peak_only)];
patch(xx, yy, 'r', 'EdgeColor', 'none')

xline(11,'--r')
xline(18,'--r')

%% FFT of all wavelengths by fly
clc
[b, a] = butter(10, 2/(Fs/2), 'high');
n_trial = size(SACCADE,1);
ALL = cell(N.fly,N.wave);
for n = 1:n_trial
    t = SACCADE.head_saccade{n}.time;
    x = SACCADE.head_saccade{n}.shift;
    if isempty(x)
        x = SACCADE.head_saccade{n}.position;
    else
        x = x.IntrpPosition;
    end
    %x = filtfilt(b, a, x);
    [Fv, Mag , ~ , ~] = FFT(t, x);
    ALL{I.fly(n),I.wave(n)}(:,end+1) = Mag;
end
emptyI = cellfun(@(x) isempty(x), ALL, 'UniformOutput', true);
ALL{emptyI} = nan(1001,1);

FLY = cellfun(@(x) nanmean(x,2), ALL, 'UniformOutput', false);
FLY_ALL = cell(1,N.wave);
for w = 1:N.wave
    FLY_ALL{w} = cat(2, FLY{:,w});
end
GRAND = cellfun(@(x) nanmean(x,2), FLY_ALL, 'UniformOutput', false);
GRAND = cat(2, GRAND{:});

fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 7 2])
ax = gobjects(N.wave,1);
pp = 1;
freq_peak = nan(1,N.wave);
mag_peak = nan(1,N.wave);
peak_range = 51:251;
for w = [2 3 4 6 1 5]
    ax(pp) = subplot(1,N.wave,pp); cla ; hold on ; title(num2str(U.wave{1}(w)))
        mag_mean = GRAND(:,w);
        mag_mean_filt = filtfilt(b, a, mag_mean);
        [pks,locs,~,~] = findpeaks(mag_mean_filt(peak_range), 'MinPeakDistance', 5, 'MinPeakHeight', 0.05, ...
                                        'MinPeakProminence', 0.01, 'MinPeakWidth', 20, 'SortStr', 'descend');
        
        %plot(Fv, ALL{w}, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.25)
        plot(Fv, FLY_ALL{w}, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.25)
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

%% Magnitude of head movements by fly
peak_width = 111:181;
Fv_peak = Fv(peak_width);
MAG = cellfun(@(x) mean(sum(x(peak_width,:),1)), ALL, 'UniformOutput', true);
MAG = nanmean(MAG,1);
MAG = mean(MAG([2:4,6])) - mean(MAG([1,5]));

fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 2 2])
clear ax
ax(1) = subplot(1,1,1); hold on ; title([num2str(MAG) '°'])
set(ax , 'LineWidth', 1, 'XLim', [-5 40], 'YLim', [-0.05 0.5])
clear mag_mean
for w = [2 3 4 6 1 5]
    %mag_mean = GRAND(:,w);
    mag_mean(:,w) = mean(FLY_ALL{w},2);
    %plot(Fv, mag_mean, 'k', 'LineWidth', 0.75)
end
mag_peak_mean = nanmean(mag_mean(:,[2:4,6]), 2);
mag_noise_mean = nanmean(mag_mean(:,[1,5]), 2);
mag_peak_only = mag_peak_mean(peak_width,:);
mag_noise_only = mag_noise_mean(peak_width,:);

plot(Fv, mag_noise_mean, 'Color', [0.5 0.5 0.5], 'LineWidth', 0.75)
plot(Fv, mag_peak_mean, 'k', 'LineWidth', 0.75)

% plot(Fv_peak, mag_peak_only, 'r', 'LineWidth', 0.75)
% plot(Fv_peak, mag_noise_only, 'r', 'LineWidth', 0.75)
xx = [Fv_peak ; flipud(Fv_peak)];
yy = [mag_noise_only ; flipud(mag_peak_only)];
patch(xx, yy, 'r', 'EdgeColor', 'none')

xline(11,'--r')
xline(18,'--r')

%% Magnitude in time domain
clc
Fs = round(SACCADE.head_saccade{1}.Fs);
Fc = [11 18];
[b, a] = butter(3, Fc/(Fs/2), 'bandpass');
n_trial = size(SACCADE,1);
ALL = cell(1,N.wave);
for n = 1:n_trial
    t = SACCADE.head_saccade{n}.time;
    x = SACCADE.head_saccade{n}.shift;
    if isempty(x)
        x = SACCADE.head_saccade{n}.position;
    else
        x = x.IntrpPosition;
    end
    x = filtfilt(b, a, x);
    %[Fv, Mag , ~ , ~] = FFT(t, x);
    ALL{I.wave(n)}(:,end+1) = x;
end
Mag_all = cellfun(@(x) 3*std(x,1)', ALL, 'UniformOutput', false);
wave_group = cellfun(@(x,y) y*ones(length(x),1), Mag_all, num2cell(1:length(Mag_all)), 'UniformOutput', false);
wave_group = cat(1, wave_group{:});
Mag_all = cat(1, Mag_all{:});

% Split into uniform & non-uniform wavelengths
wave_group_split = wave_group;
wave_group_split(any(wave_group==[1 5],2)) = 1;
wave_group_split(all(wave_group~=[1 5],2)) = 2;

fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 2 2])
clear ax
ax(1) = subplot(1,1,1); hold on
%title([num2str(MAG) '°'])
set(ax , 'LineWidth', 1, 'XLim', [-5 40], 'YLim', [-0.05 0.5])

% bx = boxplot(Mag_all, wave_group, 'Labels', U.wave, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
bx = boxplot(Mag_all, wave_group_split, 'Labels', {'Uniform', 'Grating'}, ...
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
anova1(Mag_all,wave_group_split)

%% Save magnitude data
fname = 'Static_freq_mag';
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'Fs', 'Fc', 'Mag_all', 'wave_group', ...
                                        'wave_group_split', 'U', 'N');

end