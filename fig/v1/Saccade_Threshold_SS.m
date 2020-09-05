function [] = Saccade_Threshold_SS()
%% Saccade_Threshold:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE','HEAD_SACCADE_STATS','U','N')

clearvars -except U N SACCADE HEAD_SACCADE_STATS FILE

filedata = textscan(FILE, '%s', 'delimiter', '_=');
amp = filedata{1}{6};

%% Velocity Histogram
CC = hsv(N.freq);
Freq = U.freq{1};

clear head Head
head.all = cell(N.fly,N.freq);
for n = 1:N.file
    head.all{SACCADE.fly(n),SACCADE.freq(n)}(end+1,1) = SACCADE.head_saccade{n};
end

head.fly = cellfun(@(x) cat(2,x.velocity), head.all, 'UniformOutput', false);
head.vel = cell(N.freq,1);
for v = 1:N.freq
    for f = 1:N.fly
        head.vel{v} = cat(2, head.vel{v}, head.fly{f,v});
    end
end
head.vel_all = cat(2,head.vel{:});

Head.vel = head.vel;
Head.vel_fly = head.fly;

head.vel_mean = cellfun(@(x) mean(x,'all'), head.vel, 'UniformOutput', true);
head.vel_median = cellfun(@(x) median(x,'all'), head.vel, 'UniformOutput', true);
head.vel_std = cellfun(@(x) std(x,[],'all'), head.vel, 'UniformOutput', true);

fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 7 5])
ax = gobjects(N.freq,1);
bins = -1000:5:1000;
pp = 1;
for v = 1:N.freq
    ax(v) = subplot(N.freq/2,2,pp); hold on ; title([ num2str(Freq(v)) ' Hz'])
        h = histogram(head.vel{v}, bins, 'Normalization', 'probability', ...
            'FaceColor', CC(v,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
        xlabel('Velocity (°/s)')
        ylabel('Probability')
        axis tight
	pp = pp + 1;
end
linkaxes(ax,'xy')
set(ax, 'LineWidth', 1.5, 'Box', 'off', 'XLim', 1000*[-1 1], ...
    'XTick', -1000:500:1000, 'YLim', [-0.005 0.03])

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [3 3 3 3])
clear ax
ax(1) = subplot(1,1,1); hold on
    h = histogram(head.vel_all, bins, 'Normalization', 'probability', ...
        'FaceColor', 'k', 'FaceAlpha', 1, 'EdgeColor', 'none');
    xlabel('Head Velocity (°/s)')
    ylabel('Probability')
    
thresh = 3*std(head.vel_all(:));
plot(thresh*[-1 -1],[0 ax.YLim(2)], 'm', 'LineWidth', 1)
plot(thresh*[ 1  1],[0 ax.YLim(2)], 'm', 'LineWidth', 1)
set(ax, 'LineWidth', 1, 'Box', 'off', 'XLim', 1000*[-1 1], 'YLim', [-0.005 ax.YLim(2)])

%% Position Histogram
clear head
head.all = cell(N.fly,N.freq);
for n = 1:N.file
    head.all{SACCADE.fly(n),SACCADE.freq(n)}(end+1,1) = SACCADE.head_saccade{n};
end

head.fly = cellfun(@(x) cat(2,x.position), head.all, 'UniformOutput', false);
head.pos = cell(N.freq,1);
for v = 1:N.freq
    for f = 1:N.fly
        head.pos{v} = cat(2, head.pos{v}, head.fly{f,v});
    end
end
head.pos_all = cat(2,head.pos{:});

Head.pos = head.pos;
Head.pos_fly = head.fly;

head.pos_mean = cellfun(@(x) mean(x,'all'), head.pos, 'UniformOutput', true);
head.pos_median = cellfun(@(x) median(x,'all'), head.pos, 'UniformOutput', true);
head.pos_std = cellfun(@(x) std(x,[],'all'), head.pos, 'UniformOutput', true);

fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 7 5])
ax = gobjects(N.freq,1);
bins = -25:1:25;
pp = 1;
for v = 1:N.freq
    ax(v) = subplot(N.freq/2,2,pp); hold on ; title([ num2str(Freq(v)) '(°/s)'])
        h = histogram(head.pos{v} - mean(head.pos{v}), bins, 'Normalization', 'probability', ...
            'FaceColor', CC(v,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
        xlabel('Position (°)')
        ylabel('Probability')
        axis tight
	pp = pp + 1;
end
linkaxes(ax,'xy')
set(ax, 'LineWidth', 1, 'Box', 'off', 'XLim', 25*[-1 1], 'XTick', -25:5:25, 'YLim', [-0.005 0.15])

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [3 3 3 3])
clear ax

% start_pos = -HEAD_SACCADE_STATS.StartPos.*HEAD_SACCADE_STATS.Direction;
% start_pos = start_pos(~isnan(start_pos));
% end_pos = -HEAD_SACCADE_STATS.EndPos.*HEAD_SACCADE_STATS.Direction;
% end_pos = end_pos(~isnan(end_pos));

ax(1) = subplot(1,1,1); hold on
    h(1) = histogram(head.pos_all, bins, 'Normalization', 'probability', ...
        'FaceColor', 'k', 'FaceAlpha', 1, 'EdgeColor', 'none');
%     for s = [1,6]
%         histogram(head.pos{s}, bins, 'Normalization', 'probability', ...
%             'FaceColor', CC(s,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%     h(2) = histogram(start_pos, bins, 'Normalization', 'probability', ...
%         'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
%     h(3) = histogram(end_pos, bins, 'Normalization', 'probability', ...
%         'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
    xlabel('Head Position (°)')
    ylabel('Probability')
    
set(ax, 'LineWidth', 1.5, 'FontSize', 8, 'Box', 'on', 'XLim', 25*[-1 1], 'YLim', [-0.005 ax.YLim(2)])

%% Save
amp = str2double(amp);
fname = ['SS_amp=' amp];
savedir = 'C:\Users\BC\Box\Research\Manuscripts\Head Saccade\Data';
save(fullfile(savedir, [fname '.mat']), 'Head', 'amp', 'U', 'N', 'HEAD_SACCADE_STATS');

end