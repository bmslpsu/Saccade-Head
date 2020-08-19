function [] = Head_Dist_SS()
%% Head_Dist_SS:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N')

%% Velocity Histogram
CC = hsv(N.freq);
Freq = U.freq{1};

clear head
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

head.speed = head.vel;
head.vel_all = cat(2,head.speed{:});

head.vel_mean = cellfun(@(x) mean(x,'all'), head.vel, 'UniformOutput', true);
head.vel_median = cellfun(@(x) median(x,'all'), head.vel, 'UniformOutput', true);
head.vel_std = cellfun(@(x) std(x,[],'all'), head.vel, 'UniformOutput', true);

fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 4 6])
ax = gobjects(N.freq,1);
bins = -1000:5:1000;
% bins = [-1000:50:-350 -350:5:350 350:50:1000];
pp = 1;
for v = 1:N.freq
    ax(v) = subplot(N.freq/2,2,pp); hold on ; title([ num2str(Freq(v)) '(°/s)'])
        h = histogram(head.vel{v}, bins, 'Normalization', 'probability', ...
            'FaceColor', CC(v,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
        xlabel('Velocity (°/s)')
        ylabel('Probability')
        axis tight
	pp = pp + 1;
end
linkaxes(ax,'xy')
set(ax, 'LineWidth', 1.5, 'Box', 'on', 'XLim', 1000*[-1 1], 'YLim', [-0.005 0.06])

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [3 3 3 3])
clear ax
ax(1) = subplot(1,1,1); hold on
    h = histogram(head.vel_all, bins, 'Normalization', 'probability', ...
        'FaceColor', 'k', 'FaceAlpha', 1, 'EdgeColor', 'none');
    xlabel('Head Velocity (°/s)')
    ylabel('Probability')
    
thresh = 350;
plot(thresh*[-1 -1],[0 ax.YLim(2)], 'm', 'LineWidth', 1)
plot(thresh*[ 1  1],[0 ax.YLim(2)], 'm', 'LineWidth', 1)
set(ax, 'LineWidth', 1.5, 'FontSize', 8, 'Box', 'on', 'XLim', 1000*[-1 1], 'YLim', [-0.005 ax.YLim(2)])

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

head.speed = head.pos;
head.pos_all = cat(2,head.speed{:});

head.pos_mean = cellfun(@(x) mean(x,'all'), head.pos, 'UniformOutput', true);
head.pos_median = cellfun(@(x) median(x,'all'), head.pos, 'UniformOutput', true);
head.pos_std = cellfun(@(x) std(x,[],'all'), head.pos, 'UniformOutput', true);

fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 4 6])
ax = gobjects(N.freq,1);
bins = -25:0.5:25;
pp = 1;
for v = 1:N.freq
    ax(v) = subplot(N.freq/2,2,pp); hold on ; title([ num2str(Freq(v)) '(°/s)'])
        h = histogram(head.pos{v} - 0*mean(head.pos{v}), bins, 'Normalization', 'probability', ...
            'FaceColor', CC(v,:), 'FaceAlpha', 1, 'EdgeColor', 'none');
        xlabel('Position (°)')
        ylabel('Probability')
        axis tight
	pp = pp + 1;
end
linkaxes(ax,'xy')
set(ax, 'LineWidth', 1.5, 'Box', 'on', 'XLim', 25*[-1 1], 'YLim', [-0.005 ax(1).YLim(2)])

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [3 3 3 3])
clear ax
ax(1) = subplot(1,1,1); hold on
    h = histogram(head.pos_all - 0*mean(head.pos_all), bins, 'Normalization', 'probability', ...
        'FaceColor', 'k', 'FaceAlpha', 1, 'EdgeColor', 'none');
    xlabel('Head Position (°)')
    ylabel('Probability')
    
set(ax, 'LineWidth', 1.5, 'FontSize', 8, 'Box', 'on', 'XLim', 25*[-1 1], 'YLim', [-0.005 ax.YLim(2)])
ylim([0 0.05])

end