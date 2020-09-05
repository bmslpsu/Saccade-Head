function [] = SS_Hist()
%% SS_Hist:
root = 'C:\Users\BC\Box\Research\Manuscripts\Head Saccade\Data';

[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','on');
FILE = cellstr(FILE);

n_amp = length(FILE);
amp_data = cell(n_amp,1);
Amp = nan(n_amp,1);
for a = 1:n_amp
    amp_data{a} = load(fullfile(PATH,FILE{a}));
    Amp(a) = unique(amp_data{a}.HEAD_SACCADE_STATS.amp);
end
n_freq = amp_data{1}.N.freq;
Freq = amp_data{1}.U.freq{1};
freq_color = jet(n_freq);
amp_color = parula(n_amp);

%% Position All
% fig = figure(1); clf
% set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 8 5])
% ax = gobjects(n_amp,n_freq);
% H = gobjects(n_amp,n_freq);
% pp = 1;
% edges = -25.25:0.5:25.25;
% for a = 1:n_amp
%     for f = 1:n_freq
%         ax(a,f) = subplot(n_amp, n_freq, pp); cla ; hold on
%         data = amp_data{a}.Head.pos{f};
%         data = data - mean(data);
%         H(a,f) = histogram(data, edges);
%         
%         if pp <= n_freq
%             title([num2str(Freq(f)) ' Hz'])
%         end
%         
%         if ~mod(pp-1,n_freq)
%             ylabel([num2str(Amp(a)) '°'], 'FontWeight', 'bold')
%         end
%         
%         pp = pp + 1;
%     end
% end
% set(H, 'Normalization', 'Probability')
% set(H, 'EdgeColor', 'none', 'FaceAlpha', 1, 'FaceColor', 'k')
% linkaxes(ax, 'xy')
% 
% set(ax , 'LineWidth', 1, 'FontSize', 8)
% set(ax , 'YLim', [-0.005 ax(1).YLim(2)])
% set(ax , 'XTick', -25:5:25)
% 
% set(ax(1:n_amp-1,:), 'XTickLabels', [])
% set(ax(:,2:n_freq), 'YTickLabels', [])

%% Position All
fig = figure(1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 8 7])
ax = gobjects(n_amp,n_freq);
H = gobjects(n_amp,n_freq);
pp = 1;
edges = -25.25:0.5:25.25;
for f = 1:n_freq
    for a = 1:n_amp
        ax(a,f) = subplot(n_freq, n_amp, pp); cla ; hold on
        data = amp_data{a}.Head.pos{f};
        data = data - mean(data);
        H(a,f) = histogram(data, edges);
        
        if pp <= n_amp
            title([num2str(Amp(a)) '°'])
        end
        
        if ~mod(pp-1,n_amp)
            ylabel([num2str(Freq(f)) ' Hz'], 'FontWeight', 'bold')
        end
        
        pp = pp + 1;
    end
end
set(H, 'Normalization', 'Probability')
set(H, 'EdgeColor', 'none', 'FaceAlpha', 1, 'FaceColor', 'k')
linkaxes(ax, 'xy')

set(ax , 'LineWidth', 0.75, 'FontSize', 8, 'Color', 'none')
set(ax , 'YLim', [-0.00 ax(1).YLim(2)])
set(ax , 'XLim', 25*[-1 1], 'XTick', -25:5:25)

set(ax(:,1:n_freq-1), 'XTickLabels', [], 'XColor', 'none')
set(ax(2:n_amp,:), 'YTickLabels', [],  'YColor', 'none')

%% Velocity All
fig = figure(2); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 8 7])
ax = gobjects(n_amp,n_freq);
H = gobjects(n_amp,n_freq);
pp = 1;
edges = -1000.25:5:1000.25;
for f = 1:n_freq
    for a = 1:n_amp
        ax(a,f) = subplot(n_freq, n_amp, pp); cla ; hold on
        data = amp_data{a}.Head.vel{f};
        H(a,f) = histogram(data, edges);
        
        if pp <= n_amp
            title([num2str(Amp(a)) '°'])
        end
        
        if ~mod(pp-1,n_amp)
            ylabel([num2str(Freq(f)) ' Hz'], 'FontWeight', 'bold')
        end
        
        pp = pp + 1;
    end
end
set(H, 'Normalization', 'Probability')
set(H, 'EdgeColor', 'none', 'FaceAlpha', 1, 'FaceColor', 'k')
linkaxes(ax, 'xy')

set(ax , 'LineWidth', 0.75, 'FontSize', 8, 'Color', 'none')
set(ax , 'YLim', [-0.00 ax(1).YLim(2)])
set(ax , 'XLim', 500*[-1 1], 'XTick', -1000:500:1000)

set(ax(:,1:n_freq-1), 'XTickLabels', [], 'XColor', 'none')
set(ax(2:n_amp,:), 'YTickLabels', [],  'YColor', 'none')

end