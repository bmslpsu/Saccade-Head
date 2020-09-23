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
fig = figure(1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 (9/5)*n_amp (8/6)*n_freq])
ax = gobjects(n_amp,n_freq);
H = gobjects(n_amp,n_freq);
pp = 1;
edges = -25.25:0.5:25.25;
clear data_stats
head_pos = nan(n_freq,n_amp);
scd_rate = nan(n_freq,n_amp);
for f = 1:n_freq
    for a = 1:n_amp
        ax(a,f) = subplot(n_freq, n_amp, pp); cla ; hold on
        pos = amp_data{a}.Head.pos{f};
        pos = pos - mean(pos);
        
        data_stats(f,a) = struct(basic_stats(pos(:),1));
        head_pos(f,a) = mean(abs(pos(:)));
        scd_rate(f,a) = mean(amp_data{a}.Head.count_all{f}) ./ 10;
        
        mean_speed = 4*Amp(a)*Freq(f);
        xlabel({2*data_stats(f,a).std, scd_rate(f,a)})
        
        H(a,f) = histogram(pos, edges);
        %xline( data_stats(f,a).std, '--g');
        %xline(-data_stats(f,a).std, '--g');
        xline(head_pos(f,a), '--g');
        xline(-head_pos(f,a), '--g');
        
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

set(ax(:,1:n_freq-1), 'XTick', [], 'XTickLabels', [], 'XColor', 'k')
set(ax(2:n_amp,:), 'YTickLabels', [],  'YColor', 'none')

%% Correlation OR
[P] = polyfit(2*head_pos(:), scd_rate(:), 1);
[rho,pval] = corr(2*head_pos(:), scd_rate(:));

fig = figure(10); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 3])
clear ax H
ax(1) = subplot(1,1,1) ; cla ; hold on
% xlim([2 8])
ylim([0 0.8])
cc_amp = jet(n_amp);
H = gobjects(n_freq, n_amp);
mrk = ['o','+','x','s','d','*'];
for a = 1:n_amp
    for f = 1:n_freq
        H(f,a) = scatter(2*head_pos(f,a), scd_rate(f,a), 100, '.', ...
            'MarkerEdgeColor', cc_amp(a,:), 'MarkerFaceColor', 'none');
    end
end
x_fit = ax.XLim;
y_fit = P(2) + P(1)*ax.XLim;
plot(x_fit, y_fit, 'k', 'LineWidth', 1)

legend(H(1,:) , string(Amp), 'Box', 'off')
xlabel('Operating Range (°)')
ylabel('Saccade Frequency (Hz)')
title(['R^{2} = ' num2str(rho) '   p = ' num2str(pval)])
set(ax , 'LineWidth', 0.75, 'FontSize', 8, 'Color', 'none')

%% Correlation Vel
mean_vel = (4*Amp*Freq')';
[P] = polyfit(mean_vel(:), scd_rate(:), 1);
[rho,pval] = corr(mean_vel(:), scd_rate(:));

fig = figure(10); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 3])
clear ax H
ax(1) = subplot(1,1,1) ; cla ; hold on
% xlim([2 8])
ylim([0 0.8])
cc_amp = jet(n_amp);
H = gobjects(n_freq, n_amp);
mrk = ['o','+','x','s','d','*'];
for a = 1:n_amp
    for f = 1:n_freq
        H(f,a) = scatter(mean_vel(f,a), scd_rate(f,a), 100, '.', ...
            'MarkerEdgeColor', cc_amp(a,:), 'MarkerFaceColor', 'none');
    end
end
x_fit = ax.XLim;
y_fit = P(2) + P(1)*ax.XLim;
plot(x_fit, y_fit, 'k', 'LineWidth', 1)

legend(H(1,:) , string(Amp), 'Box', 'off')
xlabel('Operating Range (°)')
ylabel('Saccade Frequency (Hz)')
title(['R^{2} = ' num2str(rho) '   p = ' num2str(pval)])
set(ax , 'LineWidth', 0.75, 'FontSize', 8, 'Color', 'none')


%% Velocity All
fig = figure(2); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 8 7])
ax = gobjects(n_amp,n_freq);
H = gobjects(n_amp,n_freq);
pp = 1;
edges = -1000.25:5:1000.25;
clear data_stats
for f = 1:n_freq
    for a = 1:n_amp
        ax(a,f) = subplot(n_freq, n_amp, pp); cla ; hold on
        pos = amp_data{a}.Head.vel{f};
        
        data_stats(f,a) = struct(basic_stats(pos(:),1));
        xlabel(2*data_stats(f,a).std)
        
        H(a,f) = histogram(pos, edges);
        xline( data_stats(f,a).std, '--g');
        xline(-data_stats(f,a).std, '--g');
        
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

set(ax(:,1:n_freq-1), 'XTick', [], 'XTickLabels', [], 'XColor', 'k')
set(ax(2:n_amp,:), 'YTickLabels', [],  'YColor', 'none')

%% Correlation by fly
n_fly = cellfun( @(x) x.N.fly, amp_data);
max_fly = max(n_fly);
stim_vel = (4*Amp*Freq')';
stim_vel = repmat(stim_vel,[1 1 max_fly]);
head_pos = nan(n_freq,n_amp,max_fly);
head_vel = nan(n_freq,n_amp,max_fly);
scd_rate = nan(n_freq,n_amp,max_fly);
for a = 1:n_amp
    for f = 1:n_freq
        n_fly = amp_data{a}.N.fly;
        for n = 1:n_fly
            pos = amp_data{a}.Head.pos_fly{n,f};
           	pos = pos - mean(pos,'all');
            vel = amp_data{a}.Head.vel_fly{n,f};
           	head_pos(f,a,n) = mean(abs(pos(:)));
            %head_pos(f,a,n) = std((pos(:)));
            head_vel(f,a,n) = mean(abs(vel(:)));
            scd_rate(f,a,n) = mean(amp_data{a}.Head.count{n,f}) ./ 10;
        end
    end
end

stim_vel_all = stim_vel(:);
head_pos_all = head_pos(~isnan(head_pos));
head_vel_all = scd_rate(~isnan(head_vel));
scd_rate_all = scd_rate(~isnan(scd_rate));

[P] = polyfit(head_pos_all(:), scd_rate_all(:), 1);
[rho,pval] = corr(head_pos_all(:), scd_rate_all(:));
%%
fig = figure(12); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 4 3])
clear ax H
ax(1) = subplot(1,1,1) ; cla ; hold on
xlim([-0.5 12])
ylim([-0.1 1])
cc_amp = prism(n_amp);
cc_freq = jet(n_freq);
H = gobjects(n_freq, n_amp);
mrk = ['o','x','p','s','d','*'];
for a = 1:n_amp
    for f = 1:n_freq
        H(f,a) = scatter(head_pos(f,a,:), scd_rate(f,a,:), 15, mrk(f), ...
            'MarkerEdgeColor', cc_amp(a,:), 'MarkerFaceColor', cc_amp(a,:));
    end
end
x_fit = 1:ax.XLim(2);
y_fit = polyval(P,x_fit);
plot(x_fit, y_fit, 'k', 'LineWidth', 1)

xlabel('Operating Range (°)')
ylabel('Saccade Frequency (Hz)')
title({['R^{2} = ' num2str(rho) '   p = ' num2str(pval)], [num2str(P(2)) ' + ' num2str(P(1)) 'x']})
set(ax , 'LineWidth', 1.5, 'FontSize', 9, 'Color', 'none', 'Box', 'off')

leg_amp = legend(H(1,:), string(Amp), 'Box', 'on', 'LineWidth', 1, ...
    'Location', 'NorthEastOutside');
leg_amp.Title.String = 'Amplitude (°)';

leg_axes = axes('position', ax.Position, 'visible', 'off');
leg_freq = legend(leg_axes, H(:,2), string(Freq),  'Box', 'on', 'LineWidth', 1, ...
    'Location','SouthEastOutside');
leg_freq.Title.String = 'Frequency (Hz)';

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


end