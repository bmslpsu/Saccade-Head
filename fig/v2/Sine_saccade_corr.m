function [] = Sine_saccade_corr()
%% Sine_saccade_corr:
root = 'H:\DATA\Rigid_Data\Saccade\processed';
[FILE,PATH] = uigetfile({'*.mat'},'Select sinusoid datasets', root, 'MultiSelect','on');
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
cc_amp = prism(n_amp);
cc_freq = prism(n_freq);

%% Position All
fig = figure(1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 (9/5)*n_amp (8/6)*n_freq])
ax = gobjects(n_amp,n_freq);
H_all = gobjects(n_amp,n_freq);
H_starts = gobjects(n_amp,n_freq);
H_ends = gobjects(n_amp,n_freq);
pp = 1;
edges = -25.25:2:25.25;
clear data_stats
head_pos = nan(n_freq,n_amp);
scd_rate = nan(n_freq,n_amp);
head_starts = nan(n_freq,n_amp);
for f = 1:n_freq
    for a = 1:n_amp
        ax(a,f) = subplot(n_freq, n_amp, pp); cla ; hold on
        pos = amp_data{a}.Head.pos{f};
        pos = pos - mean(pos);
        
        starts = amp_data{a}.Head.starts_all{f};
        ends = amp_data{a}.Head.ends_all{f};
        
        head_pos(f,a) = mean(abs(pos(:)));
        scd_rate(f,a) = mean(amp_data{a}.Head.count_all{f}) ./ 10;
        head_starts(f,a) = mean(starts);

        H_all(a,f) = histogram(pos, edges);
     	H_starts(a,f) = histogram(starts, edges);
        %H_ends(a,f) = histogram(ends, edges);
        %xline(head_pos(f,a), '--g');
        %xline(-head_pos(f,a), '--g');
        
        %xlabel({2*head_pos(f,a), scd_rate(f,a)})
        
        if pp <= n_amp
            title([num2str(Amp(a)) '°'])
        end
        
        if ~mod(pp-1,n_amp)
            ylabel([num2str(Freq(f)) ' Hz'], 'FontWeight', 'bold')
        end
        
        pp = pp + 1;
    end
end
set([H_all H_starts], 'Normalization', 'Probability', 'EdgeColor', 'none', 'FaceAlpha', 1)
set(H_all, 'FaceColor', 'k', 'FaceAlpha', 1)
set(H_starts, 'FaceColor', 'g', 'FaceAlpha', 0.5)
% set(H_ends, 'FaceColor', 'r')
linkaxes(ax, 'xy')

set(ax , 'LineWidth', 0.75, 'FontSize', 8, 'Color', 'none')
set(ax , 'YLim', [-0.00 ax(1).YLim(2)])
set(ax , 'XLim', 25*[-1 1], 'XTick', -25:5:25)

set(ax(:,1:n_freq-1), 'XTick', [], 'XTickLabels', [], 'XColor', 'k')
set(ax(2:n_amp,:), 'YTickLabels', [],  'YColor', 'none')

%% Velocity All
fig = figure (2); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 8 7])
ax = gobjects(n_amp,n_freq);
H_all = gobjects(n_amp,n_freq);
pp = 1;
edges = -1000.25:5:1000.25;
clear data_stats
head_vel = nan(n_freq,n_amp);
scd_rate = nan(n_freq,n_amp);
for f = 1:n_freq
    for a = 1:n_amp
        ax(a,f) = subplot(n_freq, n_amp, pp); cla ; hold on
        vel = amp_data{a}.Head.vel{f};
        
        head_vel(f,a) = mean(abs(vel(:)));
        scd_rate(f,a) = mean(amp_data{a}.Head.count_all{f}) ./ 10;
        
        H_all(a,f) = histogram(vel, edges);
        xline( head_vel(f,a), '--g');
        xline(-head_vel(f,a), '--g');
        
        xlabel(2*head_vel(f,a))
        
        if pp <= n_amp
            title([num2str(Amp(a)) '°'])
        end
        
        if ~mod(pp-1,n_amp)
            ylabel([num2str(Freq(f)) ' Hz'], 'FontWeight', 'bold')
        end
        
        pp = pp + 1;
    end
end
set(H_all, 'Normalization', 'Probability')
set(H_all, 'EdgeColor', 'none', 'FaceAlpha', 1, 'FaceColor', 'k')
linkaxes(ax, 'xy')

set(ax , 'LineWidth', 0.75, 'FontSize', 8, 'Color', 'none')
set(ax , 'YLim', [-0.00 ax(1).YLim(2)])
set(ax , 'XLim', 500*[-1 1], 'XTick', -1000:500:1000)

set(ax(:,1:n_freq-1), 'XTick', [], 'XTickLabels', [], 'XColor', 'k')
set(ax(2:n_amp,:), 'YTickLabels', [],  'YColor', 'none')

%% Saccade Rate by Amplitude
fI = 2;
freq_all = cell(n_freq,1);
G = cell(n_freq,1);
Count_All = cell(n_freq,n_amp);
for a = 1:n_amp
    fly_mean = cellfun(@(x) mean(x), amp_data{a}.Head.count);
    for f = 1:n_freq
        %freq_data = amp_data{a}.Head.count_all{f};
        %freq_all{f} = [freq_all{f} ; freq_data];
        %G{f} = [G{f} ; a*ones(length(freq_data),1)];
        %Count_All{f,a} = freq_data;
        Count_All{f,a} = fly_mean(:,f);
        G{f,a} = a*ones(size(Count_All{f,a}));
    end
end
Count_All  = cat(1, Count_All{fI,:});
G  = cat(1, G{fI,:});

fig = figure (10); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2 2])
clear ax H
ww = 1;
ax(ww) = subplot(1,1,1) ; cla ; hold on
    bx = boxplot(Count_All ./ 10, G, 'Labels', {Amp}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    xlabel('Amplitude (°)')
    ylabel('Saccade Frequency (Hz)')

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1}, h{kk,2}, cc_freq(kk,:),  'EdgeColor', 'none');
    end

    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);
    ylim([-0.1 1.5])
set(ax, 'LineWidth', 1, 'Box', 'off')

%% Correlation by fly
n_fly = cellfun( @(x) x.N.fly, amp_data);
max_fly = max(n_fly);
head_pos = nan(n_freq,n_amp,max_fly);
head_vel = nan(n_freq,n_amp,max_fly);
scd_rate = nan(n_freq,n_amp,max_fly);
scd_starts = nan(n_freq,n_amp,max_fly);
scd_ends = nan(n_freq,n_amp,max_fly);
for a = 1:n_amp
    for f = 1:n_freq
        n_fly = amp_data{a}.N.fly;
        for n = 1:n_fly
            pos = amp_data{a}.Head.pos_fly{n,f};
           	pos = pos - mean(pos,'all');
            vel = amp_data{a}.Head.vel_fly{n,f};
            starts = amp_data{a}.Head.starts{n,f};
            ends = amp_data{a}.Head.ends{n,f};
           	head_pos(f,a,n) = 2*mean(abs(pos(:)));
            head_vel(f,a,n) = mean(abs(vel(:)));
            scd_rate(f,a,n) = mean(amp_data{a}.Head.count{n,f}) ./ 10;
            scd_starts(f,a,n) = mean(starts(:));
            scd_ends(f,a,n) = mean(ends(:));
        end
    end
end
rmv_idx = ~isnan(head_pos) & ~isnan(scd_rate) & ~isnan(scd_starts);
head_pos_all = head_pos(rmv_idx);
head_vel_all = head_vel(rmv_idx);
scd_rate_all = scd_rate(rmv_idx);
scd_starts_all = scd_starts(rmv_idx);
scd_ends_all = scd_ends(rmv_idx);

stim_vel = (4*Amp*Freq')';
stim_vel = repmat(stim_vel,[1 1 max_fly]);
stim_vel(isnan(head_pos)) = nan;
stim_vel_all = stim_vel(~isnan(stim_vel));

%% Correlation saccade rate vs head position by fly
[P] = polyfit(head_pos_all(:), scd_rate_all(:), 1);
[rho,pval] = corr(head_pos_all(:), scd_rate_all(:));

fig = figure(11); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 2.5])
clear ax H
ax(1) = subplot(1,1,1) ; cla ; hold on
% xlim([-0.5 12])
ylim([-0.1 1])
H_all = gobjects(n_freq, n_amp);
mrk = ['<','o','p','s','d','>'];
cc_fly = hsv(max_fly);
for a = 1:n_amp
    for f = 1:n_freq
        H_all(f,a) = scatter(head_pos(f,a,:), scd_rate(f,a,:), 10, mrk(f), ...
            'MarkerEdgeColor', cc_amp(a,:), 'MarkerFaceColor', cc_amp(a,:));
    end
end
r = roots(P);
x_fit = r:ax.XLim(2);
y_fit = polyval(P,x_fit);
plot(x_fit, y_fit, 'k', 'LineWidth', 1)

xlabel('Operating Range (°)')
ylabel('Saccade Frequency (Hz)')
title({['R^{2} = ' num2str(rho) '   p = ' num2str(pval)], [num2str(P(2)) ' + ' num2str(P(1)) 'x']})
set(ax , 'LineWidth', 1, 'FontSize', 8, 'Color', 'none', 'Box', 'off')

leg_amp = legend(H_all(1,:), string(Amp), 'Box', 'on', 'LineWidth', 1, ...
    'Location', 'NorthEastOutside');
leg_amp.Title.String = 'Amplitude (°)';

leg_axes = axes('position', ax.Position, 'visible', 'off');
leg_freq = legend(leg_axes, H_all(:,end), string(Freq),  'Box', 'on', 'LineWidth', 1, ...
    'Location','SouthEastOutside');
leg_freq.Title.String = 'Frequency (Hz)';

%% Correlation saccade start vs head position by fly
[P] = polyfit(head_pos_all(:), scd_starts_all(:), 1);
[rho,pval] = corr(head_pos_all(:), scd_starts_all(:));

fig = figure(11); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 2.5])
clear ax H
ax(1) = subplot(1,1,1) ; cla ; hold on
% xlim([-0.5 12])
% ylim([-0.1 1])
H_all = gobjects(n_freq, n_amp);
mrk = ['<','o','p','s','d','>'];
cc_fly = hsv(max_fly);
for a = 1:n_amp
    for f = 1:n_freq
        H_all(f,a) = scatter(head_pos(f,a,:), scd_starts(f,a,:), 10, mrk(f), ...
            'MarkerEdgeColor', cc_amp(a,:), 'MarkerFaceColor', cc_amp(a,:));
    end
end
r = roots(P);
x_fit = r:ax.XLim(2);
y_fit = polyval(P,x_fit);
plot(x_fit, y_fit, 'k', 'LineWidth', 1)

xlabel('Operating Range (°)')
ylabel('Saccade Frequency (Hz)')
title({['R^{2} = ' num2str(rho) '   p = ' num2str(pval)], [num2str(P(2)) ' + ' num2str(P(1)) 'x']})
set(ax , 'LineWidth', 1, 'FontSize', 8, 'Color', 'none', 'Box', 'off')

leg_amp = legend(H_all(1,:), string(Amp), 'Box', 'on', 'LineWidth', 1, ...
    'Location', 'NorthEastOutside');
leg_amp.Title.String = 'Amplitude (°)';

leg_axes = axes('position', ax.Position, 'visible', 'off');
leg_freq = legend(leg_axes, H_all(:,end), string(Freq),  'Box', 'on', 'LineWidth', 1, ...
    'Location','SouthEastOutside');
leg_freq.Title.String = 'Frequency (Hz)';

%% Correlation saccade rate vs stimulus velocity by fly
[P] = polyfit(stim_vel_all(:), scd_rate_all(:), 1);
[rho,pval] = corr(stim_vel_all(:), scd_rate_all(:));

fig = figure(12); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 2.5])
clear ax H
ax(1) = subplot(1,1,1) ; cla ; hold on
% xlim([-0.5 12])
ylim([-0.1 1])
H_all = gobjects(n_freq, n_amp);
mrk = ['<','o','p','s','d','>'];
for a = 1:n_amp
    for f = 1:n_freq
        %repmat(Freq(f),[max_fly,1])
        H_all(f,a) = scatter(stim_vel(f,a,:), scd_rate(f,a,:), 10, mrk(f), ...
            'MarkerEdgeColor', cc_amp(a,:), 'MarkerFaceColor', cc_amp(a,:));
    end
end
x_fit = 1:ax.XLim(2);
y_fit = polyval(P,x_fit);
plot(x_fit, y_fit, 'k', 'LineWidth', 1)
xlim([-50 1000])
xlabel('Stimulus Speed (°/s)')
ylabel('Saccade Frequency (Hz)')
title({['R^{2} = ' num2str(rho) '   p = ' num2str(pval)], [num2str(P(2)) ' + ' num2str(P(1)) 'x']})
set(ax , 'LineWidth', 1, 'FontSize', 9, 'Color', 'none', 'Box', 'off')

leg_amp = legend(H_all(1,:), string(Amp), 'Box', 'on', 'LineWidth', 1, ...
    'Location', 'NorthEastOutside');
leg_amp.Title.String = 'Amplitude (°)';

leg_axes = axes('position', ax.Position, 'visible', 'off');
leg_freq = legend(leg_axes, H_all(:,2), string(Freq),  'Box', 'on', 'LineWidth', 1, ...
    'Location','SouthEastOutside');
leg_freq.Title.String = 'Frequency (Hz)';

end