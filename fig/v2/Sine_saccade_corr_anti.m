function [] = Sine_saccade_corr_anti()
%% Sine_saccade_corr_anti:
root = 'E:\DATA\Rigid_Data\Saccade\processed';
[FILE,PATH] = uigetfile({'*.mat'},'Select sinusoid datasets', root, 'MultiSelect','on');
FILE = cellstr(FILE);

n_amp = length(FILE);
amp_data = cell(n_amp,1);
Amp = nan(n_amp,1);
for a = 1:n_amp
    amp_data{a} = load(fullfile(PATH,FILE{a}));
    Amp(a) = unique(amp_data{a}.amp);
end
n_freq = amp_data{1}.N.freq;
Freq = amp_data{1}.U.freq{1};
cc_amp = prism(n_amp);
cc_freq = prism(n_freq);

ALL_amp = [];
ALL_amp_fly = [];
Scd_amp = [];
for a = 1:n_amp
    ALL_amp = [ALL_amp ; amp_data{a}.ALL];
    ALL_amp_fly = [ALL_amp_fly ; amp_data{a}.Fly];
    Scd_amp = [Scd_amp ; amp_data{a}.Anti_scd];
end

%% Save combined sine data
fname = 'Sine_all';
savedir = 'E:\DATA\Rigid_Data\Saccade\combined';
save(fullfile(savedir, [fname '.mat']), 'ALL_amp', 'ALL_amp_fly','Scd_amp');

%% Position All
fig = figure(1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 (9/5)*n_amp (8/6)*n_freq])
ax = gobjects(n_amp,n_freq);
H_all = gobjects(n_amp,n_freq);
H_starts = gobjects(n_amp,n_freq);
H_ends = gobjects(n_amp,n_freq);
pp = 1;
edges = -25:1:25;
clear data_stats
head_pos = nan(n_freq,n_amp);
scd_rate = nan(n_freq,n_amp);
head_starts = nan(n_freq,n_amp);
for f = 1:n_freq
    for a = 1:n_amp
        ax(a,f) = subplot(n_freq, n_amp, pp); cla ; hold on
        pos = cat(2, amp_data{a}.Head.pos{:,f});
        
        %starts = cat(2, amp_data{a}.Head_scd.start{:,f});
        %ends = cat(2, amp_data{a}.Head_scd.end{:,f});
        faI = amp_data{a}.Anti_scd.freq == f;
        starts = amp_data{a}.Anti_scd.StartPos(faI) .* amp_data{a}.Anti_scd.Direction(faI);
        
        head_pos(f,a) = mean(abs(pos(:)));
        scd_rate(f,a) = amp_data{a}.Head.vel_stats.anti_count(f).mean ./ 10;
        head_starts(f,a) = amp_data{a}.Head_scd.vel_stats.start(f).mean;

        H_all(a,f) = histogram(pos, edges);
     	%H_starts(a,f) = histogram(starts, edges);
        %H_ends(a,f) = histogram(ends, edges);
        xline(head_pos(f,a), '--g');
        xline(-head_pos(f,a), '--g');
        
        xlabel({2*head_pos(f,a), scd_rate(f,a)})
        
        if pp <= n_amp
            title([num2str(Amp(a)) '°'])
        end
        
        if ~mod(pp-1,n_amp)
            ylabel([num2str(Freq(f)) ' Hz'], 'FontWeight', 'bold')
        end
        
        pp = pp + 1;
    end
end
set([H_all ], 'Normalization', 'Probability', 'EdgeColor', 'none', 'FaceAlpha', 1)
set(H_all, 'FaceColor', 'k', 'FaceAlpha', 1)
% set(H_starts, 'FaceColor', 'g', 'FaceAlpha', 0.5)
% set(H_ends, 'FaceColor', 'r')
linkaxes(ax, 'xy')

set(ax , 'LineWidth', 0.75, 'FontSize', 8, 'Color', 'none')
set(ax , 'YLim', [-0.00 ax(1).YLim(2)])
set(ax , 'XLim', 25*[-1 1], 'XTick', -25:5:25)

set(ax(:,1:n_freq-1), 'XTick', [], 'XTickLabels', [], 'XColor', 'k')
set(ax(2:n_amp,:), 'YTickLabels', [],  'YColor', 'none')

%% Saccade Rate by Amplitude
freqI = 2; % raw freq or index depends
ALL_amp_freq = ALL_amp(ALL_amp.freq == freqI,:);

fig = figure (10); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 2 2])
clear ax H
ww = 1;
ax(ww) = subplot(1,1,1) ; cla ; hold on
    bx = boxplot(ALL_amp_freq.anti_count ./ 10, ALL_amp_freq.amp, 'Labels', {Amp}, ...
        'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    xlabel('Amplitude (°)')
    ylabel('Saccade Frequency (Hz)')

    h_sig = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h_sig,1)
       patch(h_sig{kk,1}, h_sig{kk,2}, cc_freq(kk,:),  'EdgeColor', 'none');
    end

    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);
    ylim([-0.1 1.2])
set(ax, 'LineWidth', 1, 'Box', 'off')

%% Correlation saccade rate vs head position by fly
x = ALL_amp_fly.OR;
% x = ALL_amp_fly.stim_vel;
y = ALL_amp_fly.anti_count ./ 10;
[P] = polyfit(x, y, 1);
[rho,pval] = corr(x, y);

fig = figure(11); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 2.5])
clear ax H
ax(1) = subplot(1,1,1) ; cla ; hold on
% xlim([-0.5 12])
ylim([-0.1 1])
H_all = gobjects(n_freq, n_amp);
mrk = ['<','o','p','s','d','>'];
% cc_fly = hsv(max_fly);
for a = 1:n_amp
    for f = 1:n_freq
        faI = (ALL_amp_fly.amp == Amp(a)) & (ALL_amp_fly.freq == Freq(f));
        H_all(f,a) = scatter(x(faI), y(faI), 10, mrk(f), ...
            'MarkerEdgeColor', cc_amp(a,:), 'MarkerFaceColor', cc_amp(a,:));
    end
end
r = roots(P);
x_fit = 0:ax.XLim(2);
y_fit = polyval(P,x_fit);
plot(x_fit, y_fit, 'k', 'LineWidth', 1)

xlabel('Operating Range (°)')
ylabel('Saccade Frequency (Hz)')
title({['r = ' num2str(rho) '   p = ' num2str(pval)], [num2str(P(2)) ' + ' num2str(P(1)) 'x']})
set(ax , 'LineWidth', 1, 'FontSize', 8, 'Color', 'none', 'Box', 'off')

leg_amp = legend(H_all(1,:), string(Amp), 'Box', 'on', 'LineWidth', 1, ...
    'Location', 'NorthEastOutside');
leg_amp.Title.String = 'Amplitude (°)';

leg_axes = axes('position', ax.Position, 'visible', 'off');
leg_freq = legend(leg_axes, H_all(:,end), string(Freq),  'Box', 'on', 'LineWidth', 1, ...
    'Location','SouthEastOutside');
leg_freq.Title.String = 'Frequency (Hz)';

%% Correlation saccade trigger position vs head position by fly
x = ALL_amp_fly.OR;
y = -ALL_amp_fly.start;
[P] = polyfit(x(~isnan(y)), y(~isnan(y)), 1);
[rho,pval] = corr(x(~isnan(y)), y(~isnan(y)));

fig = figure(12); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 2.5])
clear ax H
ax(1) = subplot(1,1,1) ; cla ; hold on
% xlim([-0.5 12])
ylim([-5 15])
H_all = gobjects(n_freq, n_amp);
mrk = ['<','o','p','s','d','>'];
% cc_fly = hsv(max_fly);
for a = 1:n_amp
    for f = 1:n_freq
        faI = (ALL_amp_fly.amp == Amp(a)) & (ALL_amp_fly.freq == Freq(f));
        H_all(f,a) = scatter(x(faI), y(faI), 10, mrk(f), ...
            'MarkerEdgeColor', cc_amp(a,:), 'MarkerFaceColor', cc_amp(a,:));
    end
end
r = roots(P);
x_fit = 0:ax.XLim(2);
y_fit = polyval(P,x_fit);
plot(x_fit, y_fit, 'k', 'LineWidth', 1)

xlabel('Operating range (°)')
ylabel('Saccade trigger (°)')
title({['r = ' num2str(rho) '   p = ' num2str(pval)], [num2str(P(2)) ' + ' num2str(P(1)) 'x']})
set(ax , 'LineWidth', 1, 'FontSize', 8, 'Color', 'none', 'Box', 'off')

leg_amp = legend(H_all(1,:), string(Amp), 'Box', 'on', 'LineWidth', 1, ...
    'Location', 'NorthEastOutside');
leg_amp.Title.String = 'Amplitude (°)';

leg_axes = axes('position', ax.Position, 'visible', 'off');
leg_freq = legend(leg_axes, H_all(:,end), string(Freq),  'Box', 'on', 'LineWidth', 1, ...
    'Location','SouthEastOutside');
leg_freq.Title.String = 'Frequency (Hz)';

%% Correlation saccade start vs head position by fly
ALL_amp_fly = ALL_amp_fly(~isnan(ALL_amp_fly.start),:);
[P] = polyfit(ALL_amp_fly.OR, -ALL_amp_fly.start, 1);
[rho,pval] = corr(ALL_amp_fly.OR, ALL_amp_fly.start);

fig = figure(12); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 2.5])
clear ax H
ax(1) = subplot(1,1,1) ; cla ; hold on
% xlim([-0.5 12])
ylim([-5 15])
H_all = gobjects(n_freq, n_amp);
mrk = ['<','o','p','s','d','>'];
% cc_fly = hsv(max_fly);
for a = 1:n_amp
    for f = 1:n_freq
        faI = (ALL_amp_fly.amp == Amp(a)) & (ALL_amp_fly.freq == Freq(f));
        H_all(f,a) = scatter(ALL_amp_fly.OR(faI), -ALL_amp_fly.start(faI), 10, mrk(f), ...
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

%% Trigger Comparison
clc
clear scd
[amp_freq_group, amp_group, freq_group] = findgroups(ALL_amp_fly.amp, ALL_amp_fly.freq);

scd.or = splitapply(@(x) {x}, ALL_amp_fly.OR, amp_freq_group);
scd.start = splitapply(@(x) {x}, -ALL_amp_fly.start, amp_freq_group);
scd.rate = splitapply(@(x) {x}, ALL_amp_fly.anti_count ./ 10, amp_freq_group);
scd.pos = splitapply(@(x) {cat(2,x{:})}, ALL_amp_fly.pos, amp_freq_group);
scd.pos = cellfun(@(x) x(:), scd.pos, 'UniformOutput', false);

% Anova
n_group = length(amp_group);
pval = nan(n_group,1);
alpha = 0.01;
for n = 1:n_group
    y = [scd.pos{n} ; scd.start{n}];
    g = [ones(size(scd.pos{n})) ; 2*ones(size(scd.start{n}))];
    pval(n) = anova1(y,g, 'off');
end
h_sig = pval < alpha;

start_sig = [];
for n = 1:n_group
    if h_sig(n)
        amp_freq_I = (Scd_amp.amp==amp_group(n)) & (Freq(Scd_amp.freq)==freq_group(n));
        start_sig = [start_sig ; -Scd_amp.StartPos(amp_freq_I)];
    end
end

scd_fly = structfun(@(x) splitapply(@(y) {y}, x, findgroups(amp_group)), ...
    scd, 'UniformOutput', false);
scd_fly = structfun(@(x) cat(2,x{:}), scd_fly, 'UniformOutput', false);

scd_fly.stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    scd_fly, 'UniformOutput', false);

% fnames = string(fieldnames(scd_fly.stats));
% n_field = length(fnames);
% for f = 1:n_field
%     for n = 1:n_group
%         fr = find(freq_group(n) == Freq);
%         a = find(amp_group(n) == Amp);
%         scd_fly.mean.(fnames(f)){fr,a} = cat(2, scd_fly.stats.(fnames(f))(:,n).mean);
%     end
% end
% 
% scd_fly.vel_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
%     scd_fly.mean, 'UniformOutput', false);

G = structfun(@(z) cellfun(@(x,y) y*ones(size(x)), z, num2cell((1:length(z))'), ...
    'UniformOutput', false), scd, 'UniformOutput', false);
G = structfun(@(z) cat(1,z{:}), G, 'UniformOutput', false);

CC = structfun(@(z) cellfun(@(x,y) y*ones(size(x)), z, num2cell(amp_group), ...
    'UniformOutput', false), scd, 'UniformOutput', false);
CC = structfun(@(z) cat(1,z{:}), CC, 'UniformOutput', false);

scd = structfun(@(z) cat(1,z{:}), scd, 'UniformOutput', false);

% Fly Stats
fig = figure(20); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 7 5])
clear ax b
ax(1) = subplot(4,1,1) ; cla ; hold on ; ylim(28*[-1 1])
    ylabel('Position (°)')
    b(:,1) = boxchart(G.pos, scd.pos,'GroupByColor', CC.pos, 'BoxFaceAlpha', 0, ...
        'JitterOutliers', 'on', 'MarkerStyle', 'none', 'BoxWidth', 1, 'WhiskerLineColor', 'k');
    
	b(:,2) = boxchart(G.start, scd.start,'GroupByColor', CC.start, 'BoxFaceAlpha', 1, 'BoxFaceColor', 'k', ...
        'JitterOutliers', 'on', 'MarkerStyle', 'none', 'BoxWidth', 1, 'WhiskerLineColor', 'm');
    
ax(2) = subplot(4,1,2) ; cla ; hold on ; ylim([0 16]) ; yticks([0:5:15])
    ylabel('Operating range (°)')
    b(:,3) = boxchart(G.or, scd.or,'GroupByColor', CC.or, 'BoxFaceAlpha', 0, ...
        'JitterOutliers', 'on', 'MarkerStyle', 'none', 'BoxWidth', 1, 'WhiskerLineColor', 'k');

ax(3) = subplot(4,1,3) ; cla ; hold on ; ylim([-0.1 1])
    ylabel('Saccade frequency (Hz)')
    b(:,4) = boxchart(G.rate, scd.rate, 'GroupByColor', CC.rate, 'BoxFaceAlpha', 0, ...
        'JitterOutliers', 'on', 'MarkerStyle', 'none', 'BoxWidth', 1, 'WhiskerLineColor', 'k');
    
ax(4) = subplot(4,1,4) ; cla ; hold on ; ylim([-0.3 3])
    ylabel('p-value')
    plot(pval, '-k')
    yticks([0 1])
    for n = 1:n_group
        if h_sig(n)
            plot(n, pval(n), 'r.', 'MarkerSize', 10)
        else
            plot(n, pval(n), 'k.', 'MarkerSize', 10)
        end
    end
    
amp_cc = prism(n_amp);
for a = 1:n_amp
    set(b(a,[1,3:end]), 'BoxFaceColor', amp_cc(a,:))
end
set(b,'LineWidth', 0.5)

set(ax, 'XTick', 1:30, 'XTickLabels', repmat(Freq,[n_amp,1]))
set(ax , 'LineWidth', 1, 'Color', 'none', 'Box', 'off')
set(ax, 'XLim', [0 31])
linkaxes(ax, 'x')

% set(ax, 'XColor', 'none')

%% CDF's
clear ax h H
fig = figure (1);
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 3.5])
H = gobjects(n_amp,n_freq);
starts_all = cell(n_amp,n_freq);
ax(1) = subplot(2,1,1); cla; hold on ; ylabel('Probability')
    edges = -25:25;
    for a = 1:n_amp
        for f = 1:n_freq
            amp_freq_I = (Scd_amp.amp==Amp(a)) & (Scd_amp.freq==f);
            starts_all{a,f} = -Scd_amp.StartPos(amp_freq_I);
            H(a,f) = histogram(starts_all{a,f}, edges, 'Normalization', 'pdf', ...
                'FaceAlpha', 0.3, 'EdgeColor', 'none');
            %xline(median(starts), 'k', 'LineWidth', 1)
        end
        cla
        histogram(start_sig, edges, 'Normalization', 'pdf', ...
            'FaceColor', 'k', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    
ax(2) = subplot(2,1,2); cla; hold on
    ylim([-0.1 1])
    xlabel('Trigger Position (°)')
    ylabel('Probability')
    for a = 1:n_amp
        for f = 1:n_freq
            %[h_cdf,~] = cdfplot(starts_all{a,f});
          	%set(h_cdf, 'Color', cc_amp(a,:), 'LineWidth', 1)
            title('')
        end
    end
    %[h_cdf,~] = cdfplot(cat(1,starts_all{:}));
    [h_cdf,~] = cdfplot(start_sig);
    set(h_cdf, 'Color', 'k', 'LineWidth', 1)
    title('')
    
%     [fitresult, gof] = sigmoid_fit(h_cdf.XData, h_cdf.YData);
%     
%  	xfit = -11:0.1:25;
%     yfit = 1 ./ (1 + fitresult.a*exp(fitresult.b*xfit));
%     h_fit = plot(xfit, yfit, 'r', 'LineWidth', 1);
%     
%     legend([h_cdf h_fit], 'CDF', 'Sigmoid Fit', 'Box', 'off')
    
set(ax, 'LineWidth', 1, 'Box', 'off', 'XGrid', 'off', 'YGrid', 'off')
linkaxes(ax, 'x')
set(ax, 'XTick', -25:5:25)

%% All compare start/pos
y = [scd.pos ; scd.start];
g = [ones(size(scd.pos)) ; 2*ones(size(scd.start))];
p_all = anova1(y, g, 'off');

fig = figure (1);
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 3.5])
clear ax 
ax(1) = subplot(1,1,1); cla; hold on
boxchart(g, y, 'Notch', 'on', 'BoxFaceAlpha', 0.5, ...
        'JitterOutliers', 'on', 'MarkerStyle', 'none', 'BoxWidth', 0.5, 'WhiskerLineColor', 'k');
set(ax, 'LineWidth', 1, 'Box', 'off')

%% Correlation
clc
x = ALL_amp_fly.OR;
% x = ALL_amp_fly.amp;
x = ALL_amp_fly.freq;

y = ALL_amp_fly.anti_count ./ 10;

[P] = polyfit(x(~isnan(y)), y(~isnan(y)), 1);
[rho,pval] = corr(x(~isnan(y)), y(~isnan(y)));
disp(rho)
disp(pval)

% mdl = fitlm([ALL_amp_fly.m], ALL_amp_fly.anti_count ./10);
mdl = fitlm([ALL_amp_fly.freq ALL_amp_fly.amp ALL_amp_fly.stim_vel], ALL_amp_fly.OR);

disp(sqrt(mdl.Rsquared.Ordinary))

%% Grand Mean Corr
OR_all = nan(n_freq,n_amp);
Rate_all = nan(n_freq,n_amp);
for a = 1:n_amp
    for f = 1:n_freq
       OR_all(:,a) = cat(1, amp_data{a}.Head.vel_stats.OR.mean);
       Rate_all(:,a) = cat(1, amp_data{a}.Head.vel_stats.anti_count.mean) ./ 10;
    end
end

x = OR_all(:);
y = Rate_all(:);
[P] = polyfit(x(~isnan(y)), y(~isnan(y)), 1);
[rho,pval] = corr(x(~isnan(y)), y(~isnan(y)));

fig = figure(30); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 3 3])
clear ax H
ax(1) = subplot(1,1,1) ; cla ; hold on
% xlim([-0.5 12])
ylim([-0.1 1])
H_all = gobjects(n_freq, n_amp);
mrk = ['<','o','p','s','d','>'];
% cc_fly = hsv(max_fly);
for a = 1:n_amp
    for f = 1:n_freq
        H_all(f,a) = scatter(OR_all(f,a), Rate_all(f,a), 20, mrk(f), ...
            'MarkerEdgeColor', cc_amp(a,:), 'MarkerFaceColor', cc_amp(a,:));
    end
end
r = roots(P);
x_fit = 0:ax.XLim(2);
y_fit = polyval(P,x_fit);
plot(x_fit, y_fit, 'k', 'LineWidth', 1)

xlabel('Operating range (°)')
ylabel('Saccade trigger (°)')
title({['r = ' num2str(rho) '   p = ' num2str(pval)], [num2str(P(2)) ' + ' num2str(P(1)) 'x']})
set(ax , 'LineWidth', 1, 'FontSize', 8, 'Color', 'none', 'Box', 'off')

% leg_amp = legend(H_all(1,:), string(Amp), 'Box', 'on', 'LineWidth', 1, ...
%     'Location', 'NorthEastOutside');
% leg_amp.Title.String = 'Amplitude (°)';
% 
% leg_axes = axes('position', ax.Position, 'visible', 'off');
% leg_freq = legend(leg_axes, H_all(:,end), string(Freq),  'Box', 'on', 'LineWidth', 1, ...
%     'Location','SouthEastOutside');
% leg_freq.Title.String = 'Frequency (Hz)';

end