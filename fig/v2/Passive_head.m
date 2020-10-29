function [] = Passive_head()
%% Passive_head_window:
root = 'H:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'U','N','DATA')

Ramp = load("H:\DATA\Rigid_Data\Saccade\processed\Ramp_saccade_window_wave=30.mat");

%% Get windows around saccades
clc
clearvars -except DATA U N Ramp
Data = DATA(DATA.init_pos > -12,:);
% Data = DATA(DATA.init_pos > -100,:);
fly_group = findgroups(Data.Fly);

Head.time   = splitapply(@(x) {cat(2,x{:})}, Data.time, fly_group);
Head.pos    = splitapply(@(x) {cat(2,x{:})}, Data.pos, fly_group);
Head.vel    = splitapply(@(x) {cat(2,x{:})}, Data.vel, fly_group);
Head.accel  = splitapply(@(x) {cat(2,x{:})}, Data.accel, fly_group);
Head.shift_time = splitapply(@(x) {cat(2,x{:})}, Data.shift_time, fly_group);

Head_all = structfun(@(x) cat(2,x{:}), Head, 'UniformOutput', false);

Head.fly_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    Head, 'UniformOutput', false);

fnames = string(fieldnames(Head.fly_stats));
n_field = length(fnames);
for f = 1:n_field
	Head.fly_mean.(fnames(f)) = cat(2, Head.fly_stats.(fnames(f)).mean);
end
Head.grand_stats = structfun(@(x) basic_stats(x,2), Head.fly_mean, 'UniformOutput', false);

Fs_intrp = 2000;
tintrp = (0:(1/Fs_intrp):0.1)';
Head_intrp_pass = structfun(@(x) x.mean, Head.grand_stats, 'UniformOutput', false);
fnames = fields(Head_intrp_pass);
Head_intrp_pass = rmfield(Head_intrp_pass, fnames([1,5]));
Head_intrp_pass = structfun(@(x) interp1(Head.grand_stats.time.mean, x, tintrp, 'pchip'), ...
    Head_intrp_pass, 'UniformOutput', false);

speedI = 2;
Head_intrp_actv = structfun(@(x) x.mean, Ramp.Head.vel_stats(1), 'UniformOutput', false);
fnames = fields(Head_intrp_actv);
Head_intrp_actv = rmfield(Head_intrp_actv, fnames(1));
Head_intrp_actv = structfun(@(x) interp1(Ramp.Head.vel_stats.time(speedI).mean + 0.02, x, tintrp, 'pchip'), ...
    Head_intrp_actv, 'UniformOutput', false);

%% Analytical Response
% Fit coefficents
A = mean(cellfun(@(x) x.a, Data.fit));
B = mean(cellfun(@(x) x.b, Data.fit));
C = mean(cellfun(@(x) x.b, Data.fit));

% Expressions
clear F dF d2F
syms t
F(t) = A*exp(-B*t) + C;
dF(t) = diff(F,t,1);
d2F(t) = diff(F,t,2);

% Output
tt = Head.grand_stats.time.mean;
F = double(F(tt));
F = F - F(end);
dF = double(dF(tt));
d2F = double(d2F(tt));

%% Time constant
fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', 1*[2 2 2 2])
clear ax
ax(1) = subplot(1,1,1); hold on
ylim([0 50])
ylabel('Time Constant (ms)')
bx = boxplot(Data.tau, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1}, h{kk,2}, 'k',  'EdgeColor', 'none');
end
ww = 1;
set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth', 1);
set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(ww).Children = ax(ww).Children([end 1:end-1]);
ylim([0 60])
set(ax , 'LineWidth', 1, 'Box', 'off', 'XColor', 'none')

%% Anova for fly effects on time constant
anova1(Data.tau, Data.Fly)

%% Correlations
% Initial displacement vs peak velocity
[R,P] = corr(Data.init_pos, Data.peak_vel);
[coeff,~] = polyfit(Data.init_pos, Data.peak_vel, 1);
xfit = -50:0;
yfit = polyval(coeff,xfit);

fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', 1.5*[2 2 2 3])
movegui(fig, 'center')
clear ax h
ax(1) = subplot(2,1,1); hold on ; cla
scatter(Data.init_pos, Data.peak_vel, 15, 'k', 'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'k')
plot(xfit, yfit, 'r', 'LineWidth', 0.5)
xlabel('Initial Displacement (°)')
ylabel('Peak Velocity (°/s)')
title(['R = ' num2str(R)])

set(ax(1), 'XLim', [-45 0], 'YLim', [-100 2500])

% Initial displacement vs time constant
[R,P] = corr(Data.init_pos, Data.tau);
[coeff,~] = polyfit(Data.init_pos, Data.tau, 1);
xfit = -50:0;
yfit = polyval(coeff,xfit);

ax(2) = subplot(2,1,2); hold on ; cla
scatter(Data.init_pos, Data.tau, 15, 'k', 'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'k')
plot(xfit, yfit, 'r', 'LineWidth', 0.5)
xlabel('Initial Displacement (°)')
ylabel('Time Constant (ms)')
title(['R = ' num2str(R)])
set(ax(2), 'XLim', [-45 0], 'YLim', [-2 80])

set(ax, 'LineWidth', 1)

%% Passive all trials
fig = figure (4) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', 1*[2 2 2 4])
movegui(fig, 'center')
cc_passive = 'k';
cc_active = 'r';
grand_lw = 1;
speedI = 2;

pos_passive_raw = Head.grand_stats.pos.mean;
pos_passive = pos_passive_raw - pos_passive_raw(30);
F = F - pos_passive_raw(30);

pos_active = Ramp.Head.vel_stats.pos(speedI).mean;
[~,normI] = min(abs(Ramp.Head.vel_stats.time(speedI).median  + 0.0400));
pos_active = pos_active - pos_active(normI) + pos_passive(1);
time_active = Ramp.Head.vel_stats.time(speedI).mean + 0.02;

clear ax h
ax(1) = subplot(3,1,1); hold on ; cla
    plot(Head_all.time, Head_all.pos, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.25)
    h.passive = plot(Head.grand_stats.time.mean, pos_passive, ...
        'Color', cc_passive, 'LineWidth', grand_lw);
    h.active = plot(time_active, pos_active, ...
        'Color', cc_active, 'LineWidth', grand_lw);

ax(2) = subplot(3,1,2); hold on ; cla
    plot(Head_all.time, Head_all.vel, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.25)
    h.passive = plot(Head.grand_stats.time.mean, Head.grand_stats.vel.mean, ...
        'Color', cc_passive, 'LineWidth', grand_lw);
    h.active = plot(time_active, Ramp.Head.vel_stats.vel(speedI).mean, ...
        'Color', cc_active, 'LineWidth', grand_lw);

ax(3) = subplot(3,1,3); hold on ; cla
    plot(Head_all.time, Head_all.accel, 'Color', [0.5 0.5 0.5 0.5], 'LineWidth', 0.25)
    h.passive = plot(Head.grand_stats.time.mean, Head.grand_stats.accel.mean, ...
        'Color', cc_passive, 'LineWidth', grand_lw);
    h.active = plot(time_active, Ramp.Head.vel_stats.accel(speedI).mean, ...
        'Color', cc_active, 'LineWidth', grand_lw);
        
linkaxes(ax, 'x')
set(ax, 'XLim', [0 0.05])
set(ax(1), 'YLim', [-45 5])
set(ax(2), 'YLim', [-100 2600])
set(ax(3), 'YLim', 2.6e5*[-1 1])
set(ax, 'LineWidth', 1)

YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Position (°)')
YLabelHC = get(ax(2), 'YLabel');
set([YLabelHC], 'String', 'Velocity (°/s)')
YLabelHC = get(ax(3), 'YLabel');
set([YLabelHC], 'String', 'Accleration (°/s^{2})')
XLabelHC = get(ax(3), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')

%% Active vs Passive Raw
fig = figure (13) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', 1*[2 2 2 4])
movegui(fig, 'center')
cc_passive = 'k';
cc_active = 'r';
grand_lw = 1;
speedI = 2;

pos_passive_raw = Head.grand_stats.pos.mean;
pos_passive = pos_passive_raw - pos_passive_raw(30);
F = F - pos_passive_raw(30);

pos_active = Ramp.Head.vel_stats.pos(speedI).mean;
% [~,normI] = min(abs(Ramp.Head.vel_stats.time(speedI).median  + 0.0400));
% pos_active = pos_active - pos_active(normI) + pos_passive(1);

time_active = Ramp.Head.vel_stats.time(speedI).mean + 0.02;

clear ax h
ax(1) = subplot(3,1,1); hold on ; cla
    h.passive = plot(Head.grand_stats.time.mean, pos_passive, ...
        'Color', cc_passive, 'LineWidth', grand_lw);
    h.active = plot(time_active, pos_active, ...
        'Color', cc_active, 'LineWidth', grand_lw);

ax(2) = subplot(3,1,2); hold on ; cla
    h.passive = plot(Head.grand_stats.time.mean, Head.grand_stats.vel.mean, ...
        'Color', cc_passive, 'LineWidth', grand_lw);
    h.active = plot(time_active, Ramp.Head.vel_stats.vel(speedI).mean, ...
        'Color', cc_active, 'LineWidth', grand_lw);

ax(3) = subplot(3,1,3); hold on ; cla
    h.passive = plot(Head.grand_stats.time.mean, Head.grand_stats.accel.mean, ...
        'Color', cc_passive, 'LineWidth', grand_lw);
    h.active = plot(time_active, Ramp.Head.vel_stats.accel(speedI).mean, ...
        'Color', cc_active, 'LineWidth', grand_lw);
        
linkaxes(ax, 'x')
set(ax, 'XLim', [0 0.05])
set(ax(1), 'YLim', [-10 10])
set(ax(2), 'YLim', [-100 1000])
set(ax(3), 'YLim', 1e5*[-1 1])
set(ax, 'LineWidth', 1)

YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Position (°)')
YLabelHC = get(ax(2), 'YLabel');
set([YLabelHC], 'String', 'Velocity (°/s)')
YLabelHC = get(ax(3), 'YLabel');
set([YLabelHC], 'String', 'Accleration (°/s^{2})')
XLabelHC = get(ax(3), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')

%% Active vs Passive Interpolated
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', 1*[2 2 2 4])
movegui(fig, 'center')
cc_passive = 'k';
cc_active = 'r';
cc_diff = 'g';
grand_lw = 1;

pos_passive_raw = Head_intrp_pass.pos;
normI = Fs_intrp * 0.04;
pos_passive = pos_passive_raw - pos_passive_raw(normI);
F = F - pos_passive_raw(normI);

pos_active = Head_intrp_actv.pos;
% [~,normI] = min(abs(tintrp  + 0.0400));
% pos_active = pos_active - pos_active(normI) + pos_passive(1);

clear ax h
ax(1) = subplot(3,1,1); hold on ; cla
    h.passive = plot(tintrp, pos_passive, ...
        'Color', cc_passive, 'LineWidth', grand_lw);
    h.active = plot(tintrp, pos_active, ...
        'Color', cc_active, 'LineWidth', grand_lw);
	h.diff = plot(tintrp, (pos_active - pos_passive), ...
        'Color', cc_diff, 'LineWidth', grand_lw);
    
ax(2) = subplot(3,1,2); hold on ; cla
    h.passive = plot(tintrp, Head_intrp_pass.vel, ...
        'Color', cc_passive, 'LineWidth', grand_lw);
    h.active = plot(tintrp, Head_intrp_actv.vel, ...
        'Color', cc_active, 'LineWidth', grand_lw);
	h.diff = plot(tintrp, (Head_intrp_actv.vel - Head_intrp_pass.vel), ...
        'Color', cc_diff, 'LineWidth', grand_lw);
    
ax(3) = subplot(3,1,3); hold on ; cla
    h.passive = plot(tintrp, Head_intrp_pass.accel, ...
        'Color', cc_passive, 'LineWidth', grand_lw);
    h.active = plot(tintrp, Head_intrp_actv.accel, ...
        'Color', cc_active, 'LineWidth', grand_lw);
	h.diff = plot(tintrp, (Head_intrp_actv.accel - Head_intrp_pass.accel), ...
        'Color', cc_diff, 'LineWidth', grand_lw);
    
linkaxes(ax, 'x')
set(ax, 'XLim', [0 0.05])
set(ax(1), 'YLim', [-10 10])
set(ax(2), 'YLim', [-100 700])
set(ax(3), 'YLim', 0.5e5*[-1 1])
set(ax, 'LineWidth', 1)

YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Position (°)')
YLabelHC = get(ax(2), 'YLabel');
set([YLabelHC], 'String', 'Velocity (°/s)')
YLabelHC = get(ax(3), 'YLabel');
set([YLabelHC], 'String', 'Accleration (°/s^{2})')
XLabelHC = get(ax(3), 'XLabel');
set([XLabelHC], 'String', 'Time (s)')

end