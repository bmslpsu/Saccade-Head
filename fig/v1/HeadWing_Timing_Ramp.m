function [] = HeadWing_Timing_Ramp()
%% HeadWing_Timing:

root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE','U','N')

%% Group head & wing saccades and compute stats %%
keepI = cellfun(@(x) isstruct(x) | isobject(x), SACCADE.head2wing);
Saccade = SACCADE(keepI,:);

hw_field = 'head2wing';
% hw_field = 'head2wing_align_wing';

clear head wing
head.fly = [];
wing.fly = [];
head.grand = [];
wing.grand = [];
head.fly_speed = [];
wing.fly_speed = [];
nvel = length(unique(Saccade.vel));
uvel = unique(Saccade.vel);
nfly = length(unique(Saccade.fly));
ufly = unique(Saccade.fly);
for v = 1:nvel
    velI = uvel(v) == Saccade.vel;
    for a = 1:nfly
        flyI = ufly(a) == Saccade.fly;
        sI = velI & flyI;
        hsacd = cellfun(@(x) x.hsacd, Saccade.(hw_field)(sI,:), 'UniformOutput', true);
        wsacd = cellfun(@(x) x.wsacd, Saccade.(hw_field)(sI,:), 'UniformOutput', true);
        get_fields = string(fieldnames(hsacd));
        get_fields = get_fields(1:end-1);
        stats_fields = get_fields + "_stats";
        for f = 1:length(get_fields)
            head.fly(a,v).(get_fields(f)) = cat(2, hsacd.(get_fields(f)));
            wing.fly(a,v).(get_fields(f)) = cat(2, wsacd.(get_fields(f)));
            head.fly(a,v).(stats_fields(f)) = basic_stats(head.fly(a,v).(get_fields(f)),2);
            wing.fly(a,v).(stats_fields(f)) = basic_stats(wing.fly(a,v).(get_fields(f)),2);
        end
    end
end

% Combine speeds
for v = 1:nvel/2
    ccw = v + nvel/2;
    for a = 1:nfly
        for f = 1:length(get_fields)
            if f == 1 % don't invert time
                head.fly_speed(a,v).(get_fields(f)) = cat(2, head.fly(a,v).(get_fields(f)), ...
                                                             head.fly(a,ccw).(get_fields(f)) );
                wing.fly_speed(a,v).(get_fields(f)) = cat(2, wing.fly(a,v).(get_fields(f)), ...
                                                             wing.fly(a,ccw).(get_fields(f)) );
            else
                head.fly_speed(a,v).(get_fields(f)) = cat(2, -head.fly(a,v).(get_fields(f)), ...
                                                              head.fly(a,ccw).(get_fields(f)) );
                wing.fly_speed(a,v).(get_fields(f)) = cat(2, -wing.fly(a,v).(get_fields(f)), ...
                                                              wing.fly(a,ccw).(get_fields(f)) );
            end
            head.fly_speed(a,v).(stats_fields(f)) = basic_stats(head.fly_speed(a,v).(get_fields(f)),2);
            wing.fly_speed(a,v).(stats_fields(f)) = basic_stats(wing.fly_speed(a,v).(get_fields(f)),2);
        end
    end
end

mean_win = round(0.2*200);
for v = 1:nvel/2
    for f = 1:length(stats_fields)
        temp = cat(2,head.fly_speed(:,v).(stats_fields(f)));
        head.grand(v).(get_fields(f)) = cat(2,temp.mean);
        temp = cat(2,wing.fly_speed(:,v).(stats_fields(f)));
        wing.grand(v).(get_fields(f)) = cat(2,temp.mean);
    end
    
    % Normalize dwba position
    for k = 1:size(wing.grand(v).(get_fields(f)),2)
        wing.grand(v).('pos')(:,k) = wing.grand(v).('pos')(:,k) - ...
            0*mean(wing.grand(v).('pos')(1:mean_win,k));
        wing.grand(v).('pos_desync')(:,k) = wing.grand(v).('pos_desync')(:,k) - ...
            0*mean(wing.grand(v).('pos_desync')(1:mean_win,k));
    end
    
    % Normalize head position
    for k = 1:size(head.grand(v).(get_fields(f)),2)
        head.grand(v).('pos')(:,k) = head.grand(v).('pos')(:,k) - ...
            0*mean(head.grand(v).('pos')(:,k));
        head.grand(v).('pos_desync')(:,k) = head.grand(v).('pos_desync')(:,k) - ...
            0*mean(head.grand(v).('pos_desync')(:,k));
    end
    
    for f = 1:length(stats_fields)
        wing.grand(v).(stats_fields(f)) = basic_stats(wing.grand(v).(get_fields(f)),2);
        head.grand(v).(stats_fields(f)) = basic_stats(head.grand(v).(get_fields(f)),2);
    end
end

%% Sync Stats by fly
sync_fly = nan(nfly,1);
for a = 1:nfly
    all_head = head.fly_speed(a).pos;
    sync_fly(a) = sum(isnan(all_head(1,:))) / size(all_head,2);
end

fig = figure (300) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 2 2])
movegui(fig, 'center')
clear ax
ww = 1;
ax(ww) = subplot(1,1,ww); axis tight
    bx = boxplot(100*sync_fly, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    ylabel('Simultaneous Saccade Percent')
    ylim([0 100])

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2}, [0.3 0.1 0.6]);
    end
    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);

set(ax,'LineWidth', 1, 'Box', 'off', 'XColor', 'none')

%% Synchronized head-wing saccades
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 3 3])
movegui(fig, 'center')
clear ax
head_color = [0 0 1];
wing_color = [1 0 0];
ax(1) = subplot(2,1,1); hold on ; cla
        hold on ; cla ; ylabel('Head (°)')
        ylim(20*[-1 1])
%         plot(head.grand.time, head.grand.pos,...
%             '-', 'Color', [0.7*head_color 0.3], 'LineWidth', 0.25)
        [~,~] = PlotPatch(head.grand.pos_stats.mean, head.grand.pos_stats.std, ...
            head.grand.time_stats.mean, 1, 1, head_color, 0.7*head_color, 0.3, 1);
        
%         plot(wing.grand.time, wing.grand.pos,...
%             '-', 'Color', [0.7*wing_color 0.3], 'LineWidth', 0.25)
        [~,~] = PlotPatch(wing.grand.pos_stats.mean, wing.grand.pos_stats.std, ...
            wing.grand.time_stats.mean, 1, 1, wing_color, 0.7*wing_color, 0.3, 1);
        
        
ax(2) = subplot(2,1,2); hold on ; cla ; xlabel('Time (s)')
    hold on ; cla ; ylabel('Head Velocity (°/s)')
    ylim([-200 600])
%         plot(head.grand.time, head.grand.vel,...
%             '-', 'Color', [0.7*head_color 0.3], 'LineWidth', 0.25)
        [~,~] = PlotPatch(head.grand.vel_stats.mean, head.grand.vel_stats.std, ...
            head.grand.time_stats.mean, 1, 1, head_color, 0.7*head_color, 0.3, 1);
        
%         plot(wing.grand.time, wing.grand.vel,...
%             '-', 'Color', [0.7*wing_color 0.3], 'LineWidth', 0.25)
        [~,~] = PlotPatch(wing.grand.vel_stats.mean, wing.grand.vel_stats.std, ...
            wing.grand.time_stats.mean, 1, 1, wing_color, 0.7*wing_color, 0.3, 1);
        
% ax(1).YAxis(1).Color = head_color;
% ax(1).YAxis(2).Color = wing_color;
% ax(2).YAxis(1).Color = head_color;
% ax(2).YAxis(2).Color = wing_color;
linkaxes(ax,'x')
set(ax,'XLim', 0.5*[-1 1])
% set(ax,'XTick', -0.5:0.1:0.5)
set(ax(1),'XTick',[])
set(ax,'LineWidth', 1)

%% Deynchronized head-wing saccades
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 3 3])
movegui(fig, 'center')
clear ax
head_color = [0 0 1];
wing_color = [1 0 0];
ax(1) = subplot(2,1,1); hold on ; cla
        hold on ; cla ; ylabel('Head (°)')
        ylim(20*[-1 1])
%         plot(head.grand.time, head.grand.pos_desync,...
%             '-', 'Color', [0.7*head_color 0.3], 'LineWidth', 0.25)
        [~,~] = PlotPatch(head.grand.pos_desync_stats.mean, head.grand.pos_desync_stats.std, ...
            head.grand.time_stats.mean, 1, 1, head_color, 0.7*head_color, 0.3, 1);
        
%         plot(wing.grand.time, wing.grand.pos_desync,...
%             '-', 'Color', [0.7*wing_color 0.3], 'LineWidth', 0.25)
        [~,~] = PlotPatch(wing.grand.pos_desync_stats.mean, wing.grand.pos_desync_stats.std, ...
            wing.grand.time_stats.mean, 1, 1, wing_color, 0.7*wing_color, 0.3, 1);
        
ax(2) = subplot(2,1,2); hold on ; cla ; xlabel('Time (s)')
    hold on ; cla ; ylabel('Head Velocity (°/s)')
    ylim([-200 600])
%         plot(head.grand.time, head.grand.vel_desync,...
%             '-', 'Color', [0.7*head_color 0.3], 'LineWidth', 0.25)
        [~,~] = PlotPatch(head.grand.vel_desync_stats.mean, head.grand.vel_desync_stats.std, ...
            head.grand.time_stats.mean, 1, 1, head_color, 0.7*head_color, 0.3, 1);
        
%         plot(wing.grand.time, wing.grand.vel_desync,...
%             '-', 'Color', [0.7*wing_color 0.3], 'LineWidth', 0.25)
        [~,~] = PlotPatch(wing.grand.vel_desync_stats.mean, wing.grand.vel_desync_stats.std, ...
            wing.grand.time_stats.mean, 1, 1, wing_color, 0.7*wing_color, 0.3, 1);
        
linkaxes(ax,'x')
set(ax,'XLim', 0.5*[-1 1])
set(ax,'XTick', -0.5:0.1:0.5)
set(ax(1),'XTick',[])
set(ax,'LineWidth', 1)

%% Deynchronized head-wing saccades 2
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 3 3])
movegui(fig, 'center')
clear ax
head_color = [0 0 1];
wing_color = [1 0 0];
ax(1) = subplot(2,1,1); hold on ; cla
    yyaxis left ; hold on ; cla ; ylabel('Head (°)')
        plot(head.grand.time, head.grand.pos_desync,...
            '-', 'Color', [0.7*head_color 0.3], 'LineWidth', 0.5)  
        [~,~] = PlotPatch(head.grand.pos_desync_stats.mean, head.grand.pos_desync_stats.std, ...
            head.grand.time_stats.mean, 1, 1, head_color, 0.7*head_color, 0.3, 1);
        
    yyaxis right ; hold on ; cla ; ylabel('\DeltaWBA (°)')
        plot(wing.grand.time, wing.grand.pos_desync,...
            '-', 'Color', [0.7*wing_color 0.3], 'LineWidth', 0.5)
        [~,~] = PlotPatch(wing.grand.pos_desync_stats.mean, wing.grand.pos_desync_stats.std, ...
            wing.grand.time_stats.mean, 1, 1, wing_color, 0.7*wing_color, 0.3, 1);
        
ax(2) = subplot(2,1,2); hold on ; cla ; xlabel('Time (s)')
    yyaxis left ; hold on ; cla ; ylabel('Head Velocity (°/s)')
    ylim([-200 600])
        plot(head.grand.time, head.grand.vel_desync,...
            '-', 'Color', [0.7*head_color 0.3], 'LineWidth', 0.5)
    [~,~] = PlotPatch(head.grand.vel_desync_stats.mean, head.grand.vel_desync_stats.std, ...
        head.grand.time_stats.mean, 1, 1, head_color, 0.7*head_color, 0.3, 1);
        
    yyaxis right ; hold on ; cla ; ylabel('\DeltaWBA Velocity (°/s)')
    ylim([-200 600])
        plot(wing.grand.time, wing.grand.vel_desync,...
            '-', 'Color', [0.7*wing_color 0.3], 'LineWidth', 0.5)
        [~,~] = PlotPatch(wing.grand.vel_desync_stats.mean, wing.grand.vel_desync_stats.std, ...
            wing.grand.time_stats.mean, 1, 1, wing_color, 0.7*wing_color, 0.3, 1);
        
ax(1).YAxis(1).Color = head_color;
ax(1).YAxis(2).Color = wing_color;
ax(2).YAxis(1).Color = head_color;
ax(2).YAxis(2).Color = wing_color;
linkaxes(ax,'x')
set(ax,'XLim', 0.5*[-1 1], 'XTick', -0.5:0.1:0.5)
set(ax(1),'XTick',[])
set(ax,'LineWidth', 1)

%% Time Difference
TD = cellfun(@(x) x.TimeDiff, Saccade.head2wing, 'UniformOutput', false);
fly_group = findgroups(Saccade.fly);
fly_sacd_group = cellfun(@(x,y) y*ones(size(x,1),1), TD, num2cell(fly_group), 'UniformOutput', false);
fly_sacd_group = cat(1,fly_sacd_group{:});
TD = cat(1,TD{:});
Y = splitapply(@(x) nanmean(x,1), TD, fly_sacd_group);

fig = figure (3) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 2 2])
movegui(fig, 'center')
clear ax
ww = 1;
ax(ww) = subplot(1,1,ww); axis tight
    bx = boxplot(Y, 'Labels', {'Start','Peak','Stop'}, 'Width', 0.5, 'Symbol', '', 'Whisker', 2);
    ylabel('Time Difference (ms')
    ylim(120*[-1 1])

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2}, 'k');
    end
    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);

set(ax,'LineWidth', 1, 'Box', 'on')

%% Dynamics
Amplitude = cellfun(@(x) x.Amplitude, Saccade.head2wing, 'UniformOutput', false);
Amplitude = cat(1,Amplitude{:});
PeakVel = cellfun(@(x) x.PeakVel, Saccade.head2wing, 'UniformOutput', false);
PeakVel = cat(1,PeakVel{:});
Duration = cellfun(@(x) x.Duration, Saccade.head2wing, 'UniformOutput', false);
Duration = cat(1,Duration{:});

not_nan = ~isnan(Amplitude(:,1));
Amplitude = Amplitude(not_nan,:);
PeakVel = PeakVel(not_nan,:);
Duration = Duration(not_nan,:);

amp_ratio = Amplitude(:,1) ./ Amplitude(:,2);
pkvel_ratio = PeakVel(:,1) ./ PeakVel(:,2);
dur_ratio = Duration(:,1) ./ Duration(:,2);

% amp_ratio_fly = splitapply(@(x) nanmean(x,1), amp_ratio, fly_sacd_group);
% pkvel_ratio = splitapply(@(x) nanmean(x,1), pkvel_ratio, fly_sacd_group);
% dur_ratio = splitapply(@(x) nanmean(x,1), dur_ratio, fly_sacd_group);

fig = figure (4) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 6 2])
movegui(fig, 'center')
clear ax h
ax(1) = subplot(1,3,1); hold on ; xlabel('Amplitude Ratio (°/°)')
    h(1) = histogram(amp_ratio);
    xlim([0 2.5])
    ylabel('Probability')
    
ax(2) = subplot(1,3,2); hold on ; xlabel('Peak Velocity Ratio (°s^{-1}/°s^{-1})')
    h(2) = histogram(pkvel_ratio);
    xlim([0 12])
    
ax(3) = subplot(1,3,3); hold on ; xlabel('Duration Ratio (s/s)')
    h(3) = histogram(dur_ratio);
    xlim([0 1])
    
set(h, 'FaceAlpha', 0.7, 'NumBins', 20,  'Normalization', 'probability', ...
    'FaceColor', [0.2 0.1 0.5], 'EdgeColor', 'w')
linkaxes(ax, 'y')
set(ax,'LineWidth', 1, 'Box', 'off', 'YLim', [-0.01 0.25])

%% Negative Correlation
vel_group = findgroups(Saccade.vel);
ccw = vel_group == 1;
CrossCorr = cellfun(@(x) x.int_acor, Saccade.head2wing, 'UniformOutput', false);
CrossCorr = cat(2,CrossCorr{:});

HeadWin = cellfun(@(x) x.hint, Saccade.head2wing, 'UniformOutput', false);
HeadWin(ccw) = cellfun(@(x) -x, HeadWin(ccw), 'UniformOutput', false);
HeadWin = cat(2,HeadWin{:});

WingWin = cellfun(@(x) x.wint, Saccade.head2wing, 'UniformOutput', false);
WingWin(ccw) = cellfun(@(x) -x, WingWin(ccw), 'UniformOutput', false);
WingWin = cat(2,WingWin{:});

HeadWin_vel = cellfun(@(x) x.hint_vel, Saccade.head2wing, 'UniformOutput', false);
HeadWin_vel(ccw) = cellfun(@(x) -x, HeadWin_vel(ccw), 'UniformOutput', false);
HeadWin_vel = cat(2,HeadWin_vel{:});

WingWin_vel = cellfun(@(x) x.wint_vel, Saccade.head2wing, 'UniformOutput', false);
WingWin_vel(ccw) = cellfun(@(x) -x, WingWin_vel(ccw), 'UniformOutput', false);
WingWin_vel = cat(2,WingWin_vel{:});

TimeWin = cellfun(@(x) x.tint, Saccade.head2wing, 'UniformOutput', false);
TimeWin = cat(2,TimeWin{:});

nsacd = size(CrossCorr,2);
Desync = [];
% Desync.CC = nan(size(CrossCorr,1),nsacd);
% Desync.head = nan(size(HeadWin,1),nsacd);
% Desync.wing = nan(size(WingWin,1),nsacd);
% Desync.head_vel = nan(size(HeadWin,1),nsacd);
% Desync.wing_vel = nan(size(WingWin,1),nsacd);
% Desync.time = nan(size(WingWin,1),nsacd);
cc_zr = nan(nsacd,1);
desync_thresh = 1;
pp = 1;
for s = 1:nsacd
    index_zero = ceil(size(CrossCorr,1)/2); % 0 lag index
    corr_zero = CrossCorr(index_zero,s); % cross-corr at 0 lag
    cc_zr(s) = corr_zero;
    if corr_zero < desync_thresh
        Desync.CC(:,pp) = CrossCorr(:,s);
        Desync.time(:,pp) = TimeWin(:,s);
        Desync.head(:,pp) = HeadWin(:,s) - mean(HeadWin(:,s));
        Desync.wing(:,pp) = WingWin(:,s) - mean(WingWin(:,s));
        Desync.head_vel(:,pp) = HeadWin_vel(:,s);
        Desync.wing_vel(:,pp) = WingWin_vel(:,s);
        
        figure (111) ; cla ; hold on
        plot(Desync.time, Desync.head(:,pp), 'b')
        plot(Desync.time, Desync.wing(:,pp), 'r')
        pause
        
        pp = pp + 1;
    end
end
Desync.stats = structfun(@(x) basic_stats(x,2), Desync, 'UniformOutput', false);

desync_percent = size(Desync.head,2) / sum(~isnan(HeadWin(1,:)))

%% Deynchronized (uncorrelated) head-wing saccades
fig = figure (4) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 3 3])
movegui(fig, 'center')
head_color = [0 0 1];
wing_color = [1 0 0];
clear ax
ax(1) = subplot(2,1,1); hold on ; cla
    plot(Desync.time,Desync.head, 'Color', [0.7*head_color 0.5])
    plot(Desync.time,Desync.wing, 'Color', [0.7*wing_color  0.5])
%     [~,~] = PlotPatch(Desync.stats.head.mean, Desync.stats.head.std, ...
%         Desync.stats.time.mean, 1, 1, head_color, 0.7*head_color, 0.3, 1);
%     [~,~] = PlotPatch(Desync.stats.wing.mean, Desync.stats.wing.std, ...
%         Desync.stats.time.mean, 1, 1, wing_color, 0.7*wing_color, 0.3, 1);
    
ax(2) = subplot(2,1,2); hold on ; cla
    plot(Desync.time, Desync.head_vel, 'Color', [0.7*head_color 0.5])
    plot(Desync.time, Desync.wing_vel, 'Color', [0.7*wing_color  0.5])
    [~,~] = PlotPatch(Desync.stats.head_vel.mean, Desync.stats.head_vel.std, ...
        Desync.stats.time.mean, 1, 1, head_color, 0.7*head_color, 0.3, 1);
    [~,~] = PlotPatch(Desync.stats.wing_vel.mean, Desync.stats.wing_vel.std, ...
        Desync.stats.time.mean, 1, 1, wing_color, 0.7*wing_color, 0.3, 1);
    
    xlabel('Time (s)')
    
linkaxes(ax, 'x')
set(ax,'LineWidth', 1, 'Box', 'off', 'XLim', 0.5*[-1 1])

end