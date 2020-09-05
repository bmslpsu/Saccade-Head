function [] = HeadWing_Timing_Static()
%% HeadWing_Timing_Static:

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
nwave = length(unique(Saccade.wave));
uwave = unique(Saccade.wave);
nfly = length(unique(Saccade.fly));
ufly = unique(Saccade.fly);
for v = 1:1
    %waveI = uwave(v) == Saccade.wave;
    for a = 1:nfly
        flyI = ufly(a) == Saccade.fly;
        sI = flyI;
        hsacd = cellfun(@(x) x.hsacd, Saccade.(hw_field)(sI,:), 'UniformOutput', true);
        wsacd = cellfun(@(x) x.wsacd, Saccade.(hw_field)(sI,:), 'UniformOutput', true);
        dir  = cellfun(@(x) x.SACD.Direction, Saccade.head_saccade(sI,:), 'UniformOutput', false);
        dir = cat(1,dir{:})';
        get_fields = string(fieldnames(hsacd));
        get_fields = get_fields(1:end-1);
        stats_fields = get_fields + "_stats";
        for f = 1:length(get_fields)
            head_all = cat(2, hsacd.(get_fields(f)));
            if f~=1
                dir_all = repmat(dir, size(head_all,1), 1);
            else
                dir_all = 1;
            end
            head.fly(a,v).(get_fields(f)) = head_all .* dir_all;
            wing.fly(a,v).(get_fields(f)) = cat(2, wsacd.(get_fields(f))) .* dir_all;
            head.fly(a,v).(stats_fields(f)) = basic_stats(head.fly(a,v).(get_fields(f)),2);
            wing.fly(a,v).(stats_fields(f)) = basic_stats(wing.fly(a,v).(get_fields(f)),2);
        end
    end
end

mean_win = round(0.2*200);
for v = 1:1
    for f = 1:length(stats_fields)
        temp = cat(2,head.fly(:,v).(stats_fields(f)));
        head.grand(v).(get_fields(f)) = cat(2,temp.mean);
        temp = cat(2,wing.fly(:,v).(stats_fields(f)));
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
    all_head = head.fly(a).pos;
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
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 3 3])
movegui(fig, 'center')
clear ax
head_color = [0 0 1];
wing_color = [1 0 0];
ax(1) = subplot(2,1,1); hold on ; cla ; ylim(20*[-1 1])
        hold on ; cla ; ylabel('Head (°)')
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
        
linkaxes(ax,'x')
set(ax,'XLim', 0.5*[-1 1])
% set(ax,'XTick', -0.06:0.01:0)
set(ax(1),'XTickLabels',[])
set(ax,'LineWidth', 1)

%% Deynchronized head-wing saccades
fig = figure (3) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', [2 2 3 3])
movegui(fig, 'center')
clear ax
head_color = [0 0 1];
wing_color = [1 0 0];
ax(1) = subplot(2,1,1); hold on ; cla
        hold on ; cla ; ylabel('Head (°)')
        ylim(10*[-1 1])
        plot(head.grand.time, head.grand.pos_desync,...
            '-', 'Color', [0.7*head_color 0.3], 'LineWidth', 0.25)
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

end