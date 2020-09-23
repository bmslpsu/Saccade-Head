function [] = Ramp_saccade_pos_vel()
%% Ramp_saccade_pos_vel:
root = 'H:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE', 'HEAD_SACCADE_STATS', 'FLY','GRAND','D','I','U','N')

%% Combine Speeds
clearvars -except SACCADE FLY GRAND Stim D I U N
Fs = round(SACCADE.head_saccade{1}.Fs);
n_speed = N.vel/2;
CC = repmat(hsv(n_speed),2,1);
Vel = U{1,3}{1};

norm_fields = {'time','position','velocity'};
% norm_fields_stats = string(norm_fields) + "_stats";
center = 0;
even = true;
dim = 1;
clear GRAND_comb_speed FLY_comb_speed
for jj = 1:n_speed
    speedI = [jj , jj + n_speed];
    for kk = 1:N.fly
        cw = FLY.normpeak_saccade(kk,speedI(1));
        cw.position = -cw.position;
        cw.velocity = -cw.velocity;
        ccw = GRAND.normpeak_saccade(speedI(2));
        FLY_comb_speed(kk,jj) = struct_center([cw , ccw], ...
                                            center, even, dim, norm_fields);
    end
    cw = GRAND.normpeak_saccade(speedI(1));
    cw.position = -cw.position;
    cw.velocity = -cw.velocity;
    ccw = GRAND.normpeak_saccade(speedI(2));
    GRAND_comb_speed(jj) = struct_center([cw , ccw], ...
                                        center, even, dim, norm_fields);
end
for jj = 1:n_speed
    for kk = 1:N.fly
        FLY_comb_speed(kk,jj).time_stats = basic_stats(FLY_comb_speed(kk,jj).time,2);
        FLY_comb_speed(kk,jj).position_stats = basic_stats(FLY_comb_speed(kk,jj).position,2);
        FLY_comb_speed(kk,jj).velocity_stats = basic_stats(FLY_comb_speed(kk,jj).velocity,2);
    end
    %fly_mean = [FLY_comb_speed(:,jj).position_stats]';
    
 	%GRAND_comb_speed(jj)
end

%% Saccade Position by Speed
FIG = figure (102) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 n_speed*(4/3) 1.5];
FIG.Name = 'Saccade Position';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h
ax = gobjects(n_speed,1);
med_span = Fs*0.025;
for jj = 1:n_speed
    ax(jj) = subplot(1,n_speed,jj) ; hold on
    title([num2str(Vel(jj)) ' (°/s)']);
  	
    time = 1000 * [GRAND_comb_speed(jj).time];
    med_time = nanmedian(time,2);
    pos = GRAND_comb_speed(jj).position;
    med_pos = nanmedian(pos,2);
    std_pos = nanstd(pos,[],2);
    
	cent = find(med_time==0);
    span = (cent - med_span):( cent + med_span);
    
    h.trial = plot(time, pos, 'Color', [0.5 0.5 0.5 , 0.2], 'LineWidth', 0.25);
    
%     for kk = 1:N.fly
%         plot(1000*FLY_comb_speed(kk,jj).time_stats.mean(span,:), ...
%             FLY_comb_speed(kk,jj).position_stats.mean(span,:), 'Color', 0.7*CC(jj,:))
%     end
                                        
    h.patch = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
        CC(jj,:), 0.5*CC(jj,:), 0.3, 1);
end
set(ax,'LineWidth',1,'FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.07*[-1 1],'YLim',25*[-1 1])
XLabelHC = get(ax(1), 'XLabel');
set([XLabelHC], 'String', 'Time (ms)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Head Position (°)')
set(ax(2:end), 'XTickLabels', [], 'YTickLabels', [])
linkaxes(ax,'xy')

%% Saccade Velocity by Speed
FIG = figure (103) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 n_speed*(4/3) 1.5];
FIG.Name = 'Saccade Velocity';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h
ax = gobjects(n_speed,1);
med_span = Fs*0.025;
for jj = 1:n_speed
    ax(jj) = subplot(1,n_speed,jj) ; hold on
    title([num2str(Vel(jj)) ' (°/s)']);
        
    time = 1000 * [GRAND_comb_speed(jj).time];
    med_time = nanmedian(time,2);
    vel = GRAND_comb_speed(jj).velocity;
    med_pos = nanmedian(vel,2);
    std_pos = nanstd(vel,[],2);
    
	cent = find(med_time==0);
    span = (cent - med_span):( cent + med_span);
    
    h.trial = plot(time, vel, 'Color', [0.5 0.5 0.5 , 0.2], 'LineWidth', 0.25);

%     for kk = 1:N.fly
%         plot(1000*FLY_comb_speed(kk,jj).time_stats.mean(span,:), ...
%             FLY_comb_speed(kk,jj).velocity_stats.mean(span,:), 'Color', 0.7*CC(jj,:))
%     end
                                        
    h.patch = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
        CC(jj,:),  0.5*CC(jj,:), 0.3, 1);
end
set(ax,'LineWidth',1,'FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.07*[-1 1],'YLim',[-200 1100])
XLabelHC = get(ax(1), 'XLabel');
set([XLabelHC], 'String', 'Time (ms)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Head Speed (°/s)')
set(ax(2:end), 'XTickLabels', [], 'YTickLabels', [])
linkaxes(ax,'xy')

end