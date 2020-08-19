function [] = Saccade_Pos_Vel()
%% Saccade_Pos_Vel:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE','FLY','GRAND','D','I','U','N')

n_speed = N.vel/2;
CC = repmat(hsv(n_speed),2,1);

Vel = U{1,3}{1};

% clearvars -except clms CC Vel PATH COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

%% Combine Speeds
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
med_span = 5;
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
    
    h.trial = plot(time, pos, 'Color', [0.5 0.5 0.5 , 0.2], 'LineWidth', 0.5);
    
%     for kk = 1:N.fly
%         plot(1000*FLY_comb_speed(kk,jj).time_stats.mean(span,:), ...
%             FLY_comb_speed(kk,jj).position_stats.mean(span,:), 'Color', 0.7*CC(jj,:))
%     end
                                        
    h.patch = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
        CC(jj,:),  0.5*CC(jj,:), 0.3, 2);
end
set(ax,'LineWidth',1,'FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.05*[-1 1],'YLim',25*[-1 1])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Head Position (°)')
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
med_span = 5;
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
    
    h.trial = plot(time, vel, 'Color', [0.5 0.5 0.5 , 0.2], 'LineWidth', 0.5);

%     for kk = 1:N.fly
%         plot(1000*FLY_comb_speed(kk,jj).time_stats.mean(span,:), ...
%             FLY_comb_speed(kk,jj).velocity_stats.mean(span,:), 'Color', 0.7*CC(jj,:))
%     end
                                        
    h.patch = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
        CC(jj,:),  0.5*CC(jj,:), 0.3, 2);
end
set(ax,'LineWidth',1,'FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.05*[-1 1],'YLim',[-500 1500])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Head Speed (°/s)')
linkaxes(ax,'xy')

%% Saccade Position %%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 clms*(4/3) 3];
FIG.Name = 'Saccade Position';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h
ax = gobjects(N{1,3},1);
med_span = 5;
for jj = 1:N.vel
    ax(jj) = subplot(ceil(N.vel/clms),clms,jj) ; hold on
    title([num2str(Vel(jj)) ' (°/s)']);
       
    time = 1000 * [GRAND.normpeak_saccade(jj).time];
    med_time = nanmedian(time,2);
    pos = GRAND.normpeak_saccade(jj).position;
    med_pos = nanmedian(pos,2);
    std_pos = nanstd(pos,[],2);
    
	cent = find(med_time==0);
    span = (cent - med_span):( cent + med_span);
    
    h.trial = plot(time, pos, 'Color', [0.5 0.5 0.5 , 0.2], 'LineWidth', 0.5);
                                        
    h.patch = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
        CC(jj,:),  0.5*CC(jj,:), 0.3, 2);
    %h.patch.EdgeColor = CC(jj,:);
end
set(ax,'LineWidth',1,'FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.05*[-1 1],'YLim',20*[-1 1])
XLabelHC = get(ax(6:10), 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax([1,clms+1]), 'YLabel');
set([YLabelHC{:}], 'String', 'Head Position (°)')
linkaxes(ax(1:clms),'y')
linkaxes(ax((clms+1):N.vel),'y')
linkaxes(ax,'x')

%% Saccade Velocity %%
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 clms*(4/3) 3];
FIG.Name = 'Saccade Velocity';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h
ax = gobjects(N{1,3},1);
med_span = 5;
for jj = 1:N.vel
    ax(jj) = subplot(ceil(N{1,3}/clms),clms,jj) ; hold on
    title([num2str(Vel(jj)) ' (°/s)'])
    
    time = 1000*GRAND.normpeak_saccade(jj).time;
    med_time = nanmedian(time,2);
    vel = GRAND.normpeak_saccade(jj).velocity;
    med_vel = nanmedian(vel,2);
    std_vel = nanstd(vel,[],2);
%     med_vel = GRAND.normpeak_saccade(jj).velocity_stats.mean;
%     std_vel = GRAND.normpeak_saccade(jj).velocity_stats.std;
    
	cent = find(med_time==0);
    span = (cent - med_span):( cent + med_span);
    
    h.trial = plot(time, vel, 'Color', [0.5 0.5 0.5 , 0.2], 'LineWidth', 0.5);
                                        
    h.patch = PlotPatch(med_vel(span), std_vel(span), med_time(span), 1, 1, ...
        CC(jj,:), 0.5*CC(jj,:), 0.3, 2);
    %h.patch.EdgeColor = CC(jj,:);
    
    if sign(Vel(jj))==1
        ax(jj).YLim = 1000*[-1.2 0.4];
    else
        ax(jj).YLim =  1000*[-0.4 1.2];
    end
end
set(ax,'LineWidth',1,'FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.05*[-1 1])
XLabelHC = get(ax(6:10), 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax([1,clms+1]), 'YLabel');
set([YLabelHC{:}], 'String', 'Head Velocity (°)')
linkaxes(ax(1:clms),'y')
linkaxes(ax((clms+1):N{1,3}),'y')
linkaxes(ax,'x')

%% Saccade Position: Both Directions %%
FIG = figure (3) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 clms*(4/3) 3/2];
FIG.Name = 'Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel/2,1);
clear ax h

med_span = 5;
for jj = 1:N.vel/2
    ax(jj) = subplot(1,clms,jj) ; hold on
    title([num2str(Vel(jj)) ' (°/s)'])
    
    time = 1000*GRAND.normpeak_saccade(jj).time;
    med_time = nanmedian(time,2);
    vel = GRAND.normpeak_saccade(jj).position;
    med_vel = nanmedian(vel,2);
    std_vel = nanstd(vel,[],2);
    
	cent = find(med_time==0);
    span = (cent - med_span):( cent + med_span);
    
    h.trial = plot(time, vel, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
                                        
    h.patch = PlotPatch(med_vel(span), std_vel(span), med_time(span), 1, 1, ...
        CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    h.patch.EdgeColor = CC(jj,:);
    
    njj = jj + N.vel/2;
    time = 1000*GRAND.normpeak_saccade(njj).time;
    med_time = nanmedian(time,2);
    vel = GRAND.normpeak_saccade(njj).position;
    med_vel = nanmedian(vel,2);
    std_vel = nanstd(vel,[],2);
    
	cent = find(med_time==0);
    span = (cent - med_span):( cent + med_span);
    
    h.trial = plot(time, vel, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
                                        
    h.patch = PlotPatch(med_vel(span), std_vel(span), med_time(span), 1, 1, ...
        CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    h.patch.EdgeColor = CC(jj,:);
end
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.05*[-1 1],'YLim',20*[-1 1])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Head Velocity (°/s)')
linkaxes(ax,'xy')
set(ax(2:end),'YTickLabels',[])

%% Saccade Velocity: Both Directions %%
FIG = figure (4) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 clms*(4/3) 3/2];
FIG.Name = 'Saccade Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.vel/2,1);
clear ax h

med_span = 5;
for jj = 1:N.vel/2
    ax(jj) = subplot(1,clms,jj) ; hold on
    title([num2str(Vel(jj)) ' (°/s)'])
    
    time = 1000*GRAND.normpeak_saccade(jj).time;
    med_time = nanmedian(time,2);
    vel = GRAND.normpeak_saccade(jj).velocity;
    med_vel = nanmedian(vel,2);
    std_vel = nanstd(vel,[],2);
    
	cent = find(med_time==0);
    span = (cent - med_span):( cent + med_span);
    
    h.trial = plot(time, vel, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
                                        
    h.patch = PlotPatch(med_vel(span), std_vel(span), med_time(span), 1, 1, ...
        CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    h.patch.EdgeColor = CC(jj,:);
    
    njj = jj + N.vel/2;
    time = 1000*GRAND.normpeak_saccade(njj).time;
    med_time = nanmedian(time,2);
    vel = GRAND.normpeak_saccade(njj).velocity;
    med_vel = nanmedian(vel,2);
    std_vel = nanstd(vel,[],2);
    
	cent = find(med_time==0);
    span = (cent - med_span):( cent + med_span);
    
    h.trial = plot(time, vel, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
                                        
    h.patch = PlotPatch(med_vel(span), std_vel(span), med_time(span), 1, 1, ...
        CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    h.patch.EdgeColor = CC(jj,:);
    
    
end
set(ax,'LineWidth',1,'FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.05*[-1 1], 'YLim', 1100*[-1 1])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Head Velocity (°/s)')
linkaxes(ax,'xy')
set(ax(2:end),'YTickLabels',[])

end