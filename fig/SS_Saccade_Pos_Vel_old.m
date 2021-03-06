function [] = SS_Saccade_Pos_Vel()
%% SS_Saccade_Pos_Vel:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'PATH','COUNT','SACCADE','GRAND','D','I','U','N')

clearvars -except clms CC Freq PATH COUNT SACCADE SACCADE_STATS FLY GRAND D I U N

CC = hsv(N.freq);

Freq = U.freq{1};

%% Saccade Position %%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 N.freq*(4/3) 3/2];
FIG.Name = 'Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.freq,1);
clear ax h

med_span = 5;
for jj = 1:N.freq
    ax(jj) = subplot(1,N.freq,jj) ; hold on
    title([num2str(Freq(jj)) ' Hz']);
    
    time = 1000*GRAND.normpeak_saccade(jj).time;
    med_time = nanmedian(time,2);
    pos = GRAND.normpeak_saccade(jj).position;
    med_pos = nanmedian(pos,2);
    std_pos = nanstd(pos,[],2);
    
	cent = find(med_time==0);
    span = (cent - med_span):( cent + med_span);
    
    h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
                                        
    h.patch = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
        CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    h.patch.EdgeColor = CC(jj,:);
end
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.05*[-1 1],'YLim',20*[-1 1])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Head Position (�)')
linkaxes(ax,'xy')

%% Saccade Velocity %%
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 N.freq*(4/3) 3/2];
FIG.Name = 'Saccade Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(N.freq,1);
clear ax h

med_span = 5;
for jj = 1:N.freq
    ax(jj) = subplot(1,N.freq,jj) ; hold on
    title([num2str(Freq(jj)) ' Hz'])
    
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
    
    if sign(Freq(jj))==1
        ax(jj).YLim = 1000*[-1 0.2];
    else
        ax(jj).YLim =  1000*[-0.2 1];
    end
end
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.05*[-1 1])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Head Velocity (�/s)')
linkaxes(ax,'xy')

%% Saccade Velocity: Both Directions %%
FIG = figure (3) ; clf
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
    title([num2str(Freq(jj)) ' (�/s)'])
    
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
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.05*[-1 1], 'YLim', 1000*[-1 1])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Head Velocity (�)')
linkaxes(ax,'xy')
set(ax(2:end),'YTickLabels',[])

end