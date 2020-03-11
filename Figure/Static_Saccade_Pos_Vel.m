function [] = Static_Saccade_Pos_Vel()
%% Static_Saccade_Pos_Vel:
root = 'H:\DATA\Rigid_Data\';

[CW,PATH_CW] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select CW data', root, 'MultiSelect','off');

[CCW,PATH_CCW] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select CCW data', root, 'MultiSelect','off');

CW = load(fullfile(PATH_CW,CW),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','D','I','U','N');
CCW = load(fullfile(PATH_CCW,CCW),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','D','I','U','N');

clearvars -except CW CCW

clms = CW.N.wave;
CC = repmat(jet(clms),2,1);

Wave = CW.U.wave{1};
WaveLabel = string(Wave);
WaveLabel(1) = 'All Off';
WaveLabel(5) = 'All On';
WaveLabel(6) = 'Random';
WaveLabel([2 3 4]) = WaveLabel([2 3 4]) + '°';
wave_order = [2 3 4 6 5 1];

%% Saccade Position %%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 1.1384*clms*(4/3) 3/2];
FIG.Name = 'Static Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h
ax = gobjects(CW.N.wave,1);

med_span = 5;
pp = 1;
for jj = wave_order
    ax(pp) = subplot(ceil(CW.N.wave/clms),clms,pp) ; hold on
    title(WaveLabel(jj))
    
    % CW
    time = 1000*CW.GRAND.normpeak_saccade(jj).time;
    med_time = nanmedian(time,2);
    pos = CW.GRAND.normpeak_saccade(jj).position;
    med_pos = nanmedian(pos,2);
    std_pos = nanstd(pos,[],2);
    
	cent = find(med_time==0);
    span = (cent - med_span):( cent + med_span);
    
    h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
                                        
    h.patch = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
        CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    h.patch.EdgeColor = CC(jj,:);
    
    % CCW
   	time = 1000*CCW.GRAND.normpeak_saccade(jj).time;
    med_time = nanmedian(time,2);
    pos = CCW.GRAND.normpeak_saccade(jj).position;
    med_pos = nanmedian(pos,2);
    std_pos = nanstd(pos,[],2);
    
	cent = find(med_time==0);
    span = (cent - med_span):( cent + med_span);
    
    h.trial = plot(time, pos, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
                                        
    h.patch = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
        CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    h.patch.EdgeColor = CC(jj,:);
    
    pp = pp +1;
end
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.05*[-1 1],'YLim',20*[-1 1])
set(ax(2:end),'YTickLabels',[])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC{:}], 'String', 'Head Position (°)')
linkaxes(ax,'xy')


%% Saccade Velocity %%
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 1.1384*clms*(4/3) 3/2];
FIG.Name = 'Static Saccade Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
ax = gobjects(CW.N{1,3},1);
clear ax h

med_span = 5;
pp = 1;
for jj = wave_order
    ax(pp) = subplot(ceil(CW.N.wave/clms),clms,pp) ; hold on
    title(WaveLabel(jj))
    
    % CW
    time = 1000*CW.GRAND.normpeak_saccade(jj).time;
    med_time = nanmedian(time,2);
    vel = CW.GRAND.normpeak_saccade(jj).velocity;
    med_vel = nanmedian(vel,2);
    std_vel = nanstd(vel,[],2);
    
	cent = find(med_time==0);
    span = (cent - med_span):( cent + med_span);
    
    h.trial = plot(time, vel, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
                                        
    h.patch = PlotPatch(med_vel(span), std_vel(span), med_time(span), 1, 1, ...
        CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    h.patch.EdgeColor = CC(jj,:);
    
    % CCW
    time = 1000*CCW.GRAND.normpeak_saccade(jj).time;
    med_time = nanmedian(time,2);
    vel = CCW.GRAND.normpeak_saccade(jj).velocity;
    med_vel = nanmedian(vel,2);
    std_vel = nanstd(vel,[],2);
    
	cent = find(med_time==0);
    span = (cent - med_span):( cent + med_span);
    
    h.trial = plot(time, vel, 'Color', [0.7*CC(jj,:) , 0.2], 'LineWidth', 0.5);
                                        
    h.patch = PlotPatch(med_vel(span), std_vel(span), med_time(span), 1, 1, ...
        CC(jj,:), [0.7 0.7 0.7], 0.4, 3);
    h.patch.EdgeColor = CC(jj,:);
    
    pp = pp +1;    
end
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.05*[-1 1],'YLim',1000*[-1 1])
set(ax(2:end),'YTickLabels',[])
XLabelHC = get(ax, 'XLabel');
set([XLabelHC{:}], 'String', 'Time (ms)')
YLabelHC = get(ax(1), 'YLabel');
set([YLabelHC], 'String', 'Head Velocity (°)')
linkaxes(ax,'xy')

end