function [] = Static_Saccade_TriggerPos()
%% Static_Saccade_TriggerPos:
root = 'H:\DATA\Rigid_Data\';

[CW,PATH_CW] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select CW data', root, 'MultiSelect','off');

[CCW,PATH_CCW] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select CCW data', root, 'MultiSelect','off');

CW = load(fullfile(PATH_CW,CW),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','D','I','U','N');
CCW = load(fullfile(PATH_CCW,CCW),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','D','I','U','N');

clearvars -except CW CCW

%% Saccade Trigger Polar Plot
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Saccade Trigger';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax

ax(1) = subplot(1,2,1,polaraxes); grid off ; axis tight
    h(1) = polarhistogram(deg2rad(CW.SACCADE_STATS.StartPos),75,...
        'FaceColor','g','FaceAlpha',0.9,'Normalization','Probability'); hold on
    h(2) = polarhistogram(deg2rad(CW.SACCADE_STATS.EndPos),'BinEdges', h(1).BinEdges, ...
        'FaceColor','r','FaceAlpha',0.9,'Normalization','Probability');
    ax(1).ThetaAxis.Label.String = 'Head Position (°)';
    title('CW')
    
ax(2) = subplot(1,2,2,polaraxes); grid off ; axis tight
    h(3) = polarhistogram(deg2rad(CCW.SACCADE_STATS.StartPos),75,...
        'FaceColor','g','FaceAlpha',0.9,'Normalization','Probability'); hold on
    h(4) = polarhistogram(deg2rad(CCW.SACCADE_STATS.EndPos),'BinEdges', h(1).BinEdges, ...
        'FaceColor','r','FaceAlpha',0.9,'Normalization','Probability');
    ax(2).ThetaAxis.Label.String = 'Head Position (°)';
    title('CCW')

set(h,'EdgeColor','none')
set(ax,'FontSize',8);
set(ax,'Color','w');
set(ax,'ThetaLim',[-20 20]);
set(ax,'RLim',[0 0.11]);
set(ax,'ThetaDir','clockwise')
set(ax,'ThetaTick',-20:10:20);
set(ax,'ThetaZeroLocation','top');
% linkaxes(ax)
% leg = legend('Start','End','Location','North');
% leg.Location = 'northwest';
% leg.Box = 'off';

set(ax,'LineWidth',1,'FontWeight','bold')

%% Saccade Trigger Polar Plot Normalized
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Saccade Trigger';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h

ax(1) = subplot(1,1,1,polaraxes); grid off ; axis tight
    h(1) = polarhistogram(deg2rad([CW.SACCADE_STATS.StartPos;CCW.SACCADE_STATS.StartPos]),100,...
        'FaceColor','g','FaceAlpha',0.9,'Normalization','Probability'); hold on
    h(2) = polarhistogram(deg2rad([CW.SACCADE_STATS.EndPos;CCW.SACCADE_STATS.EndPos]),'BinEdges', h(1).BinEdges, ...
        'FaceColor','r','FaceAlpha',0.9,'Normalization','Probability');
    ax(1).ThetaAxis.Label.String = 'Head Position (°)';
    
set(h,'EdgeColor','none')
set(ax,'FontSize',8);
set(ax,'Color','w');
set(ax,'ThetaLim',[-20 20]);
% set(ax,'RLim',[0 300]);
set(ax,'ThetaDir','clockwise')
set(ax,'ThetaTick',-20:10:20);
set(ax,'ThetaZeroLocation','top');
% leg = legend('Start','End','Location','North');
% leg.Location = 'northwest';
% leg.Box = 'off';

set(ax,'LineWidth',1,'FontWeight','bold')

end