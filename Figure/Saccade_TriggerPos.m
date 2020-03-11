function [] = Saccade_TriggerPos()
%% Saccade_TriggerPos:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N')

%% Saccade Trigger Polar Plot
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 4];
FIG.Name = 'Saccade Trigger';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax

stim_dir = sign(SACCADE_STATS.Vel);
CW  = stim_dir ==  1;
CCW = stim_dir == -1;

ax(1) = subplot(1,2,1,polaraxes); grid off ; axis tight
    h(1) = polarhistogram(deg2rad(SACCADE_STATS.StartPos(CW)),100,...
        'FaceColor','g','FaceAlpha',0.9,'Normalization','Probability'); hold on
    h(2) = polarhistogram(deg2rad(SACCADE_STATS.EndPos(CW)),'BinEdges', h(1).BinEdges, ...
        'FaceColor','r','FaceAlpha',0.9,'Normalization','Probability');
    title('CW Stimulus Direction   ===>')
    ax(1).ThetaAxis.Label.String = 'Head Position (°)';
    
ax(2) = subplot(1,2,2,polaraxes); grid off ; axis tight
    h(3) = polarhistogram(deg2rad(SACCADE_STATS.StartPos(CCW)),'BinEdges', h(1).BinEdges, ...
        'FaceColor','g','FaceAlpha',0.9,'Normalization','Probability'); hold on
    h(4) = polarhistogram(deg2rad(SACCADE_STATS.EndPos(CCW)),'BinEdges', h(1).BinEdges, ...
        'FaceColor','r','FaceAlpha',0.9,'Normalization','Probability');
    title('CCW Stimulus Direction  <===')
    ax(2).ThetaAxis.Label.String = 'Head Position (°)';

set(h,'EdgeColor','none')
set(ax,'FontSize',8);
set(ax,'Color','w');
set(ax,'ThetaLim',[-20 20]);
% set(ax,'RLim',[0 300]);
set(ax,'ThetaDir','clockwise')
set(ax,'ThetaTick',-20:10:20);
set(ax,'ThetaZeroLocation','top');
leg = legend('Start','End','Location','North');
leg.Location = 'northwest';
leg.Box = 'off';

set(ax,'LineWidth',1,'FontWeight','bold')

% set(ax,'LineWidth',1,'FontWeight','bold','Color','w','GridColor','k','GridAlpha',1,...
%     'MinorGridAlpha',1,'MinorGridColor','k')

%% Saccade Trigger Polar Plot, normalized to stimulus direction
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 3 4];
FIG.Name = 'Normalized Saccade Trigger';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h

stim_dir = sign(SACCADE_STATS.Vel);
CW  = stim_dir ==  1;
CCW = stim_dir == -1;

ax(1) = subplot(1,1,1,polaraxes); grid off ; axis tight
    Trig_Positive = [ SACCADE_STATS.StartPos(CW); -SACCADE_STATS.StartPos(CCW)];
    Trig_Negative = [ SACCADE_STATS.EndPos(CW);  -SACCADE_STATS.EndPos(CCW)];
    h(1) = polarhistogram(deg2rad(Trig_Positive),100,...
        'FaceColor','g','FaceAlpha',0.9,'Normalization','Probability'); hold on
    h(2) = polarhistogram(deg2rad(Trig_Negative),'BinEdges', h(1).BinEdges, ...
        'FaceColor','r','FaceAlpha',0.9,'Normalization','Probability');
    title('CW Stimulus Direction   ===>')
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