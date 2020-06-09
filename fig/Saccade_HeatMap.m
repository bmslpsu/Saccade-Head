function [] = Saccade_HeatMap()
%% Saccade_HeatMap:

root = 'H:\DATA\Rigid_Data\';

[Move,MovePATH] = uigetfile({'*.mat'},'Select moving data', root, 'MultiSelect','off');
[Static,StaticPATH] = uigetfile({'*.mat'},'Select static data', root, 'MultiSelect','off');

Move = load(fullfile(MovePATH,Move),'PATH','COUNT','SACCADE_STATS','D','I','U','N');
Static = load(fullfile(StaticPATH,Static),'PATH','COUNT','SACCADE_STATS','D','I','U','N');

%% Heat Maps
FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 6 6];
movegui(FIG,'center')
nbins = 50;
colormap(jet)

clear h
ax(1) = subplot(3,2,1);
    X = abs(Move.SACCADE_STATS.Amplitude);
    Y = abs(Move.SACCADE_STATS.PeakVel);
    h(1) = histogram2(X, Y, [nbins,nbins], 'Normalization','Probability','FaceColor','flat',...
        'ShowEmptyBins','on','DisplayStyle','tile');
    grid off
    xlabel('Amplitude (°)')
    ylabel('Peak Velocity (°/s)')
    zlabel('Probability')
	title('Moving')

ax(2) = subplot(3,2,2);
    X = abs(Static.SACCADE_STATS.Amplitude);
    Y = abs(Static.SACCADE_STATS.PeakVel);
    h(2) = histogram2(X, Y, h(1).XBinEdges, h(1).YBinEdges, 'Normalization','Probability','FaceColor','flat',...
        'ShowEmptyBins','on','DisplayStyle','tile');
    grid off
    xlabel('Amplitude (°)')
    ylabel('Peak Velocity (°/s)')
    zlabel('Probability')
    title('Static')
    
ax(3) = subplot(3,2,3);
    X = abs(Move.SACCADE_STATS.Amplitude);
    Y = abs(Move.SACCADE_STATS.Duration);
    dur_bins = 0:0.005:0.1;
    h(3) = histogram2(X, Y, h(1).XBinEdges, dur_bins, 'Normalization','Probability','FaceColor','flat',...
        'ShowEmptyBins','on','DisplayStyle','tile');
    grid off
    xlabel('Amplitude (°)')
    ylabel('Duration (s)')
    zlabel('Probability')
    
ax(4) = subplot(3,2,4);
    X = abs(Static.SACCADE_STATS.Amplitude);
    Y = abs(Static.SACCADE_STATS.Duration);
    dur_bins = 0:0.005:0.1;
    h(4) = histogram2(X, Y, h(1).XBinEdges, dur_bins, 'Normalization','Probability','FaceColor','flat',...
        'ShowEmptyBins','on','DisplayStyle','tile');
    grid off
    xlabel('Amplitude (°)')
    ylabel('Duration (s)')
    zlabel('Probability')
    
ax(5) = subplot(3,2,5);
    X = abs(Move.SACCADE_STATS.PeakVel);
    Y = abs(Move.SACCADE_STATS.Duration);
    dur_bins = 0:0.005:0.1;
    h(5) = histogram2(X, Y, h(1).YBinEdges, dur_bins, 'Normalization','Probability','FaceColor','flat',...
        'ShowEmptyBins','on','DisplayStyle','tile');
    grid off
    xlabel('Peak Velocity (°/s)')
    ylabel('Duration (s)')
    zlabel('Probability')
    
ax(6) = subplot(3,2,6);
    X = abs(Static.SACCADE_STATS.PeakVel);
    Y = abs(Static.SACCADE_STATS.Duration);
    dur_bins = 0:0.005:0.1;
    h(6) = histogram2(X, Y, h(1).YBinEdges, dur_bins, 'Normalization','Probability','FaceColor','flat',...
        'ShowEmptyBins','on','DisplayStyle','tile');
    grid off
    xlabel('Peak Velocity (°/s)')
    ylabel('Duration (s)')
    zlabel('Probability')
    
    
set(ax,'LineWidth',1,'FontWeight','bold')
set(h,'EdgeColor','none')

hp = get(subplot(3,2,6),'Position');
cbar = colorbar('Position', [hp(1)+hp(3)+0.02  0.15+hp(2)  0.03  hp(2)+hp(3)*1]);
cbar.Label.String = 'Probability';

end