function [] = Saccade_2DHist()
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
FIG.Position = [2 2 5 6];
movegui(FIG,'center')
colormap(jet)

edges.dur = (0:5:100)';
edges.amp = (7:0.5:30)';
edges.pkv = (300:10:1000)';

h = gobjects(6,1);
ax(1) = subplot(3,2,1);
    X = abs(Move.SACCADE_STATS.Amplitude);
    Y = abs(Move.SACCADE_STATS.PeakVel);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    [rho_move.amp_pkv,pval_move.amp_pkv] = corr(X,Y);
    h(1) = histogram2(X, Y, edges.amp, edges.pkv, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'ShowEmptyBins','on', 'Normalization','Probability','DisplayStyle','tile');
    axis tight
    grid off
    xlabel('Amplitude (°)')
    ylabel('Peak Velocity (°/s)')
    zlabel('Probability')
	title('Moving')

ax(2) = subplot(3,2,2);
    X = abs(Static.SACCADE_STATS.Amplitude);
    Y = abs(Static.SACCADE_STATS.PeakVel);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    [rho_static.amp_pkv,pval_static.amp_pkv] = corr(X,Y);
    h(2) = histogram2(X, Y, edges.amp, edges.pkv, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'ShowEmptyBins','on', 'Normalization','Probability','DisplayStyle','tile');
    grid off
    xlabel('Amplitude (°)')
    ylabel('Peak Velocity (°/s)')
    zlabel('Probability')
    title('Static')
   
ax(3) = subplot(3,2,3);
    X = abs(Move.SACCADE_STATS.Amplitude);
    Y = 1000*abs(Move.SACCADE_STATS.Duration);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    [rho_move.amp_dur,pval_move.amp_dur] = corr(X,Y);
    h(3) = histogram2(X, Y, edges.amp, edges.dur, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'ShowEmptyBins','on', 'Normalization','Probability','DisplayStyle','tile');
    grid off
    xlabel('Amplitude (°)')
    ylabel('Duration (s)')
    zlabel('Probability')
    
ax(4) = subplot(3,2,4);
    X = abs(Static.SACCADE_STATS.Amplitude);
    Y = 1000*abs(Static.SACCADE_STATS.Duration);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    [rho_static.amp_dur,pval_static.amp_dur] = corr(X,Y);
    h(4) = histogram2(X, Y, edges.amp, edges.dur, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'ShowEmptyBins','on', 'Normalization','Probability','DisplayStyle','tile');
    grid off
    xlabel('Amplitude (°)')
    ylabel('Duration (s)')
    zlabel('Probability')
    
ax(5) = subplot(3,2,5);
    X = abs(Move.SACCADE_STATS.PeakVel);
    Y = 1000*abs(Move.SACCADE_STATS.Duration);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    [rho_move.pkv_dur,pval_move.pkv_dur] = corr(X,Y);
    h(5) = histogram2(X, Y, edges.pkv, edges.dur, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'ShowEmptyBins','on', 'Normalization','Probability','DisplayStyle','tile');
    grid off
    xlabel('Peak Velocity (°/s)')
    ylabel('Duration (s)')
    zlabel('Probability')
    
ax(6) = subplot(3,2,6);
    X = abs(Static.SACCADE_STATS.PeakVel);
    Y = 1000*abs(Static.SACCADE_STATS.Duration);
    X = X(~isnan(X));
    Y = Y(~isnan(Y));
    [rho_static.pkv_dur,pval_static.pkv_dur] = corr(X,Y);
    h(6) = histogram2(X, Y, edges.pkv, edges.dur, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'ShowEmptyBins','on', 'Normalization','Probability','DisplayStyle','tile');
    grid off
    xlabel('Peak Velocity (°/s)')
    ylabel('Duration (s)')
    zlabel('Probability')
    
set(ax, 'LineWidth', 1, 'FontWeight', 'bold')
% set(h, 'ShowEmptyBins', 'off')
% set(h,'EdgeColor','none')

hp = get(subplot(3,2,6),'Position');
cbar = colorbar('Position', [hp(1)+hp(3)+0.02  0.15+hp(2)  0.03  hp(2)+hp(3)*1]);
cbar.Label.String = 'Probability';

%% Surface
x = 1000*abs(Move.SACCADE_STATS.Duration);
y = abs(Move.SACCADE_STATS.Amplitude);
z = abs(Move.SACCADE_STATS.PeakVel);

nanI = isnan(x);
x = x(~nanI);
y = y(~nanI);
z = z(~nanI);

n = 200;
edges.dur = linspace(0,100,n);
edges.amp = linspace(7,31,n);
edges.pkv = linspace(300,1000,n);

% vq = griddata(x, y, z, edges.dur, edges.amp, edges.pkv);
clc
[xq,yq] = meshgrid(edges.dur, edges.amp);
vq = griddata(x, y, z, edges.dur, edges.amp');

surface(xq,yq,vq, 'EdgeColor', 'none')

hold on
view(3)
xlabel('Duration (s)')
ylabel('Ampltitude (°)')
zlabel('Peak Velocity (°/s)')

end