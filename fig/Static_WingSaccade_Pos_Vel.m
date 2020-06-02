function [] = Static_WingSaccade_Pos_Vel()
%% Static_WingSaccade_Pos_Vel:
root = 'H:\DATA\Rigid_Data\';

[CW,PATH_CW] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select CW data', root, 'MultiSelect','off');

[CCW,PATH_CCW] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select CCW data', root, 'MultiSelect','off');

CW = load(fullfile(PATH_CW,CW),'PATH','HEAD_COUNT','WING_COUNT','SACCADE','HEAD_STATS',...
      'WING_STATS','HEAD_FLY','HEAD_GRAND','WING_FLY','WING_GRAND','D','I','U','N','T','direction');
CCW = load(fullfile(PATH_CCW,CCW),'PATH','HEAD_COUNT','WING_COUNT','SACCADE','HEAD_STATS',...
      'WING_STATS','HEAD_FLY','HEAD_GRAND','WING_FLY','WING_GRAND','D','I','U','N','T','direction');

clearvars -except CW CCW

%% Saccade Position %%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = 2*[2 2 2 2];
FIG.Name = 'Static Saccade Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h

med_span = [80 80];
cw_color = [1 0 0];
ccw_color = [0 0 1];

ax = subplot(1,1,1); hold on

% CW
time = 1000*CW.WING_GRAND.normpeak_saccade.time;
med_time = nanmedian(time,2);
pos = CW.WING_GRAND.normpeak_saccade.position;
for n = 1:size(pos,2)
    sIdx = find(time(:,n) == 0, 1, 'first');
    pos(:,n) =  pos(:,n) - pos(sIdx,n);
end
med_pos = nanmedian(pos,2);
std_vel = nanstd(pos,[],2);

cent = find(med_time==0);
span = (cent - med_span(1)):( cent + med_span(2));

h.trial = plot(time, pos, 'Color', [0.7*cw_color , 0.2], 'LineWidth', 0.5);

[h.patch(1), h.mean(1)] = PlotPatch(med_pos(span), std_vel(span), med_time(span), 1, 1, ...
    cw_color, [0.7 0.7 0.7], 0.4, 3);
h.patch(1).EdgeColor = cw_color;

% CCW
time = 1000*CCW.WING_GRAND.normpeak_saccade.time;
med_time = nanmedian(time,2);
pos = CCW.WING_GRAND.normpeak_saccade.position;
for n = 1:size(pos,2)
    sIdx = find(time(:,n) == 0, 1, 'first');
    pos(:,n) =  pos(:,n) - pos(sIdx,n);
end
med_pos = nanmedian(pos,2);
std_vel = nanstd(pos,[],2);

cent = find(med_time==0);
span = (cent - med_span(1)):( cent + med_span(2));

h.trial = plot(time, pos, 'Color', [0.7*ccw_color , 0.2], 'LineWidth', 0.5);

[h.patch(2), h.mean(2)] = PlotPatch(med_pos(span), std_vel(span), med_time(span), 1, 1, ...
    ccw_color, [0.7 0.7 0.7], 0.4, 3);
h.patch(2).EdgeColor = ccw_color;

uistack(h.patch,'top')
uistack(h.mean,'top')

plot(1000*[-1 1], [0 0], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2)

set(ax,'LineWidth',1.5,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.4*[-1 1],'YLim',40*[-1 1])
xlabel('Time (ms)')
ylabel('\DeltaWBA (°)')

legend(h.mean, 'CW', 'CCW', 'Box', 'off')

%% Saccade Velocity %%
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = 2*[2 2 2 2];
FIG.Name = 'Static Saccade Velocity';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h

med_span = 5;
cw_color = [1 0 0];
ccw_color = [0 0 1];

ax = subplot(1,1,1); hold on

% CW
time = 1000*CW.WING_GRAND.normpeak_saccade.time;
med_time = nanmedian(time,2);
vel = CW.WING_GRAND.normpeak_saccade.velocity;
med_vel = nanmedian(vel,2);
std_vel = nanstd(vel,[],2);

cent = find(med_time==0);
span = (cent - med_span):( cent + med_span);

h.trial = plot(time, vel, 'Color', [0.7*cw_color , 0.2], 'LineWidth', 0.5);

[h.patch, h.mean(1)] = PlotPatch(med_vel(span), std_vel(span), med_time(span), 1, 1, ...
    cw_color, [0.7 0.7 0.7], 0.4, 3);
h.patch.EdgeColor = cw_color;

% CCW
time = 1000*CCW.WING_GRAND.normpeak_saccade.time;
med_time = nanmedian(time,2);
vel = CCW.WING_GRAND.normpeak_saccade.velocity;
med_vel = nanmedian(vel,2);
std_vel = nanstd(vel,[],2);

cent = find(med_time==0);
span = (cent - med_span):( cent + med_span);

h.trial = plot(time, vel, 'Color', [0.7*ccw_color , 0.2], 'LineWidth', 0.5);

[h.patch, h.mean(2)] = PlotPatch(med_vel(span), std_vel(span), med_time(span), 1, 1, ...
    ccw_color, [0.7 0.7 0.7], 0.4, 3);
h.patch.EdgeColor = ccw_color;

set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.05*[-1 1],'YLim',1500*[-1 1])
xlabel('Time (ms)')
ylabel('Wing Velocity (°/s)')

legend(h.mean, 'CW', 'CCW', 'Box', 'off')

end