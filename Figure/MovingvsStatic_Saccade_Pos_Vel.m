function [] = MovingvsStatic_Saccade_Pos_Vel()
%% Static_Saccade_Pos_Vel:
root = 'H:\DATA\Rigid_Data\';

[CW,PATH_CW] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select CW data', root, 'MultiSelect','off');

[CCW,PATH_CCW] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select CCW data', root, 'MultiSelect','off');

[MOVE,PATH_Move] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select Moving data', root, 'MultiSelect','off');

CW = load(fullfile(PATH_CW,CW),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','D','I','U','N');
CCW = load(fullfile(PATH_CCW,CCW),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','D','I','U','N');
MOVE = load(fullfile(PATH_Move,MOVE),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','D','I','U','N');

clearvars -except CW CCW MOVE

%% Get all saccades for Moving, CW, CCW data
center = 0;
even = true;
dim = 1;

ALL = [];

ALL.MovingCCW = struct_center(MOVE.GRAND.normpeak_saccade(1:5), center, even, dim);
ALL.MovingCCW_time = basic_stats(ALL.MovingCCW.time,2);
ALL.MovingCCW_vel = basic_stats(ALL.MovingCCW.velocity,2);
ALL.MovingCCW_pos = basic_stats(ALL.MovingCCW.position,2);

ALL.MovingCW = struct_center(MOVE.GRAND.normpeak_saccade(6:10), center, even, dim);
ALL.MovingCW_time = basic_stats(ALL.MovingCW.time,2);
ALL.MovingCW_vel = basic_stats(ALL.MovingCW.velocity,2);
ALL.MovingCW_pos = basic_stats(ALL.MovingCW.position,2);

ALL.CW = struct_center(CW.GRAND.normpeak_saccade, center, even, dim);
ALL.CW_time = basic_stats(ALL.CW.time,2);
ALL.CW_vel = basic_stats(ALL.CW.velocity,2);
ALL.CW_pos = basic_stats(ALL.CW.position,2);

ALL.CCW = struct_center(CCW.GRAND.normpeak_saccade, center, even, dim);
ALL.CCW_time = basic_stats(ALL.CCW.time,2);
ALL.CCW_vel = basic_stats(ALL.CCW.velocity,2);
ALL.CCW_pos = basic_stats(ALL.CCW.position,2);


%% Saccade Position %%
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = 2*[2 2 1.1384*(4/3) 3/2];
FIG.Name = 'Static Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';

move_color = [0.8 0 0];
static_color = [0 0 0.8];
med_span = 5;

ax = subplot(1,1,1); hold on
    % Moving CCW
    cent = find(ALL.MovingCCW_time.median==0);
    span = (cent - med_span):( cent + med_span);

    % plot(1000*ALL.MovingCCW.time, ALL.MovingCCW.position, 'Color', [0.5 0.5 0.5 0.2])
    [h.patch,h.move] = PlotPatch(ALL.MovingCCW_pos.median(span), ALL.MovingCCW_pos.std(span), ...
                1000*ALL.MovingCCW_time.median(span), 1, 1, move_color, move_color, 0.2, 3);
    h.patch.EdgeColor = move_color;

    % Moving CW
    cent = find(ALL.MovingCW_time.median==0);
    span = (cent - med_span):( cent + med_span);

    % plot(1000*ALL.MovingCW.time, ALL.MovingCW.position, 'Color', [0.5 0.5 0.5 0.2])
    [h.patch,h.move] = PlotPatch(ALL.MovingCW_pos.median(span), ALL.MovingCW_pos.std(span), ...
                1000*ALL.MovingCW_time.median(span), 1, 1, move_color, move_color, 0.2, 3);
    h.patch.EdgeColor = move_color;
    
    % Static CCW
    cent = find(ALL.CCW_time.median==0);
    span = (cent - med_span):( cent + med_span);

    % plot(1000*ALL.CCW.time, ALL.CCW.position, 'Color', [0.5 0.5 0.5 0.2])
    [h.patch,h.stat] = PlotPatch(ALL.CCW_pos.median(span), ALL.CCW_pos.std(span), ...
                1000*ALL.CCW_time.median(span), 1, 1, static_color, static_color, 0.2, 3);
    h.patch.EdgeColor = move_color;

    % Static CW
    cent = find(ALL.CW_time.median==0);
    span = (cent - med_span):( cent + med_span);

    % plot(1000*ALL.CW.time, ALL.CW.position, 'Color', [0.5 0.5 0.5 0.2])
    [h.patch,h.stat] = PlotPatch(ALL.CW_pos.median(span), ALL.CW_pos.std(span), ...
                1000*ALL.CW_time.median(span), 1, 1, static_color, static_color, 0.2, 3);
    h.patch.EdgeColor = move_color;
    
xlabel('Time (ms)')
ylabel('Head Position (°)')

leg = legend([h.move,h.stat],'Moving','Static');
leg.Box = 'off';

set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.03*[-1 1],'YLim',20*[-1 1])

%% Saccade Velocity %%
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = 2*[2 2 1.1384*(4/3) 3/2];
FIG.Name = 'Static Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';

move_color = [0.8 0 0];
static_color = [0 0 0.8];
med_span = 5;

ax = subplot(1,1,1); hold on
    % Moving CCW
    cent = find(ALL.MovingCCW_time.median==0);
    span = (cent - med_span):( cent + med_span);

    % plot(1000*ALL.MovingCCW.time, ALL.MovingCCW.velocity, 'Color', [0.5 0.5 0.5 0.2])
    [h.patch,h.move] = PlotPatch(ALL.MovingCCW_vel.median(span), ALL.MovingCCW_vel.std(span), ...
                1000*ALL.MovingCCW_time.median(span), 1, 1, move_color, move_color, 0.2, 3);
    h.patch.EdgeColor = move_color;

    % Moving CW
    cent = find(ALL.MovingCW_time.median==0);
    span = (cent - med_span):( cent + med_span);

    % plot(1000*ALL.MovingCW.time, ALL.MovingCW.velocity, 'Color', [0.5 0.5 0.5 0.2])
    [h.patch,h.move] = PlotPatch(ALL.MovingCW_vel.median(span), ALL.MovingCW_vel.std(span), ...
                1000*ALL.MovingCW_time.median(span), 1, 1, move_color, move_color, 0.2, 3);
    h.patch.EdgeColor = move_color;
    
    % Static CCW
    cent = find(ALL.CCW_time.median==0);
    span = (cent - med_span):( cent + med_span);

    % plot(1000*ALL.CCW.time, ALL.CCW.velocity, 'Color', [0.5 0.5 0.5 0.2])
    [h.patch,h.stat] = PlotPatch(ALL.CCW_vel.median(span), ALL.CCW_vel.std(span), ...
                1000*ALL.CCW_time.median(span), 1, 1, static_color, static_color, 0.2, 3);
    h.patch.EdgeColor = move_color;

    % Static CW
    cent = find(ALL.CW_time.median==0);
    span = (cent - med_span):( cent + med_span);

    % plot(1000*ALL.CW.time, ALL.CW.velocity, 'Color', [0.5 0.5 0.5 0.2])
    [h.patch,h.stat] = PlotPatch(ALL.CW_vel.median(span), ALL.CW_vel.std(span), ...
                1000*ALL.CW_time.median(span), 1, 1, static_color, static_color, 0.2, 3);
    h.patch.EdgeColor = move_color;
    
xlabel('Time (ms)')
ylabel('Head Velocity (°)')

leg = legend([h.move,h.stat],'Moving','Static');
leg.Box = 'off';

set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.05*[-1 1],'YLim',1000*[-1 1])

end