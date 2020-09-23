function [] = Static_saccade_pos_vel()
%% Saccade_Pos_Vel:
root = 'H:\DATA\Rigid_Data\Saccade';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'HEAD_SACCADE','FLY','GRAND','D','I','U','N')

%% Get saccade at given wavelengths & normalize direction
% wave_group = [1 2 3 4 5 6];
% wave_group = [2 3 4];
wave_group = 3;

Sacd_All = GRAND.normpeak_saccade(wave_group);
norm_fields = {'time','position','velocity'};
Sacd_All = struct_center(Sacd_All, 0, true , 1, norm_fields);

% Normalize saccades
[~,max_vel] = max( abs(Sacd_All.velocity), [], 1);
dir = sign( Sacd_All.velocity(max_vel(1),:) );
dir = repmat(dir, size(Sacd_All.time,1), 1);
Sacd_All.position = dir .* Sacd_All.position;
Sacd_All.velocity = dir .* Sacd_All.velocity;

CC = [0 0 0]; % color

%% Saccade Position by Speed
FIG = figure (102) ; clf
FIG.Units = 'inches';
FIG.Position = 1*[2 2 1*(4/3) 1.5];
FIG.Name = 'Saccade Position';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h
med_span = 5;

ax(1) = subplot(1,1,1) ; hold on
title([num2str(0) ' (°/s)']);

time = 1000 * [Sacd_All.time];
med_time = nanmean(time,2);
pos = Sacd_All.position;
med_pos = nanmean(pos,2);
std_pos = nanstd(pos,[],2);

cent = find(med_time==0);
span = (cent - med_span):( cent + med_span);

h.trial = plot(time, pos, 'Color', [0.5 0.5 0.5 , 0.2], 'LineWidth', 0.25);

h.patch = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
    CC,  0.5*CC, 0.3, 1);
    
set(ax,'LineWidth',1,'FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.07*[-1 1],'YLim',25*[-1 1])
xlabel('Time (ms)')
ylabel('Head Position (°)')

%% Saccade Velocity by Speed
FIG = figure (103) ; clf
FIG.Units = 'inches';
FIG.Position = 1*[2 2 1*(4/3) 1.5];
FIG.Name = 'Saccade Velocity';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h
med_span = 5;

ax(1) = subplot(1,1,1) ; hold on
title([num2str(0) ' (°/s)']);

time = 1000 * [Sacd_All.time];
med_time = nanmean(time,2);
pos = Sacd_All.velocity;
med_pos = nanmean(pos,2);
std_pos = nanstd(pos,[],2);

cent = find(med_time==0);
span = (cent - med_span):( cent + med_span);

h.trial = plot(time, pos, 'Color', [0.5 0.5 0.5 , 0.2], 'LineWidth', 0.25);

h.patch = PlotPatch(med_pos(span), std_pos(span), med_time(span), 1, 1, ...
    CC,  0.5*CC, 0.3, 1);
    
set(ax,'LineWidth',1,'FontSize',8,'Color','w',...
    'YColor','k','XColor','k','XLim',1000*0.07*[-1 1],'YLim',[-200 1100])
xlabel('Time (ms)')
ylabel('Head Velocity (°/s)')

end