function [] = Wing_Trigger_left_right()
%% Wing_Trigger_left_right:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'U','N','SACCADE')

%% Get head and wing intervals
clearvars -except SACCADE U N

keepI = cellfun(@(x) isstruct(x) | isobject(x), SACCADE.head2wing);
Saccade = SACCADE(keepI,:);
n_file = size(Saccade,1);
n_speed = N.vel/2;
Vel = U.vel{1};

Fs = round(Saccade.wing_saccade{1}.Fs);
Fc = 7;
[b,a] = butter(3, Fc/(Fs/2), 'low');

[vel_idx,vel] = findgroups(Saccade.vel);
n_vel = length(vel);

% vel_idx(vel_idx >= n_vel) = vel_idx(vel_idx >= n_vel) - n_vel + 1;
% n_vel = 1;

[fly_idx,fly] = findgroups(Saccade.fly);
n_fly = length(fly);

LWing = cell(n_fly, n_vel);
RWing = cell(n_fly, n_vel);
WBA_sacd = cell(n_fly, n_vel);

Left = cell(n_fly, n_vel);
Right = cell(n_fly, n_vel);
WBA = cell(n_fly, n_vel);
for n = 1:n_file
  	vel = vel_idx(n);
    fly = fly_idx(n);
        
    head_saccade = Saccade.head_saccade{n};
    lwing = Saccade.wing_saccade{n}.extra.lwing;
    rwing = Saccade.wing_saccade{n}.extra.rwing;
    dwba = Saccade.wing_saccade{n}.extra.dwba;
    
    lwing_filt = filtfilt(b, a, lwing);
    rwing_filt = filtfilt(b, a, rwing);
    dwba_filt = filtfilt(b, a, dwba);
    
    Left{fly,vel}(:,end+1) = lwing_filt;
    Right{fly,vel}(:,end+1) = rwing_filt;
    WBA{fly,vel}(:,end+1) = dwba_filt;
    
    [lwing_sacd,~] = getSaccade(head_saccade, lwing_filt);
    [rwing_sacd,~] = getSaccade(head_saccade, rwing_filt);
    [dwba_sacd,~] = getSaccade(head_saccade, dwba_filt);
    
    n_sacd = length(lwing_sacd);
    sz = size(LWing{fly,vel},1);
    span = sz+1:sz+n_sacd;
    
    lwing_start_end = cellfun(@(x) [x(1), x(end)], lwing_sacd, 'UniformOutput', false);
    rwing_start_end = cellfun(@(x) [x(1), x(end)], rwing_sacd, 'UniformOutput', false);
    dwba_start_end = cellfun(@(x) [x(1), x(end)], dwba_sacd, 'UniformOutput', false);
    
    LWing{fly,vel}(span,1:2) = cat(1, lwing_start_end{:});
    RWing{fly,vel}(span,1:2) = cat(1, rwing_start_end{:});
    WBA_sacd{fly,vel}(span,1:2) = cat(1, dwba_start_end{:});
end

Left_all = cell(1, n_vel);
Right_all = cell(1, n_vel);
WBA_all = cell(1, n_vel);

LWing_all = cell(1, n_vel);
RWing_all = cell(1, n_vel);
WBA_sacd_all = cell(1, n_vel);
for v = 1:n_vel
    Left_all{v} = cat(2, Left{:,v});
    Right_all{v} = cat(2, Right{:,v});
    WBA_all{v} = cat(2, WBA{:,v});
    
    LWing_all{v} = cat(1, LWing{:,v});
    RWing_all{v} = cat(1, RWing{:,v});
    WBA_sacd_all{v} = cat(1, WBA_sacd{:,v});
end

%% Near vs Far side wings plot
FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 5 7];
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h

ax(1) = subplot(2,1,1); cla ; hold on ; grid off ; axis tight ; title('Near vs Far Wing')
    edges = -10:2:90;
    Near = [Right_all{1},Left_all{2}];
    Far = [Right_all{2},Left_all{1}];
    h(1) = histogram(Near, edges, 'Normalization', 'Probability');
    h(2) = histogram(Far, edges, 'Normalization', 'Probability');
    legend(h, 'Near', 'Far', 'Box', 'off')
    ax(1).YLim(1) = -0.001;

ax(2) = subplot(2,1,2); cla ; hold on ; grid off ; axis tight ; title('\DeltaWBA')
    edges = -70:2:70;
    All = [WBA_all{1},WBA_all{2}];
    h(3) = histogram(All, edges, 'Normalization', 'Probability');
    legend(h, 'Near', 'Far', 'Box', 'off')
    ax(1).YLim(1) = -0.001;
    
% set(h, 'EdgeColor', 'k')
set(ax, 'LineWidth', 1)

%% Saccade Left-Right Trigger Polar Plot, normalized to stimulus direction
FIG = figure (2) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 5 2];
FIG.Name = 'Normalized Saccade Wing Trigger';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h

edges = deg2rad(0:5:90);

ax(1) = subplot(1,2,1,polaraxes); grid off ; axis tight
    All_Near = [ Right_all{1}(:) ; Left_all{2}(:) ];
    Start_Near = [ RWing_all{1}(:,1) ; LWing_all{2}(:,1) ];
    End_Near = [ RWing_all{1}(:,2) ; LWing_all{2}(:,2) ];
    
%     h(1,1) = polarhistogram(deg2rad(All_Near), 'BinEdges', edges, ...
%         'FaceColor','k','FaceAlpha',0.5,'Normalization','Probability'); hold on
    h(2,1) = polarhistogram(deg2rad(End_Near), 'BinEdges', edges, ...
        'FaceColor','r','FaceAlpha',0.5,'Normalization','Probability'); hold on
    h(3,1) = polarhistogram(deg2rad(Start_Near), 'BinEdges', edges, ...
        'FaceColor','g','FaceAlpha',0.5, 'Normalization','Probability');
    ax(1).ThetaAxis.Label.String = 'Left Wing (°)';
    
ax(2) = subplot(1,2,2,polaraxes); grid off ; axis tight
    All_Far = [ Right_all{2}(:) ; Left_all{1}(:) ];
    Start_Far = [ RWing_all{2}(:,1) ; LWing_all{1}(:,1) ];
    End_Far = [ RWing_all{2}(:,2) ; LWing_all{1}(:,2) ];
    
%     h(1,2) = polarhistogram(deg2rad(All_Far), 'BinEdges', edges, ...
%         'FaceColor','k','FaceAlpha',0.5,'Normalization','Probability'); hold on
    h(2,2) = polarhistogram(deg2rad(End_Far), 'BinEdges', edges, ...
        'FaceColor','r','FaceAlpha',0.5, 'Normalization','Probability'); hold on
    h(3,2) = polarhistogram(deg2rad(Start_Far), 'BinEdges', edges, ...
        'FaceColor','g','FaceAlpha',0.5,'Normalization','Probability');
    ax(2).ThetaAxis.Label.String = 'Right Wing (°)';
    
set(h(2:3,:),'EdgeColor','none')
set(ax,'FontSize',8)
set(ax,'Color','w')
set(ax,'RLim',[0 0.21])
set(ax,'ThetaLim',[0 90])
set(ax,'ThetaTick',0:10:90);

set(ax(1),'ThetaZeroLocation','left')
set(ax(2),'ThetaZeroLocation','right')

set(ax(1),'ThetaDir','clockwise')
set(ax(2),'ThetaDir','counterclockwise')

%% Saccade WBA Trigger Polar Plot, normalized to stimulus direction
FIG = figure (3) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 2.5 2];
FIG.Name = 'Normalized Saccade Wing Trigger';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h

edges = deg2rad(-60:5:60);

ax(1) = subplot(1,1,1,polaraxes); grid off ; axis tight
    All = [ -WBA_all{1}(:) ; WBA_all{2}(:) ];
    Start = [ -WBA_sacd_all{1}(:,1) ; WBA_sacd_all{2}(:,1) ];
    End = [ -WBA_sacd_all{1}(:,2) ; WBA_sacd_all{2}(:,2) ];
    
%     h(1,1) = polarhistogram(deg2rad(All), 'BinEdges', edges, ...
%         'FaceColor','k','FaceAlpha',0.5,'Normalization','Probability'); hold on
    h(2,1) = polarhistogram(deg2rad(End), 'BinEdges', edges, ...
        'FaceColor','r','FaceAlpha',0.5,'Normalization','Probability'); hold on
    h(3,1) = polarhistogram(deg2rad(Start), 'BinEdges', edges, ...
        'FaceColor','g','FaceAlpha',0.5, 'Normalization','Probability');
    
    ax(1).ThetaAxis.Label.String = '\DeltaWBA (°)';
    
set(h(2:3),'EdgeColor','none')
set(ax,'FontSize',8)
set(ax,'Color','w')
% set(ax,'RLim',[0 0.21])
set(ax,'ThetaLim',60*[-1 1])
set(ax,'ThetaTick',-70:10:70);
set(ax(1),'ThetaZeroLocation','top')
set(ax(1),'ThetaDir','clockwise')

%% ANOVA
D = [All_Near ; Start_Near];
G = [ones(length(All_Near),1) ; 2*ones(length(Start_Near),1)];
p = anovan(D,{G});

end