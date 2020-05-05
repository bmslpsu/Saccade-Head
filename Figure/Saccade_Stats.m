function [] = Saccade_Stats()
%% Saccade_Stats:

root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'COUNT','SACCADE_STATS','U','N')

clearvars -except clms CC Vel PATH COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

clms = N.vel/2;
CC = repmat(hsv(clms),2,1);
Vel = U.vel{1};
Speed = Vel(1:N.vel/2);

%% Saccade Statistics %%
absflag = true;
vel = U.vel{1};
vel = sort(unique(abs(vel)),'ascend');
saccade_vel = SACCADE_STATS.Vel;
if ~absflag
    vel = [vel ; -vel];
else
    saccade_vel = abs(saccade_vel);
end
n_vel = length(vel);

n_sacd = size(SACCADE_STATS,1);
G = nan(n_sacd,1);
for jj = 1:n_vel
    [rr,~] = find( saccade_vel == vel(jj) );
    G(rr,1) = jj;
end
w_scale = 2*n_vel/10;

SS = [6,7,19,15,17,21,22,24,23,25];
ylim_list = [120 35 1200];
% n_plot = length(SS);
n_plot = 3;
clms = 4;
rows = ceil(n_plot/clms);

YY = [  "Duration (ms)",...
        "Amplitude (°)",...
        "Peak Velocity (°/s)",...
      	"Trigger Position (°)",...
        "End Position (°)",...
        "Interval Duration", ...
        "Error (°)",...
        "Integrated Error (°*s)",...
      	"Velocity Err (°/s)",...
        "Integrated Velocity Err (°)" ];

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 2*clms*w_scale 2*rows];
movegui(FIG,'center')
ax = gobjects(n_plot,1);
for ww = 1:n_plot
    ax(ww) = subplot(rows,clms,ww); axis tight
    data = SACCADE_STATS{:,SS(ww)};
    
    if any(SS(ww)==[6,7,19,15,17,21,22,24,23,25])
        data = abs(data);
    end
    if ww==1
        data = 1000*data;
    end
    
    bx = boxplot(data, G, 'Labels', {vel}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
    xlabel('Stimulus Speed (°/s)')
    ylabel(YY(ww))

    h = get(bx(5,:),{'XData','YData'});
    for kk = 1:size(h,1)
       patch(h{kk,1},h{kk,2},CC(kk,:));
    end

    set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth',1.5);
    set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
    set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
    set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
    ax(ww).Children = ax(ww).Children([end 1:end-1]);
    ax(ww).YLim(1) = 0;
    ax(ww).YLim(2) = ylim_list(ww);
end

set(ax,'LineWidth',1,'FontWeight','bold', 'Box', 'off')

%% Saccade Count/Rate
count.stats = cellfun(@(x) basic_stats(x,1), COUNT, 'UniformOutput', true);
for v = 1:N.vel
    count.all{v,1} = cat(1,COUNT{:,v});
    count.med(:,v)  = cat(1,count.stats(:,v).median);
    count.mean(:,v) = cat(1,count.stats(:,v).mean);
end
n_length = cellfun(@(x) length(x), count.all, 'UniformOutput', true);
n_length = sum(reshape(n_length,N.vel/2,2),2);
G = [];
for v = 1:N.vel/2
   G = [G ; v*ones(n_length(v),1)]; 
end
count_all = cat(1,count.all{:});


% FIG = figure (3) ; clf
% FIG.Color = 'w';
% FIG.Units = 'inches';
% FIG.Position = [2 2 2 2];
% movegui(FIG,'center')
ax(4) = subplot(1,4,4); hold on

bx = boxplot(count_all./10, G, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
xlabel('Stimulus Speed (°/s)')
ylabel('Rate (#/s)')
ylim([-0.1 3])

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1},h{kk,2},CC(kk,:));
end

set(findobj(ax(4),'tag','Median'), 'Color', 'w','LineWidth',1.5);
set(findobj(ax(4),'tag','Box'), 'Color', 'none');
set(findobj(ax(4),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(4),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(4).Children = ax(4).Children([end 1:end-1]);

set(ax,'LineWidth',1,'FontWeight','bold', 'Box', 'off')

end