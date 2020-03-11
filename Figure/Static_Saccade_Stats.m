function [] = Static_Saccade_Stats()
%% Static_Saccade_Stats:

root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','D','I','U','N')

clms = N.wave;
CC = repmat(jet(clms),2,1);

%% Saccade Statistics %%
absflag = true;

wave = U.wave{1};
wave = sort(unique(abs(wave)),'ascend');
saccade_wave = SACCADE_STATS.Wave;
if ~absflag
    wave = [wave ; -wave];
else
    saccade_wave = abs(saccade_wave);
end
n_wave = length(wave);

n_sacd = size(SACCADE_STATS,1);
G = nan(n_sacd,1);
for jj = 1:n_wave
    if ~isnan(wave(jj))
        [rr,~] = find( saccade_wave == wave(jj) );
    else
       [rr,~] = find(isnan(saccade_wave));
    end
    G(rr,1) = jj;
end
w_scale = 2*n_wave/10;

SS = [6,7,19,15,17,21];
ylim_list = [100 30 1050];

% n_plot = length(SS);
n_plot = 3;
clms = 3;
rows = ceil(n_plot/clms);

YY = [  "Duration (ms)",...
        "Amplitude (°)",...
        "Peak Velocity (°/s)",...
      	"Trigger Position (°)",...
        "End Position (°)",...
        "Interval Duration" ];

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
    
    bx = boxplot(data,G,'Labels',{wave},'Width',0.5,'Symbol','.','Whisker',2);
    xlabel('Stimulus Speed (°/s)')
    ylabel(YY(ww))
    box off

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

set(ax,'LineWidth',1,'FontWeight','bold')
end