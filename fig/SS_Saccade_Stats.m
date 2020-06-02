function [] = SS_Saccade_Stats()
%% SS_Saccade_Stats:

root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','D','I','U','N')

clms = N.freq;
CC = hsv(clms);

%% Saccade Statistics %%
absflag = true;

Freq = U.freq{1};
Freq = sort(unique(abs(Freq)),'ascend');
n_freq = length(Freq);
freq = 1:n_freq;

saccade_freq = SACCADE_STATS.freq;
if ~absflag
    freq = [freq ; -freq];
else
    saccade_freq = abs(saccade_freq);
end

n_sacd = size(SACCADE_STATS,1);
G = nan(n_sacd,1);
for jj = 1:n_freq
    if ~isnan(freq(jj))
        [rr,~] = find( saccade_freq == freq(jj) );
    else
       [rr,~] = find(isnan(saccade_freq));
    end
    G(rr,1) = jj;
end
w_scale = 2*n_freq/10;

SS = [5,6,18,14,15,20,21];
ylim_list = [100 30 1050];

n_plot = length(SS);
% n_plot = 3;
clms = 3;
rows = ceil(n_plot/clms);

YY = [  "Duration (ms)",...
        "Amplitude (°)",...
        "Peak Velocity (°/s)",...
      	"Trigger Position (°)",...
        "End Position (°)",...
        "Interval Duration",...
        "Error (°)"];

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 2*clms*w_scale 2*rows];
movegui(FIG,'center')
ax = gobjects(n_plot,1);
for ww = 1:n_plot
    ax(ww) = subplot(rows,clms,ww); axis tight
    data = SACCADE_STATS{:,SS(ww)};
    
    if any(SS(ww)==[6,7,18,14,15,20,21])
        data = abs(data);
    end
    if ww==1
        data = 1000*data;
    end
    
    bx = boxplot(data,G,'Labels',{Freq},'Width',0.5,'Symbol','.','Whisker',2);
    xlabel('Stimulus Frequency (Hz)')
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
    try
        ax(ww).YLim(2) = ylim_list(ww);
    end
end

set(ax,'LineWidth',1,'FontWeight','bold')

%% Saccade Count
count.stats = cellfun(@(x) basic_stats(x,1), COUNT, 'UniformOutput', true);
for f  = 1:N.freq
    count.all{f,1} = cat(1,COUNT{:,f});
    count.med(:,f)  = cat(1,count.stats(:,f).median);
    count.mean(:,f) = cat(1,count.stats(:,f).mean);
end
n_length = cellfun(@(x) length(x), count.all, 'UniformOutput', true);
G = [];
for f = 1:N.freq
   G = [G ; f*ones(n_length(f),1)]; 
end
count_all = cat(1,count.all{:});


FIG = figure (3) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 2 2];
movegui(FIG,'center')
ax = subplot(1,1,1); hold on

bx = boxplot(count_all./10, G, 'Labels', {Freq}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
xlabel('Stimulus Frequency (Hz)')
ylabel('Rate (#/s')

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1},h{kk,2},CC(kk,:));
end

set(findobj(ax,'tag','Median'), 'Color', 'w','LineWidth',1.5);
set(findobj(ax,'tag','Box'), 'Color', 'none');
set(findobj(ax,'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax,'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax.Children = ax.Children([end 1:end-1]);

set(ax,'LineWidth',1,'FontWeight','bold')


end