function [] = Saccade_Norm_Rate()
%% Saccade_Norm_Rate:

% Main saccade stats
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat',}, 'Select data', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'COUNT','SACCADE_STATS','U','N')

clearvars -except clms CC Vel PATH COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

clms = N.vel/2;
CC = repmat(hsv(clms),2,1);
Vel = U.vel{1};
Speed = Vel(1:N.vel/2);

% Gain stats
root = 'C:\Users\BC\Box\Research\Manuscripts\Head Saccade\Data';

[FILE,PATH] = uigetfile({'*.mat'},'Select data files', root, 'MultiSelect','off');
FILE = cellstr(FILE);

n_wave = length(FILE);
wave_data = cell(n_wave,1);
for n = 1:n_wave
    wave_data{n} = load(fullfile(PATH,FILE{n}));
end
wave_color = parula(n_wave);


%% Normalized Saccade Count/Rate
range = 15;
gain_mean = [wave_data{1}.Mean_Gain_Stats.mean]';
optm_rate = range ./ ( abs(Vel) .* repmat(gain_mean,2,1) );

count.stats = cellfun(@(x) basic_stats(x,1), COUNT, 'UniformOutput', true);
for v = 1:N.vel
    count.all{v,1} = cat(1,COUNT{:,v}) / optm_rate(v);
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

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 2 2];
movegui(FIG,'center')
clear ax h
ax(1) = subplot(1,1,1); hold on

bx = boxplot(count_all./10, G, 'Labels', {Speed}, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
xlabel('Stimulus Speed (°/s)')
ylabel('Rate (#/s)')
ylim([-0.1 5])

h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1},h{kk,2},CC(kk,:));
end

set(findobj(ax(1),'tag','Median'), 'Color', 'w','LineWidth',1.5);
set(findobj(ax(1),'tag','Box'), 'Color', 'none');
set(findobj(ax(1),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(1),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(1).Children = ax(1).Children([end 1:end-1]);

set(ax, 'LineWidth', 1, 'Box', 'off')

end