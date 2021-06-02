function [] = All_stats_wing_new()
%% All_stats_wing_new:
root = 'E:\DATA\Rigid_Data\Saccade';

[Free,FreePath] = uigetfile({'*.mat'},'Select free data', root, 'MultiSelect','off');
[Fixed30,FixedPath30] = uigetfile({'*.mat'},'Select fixed 30 data', root, 'MultiSelect','off');
[Fixed60,FixedPath60] = uigetfile({'*.mat'},'Select fixed 60 data', root, 'MultiSelect','off');
[Static,StaticPath] = uigetfile({'*.mat'},'Select static data', root, 'MultiSelect','off');

Free = load(fullfile(FreePath,Free),'SACCADE','U','N','D');
Fixed30 = load(fullfile(FixedPath30,Fixed30),'SACCADE','U','N','D');
Fixed60 = load(fullfile(FixedPath60,Fixed60),'SACCADE','U','N','D');
Static = load(fullfile(StaticPath,Static),'SACCADE','U','N','D');

%% Saccade Frequency
clearvars -except Free Static Fixed30 Fixed60 root
clc

% Get saccade frequencies
all.head_free_30    = get_count(Free.SACCADE, [], [1 3], 1);
all.head_free_60    = get_count(Free.SACCADE, [], [2 4], 1);
all.wing_free_30    = get_count(Free.SACCADE, [], [1 3], 2);
all.wing_free_60    = get_count(Free.SACCADE, [], [2 4], 2);
all.wing_fixed_30   = get_count(Fixed30.SACCADE, [], [], 2);
all.wing_fixed_60   = get_count(Fixed60.SACCADE, [], [], 2);
all.head_static     = get_count(Static.SACCADE, [2 3 4], [], 1);

all_count = struct2cell( structfun(@(x) x.count(:,1), all, 'UniformOutput', false) );
% all_count = struct2cell( structfun(@(x) x.count_vel_fly_mean(:,1), all, 'UniformOutput', false) );
all_group = cellfun(@(x,y) y*ones(size(x)), all_count, num2cell((1:length(all_count))'),  ...
    'UniformOutput', false);

all_count = cat(1, all_count{:});
all_group = cat(1, all_group{:});

labels = fieldnames(all);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 2 2])
movegui(fig, 'center')
clear ax
ax(1) = subplot(1,1,1); hold on ; ylim([-0.1 2])
    b = boxchart(all_count./10, 'GroupByColor', all_group, ...
        'LineWidth', 0.5, 'MarkerStyle', '.', 'JitterOutliers', 'on', 'Notch', 'off');
    
    legend(labels, 'Interpreter', 'none', 'Box', 'off')
    
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'XColor', 'none')

%% Stats
[p,tbl,stats] = anova1(all_count, all_group);
% [p,tbl,stats] = kruskalwallis(all_count, ALL.group);
c = multcompare(stats, 'alpha', 0.001);

%% T-test
% close all
clc
g1 = 6;
g2 = 7;

G = all_group;
% G = findgroups(G_count);

G1 = G == g1;
Y1 = all_count(G1);
G1 = G(G1);
G2 = G == g2;
Y2 = all_count(G2);
G2 = G(G2);

figure (100)
boxplot([Y1;Y2], [G1;G2])

[~,p] = ttest2(Y1,Y2)


%% Head-wing sync
keepI = cellfun(@(x) isstruct(x) | isobject(x), Free.SACCADE.head2wing);
Saccade = Free.SACCADE(keepI,:);
head_sync = cellfun(@(x) x.sync_head_rate, Saccade.head2wing);
wing_sync = cellfun(@(x) x.sync_wing_rate, Saccade.head2wing);

close all
figure (10) ; clf
subplot(1,1,1)
    boxplot([head_sync wing_sync])

%% Save
fname = 'Wing_count';
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'all', 'all_count', 'all_group', 'labels');

end

function data = get_count(saccade_struct, wave, vel, hw)
if hw == 1
    fname = 'head_saccade';
elseif hw == 2
    fname = 'wing_saccade';
end

keepI = cellfun(@(x) isstruct(x) | isobject(x), saccade_struct.(fname));
Saccade = saccade_struct(keepI,:);

if any(wave)
    waveI = any(Saccade.wave == wave, 2);
    Saccade = Saccade(waveI,:);
end

if any(vel)
    velI = any(Saccade.vel == vel, 2);
    Saccade = Saccade(velI,:);
end
n_speed = length(unique(Saccade.vel))/2;

fly_group = Saccade.fly;
vel_group_all = Saccade.vel;

if length(unique(vel_group_all)) > 1
    vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;
end
[fly_vel_group,fly_group,vel_group] = findgroups(fly_group, findgroups(vel_group_all));

count = cellfun(@(x) x.count, Saccade.(fname), 'UniformOutput', true);
count_vel_fly_mean = splitapply(@(x) nanmean(x,1), count, fly_vel_group);
count_vel_fly_mean = splitapply(@(x) {x}, count_vel_fly_mean, vel_group);
count_vel_fly_mean = cat(2, count_vel_fly_mean{:});

data.vel_group_all = vel_group_all;
data.count = count;
data.count_vel_fly_mean = count_vel_fly_mean;
end