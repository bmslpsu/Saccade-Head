function [] = Saccade_WingMean()
%% Saccade_WingMean:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat'},'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE','U','N','D')

%% Active Tracking vs Landing Behavior %%
clearvars -except SACCADE U N D
n_speed = N.vel/2;
keepI = cellfun(@(x) isstruct(x) | isobject(x), SACCADE.wing_saccade);
Saccade = SACCADE(keepI,:);
% Saccade = Saccade(98:end,:);

tt = Saccade.wing_saccade{1}.time;
wba = cellfun(@(x) x.position - x.position(1), Saccade.wing_saccade, 'UniformOutput', false);

flip_vel = Saccade.vel > n_speed;
wba(flip_vel) = cellfun(@(x) -x, wba(flip_vel), 'UniformOutput', false);
Saccade.vel(Saccade.vel > n_speed) = Saccade.vel(Saccade.vel > n_speed) - n_speed;

[vel_group,Vel] = findgroups(Saccade.vel);
n_vel = length(Vel);

[fly_group,Fly] = findgroups(Saccade.fly);
n_fly = length(Fly);

WBA = cell(n_fly,n_vel);
for n = 1:size(Saccade,1)
    v = vel_group(n);
    f = fly_group(n);
    if Saccade.vel  > 5
        flip = -1;
    else
        flip = 1;
    end
    WBA{f,v} = flip * [WBA{f,v} , wba{n}];
end

WBA_fly = cellfun(@(x) mean(x,2), WBA, 'UniformOutput', false);
WBA_all = cell(1,n_vel);
for v = 1:n_vel
    WBA_all{v} = cat(2, WBA{:,v});
    WBA_fly_all{v} = cat(2, WBA_fly{:,v});
end

WBA_all_mean = cellfun(@(x) mean(x,2), WBA_all, 'UniformOutput', false);
WBA_grand_mean = cellfun(@(x) mean(x,2), WBA_fly_all, 'UniformOutput', false);
WBA_grand_std = cellfun(@(x) std(x,[],2), WBA_fly_all, 'UniformOutput', false);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 8 5])
movegui(fig, 'center')
ax(1) = subplot(1,1,1) ; cla ; hold on
    %plot(tt, WBA_all{1}, 'Color', [0.5 0.5 0.5 0.3], 'LineWidth', 0.25)
    %plot(tt, WBA_all{2}, 'Color', [0 0 0.7 0.3], 'LineWidth', 0.25)
    plot(tt, WBA_fly_all{1}, 'r', 'LineWidth', 1)
    %plot(tt, WBA_fly_all{2},'b', 'LineWidth', 1)
    [~] = PlotPatch(WBA_grand_mean{1}, WBA_grand_std{1}, tt, 1, 1, 'k', 'r', 0.3, 3);
    %[~] = PlotPatch(WBA_grand_mean{2}, WBA_grand_std{2}, tt, 1, 1, 'k', 'b', 0.3, 3);   
    %plot(tt, WBA_all_mean{1}, 'g', 'LineWidth', 2)
    %plot(tt, WBA_all_mean{2}, 'g', 'LineWidth', 2)
    %plot(WBA_grand{1}+WBA_grand{2},'g','LineiIdth',3)
    xlabel('Time (s)')
    ylabel('\DeltaWBA (°)')
    ylim(30*[-1 1])
    
set(ax, 'LineWidth', 2)

end