function [] = Ramp_roll_intervals()
%% Ramp_saccade_window_roll:
root = 'E:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'U','N','SACCADE')

%%
clc
clearvars -except SACCADE U N
keepI = cellfun(@(x) length(x) > 1, SACCADE.yaw_scd_pos);
Saccade = SACCADE(keepI,:);
if N.vel > 1
    n_speed = N.vel / 2;
else
    n_speed = N.vel;
end
velI = Saccade.vel;
velI(velI > n_speed) = velI(velI > n_speed) - n_speed;

Fs = Saccade.head_yaw{1}.Fs;
[b,a] = butter(5, 10 / (Fs/2), 'low');
n_file = size(Saccade,1);
Roll = cell(N.fly,1);
Yaw = cell(N.fly,1);
for n = 1:n_file
    [~,ints,~,~] = getSaccade(Saccade.head_yaw{n}, Saccade.head_roll{n}, [], [], false);
    dir = sign(U.vel{1}(Saccade.vel(n)));
    ints = ints(2:end);
    ints = cellfun(@(x) x*dir, ints, 'UniformOutput', false);
    try
        ints = cellfun(@(x) filtfilt(b, a, x), ints, 'UniformOutput', false);
    catch
        disp('too short')
    end
    ints = cellfun(@(x) x - x(1), ints, 'UniformOutput', false);
    Roll{Saccade.fly(n)} = [Roll{Saccade.fly(n)} ; ints];
    
    [~,ints,~,~] = getSaccade(Saccade.head_yaw{n}, Saccade.head_yaw{n}.position, [], [], false);
    dir = sign(U.vel{1}(Saccade.vel(n)));
    ints = ints(2:end);
    ints = cellfun(@(x) x*dir, ints, 'UniformOutput', false);
    try
        ints = cellfun(@(x) filtfilt(b, a, x), ints, 'UniformOutput', false);
    catch
        disp('too short')
    end
    ints = cellfun(@(x) x - x(1), ints, 'UniformOutput', false);
    Yaw{Saccade.fly(n)} = [Yaw{Saccade.fly(n)} ; ints];
end

for f = 1:N.fly
    Roll{f} = nancat(2, Roll{f}{:});
    nanI = isnan(Roll{f}(1,:));
    Roll{f} = Roll{f}(:,~nanI);
    
    Yaw{f} = nancat(2, Yaw{f}{:});
    Yaw{f} = Yaw{f}(:,~nanI);
end
Roll_fly = cellfun(@(x) nanmean(x,2), Roll, 'UniformOutput', false);
Roll_fly = nancat(2, Roll_fly{:});
Roll_grand_mean = nanmean(Roll_fly,2);
Roll_grand_std = nanstd(Roll_fly, [], 2);

Yaw_fly = cellfun(@(x) nanmean(x,2), Yaw, 'UniformOutput', false);
Yaw_fly = nancat(2, Yaw_fly{:});
Yaw_grand_mean = nanmean(Yaw_fly,2);
Yaw_grand_std = nanstd(Yaw_fly, [], 2);

% int_time = Saccade.head_yaw{n}.Ts * (1:size(Roll_all,1))';
int_time_fly = Saccade.head_yaw{n}.Ts * (1:size(Roll_fly,1))';

yaw_color = [0 0 1];
roll_color = [0.7 0.4 0.1];

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', 2*[2 2 2.5 1.5])
movegui(fig, 'center')
clear ax H
ax(1,1) = subplot(1,1,1); hold on ; cla
plot(int_time_fly, Yaw_fly, 'Color', 0.7*yaw_color)
PlotPatch(Yaw_grand_mean, Yaw_grand_std, int_time_fly, 1, 1, yaw_color, 0.7*yaw_color, 0.2, 1.5);

plot(int_time_fly, Roll_fly, 'Color', 0.7*roll_color)
PlotPatch(Roll_grand_mean, Roll_grand_std, int_time_fly, 1, 1, roll_color , 0.7*roll_color, 0.2, 1.5);

yline(0, '--k')

xlim([-0.05 2])
ylim([-20 20])
xlabel('time (s)')
ylabel('roll (°)')
set(ax, 'LineWidth', 1, 'Color', 'none')

end