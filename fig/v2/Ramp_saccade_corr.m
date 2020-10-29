function [] = Ramp_saccade_corr()
%% Ramp_saccade_corr:
load('H:\DATA\Rigid_Data\Saccade\combined\Ramp_All_Stats.mat')

%% Saccade Statistics by velocity
wave_group_all = All_Stats.wave;
wave = unique(wave_group_all);
% All_Stats = All_Stats(wave_group_all==1,:);

vel_group_all = All_Stats.vel;
fly_group_all = All_Stats.fly;

Vel = unique(All_Stats.Vel);
n_speed = length(Vel)/2;
Speed = Vel(n_speed+1:end);
CC = repmat(hsv(n_speed),2,1);

vel_group_all(vel_group_all > n_speed) = vel_group_all(vel_group_all > n_speed) - n_speed;

[fly_vel_group, vel_group, fly_group] = findgroups(vel_group_all, fly_group_all);
[wave_vel_group, ~, wave_group] = findgroups(vel_group_all, wave_group_all);


%% Start position vs peak velocity
All_Stats = All_Stats(~isnan(All_Stats.StartPos),:);
x = All_Stats.StartPos .* All_Stats.Direction;
y = All_Stats.PeakVel .* All_Stats.Direction;
% x = (All_Stats.StartPos - 1*mean(All_Stats.StartPos));
% y = All_Stats.PeakVel;

[R,P] = corr(x, y);
[coeff,~] = polyfit(x, y, 1);
xfit = -20:20;
yfit = polyval(coeff,xfit);

fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches')
set(fig, 'Position', 1.5*[2 2 2 2])
movegui(fig, 'center')
clear ax h
ax(1) = subplot(1,1,1); cla ; hold on
    h(1) = histogram2(x, y, -30:1:30, -1200:30:1200 , 'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'ShowEmptyBins','on', 'Normalization','Probability','DisplayStyle','tile');
    plot(xfit, yfit, 'r')
    axis tight
    grid off
    xlabel('Start Position (°)')
    ylabel('Peak Velocity (°/s)')
    zlabel('Probability')
    cbar = colorbar('Location', 'eastoutside');
    cbar.Ticks = [0 max(h(1).Values,[],'all')];
	title('Moving')

end