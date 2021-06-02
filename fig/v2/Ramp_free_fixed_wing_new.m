function [] = Ramp_free_fixed_wing_new()
%% Ramp_free_fixed_wing_new:
root = 'E:\DATA\Rigid_Data\Saccade';

[Free,FreePath] = uigetfile({'*.mat'},'Select free data', root, 'MultiSelect','off');
[Fixed,FixedPath] = uigetfile({'*.mat'},'Select free data', root, 'MultiSelect','off');

Free = load(fullfile(FreePath,Free),'SACCADE','U','N','D');
Fixed = load(fullfile(FixedPath,Fixed),'SACCADE','U','N','D');

%% Head-free vs Head-fixed dwba response
clearvars -except SACCADE U N D Free Fixed root

head_color = [0 0 1];
free_color_30 = [1 0 0];
free_color_60 = [0.1 0.8 0.3];
fixed_color = [0 0 0];

norm = true;
head_prop = 'head_saccade';
wing_prop = 'wing_saccade';
Fc = []; % high-pass filter cutoff frequency [Hz]
detrend_n = [];
[Head,~] = get_data(Free, head_prop, norm, [], [], [1 3]);
[WBA_Free_30,time_free_30] = get_data(Free, wing_prop, norm, Fc, detrend_n, [1 3]);
[WBA_Free_60,time_free_60] = get_data(Free, wing_prop, norm, Fc, detrend_n, [2 4]);
[WBA_Fixed,time_fixed] = get_data(Fixed, wing_prop, norm, Fc, detrend_n, [1 2]);

%% Free vs Fixed time series
fig = figure (10) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 5 3])
movegui(fig, 'center')
clear ax h
ax(1) = subplot(1,1,1) ; cla ; hold on
    plot([0 10], [0 0], 'k', 'LineWidth', 1)
    
%     plot(time_free, WBA_Free.all{1}(:,:), 'Color', [0.7*free_color 0.3], 'LineWidth', 0.5)
%     plot(time_fixed, WBA_Fixed.all{1}(:,:), 'Color', [0.7*fixed_color 0.3], 'LineWidth', 0.5)
%     
%     [~] = PlotPatch(Head.all_mean{1}, Head.all_std{1}, time_free, 1, 1, head_color, head_color, 0.3, 3);
%     [~] = PlotPatch(WBA_Free.all_mean{1}, WBA_Free.all_std{1}, time_free, 1, 1, free_color, free_color, 0.3, 3);
%     [~] = PlotPatch(WBA_Fixed.all_mean{1}, WBA_Fixed.all_std{1}, time_fixed, 1, 1, fixed_color, fixed_color, 0.3, 3);
    
    %plot(time_free_30, WBA_Free_30.fly_mean_all {1}, 'Color', [free_color_30 0.7], 'LineWidth', 0.5)
  	%plot(time_free_60, WBA_Free_60.fly_mean_all {1}, 'Color', [free_color_60 0.7], 'LineWidth', 0.5)
    %plot(time_fixed, WBA_Fixed.fly_mean_all {1}, 'Color', [fixed_color 0.7], 'LineWidth', 0.5)
    %plot(time_free, Head.fly_mean_all {1}, 'Color', [head_color 0.5], 'LineWidth', 1)
   
 	[~,h.mean(1)] = PlotPatch(WBA_Free_30.grand_mean{1}, WBA_Free_30.grand_std{1}, time_free_30, ...
        1, 1, free_color_30, free_color_30, 0.3, 2);
 	[~,h.mean(2)] = PlotPatch(WBA_Free_60.grand_mean{1}, WBA_Free_60.grand_std{1}, time_free_30, ...
        1, 1, free_color_60, free_color_60, 0.3, 2);
    [~,h.mean(3)] = PlotPatch(WBA_Fixed.grand_mean{1}, WBA_Fixed.grand_std{1}, time_fixed, ...
        1, 1, fixed_color, fixed_color, 0.3, 2);
%   	[~,h.mean(3)] = PlotPatch(Head.grand_mean{1}, Head.grand_std{1}, time_free_30, ...
%         1, 1, head_color, head_color, 0.3, 2);
    
    xlabel('Time (s)')
    ylabel('\DeltaWBA (°)')
    ylim(30*[-1 1])
    xlim([-0.2 10])
    xticks(0:10)
    
    uistack(h.mean, 'top')
    
set(ax, 'LineWidth', 1)

%% Mean & variance distribution in time
ALL.mean = [ cat(2,Head.fly_mean_time {:})' ; cat(2,WBA_Free_30.fly_mean_time{:})' ; ...
                    cat(2,WBA_Fixed.fly_mean_time{:})' ];
                
ALL.std = [ cat(2,Head.fly_std_time {:})' ; cat(2,WBA_Free_30.fly_std_time{:})' ; ...
                    cat(2,WBA_Fixed.fly_std_time{:})' ];

ALL.range = [ cat(2,Head.fly_range_time{:})' ; cat(2,WBA_Free_30.fly_range_time{:})' ; ...
                    cat(2,WBA_Fixed.fly_range_time{:})' ];
                            
ALL.group = [ ones(size(cat(2,WBA_Free_30.fly_mean_time{:})')) ; ...
                    2*ones(size(cat(2,WBA_Free_30.fly_mean_time{:})')) ; ...
                    3*ones(size(cat(2,WBA_Fixed.fly_mean_time{:})')) ];
                
fig = figure (12) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 6 1.5])
movegui(fig, 'center')
clear ax b
ax(1) = subplot(1,3,1) ; cla ; hold on ; ylabel('Mean (°)')
    ylim(50*[-1 1])
    yticks(-50:25:50)
    b(1) = boxchart(ALL.group, ALL.mean);

ax(2) = subplot(1,3,2) ; cla ; hold on ; ylabel('STD (°)') ; 
    ylim([-1 25])
    yticks(0:5:25)
    b(2) = boxchart(ALL.group, ALL.std);
    
ax(3) = subplot(1,3,3) ; cla ; hold on ; ylabel('Range (°)') ; 
%     ylim([-1 50])
%     yticks(0:10:50)
    b(3) = boxchart(ALL.group, ALL.range);
    
set(ax, 'LineWidth', 1, 'Box', 'off', 'Color', 'none')
set(ax(1:end), 'XColor', 'none')
set(b, 'MarkerStyle', 'none', 'JitterOutliers', 'on', 'BoxFaceAlpha', 0.7)
linkaxes(ax,'x')

%% Stats
% data = ALL.mean;
data = ALL.std;
[p,tbl,stats] = anova1(data, ALL.group);
% [p,tbl,stats] = kruskalwallis(data, ALL.group);
c = multcompare(stats, 'alpha', 0.001);

%% Save
fname = 'Wing_response';
savedir = fullfile(root,'processed');
mkdir(savedir)
save(fullfile(savedir, [fname '.mat']), 'Fc', 'detrend_n', 'WBA_Free', 'WBA_Fixed', 'Head', 'ALL');

end

%% Function to get wing data from stucture
function [DATA,tt] = get_data(WingStruct, prop, norm, Fc, detrend_n, velI)
    % Get wing saccade objects
    n_speed = WingStruct.N.vel/2;
    keepI = cellfun(@(x) isstruct(x) | isobject(x), WingStruct.SACCADE.wing_saccade);
    Saccade = WingStruct.SACCADE(keepI,:);
    %Saccade = Saccade(1:102,:);
    
    if any(velI)
        keepI = any(Saccade.vel == velI, 2);
        Saccade = WingStruct.SACCADE(keepI,:);
    end

    % Pull out wba and time
  	tt = Saccade.(prop){1}.time;
    Fs = round(Saccade.(prop){1}.Fs);
    N = round(Saccade.(prop){1}.n);
    if strcmp(prop, '')
        data = cellfun(@(x) x.extra.dwba, Saccade.(prop), 'UniformOutput', false);
    else
        data = cellfun(@(x) x.position - 0*x.position(1), Saccade.(prop), 'UniformOutput', false);
    end
    %wba = cellfun(@(x) x.position - mean(x.position), Saccade.(prop), 'UniformOutput', false);
    
    if ~isempty(Fc)
       [b,a] = butter(3, Fc / (Fs/2), 'low');
       data = cellfun(@(x) filtfilt(b, a, x), data, 'UniformOutput', false);
    end
    
    if ~isempty(detrend_n)
        data = cellfun(@(x) detrend(x, detrend_n), data, 'UniformOutput', false);
    end
    
    % Flip CCW to CW if specified
    if norm
        flip_vel = Saccade.vel > n_speed;
        data(flip_vel) = cellfun(@(x) -x, data(flip_vel), 'UniformOutput', false);
        Saccade.vel(Saccade.vel > n_speed) = Saccade.vel(Saccade.vel > n_speed) - n_speed;
    end

    % Group by fly and speed
    [vel_group,Vel] = findgroups(Saccade.vel);
    n_vel = length(Vel);

    [fly_group,Fly] = findgroups(Saccade.fly);
    n_fly = length(Fly);
    
    % Calculate auto-correlation
   	%[acor,lags] = cellfun(@(x) autocorr(x, 'NumLags', N-1), data, 'UniformOutput', false);

    % Group wba signals by fly and speed
    DATA.fly = cell(n_fly,n_vel);
    for n = 1:size(Saccade,1)
        v = vel_group(n);
        f = fly_group(n);
        if Saccade.vel  > 5
            flip = -1;
        else
            flip = 1;
        end
        DATA.fly{f,v} = flip * [DATA.fly{f,v} , data{n}];
    end
    
    % Fly stats
    DATA.fly_mean    	= cellfun(@(x) mean(x,2), DATA.fly, 'UniformOutput', false);
    DATA.fly_median   	= cellfun(@(x) mean(x,2), DATA.fly, 'UniformOutput', false);
    DATA.fly_std      	= cellfun(@(x) std(x,[],2), DATA.fly, 'UniformOutput', false);
    
    DATA.fly_mean_time 	= cellfun(@(x) mean(x,1), DATA.fly, 'UniformOutput', false);
    DATA.fly_std_time 	= cellfun(@(x) std(x,[],1), DATA.fly, 'UniformOutput', false);
    DATA.fly_range_time	= cellfun(@(x) range(x,1), DATA.fly, 'UniformOutput', false);
    
  	DATA.fly_std_time_mean	= cellfun(@(x) mean(x), DATA.fly_std_time, 'UniformOutput', true);
    
    % Group wba signals by speed
    DATA.all = cell(1,n_vel);
    for v = 1:n_vel
        DATA.all{v} = cat(2, DATA.fly{:,v});
        DATA.fly_mean_all{v} = cat(2, DATA.fly_mean{:,v});
    end

    % Grand stats
    DATA.all_mean = cellfun(@(x) mean(x,2), DATA.all, 'UniformOutput', false);
    DATA.all_std = cellfun(@(x) std(x,[],2), DATA.all, 'UniformOutput', false);
    
    DATA.grand_mean = cellfun(@(x) mean(x,2), DATA.fly_mean_all, 'UniformOutput', false);
    DATA.grand_std = cellfun(@(x) std(x,[],2), DATA.fly_mean_all, 'UniformOutput', false);  
end
