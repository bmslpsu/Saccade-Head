function [] = Saccade_ExampleTrial_Leg()
%% Saccade_ExampleTrials:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat'},'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'SACCADE','U','N','D')

%% Active Tracking vs Landing Behavior %%
clearvars -except SACCADE U N D

FIG = figure (1) ; clf
FIG.Units = 'inches';
FIG.Position = [2 2 6 5];
FIG.Name = 'Saccade Position';
FIG.PaperPositionMode = 'auto';
movegui(FIG,'center')
FIG.Color = 'w';
clear ax h

idx = 153;
idx = 158;

trial = SACCADE.head_saccade{idx};
leg_data = SACCADE.leg{idx};
leg_class = max(leg_data,[],2) > 0.9;
leg_class_filt = logical(ceil(medfilt1(double(leg_class),40)));
leg_start = trial.Ts*find(leg_class_filt, 1, 'first');

ax(1) = subplot(4,1,1) ; hold on ; title(['Stimulus: ' num2str(U.vel{1}(SACCADE.vel(idx))) ' (°/s)'])
    ylabel('Position (°)')
    h(1) = plot(trial.time,trial.position,'k');
    plot(trial.time, zeros(trial.n,1),'--','Color',[0.5 0.5 0.5])
    if trial.count~=0
        for ww = 1:trial.count
           plot(trial.saccades{ww}.Time,trial.saccades{ww}.Position,...
               'LineWidth', 1, 'Color', 'c')
%            plot(trial.intervals{ww}.Time,trial.intervals{ww}.Position,...
%                'LineWidth', 1, 'Color', 0.7*trial.cmap(ww,:))
        end
%         plot(trial.starts.time , trial.starts.position , '*g')
%         plot(trial.peaks.time  , trial.peaks.position  , '*b')
%         plot(trial.ends.time   , trial.ends.position   , '*r')
    end

    ax(1).YLim = 20*[-1 1];
  	%plot([leg_time leg_time], ax(1).YLim, 'r', 'LineWidth', 2)
    
    x_leg = [leg_start leg_start ax(1).XLim(2) ax(1).XLim(2)];
    y_leg = [ax(1).YLim(1) ax(1).YLim(2) ax(1).YLim(2) ax(1).YLim(1)];
    hp = patch(x_leg,y_leg,'r','FaceAlpha',0.3);
    uistack(hp,'bottom')

ax(2) = subplot(4,1,2) ; hold on
    ylabel('Velocity (°/s)')
    h(2) = plot(trial.time,trial.velocity,'k'); 
    if trial.count~=0
        for ww = 1:trial.count
           plot(trial.saccades{ww}.Time,trial.saccades{ww}.Velocity,...
               'LineWidth', 1, 'Color', 'c')
%            plot(trial.intervals{ww}.Time,trial.intervals{ww}.Velocity,...
%                'LineWidth', 1, 'Color', 0.7*trial.cmap(ww,:))
        end
%         plot(trial.starts.time , trial.starts.velocity  , '*g')
%         plot(trial.peaks.time  , trial.peaks.velocity   , '*b')
%         plot(trial.ends.time   , trial.ends.velocity    , '*r')
    end
%     plot(trial.time, -trial.threshold*ones(trial.n,1) , '--', 'Color', [0.5 0.5 0.5])
%     plot(trial.time,  trial.threshold*ones(trial.n,1) , '--', 'Color', [0.5 0.5 0.5])
%     ax(2).YLim = 700*[-1 1];

ax(3) = subplot(4,1,3) ; hold on
    if trial.count~=0
        ylabel({'Removed', 'Position (°)'})                
        plot(trial.shift.Time, trial.shift.IntrpPosition,      	 'k', 'LineWidth', 1);
        % plot(trial.shift.Time, trial.shift.Position,             'k', 'LineWidth', 1);
        % plot(trial.removed_all.Time, trial.removed_all.Position, 'c', 'LineWidth', 1);
    end
    %plot(trial.time, trial.stimlus_position, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    
ax(4) = subplot(4,1,4) ; hold on
    if trial.count~=0
        ylabel({'Removed', 'Velocity (°/s)'})
        xlabel('Time (s)')
        [b,a] = butter(2,2/(100/2),'low');
        filtvel = filtfilt(b, a, trial.shift.IntrpVelocity);
        plot(trial.shift.Time, trial.shift.IntrpVelocity,   'r', 'LineWidth', 0.5);
        plot(trial.shift.Time, trial.shift.Velocity,        'k', 'LineWidth', 0.5);
        plot(trial.shift.Time, filtvel,                     'c', 'LineWidth', 1);
    end
    % ax(4).YLim = ax(2).YLim;
    ylim(400*[-1 1])
    
set(ax,'LineWidth',1,'FontWeight','bold','FontSize',8)
linkaxes(ax,'x')
set(h,'LineWidth',0.5)
align_Ylabels(FIG)

end