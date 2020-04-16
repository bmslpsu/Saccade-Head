function [] = Saccade_Wing()
%% Saccade_Wing:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N')

clms = N.vel/2;
CC = repmat(hsv(clms),2,1);

Vel = U{1,3}{1};

clearvars -except clms CC Vel PATH COUNT SACCADE SACCADE_STATS FLY GRAND Stim D I U N

%% Saccade Position %%
wing_all = cell(N.fly,N.vel);
wing_interval = cell(N.fly,N.vel);
P = cell(N.fly,N.vel);
head_all = cell(N.fly,N.vel);
for n = 1:N.file
	wing_all{SACCADE.fly(n),SACCADE.vel(n)}(:,end+1) = SACCADE.sWBA{n};
    head_all{SACCADE.fly(n),SACCADE.vel(n)}(end+1,1) = SACCADE.head_saccade{n};
end

bins = -10.5:0.25:10.5;
window = 0.04; % windows size [s]
alpha = 0.01;
showplot = false;

if showplot
    fig = figure (1); clf
    set(fig, 'Color', 'w')
    clear ax h
end

for f = 1:N.fly
    for v = 1:N.vel
        wing_data = wing_all{f,v};
        head_data = head_all{f,v};
        if showplot
            ax(1) = subplot(4,2,1:2); hold on ; ylabel('Wing Probability')
                cla
                histogram(wing_data(:), bins, 'Normalization', 'probability')
        end
        si = 1;
        for t = 1:length(head_data)
            head_saccade = head_data(t);
            span = round((window / head_saccade.Ts) / 2);
            wing_pos = wing_data(:,t);
            wing_vel = diff(wing_pos) / mean(diff(head_saccade.time));
            wing_vel = [wing_vel(1) ; wing_vel];
            if showplot
                ax(2) = subplot(4,2,3:4); hold on ; cla
                    yyaxis left ; hold on ; cla ; ylabel('Position')
                        plot(head_saccade.time, head_saccade.position)
                        plot(head_saccade.SACD.PeakTime, head_saccade.SACD.PeakPos, 'g*')
                    yyaxis right
                        plot(head_saccade.time, wing_pos)

                ax(3) = subplot(4,2,5:6); hold on ; cla
                    yyaxis left ; hold on ; cla
                        plot(head_saccade.time, head_saccade.velocity)
                        plot(head_saccade.SACD.PeakTime, head_saccade.SACD.PeakVel, 'g*')
                    yyaxis right
                        plot(head_saccade.time, wing_vel)
            end
        	
            for s = 1:head_saccade.count
                peak = head_saccade.SACD.PeakIdx(s);
                % peak = randi([21 2000-21],1,1);
                
                interval = (peak - span):(peak + span);
                time = head_saccade.time(interval);
                wing_saccade = wing_data(interval);
                wing_interval{f,v}(:,end+1) = wing_saccade;
                                
                [h,p,~,~] = ttest2(wing_data(:), wing_saccade, ...
                    'Vartype','unequal','Alpha', alpha);
                [rho,p_r] = corr(head_saccade.position(interval), wing_pos(interval));
                
                P{f,v}(si,1) = p;
                P{f,v}(si,2) = h;
                P{f,v}(si,3) = rho;
                P{f,v}(si,4) = p_r;
                si = si + 1;
                
                if showplot
                    x = [time(1), time(1), time(end), time(end)];
                    pos_y = [ax(2).YLim(1), ax(2).YLim(2), ax(2).YLim(2), ax(2).YLim(1)];
                    vel_y = [ax(3).YLim(1), ax(3).YLim(2), ax(3).YLim(2), ax(3).YLim(1)];

                    subplot(4,2,3:4); hold on
                        h1 = patch(x, pos_y, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

                    subplot(4,2,5:6); hold on
                        h2 = patch(x, vel_y, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

                    subplot(4,2,1:2); hold on
                        title([num2str(Vel(v)) ' (°/s)     P = ' num2str(p)], 'FontSize', 12)
                        h3 = histogram(wing_saccade, bins, 'Normalization', 'probability', 'FaceColor', 'g');

                    ax(4) = subplot(4,2,7); cla
                        yyaxis left ; hold on ; cla ; ylabel('Position')
                            plot(head_saccade.time(interval), head_saccade.position(interval))
                            % plot(head_saccade.SACD.PeakTime(s), head_saccade.SACD.PeakPos(s), 'g*')
                        yyaxis right ; hold on ; cla
                            plot(head_saccade.time(interval), wing_pos(interval))

                    ax(5) = subplot(4,2,8); cla
                        yyaxis left ; hold on ; cla ; ylabel('Velocity')
                            plot(head_saccade.time(interval), head_saccade.velocity(interval))
                            % plot(head_saccade.SACD.PeakTime(s), head_saccade.SACD.PeakVel(s), 'g*')
                        yyaxis right ; hold on ; cla
                            plot(head_saccade.time(interval), wing_vel(interval))

                    set(ax, 'LineWidth', 1.5, 'Box', 'on')
                    linkaxes(ax(2:3),'x')
                    linkaxes(ax(4:5),'x')
                    axis tight

                    if rho < 0
                       pause
                    end

                    % pause(0.05)
                    delete([h1 h2])
                    delete(h3)
                end
            end
        end
    end
end

figure (2) ; clf
p_bins = -1:0.1:1;
for v = 1:N.vel
    test = cat(1, P{:,v});
    subplot(2,5,v) ; hold on
    histogram(test(:,3),p_bins)
end


%%

bins = -10.5:0.5:10.5;
window = 0.05; % windows size [s]
alpha = 0.01;
fly = 1:N.fly;

fig = figure (1); clf
set(fig, 'Color', 'w')
ax = gobjects(N.vel,1);
for f = fly
    for v = 1:N.vel
        wing_data = wing_all{f,v}(:);
        head_data = head_all{f,v};
        ax(v) = subplot(2,clms,v); hold on
            cla
            title([num2str(Vel(v)) ' (°/s)'])
            histogram(wing_data, bins, 'Normalization', 'pdf')
        for t = 1:length(head_data)
            head_saccade = head_data(t);
            span = (window / head_saccade.Ts) / 2;
            for s = 1:head_saccade.count
                peak = head_saccade.SACD.PeakIdx(s);
                interval = (peak - span):(peak + span);
                wing_saccade = wing_data(interval);
                wing_interval{f,v}(:,end+1) = wing_saccade;
            end
        end
        histogram(wing_interval{f,v}, bins, 'Normalization', 'pdf')
    end
    % pause
end
set(ax, 'LineWidth', 1.5, 'Box', 'off')
linkaxes(ax,'xy')


end