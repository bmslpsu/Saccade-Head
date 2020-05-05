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
wing_interval.pos = cell(N.fly,N.vel);
wing_interval.vel = cell(N.fly,N.vel);
P = cell(N.fly,N.vel);
head_all = cell(N.fly,N.vel);
for n = 1:N.file
	wing_all{SACCADE.fly(n),SACCADE.vel(n)}(:,end+1) = SACCADE.dWBA{n};
    head_all{SACCADE.fly(n),SACCADE.vel(n)}(end+1,1) = SACCADE.head_saccade{n};
end

bins = -10.5:0.25:10.5;
window = 0.1; % windows size [s]
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
            ax(1) = subplot(4,2,1:2); hold on ; ylabel('Wing Probability') ; xlabel('Wing')
                cla
                histogram(wing_data(:), bins, 'Normalization', 'probability', 'FaceColor', 'r')
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
                    yyaxis left ; hold on ; cla ; ylabel('Head Position')
                        plot(head_saccade.time, head_saccade.position)
                        plot(head_saccade.SACD.PeakTime, head_saccade.SACD.PeakPos, 'g*')
                    yyaxis right ; hold on ; cla ; ylabel('Wing Position')
                        plot(head_saccade.time, wing_pos)

                ax(3) = subplot(4,2,5:6); hold on ; cla ; ylabel('Wing Velocity')
                    yyaxis left ; hold on ; cla ; ylabel('Wing Position')
                        plot(head_saccade.time, head_saccade.velocity)
                        plot(head_saccade.SACD.PeakTime, head_saccade.SACD.PeakVel, 'g*')
                    yyaxis right ; hold on ; cla ; ylabel('Wing Velocity')
                        plot(head_saccade.time, wing_vel)
            end
        	
            for s = 1:head_saccade.count
                peak = head_saccade.SACD.PeakIdx(s);
                % peak = randi([1+span 2000-span-1],1,1);
                
                interval = (peak - span):(peak + span);
                time = head_saccade.time(interval);
                
                wing_saccade_pos = wing_pos(interval);
                wing_saccade_vel = wing_vel(interval);
                
                wing_interval.pos{f,v}(:,end+1) = wing_saccade_pos;
                wing_interval.vel{f,v}(:,end+1) = wing_saccade_vel;
                                
                [h,p,~,~] = ttest2(wing_data(:), wing_saccade_pos, ...
                    'Vartype','unequal','Alpha', alpha);
                [rho_pos,p_r_pos] = corr(head_saccade.position(interval), wing_pos(interval));
                [rho_vel,p_r_vel] = corr(head_saccade.velocity(interval), wing_vel(interval));
                
                P{f,v}(si,1) = p;
                P{f,v}(si,2) = h;
                P{f,v}(si,3) = rho_pos;
                P{f,v}(si,4) = p_r_pos;
              	P{f,v}(si,5) = rho_vel;
                P{f,v}(si,6) = p_r_vel;
                si = si + 1;
                
                if showplot
                    x = [time(1), time(1), time(end), time(end)];
                    
                    vel_y = [ax(3).YLim(1), ax(3).YLim(2), ax(3).YLim(2), ax(3).YLim(1)];

                    subplot(4,2,3:4); hold on
                     pos_y = [ax(2).YLim(1), ax(2).YLim(2), ax(2).YLim(2), ax(2).YLim(1)];
                        h1 = patch(x, pos_y, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

                    subplot(4,2,5:6); hold on
                        
                        h2 = patch(x, vel_y, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');

                    subplot(4,2,1:2); hold on
                        title([num2str(Vel(v)) ' (°/s)     P = ' num2str(p)], 'FontSize', 12)
                        h3 = histogram(wing_saccade_pos, bins, 'Normalization', 'probability', 'FaceColor', 'g');

                    ax(4) = subplot(4,2,7); cla ; title(['R^{2} = ' num2str(rho_pos) ...
                        '   P = ' num2str(p_r_pos)])
                        yyaxis left ; hold on ; cla ; ylabel('Head Position')
                            plot(head_saccade.time(interval), head_saccade.position(interval))
                            % plot(head_saccade.SACD.PeakTime(s), head_saccade.SACD.PeakPos(s), 'g*')
                        yyaxis right ; hold on ; cla ; ylabel('Wing Velocity')
                            plot(head_saccade.time(interval), wing_pos(interval))

                    ax(5) = subplot(4,2,8); cla ; title(['R^{2} = ' num2str(rho_vel) ...
                        '   P = ' num2str(p_r_vel)])
                        yyaxis left ; hold on ; cla ; ylabel('Head Velocity')
                            plot(head_saccade.time(interval), head_saccade.velocity(interval))
                            % plot(head_saccade.SACD.PeakTime(s), head_saccade.SACD.PeakVel(s), 'g*')
                        yyaxis right ; hold on ; cla ; ylabel('Wing Velocity')
                            plot(head_saccade.time(interval), wing_vel(interval))

                    set(ax, 'LineWidth', 1.5, 'Box', 'on')
                    linkaxes(ax(2:3),'x')
                    linkaxes(ax(4:5),'x')
                    axis tight
                    
                    %pause()
%                  %   if rho_pos < -0.5
%                        pause
%                     end

                    delete([h1 h2])
                    delete(h3)
                end
            end
        end
    end
end

%% Correlation R^2 & P values
all = reshape(P,N.vel*N.fly,1);
all = cat(1, all{:});
r_bins = -1:0.1:1;
p_bins = 0:0.05:1;
clear ax
fig = figure (2) ; clf
set(fig, 'Color', 'w')
ax(1) = subplot(2,2,1); hold on ; title('Position')
    histogram(all(:,3), r_bins, 'Normalization', 'probability')
    xlabel('R^{2}')
    ylabel('Probability')
ax(2) = subplot(2,2,3); hold on
    histogram(all(:,4), p_bins, 'Normalization', 'probability')
    xlabel('p-value')
    ylabel('Probability')
ax(3) = subplot(2,2,2); hold on ; title('Velocity')
    histogram(all(:,5), r_bins, 'Normalization', 'probability')
    xlabel('R^{2}')
ax(4) = subplot(2,2,4); hold on
    histogram(all(:,6), p_bins, 'Normalization', 'probability')
    xlabel('p-value')
set(ax, 'LineWidth', 1.5, 'Box', 'on')

%% Correlation R^2 & P values
fig = figure (2) ; clf
set(fig, 'Color', 'w')
ax(1) = subplot(2,2,1); hold on ; title('Position')
    histogram(all(:,3), r_bins, 'Normalization', 'probability')
    xlabel('R^{2}')
    ylabel('Probability')
end