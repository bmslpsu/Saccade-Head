function [] = Saccade_time_constant()
%% Passive_head_window:
root = 'H:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat'},'Select data file', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'U','N','SACCADE')

%%
clc
clearvars -except SACCADE U N
keepI = cellfun(@(x) length(x) > 1, SACCADE.head_scd_pos);
Saccade = SACCADE(keepI,:);
n_speed = N.vel / 2;
n_trial = size(Saccade,1);

tt = Saccade.scd_time{1,1};
Fs = round(1 / mean(diff(tt)));
[~,zI] = min(abs(tt));
t_win = [-0.02 0.08];
fit_win = round( zI + (t_win * Fs) );
fit_win = fit_win(1):fit_win(2);
tt = tt(fit_win);
time_active = tt - tt(1);

n_scd = cellfun(@(x) size(x,2), Saccade.head_scd_pos, 'UniformOutput', false);
n_scd_all = sum(cat(1,n_scd{:}));

Saccade_Table = splitvars(table(nan(n_scd_all,5)));
Saccade_Table.Properties.VariableNames = {'fly', 'vel', 'wave', 'tau', 'r2'};
pp = 1;
bad_count = 0;
for n = 1:n_trial
    scd_all = Saccade.head_scd_pos{n}(fit_win,:);
    fly = Saccade.fly(n);
    vel = Saccade.vel(n);
    for s = 1:n_scd{n}
        xx = scd_all(:,s);
        [fitresult, gof] = exp_fit(time_active, xx, false);
        tau = 1000 * (1 / fitresult.b);
        r2 = gof.rsquare;
        Saccade_Table{pp,:} = [fly vel U.wave tau r2];
        if r2 < 0.8
            bad_count = bad_count + 1;
        end
        pp = pp + 1;

%         title(tau)
%         pause
%         close all
    end
end
disp('Done')

%% Time constant
fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', 1*[2 2 2 2])
clear ax
ax(1) = subplot(1,1,1); hold on
ylabel('Time Constant (ms)')

Saccade_Table_cut = Saccade_Table(Saccade_Table.r2 > 0.8,:);
velI = Saccade_Table_cut.vel;
velI(velI > n_speed) = velI(velI > n_speed) - n_speed;

bx = boxplot(Saccade_Table_cut.tau, 'Width', 0.5, 'Symbol', '.', 'Whisker', 2);
h = get(bx(5,:),{'XData','YData'});
for kk = 1:size(h,1)
   patch(h{kk,1}, h{kk,2}, 'k',  'EdgeColor', 'none');
end
ww = 1;
set(findobj(ax(ww),'tag','Median'), 'Color', 'w','LineWidth', 1);
set(findobj(ax(ww),'tag','Box'), 'Color', 'none');
set(findobj(ax(ww),'tag','Upper Whisker'), 'Color', 'k','LineStyle','-');
set(findobj(ax(ww),'tag','Lower Whisker'), 'Color', 'k','LineStyle','-');
ax(ww).Children = ax(ww).Children([end 1:end-1]);
ylim([0 60])

set(ax , 'LineWidth', 1, 'Box', 'off', 'XColor', 'none')

%% Anova for speed & fly effects
[p,tbl,stats] = anovan(Saccade_Table_cut.tau, {Saccade_Table_cut.vel Saccade_Table_cut.fly}, ...
    'model','interaction', 'varnames', {'vel','Fly'});

%% Save
savedir = 'H:\DATA\Rigid_Data\Saccade\processed';
filename = ['Ramp_time_constant_wave=' num2str(U.wave) '.mat'];
save(fullfile(savedir, filename), 'Saccade_Table','U','N','-v7.3')

end