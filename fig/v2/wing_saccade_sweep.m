function [] = wing_saccade_sweep()
%% wing_saccade_sweep:
%% Head free
root = 'E:\DATA\Rigid_Data\Saccade\sweep';
[FILE,PATH] = uigetfile({'*.mat'},'Select files', root, 'MultiSelect','on');
FILE = cellstr(FILE);

n_sweep = length(FILE);
Free = splitvars(table(nan(n_sweep,2)));
Free.Properties.VariableNames = {'stdThresh','ampCut'};
Free = [Free splitvars(table(num2cell(nan(n_sweep,6))))];
Free.Properties.VariableNames(3:end) = ...
    {'vel','head_sync','wing_sync','wing_rate','wing_rate_30', 'wing_rate_60'};
for n = 1:n_sweep
    data = load(fullfile(PATH,FILE{n}),'wing','SACCADE','D','U','N');
    Free.stdThresh(n) = data.wing.thresh(3);
    Free.ampCut(n) = data.wing.amp_cut;
    
    clear Saccade
    keepI = cellfun(@(x) isstruct(x) | isobject(x), data.SACCADE.head2wing);
    Saccade = data.SACCADE(keepI,:);
    Free.vel{n} = data.D.vel(keepI);
    Free.head_sync{n} = cellfun(@(x) x.sync_head_rate, Saccade.head2wing);
    Free.wing_sync{n} = cellfun(@(x) x.sync_wing_rate, Saccade.head2wing);
    
    Free.wing_rate{n} = cellfun(@(x) x.rate, Saccade.wing_saccade);
    Free.wing_rate_30{n} = Free.wing_rate{n}(abs(Free.vel{n}) == 30);
    Free.wing_rate_60{n} = Free.wing_rate{n}(abs(Free.vel{n}) == 60);
end
Free = sortrows(Free, 'stdThresh');
Free = sortrows(Free, 'ampCut');

% Collect data
Free_all = [];
Free_all.stdThresh_all = cellfun(@(x,y) y*ones(size(x)), Free.wing_rate, num2cell(Free.stdThresh), ...
    'UniformOutput', false);
Free_all.stdThresh_all = cat(1,Free_all.stdThresh_all{:});

Free_all.ampCut_all = cellfun(@(x,y) y*ones(size(x)), Free.wing_rate, num2cell(Free.ampCut), ...
    'UniformOutput', false);
Free_all.ampCut_all = cat(1,Free_all.ampCut_all{:});

Free_all.stdThresh_wing_rate_30 = cellfun(@(x,y) y*ones(size(x)), Free.wing_rate_30, num2cell(Free.stdThresh), ...
    'UniformOutput', false);
Free_all.stdThresh_wing_rate_30 = cat(1,Free_all.stdThresh_wing_rate_30{:});

Free_all.stdThresh_wing_rate_60 = cellfun(@(x,y) y*ones(size(x)), Free.wing_rate_60, num2cell(Free.stdThresh), ...
    'UniformOutput', false);
Free_all.stdThresh_wing_rate_60 = cat(1,Free_all.stdThresh_wing_rate_60{:});

Free_all.ampCut_wing_rate_30 = cellfun(@(x,y) y*ones(size(x)), Free.wing_rate_30, num2cell(Free.ampCut), ...
    'UniformOutput', false);
Free_all.ampCut_wing_rate_30 = cat(1,Free_all.ampCut_wing_rate_30{:});

Free_all.ampCut_wing_rate_60 = cellfun(@(x,y) y*ones(size(x)), Free.wing_rate_60, num2cell(Free.ampCut), ...
    'UniformOutput', false);
Free_all.ampCut_wing_rate_60 = cat(1,Free_all.ampCut_wing_rate_60{:});

Free_all.head_sync = cat(1, Free.head_sync{:});
Free_all.wing_sync = cat(1, Free.wing_sync{:});
Free_all.wing_rate_30_all = cat(1, Free.wing_rate_30{:});
Free_all.wing_rate_60_all = cat(1, Free.wing_rate_60{:});

%% Head fixed
root = 'E:\DATA\Rigid_Data\Saccade\sweep';
[FILE,PATH] = uigetfile({'*.mat'},'Select files', root, 'MultiSelect','on');
FILE = cellstr(FILE);

n_sweep = length(FILE);
Fixed = splitvars(table(nan(n_sweep,2)));
Fixed.Properties.VariableNames = {'stdThresh','ampCut'};
Fixed = [Fixed splitvars(table(num2cell(nan(n_sweep,2))))];
Fixed.Properties.VariableNames(3:end) = {'vel','wing_rate_30'};
for n = 1:n_sweep
    data = load(fullfile(PATH,FILE{n}),'wing','SACCADE','D','U','N');
    Fixed.stdThresh(n) = data.wing.thresh(3);
    Fixed.ampCut(n) = data.wing.amp_cut;
    
    clear Saccade
    %keepI = cellfun(@(x) isstruct(x) | isobject(x), data.SACCADE.head2wing);
    Saccade = data.SACCADE;
    Fixed.vel{n} = data.D.vel;    
    Fixed.wing_rate{n} = cellfun(@(x) x.rate, Saccade.wing_saccade);
    Fixed.wing_rate_30{n} = Fixed.wing_rate{n}(abs(Fixed.vel{n}) == 30);
end
Fixed = sortrows(Fixed, 'stdThresh');
Fixed = sortrows(Fixed, 'ampCut');

Fixed_all = [];
Fixed_all.stdThresh_all = cellfun(@(x,y) y*ones(size(x)), Fixed.wing_rate, num2cell(Fixed.stdThresh), ...
    'UniformOutput', false);
Fixed_all.stdThresh_all = cat(1,Fixed_all.stdThresh_all{:});

Fixed_all.ampCut_all = cellfun(@(x,y) y*ones(size(x)), Fixed.wing_rate, num2cell(Fixed.ampCut), ...
    'UniformOutput', false);
Fixed_all.ampCut_all = cat(1,Fixed_all.ampCut_all{:});

Fixed_all.wing_rate_30 = cat(1, Fixed.wing_rate_30{:});

%% Compare free & fixed
stim_speed = 30;
ampCut = unique(Free.ampCut);
stdThresh = unique(Free.stdThresh);
P = [Free(:,1:2) , table(nan(size(Free,1),1), 'VariableNames', {'pval'})];
for a = 1:length(ampCut)
   for s = 1:length(stdThresh)
       freeI = (Free.ampCut == ampCut(a)) & (Free.stdThresh == stdThresh(s));
       free_rate = Free.(['wing_rate_' num2str(stim_speed)]){freeI};
       fixedI = (Fixed.ampCut == ampCut(a)) & (Fixed.stdThresh == stdThresh(s));
       fixed_rate = Fixed.wing_rate_30{fixedI};
       
       Y = [free_rate ; fixed_rate];
       G = [ones(size(free_rate)) ; 2*ones(size(fixed_rate)) ];
       
       %[p,tbl,stats] = anova1(Y, G, 'off');
       [p,tbl,stats] = kruskalwallis(Y, G, 'off');
       
       P.pval(freeI) = p;
   end
end

%% Head-wing synchronization
fig = figure (1) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 4 3])
movegui(fig, 'center')
clear ax bx
ax(1) = subplot(2,1,1); cla ; hold on ; title('Head')
    ylabel('True positive')
    bx(1,:) = boxchart(Free_all.ampCut_all, Free_all.head_sync, ...
        'GroupByColor', Free_all.stdThresh_all);
ax(2) = subplot(2,1,2); cla ; hold on ; title('Wing')
    ylabel('False positive')
    bx(2,:) = boxchart(Free_all.ampCut_all, Free_all.wing_sync, ...
    'GroupByColor', Free_all.stdThresh_all);
 
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'YLim', [0 1.5], 'XColor', 'none')
set(bx, 'MarkerStyle', 'none')

%% Saccade rate free
fig = figure (2) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 4 3])
movegui(fig, 'center')
clear ax bx
ax(1) = subplot(2,1,1); cla ; hold on ; title('30')
    ylabel('Saccade frequency (Hz)')
    g = findgroups(Free_all.ampCut_wing_rate_30, Free_all.stdThresh_wing_rate_30);
    bx(1,:) = boxchart(g, Free_all.wing_rate_30_all, ...
        'GroupByColor', Free_all.ampCut_wing_rate_30);
    
ax(2) = subplot(2,1,2); cla ; hold on ; title('60')
    ylabel('Saccade frequency (Hz)')
    g = findgroups(Free_all.ampCut_wing_rate_60, Free_all.stdThresh_wing_rate_60);
    bx(2,:) = boxchart(g, Free_all.wing_rate_60_all, ...
    'GroupByColor', Free_all.ampCut_wing_rate_60);
 
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'YLim', [0 1.5], ...
    'XColor', 'none', 'XLim', [0 length(unique(g))+1])
set(bx, 'MarkerStyle', '.', 'BoxWidth', 1)

%% Saccade rate fixed
fig = figure (3) ; clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', 1*[2 2 4 1.5])
movegui(fig, 'center')
clear ax bx
ax(1) = subplot(1,1,1); cla ; hold on ; title('30')
    ylabel('Saccade frequency (Hz)')
    g = findgroups(Fixed_all.ampCut_all, Fixed_all.stdThresh_all);
    bx(1,:) = boxchart(g, Fixed_all.wing_rate_30, ...
        'GroupByColor', Fixed_all.ampCut_all);
 
set(ax, 'Color', 'none', 'LineWidth', 0.75, 'YLim', [0 1.5], ...
    'XColor', 'none', 'XLim', [0 length(unique(g))+1])
set(bx, 'MarkerStyle', 'none', 'BoxWidth', 1)

%% Save
savedir = 'E:\DATA\Rigid_Data\Saccade\sweep\combined';
filename = 'Ramp_wing_saccade_sweep_free_fixed.mat';
save(fullfile(savedir, filename), 'Free','Fixed','Free_all','Fixed_all','-v7.3')
end