function [] = Process_Saccade_Sine_anti(showplot)
%% Process_Saccade_Sine_anti:
root = 'H:\DATA\Rigid_Data\Saccade';
[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');
load(fullfile(PATH,FILE),'SACCADE','HEAD_SACCADE_STATS','U','N','I')

filedata = textscan(FILE, '%s', 'delimiter', '_=');
amp = str2double(filedata{1}{5});

%%
clearvars -except amp U N I SACCADE showplot

showplot = false;
Fs = round(SACCADE.head_saccade{1}.Fs);
Freq = U.freq{1};
% vel_win = round(0.5 *Fs);
I_new = [I(:,1:2), table(repmat(amp, [N.file, 1]), 'VariableNames', {'amp'}) , I(:,3)];
ALL = [I_new , splitvars(table(num2cell(zeros(N.file,2))))];
ALL.Properties.VariableNames(5:end) = {'pos', 'head_saccade'};
ALL = [ALL , splitvars(table(zeros(N.file,5)))];
ALL.Properties.VariableNames(7:end) = {'anti_count', 'co_count','OR','mag','start'};

if showplot
    fig = figure (1);
    set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [3 4 13 5])
end
for n = 1:N.file
    time = SACCADE.head_saccade{n}.time;
    pos = SACCADE.head_saccade{n}.position;
    vel = SACCADE.head_saccade{n}.velocity;
    stim_pos = SACCADE.head_saccade{n}.stimlus_position;
    stim_vel = SACCADE.head_saccade{n}.stimlus_velocity;
    
    ALL.pos{n} = pos;
    
    %amp = SACCADE.amp(n);
    freq = SACCADE.freq(n);
    vel_win = round(0.3 * (1/freq) * Fs);
    stim_mean_vel = 4*amp*freq;
    
    pos_norm = pos - median(pos);
    ALL.OR(n) = 2*median(abs(pos_norm));
    
    [Fv, Mag , Phs , FREQ] = FFT(time, pos_norm);
    [~,freqI] = min(abs(Fv - U.freq{1}(freq)));
    search_win = freqI - 3 : freqI + 3;
    mag_peak = max(Mag(search_win));
    ALL.mag(n) = 2*mag_peak;
    
    n_scd = SACCADE.head_saccade{n}.count;
    if n_scd > 0
        all_tbl = splitvars(table(nan(n_scd, 2)));
        all_tbl.Properties.VariableNames = {'vel_pre','co_anti'};
        stim_tbl = splitvars(table( repmat([SACCADE.fly(n) amp freq stim_mean_vel], [n_scd 1])));
        stim_tbl.Properties.VariableNames = {'fly','amp','freq','stim_mean_vel'};
        startI = SACCADE.head_saccade{n}.SACD.StartIdx;
        stim_win = cell(n_scd,1);
        for s = 1:n_scd
            winI = ( startI(s) - vel_win ) : startI(s);
            winI(winI < 1) = [];
            winI(winI > SACCADE.head_saccade{n}.n) = [];
            vel_pre = median(stim_vel(winI));
            co_anti = sign(vel_pre) * SACCADE.head_saccade{n}.SACD.Direction(s);
            all_tbl.vel_pre(s) = vel_pre;
            all_tbl.co_anti(s) = co_anti;
            stim_win{s} = winI;
        end
        all_tbl = [stim_tbl, all_tbl , SACCADE.head_saccade{n}.SACD(:,[1:3,10,12,14]) ];
        all_tbl{:,[8,10:12]} = all_tbl{:,[8,10:12]} .* repmat(all_tbl.Direction,[1,4]);
        anti_count = sum( all_tbl.co_anti == -1 );
        co_count = sum( all_tbl.co_anti == 1 );
        start_fly = mean(all_tbl.StartPos);
    else
        all_tbl = [];
        anti_count = 0;
        co_count = 0;
        start_fly = nan;
    end
    ALL.head_saccade{n} = all_tbl;
    ALL.anti_count(n) = anti_count;
	ALL.co_count(n) = co_count;
    ALL.start(n) = start_fly;
    
    if showplot
    subplot(2,1,1) ; cla ; hold on
        plot(time, stim_pos, 'Color', [0.5 0.5 0.5 0.7], 'LineWidth', 1)
        plot(time, pos, 'k', 'LineWidth', 1)
        if n_scd > 0
            for s = 1:n_scd
                if all_tbl.co_anti(s) == -1
                    cc = 'r';
                else
                    cc = 'b';
                end
                plot(time(stim_win{s}), stim_pos(stim_win{s}), 'g', 'LineWidth', 1)
                plot(SACCADE.head_saccade{n}.saccades{s}.Time, ...
                    SACCADE.head_saccade{n}.saccades{s}.Position, 'Color', cc, 'LineWidth', 1.5)
            end
        end
        
    subplot(2,1,2) ; cla ; hold on
        plot(time, stim_vel, 'Color', [0.5 0.5 0.5 0.7], 'LineWidth', 1)
        plot(time, vel, 'k', 'LineWidth', 1)
        if n_scd > 0
            for s = 1:n_scd
                if all_tbl.co_anti(s) == -1
                    cc = 'r';
                else
                    cc = 'b';
                end
                plot(time(stim_win{s}), stim_vel(stim_win{s}), 'g', 'LineWidth', 1)
                plot(SACCADE.head_saccade{n}.saccades{s}.Time, ...
                    SACCADE.head_saccade{n}.saccades{s}.Velocity, 'Color', cc, 'LineWidth', 1.5)
            end
        end
	pause
    end
end

All_scd = cat(1, ALL.head_saccade{:});
Anti_scd = All_scd(All_scd.co_anti == -1,:);
Co_scd = All_scd(All_scd.co_anti == 1,:);

%% Combine data for each fly/frequency combo
clear Head

[freq_fly_group, freq_group, fly_group] = findgroups(ALL.freq, ALL.fly);
Head.pos = splitapply(@(x) {cat(2,x{:})}, ALL.pos, freq_fly_group);
Head.OR = splitapply(@(x) {x'}, ALL.OR, freq_fly_group);
Head.mag = splitapply(@(x) {x'}, ALL.mag, freq_fly_group);
Head.anti_count = splitapply(@(x) {x'}, ALL.anti_count, freq_fly_group);
Head.co_count = splitapply(@(x) {x'}, ALL.anti_count, freq_fly_group);

Head = structfun(@(x) splitapply(@(y) {y}, x, freq_group), Head, 'UniformOutput', false);
Head = structfun(@(x) cat(2,x{:}), Head, 'UniformOutput', false);
Head.pos = cellfun(@(x) x - median(x(:)), Head.pos, 'UniformOutput', false);

Head.fly_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    Head, 'UniformOutput', false);

fnames = string(fieldnames(Head.fly_stats));
n_field = length(fnames);
for f = 1:n_field
    for v = 1:N.freq
        Head.fly_mean.(fnames(f)){v} = cat(2, Head.fly_stats.(fnames(f))(:,v).mean);
    end
end

Head.vel_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    Head.fly_mean, 'UniformOutput', false);

%% Head saccades starts and ends
clc
clear Head_scd
Head_scd.start = cell(N.fly, N.freq);
Head_scd.end = cell(N.fly, N.freq);
for n = 1:size(Anti_scd, 1)
    v = Anti_scd.freq(n);
    f = Anti_scd.fly(n);
    Head_scd.start{f,v}(1,end+1) = Anti_scd.StartPos(n);
 	Head_scd.end{f,v}(1,end+1) = Anti_scd.EndPos(n);
end

Head_scd.fly_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    Head_scd, 'UniformOutput', false);

fnames = string(fieldnames(Head_scd.fly_stats));
n_field = length(fnames);
for f = 1:n_field
    for v = 1:N.freq
        Head_scd.fly_mean.(fnames(f)){v} = nancat(2, Head_scd.fly_stats.(fnames(f))(:,v).mean);
        Head_scd.fly_mean.(fnames(f)){v}(isnan( Head_scd.fly_mean.(fnames(f)){v})) = [];
    end
end

Head_scd.vel_stats = structfun(@(x) cellfun(@(y) basic_stats(y,2), x, 'UniformOutput', true), ...
    Head_scd.fly_mean, 'UniformOutput', false);

%%  Store fly data in table
Fly = ALL(1,:); % initialiaze fly table
Fly.trial = [];
Fly.head_saccade = [];
warning('off', 'MATLAB:table:RowsAddedExistingVars');
pp = 1;
for n = 1:N.fly
    for v = 1:N.freq
        Fly.fly(pp) = n;
        Fly.amp(pp) = amp;
        Fly.freq(pp) = Freq(v);
        Fly.stim_vel(pp) = 4 * amp * Freq(v);
        Fly.pos(pp) = Head.pos(n,v);
        Fly.OR(pp) = Head.fly_stats.OR(n,v).mean;
        Fly.mag(pp) = Head.fly_stats.mag(n,v).mean;
        Fly.anti_count(pp) = Head.fly_stats.anti_count(n,v).mean;
        Fly.co_count(pp) = Head.fly_stats.co_count(n,v).mean;
        Fly.co_count(pp) = Head.fly_stats.co_count(n,v).mean;
        
        start = Head_scd.fly_stats.start(n,v).mean;
        if ~isempty(start)
            Fly.start(pp) = start;
            Fly.end(pp) = Head_scd.fly_stats.end(n,v).mean;
        else
            Fly.start(pp) = nan;
        end
        
        pp = pp + 1;
    end
end
Fly = movevars(Fly, 'stim_vel', 'After', 'freq');

%% Position Histogram
fig = figure (1) ; clf
set(fig, 'Color', 'w','Units', 'inches', 'Position', [2 2 5 7])
ax = gobjects(N.freq,1);
bins = -25:2:25;
pp = 1;
for v = 1:N.freq
    ax(v) = subplot(N.freq/2,2,pp); hold on ; title([ num2str(U.freq{1}(v)) '(Hz)'])
        h = histogram(cat(2,Head.pos{:,v}), bins, 'Normalization', 'probability', ...
            'FaceColor', 'k', 'FaceAlpha', 1, 'EdgeColor', 'none');
        h = histogram(cat(2,Head_scd.start{:,v}), bins, 'Normalization', 'probability', ...
            'FaceColor', 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
        xlabel('Position (°)')
        ylabel('Probability')
        axis tight
	pp = pp + 1;
end
linkaxes(ax,'xy')
set(ax, 'LineWidth', 1, 'Box', 'off', 'XLim', 25*[-1 1], 'XTick', -25:5:25, 'YLim', [-0.005 0.28])

%% Save
fname = ['Sine_amp=' num2str(amp)];
savedir = 'H:\DATA\Rigid_Data\Saccade\processed';
save(fullfile(savedir, [fname '.mat']), 'amp', 'U', 'N', 'I', 'ALL', 'All_scd', 'Anti_scd', ...
                                            'Head', 'Head_scd', 'Fly');

end