function [] = Combine_Stats_Table()
%% Combine_Stats_Table:

root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','on');
FILE = cellstr(FILE);
n_file = length(FILE);

%% Get tables
T = cell(n_file,1);
for n = 1:n_file
    T{n} = load(fullfile(PATH,FILE{n}),'SACCADE','HEAD_SACCADE_STATS','U','N');
end

%% Combine tables
All_Stats = [];
Count_Stats = [];
flyI = 0;
for n = 1:n_file
    tbl = T{n}.HEAD_SACCADE_STATS;
    tbl.wave = tbl.wave + n - 1;
    tbl.fly = tbl.fly + flyI;
    
    All_Stats = cat(1, All_Stats, tbl);
    
    count_table = T{n}.SACCADE(:,1:4);
    count_table.wave = count_table.wave + n - 1;
    count_table.fly = count_table.fly + flyI;
    
    count = table(cellfun(@(x) x.count, T{n}.SACCADE.head_saccade),'VariableNames', {'count'});
    count_table = [count_table, count];
    Count_Stats = cat(1, Count_Stats, count_table);
    
    flyI = tbl.fly(end);
end

%% Save
fname = 'Ramp_All_Stats';
savedir = 'C:\Users\BC\Box\Research\Manuscripts\Head Saccade\Data';
save(fullfile(savedir, [fname '.mat']), 'All_Stats', 'Count_Stats');

end