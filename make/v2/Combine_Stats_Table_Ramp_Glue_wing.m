function [] = Combine_Stats_Table_Ramp_Glue_wing(root)
%% Combine_Stats_Table_Ramp_Glue_wing:
root = 'E:\DATA\Rigid_Data\Saccade';
cmb_dir = fullfile(root,'processed');
mkdir(cmb_dir)

[FILE,PATH] = uigetfile({'*.mat'}, root, 'MultiSelect','on');
FILE = cellstr(FILE);
n_file = length(FILE);

%% Get tables
T = cell(n_file,1);
for n = 1:n_file
    T{n} = load(fullfile(PATH,FILE{n}),'SACCADE','WING_SACCADE_STATS','U','N');
end

%% Combine tables & get position and velocity structures
All_Stats = [];
Count_Stats = [];
Data = [];
flyI = 0;
for n = 1:n_file
    tbl = T{n}.WING_SACCADE_STATS;
    tbl.wave = tbl.wave + n - 1;
    tbl.fly = tbl.fly + flyI;
    
    All_Stats = cat(1, All_Stats, tbl);
    
    keepI = cellfun(@(x) isstruct(x) || isobject(x), T{n}.SACCADE.wing_saccade);
    Saccade =  T{n}.SACCADE(keepI,:);
    
    count_table = Saccade(:,1:4);
    count_table.wave = count_table.wave + n - 1;
    count_table.fly = count_table.fly + flyI;
    
    count = table(cellfun(@(x) x.count, Saccade.wing_saccade),'VariableNames', {'count'});
    count_table = [count_table, count];
    Count_Stats = cat(1, Count_Stats, count_table);
    
    flyI = tbl.fly(end);
        
    n_speed = T{n}.N.vel/2;
    velI = findgroups(Saccade.vel);
    velI(velI > n_speed) = velI(velI > n_speed) - n_speed;
    pos = cellfun(@(x,y) x.position .* y, ...
        Saccade.wing_saccade, num2cell(sign(Saccade.vel)), 'UniformOutput', false);
    vel = cellfun(@(x,y) x.velocity .* y, ...
        Saccade.wing_saccade, num2cell(sign(Saccade.vel)), 'UniformOutput', false);
    pos = splitapply(@(x) {cat(2,x{:})}, pos, velI);
    vel = splitapply(@(x) {cat(2,x{:})}, vel, velI);
    Data(n).wave = T{n}.U.wave;
    Data(n).pos = pos;
    Data(n).vel = vel;
end

%% Save
fname = 'Ramp_Glue_All_Stats_wing';
save(fullfile(cmb_dir, [fname '.mat']), 'All_Stats', 'Count_Stats', 'Data');

end