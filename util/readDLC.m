function DLC_table = readDLC(filename, dataLines)
% readDLC: Import Deep-Lab-Cut data from .csv
%  
%   Returns the data as a table.
%

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [4, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 22);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["frame", "antenna_left_x", "antenna_left_y", "antenna_left_prob", "antenna_right_x", ...
                        "antenna_right_y", "antenna_right_prob", "neck_joint_x", "neck_joint_y", "neck_joint_prob", ...
                        "wing_joint_left_x", "wing_joint_left_y", "wing_joint_left_prob", "wing_joint_right_x", ...
                        "wing_joint_right_y", "wing_joint_right_prob", "front_leg_left_x", "front_leg_left_y", ...
                        "front_leg_left_prob", "front_leg_right_x", "front_leg_right_y", "front_leg_right_prob"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", ...
                        "double", "double", "double", "double", "double", "double", "double", "double", ...
                        "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
DLC_table = readtable(filename, opts);

end