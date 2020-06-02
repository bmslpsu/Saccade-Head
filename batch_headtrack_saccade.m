function [] = batch_headtrack_saccade(root, npoints, center_file, playback, showpoint)
%% batch_headtrack_saccade: runs head tracker for user selected video files
%
%   INPUT:
%       root        :   root directory
%       npoint      :   # of points for tracker
%       playback    :   playback rate (show a frame in increments of "playback")
%                       If false, then don't show anything (default = 1)
%       showpoint 	:  	logical >>> true = debug mode
%
%   OUTPUT:
%       -
%

showpoint = true;
npoints = 1;
playback = 5;
root = 'H:\EXPERIMENTS\RIGID\Experiment_Asymmetry_Control_Verification\HighContrast\30\Vid';
center_file = ['H:\EXPERIMENTS\RIGID\Experiment_Asymmetry_Control_Verification\HighContrast\30\Vid' ...
    '\tracked_head\Fly_2_Trial_1_Vel_-120_SpatFreq_30.mat'];

[FILES, PATH] = uigetfile({'*.mat', 'MAT-files'},'Select videos', root, 'MultiSelect','on');
FILES = string(FILES);
nfile = length(FILES);

headdir = fullfile(PATH,'tracked_head');

if ~isempty(center_file)
    load(center_file,'cPoint')
else
    cPoint = [];
end

for file = 1:nfile
    disp(FILES(file))
    disp('---------------------------------------')
    load(fullfile(PATH,FILES(file)),'vidData','t_v')
    [hAngles,cPoint,validity,ROI,initframe,finalframe] = headtracker(vidData, npoints, cPoint, playback, showpoint);
    
    figure
    imagesc(validity)
    
    pause
    
    close all
    
 	save(fullfile(headdir,FILES{file}),'-v7.3','hAngles','cPoint','validity',...
                                            'ROI','initframe','finalframe','t_v')
                                                                       
end
disp('ALL DONE')
end