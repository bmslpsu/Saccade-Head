function [MOV] = make_montage_rigid_head_wing_fixed(root_free,rootpat,vidFs,export)
%% make_montage_rigid_head_wing_fixed: makes movie for fly in rigid tether
%
% 	Includes fly video, head tracking, wing tracking, leg tracking 
%   & pattern position
%
%   INPUT:
%       rootdir     : directory containing .mat DAQ & VIDEO files
%       rootpat     : directory containing PATTERN file
%       rootleg     : directory containing DLC tracked leg files (.csv)
%       vidFs       : video display FPS
%       export      : boolean (1=export video to images)
%
%   OUTPUT:
%       MOV         : structure containing movie 
%

% Example Input %
clear ; clc ; close all
export = true;
vidFs = 50;
root_free = 'H:\EXPERIMENTS\RIGID\Experiment_Asymmetry_Control_Verification\HighContrast\30';
root_fixed = 'H:\EXPERIMENTS\RIGID\Experiment_Ramp_30_HeadFixed';
rootpat = 'C:\Users\BC\Box\Git\Arena\Patterns';
pat_ypos = 5;

% Select pattern file
[FILE.pat, PATH.pat] = uigetfile({'*.mat'}, ...
    'Select pattern file', rootpat, 'MultiSelect','off');

% Create data paths
PATH.raw_free      	= root_free;
PATH.vid_free      	= fullfile(PATH.raw_free,'Vid');
PATH.head_free      = fullfile(PATH.vid_free,'tracked_head');
PATH.beninfly_free	= fullfile(PATH.vid_free,'wing_filt','tracked_head_wing');
PATH.raw_fixed    	= root_fixed;
PATH.vid_fixed    	= fullfile(PATH.raw_fixed,'');
PATH.beninfly_fixed	= fullfile(PATH.vid_fixed,'wing_filt','tracked_wing');

% Select tracked angle file (use head tracked files to select)
[FILE.raw_free, PATH.beninfly_free] = uigetfile({'*.csv'}, ...
    'Select fly file', PATH.beninfly_free, 'MultiSelect','off');

[FILE.raw_fixed, PATH.beninfly_fixed] = uigetfile({'*.csv'}, ...
    'Select fly file', PATH.beninfly_fixed, 'MultiSelect','off');

% Create movie output directory
PATH.mov = fullfile(PATH.raw_fixed,'movie'); % movie directory
mkdir(PATH.mov) % create directory for export images

% Set file names
[~,FILE.basename_free,~] = fileparts(FILE.raw_free);
[~,FILE.basename_fixed,~] = fileparts(FILE.raw_fixed);
FILE.mat_free = [FILE.basename_free '.mat'];
FILE.mat_fixed = [FILE.basename_fixed '.mat'];
FILE.benifly_free = [FILE.basename_free '.csv'];
FILE.benifly_fixed = [FILE.basename_fixed '.csv'];
FILE.mask_free 	= [FILE.basename_free '.json'];
FILE.mask_fixed = [FILE.basename_fixed '.json'];

FILE.montage  = ['Free=' FILE.basename_free '_Fxed=' FILE.basename_fixed '_Montage.mp4'];

% Load data
disp('Loading Data ...')
pattern_data    = load(fullfile(PATH.pat,FILE.pat),'pattern'); % load pattern

benifly_free    = ImportBenifly(fullfile(PATH.beninfly_free,FILE.benifly_free)); % load Benifly tracked kinematics
raw_free        = load(fullfile(PATH.raw_free,FILE.mat_free),'data','t_p'); % load DAQ pattern positions
vid_free        = load(fullfile(PATH.vid_free,FILE.mat_free),'vidData'); % load raw video
head_free    	= load(fullfile(PATH.head_free,FILE.mat_free),'hAngles','cPoint'); % load head angles

benifly_fixed    = ImportBenifly(fullfile(PATH.beninfly_fixed,FILE.benifly_fixed)); % load Benifly tracked kinematics
raw_fixed        = load(fullfile(PATH.raw_fixed,FILE.mat_fixed),'data','t_p','vidData'); % load DAQ pattern position & raw video

disp('DONE')

%% Get pattern data & sync with start of visual stimulus

pattern_total_time = 10; % [s]
reg = true; % use interp times
debug = false;

[TRIG_free,PAT_free] = sync_pattern_trigger(raw_free.t_p, raw_free.data(:,2), pattern_total_time, ... 
                                raw_free.data(:,1), reg, nan, true, false);

[TRIG_fixed,PAT_fixed] = sync_pattern_trigger(raw_fixed.t_p, raw_fixed.data(:,2), pattern_total_time, ... 
                                raw_fixed.data(:,1), reg, [], false, false);
                            
%% Get kinematic data
FLY.Fc = 20; % cut off frequency for lpf

% Free
FLY.time_free 	= TRIG_free.time_sync; % video time
FLY.Fs_free   	= round(1/mean(diff(FLY.time_free))); % video sampling rate
[b_free,a_free]	= butter(2,FLY.Fc/(FLY.Fs_free/2),'low'); % make lpf
FLY.head_free  	= filtfilt(b_free,a_free,head_free.hAngles); % head angles [°]
FLY.lwing_free	= rad2deg(hampel(FLY.time_free,benifly_free.LWing)); % left wing angles [°]
FLY.rwing_free 	= rad2deg(hampel(FLY.time_free,benifly_free.RWing)); % right wing angles [°]
FLY.lwing_free 	= filtfilt(b_free,a_free,FLY.lwing_free); % left wing angles [°]
FLY.rwing_free 	= filtfilt(b_free,a_free,FLY.rwing_free); % right wing angles [°]
FLY.wba_free   	= FLY.lwing_free - FLY.rwing_free; % delta wing-beat-amplitude [°]
FLY.wba_free  	= filtfilt(b_free,a_free,FLY.wba_free); % delta wing-beat-amplitude [°]
% FLY.wba_free 	= FLY.wba_free - mean(FLY.wba_free); % delta wing-beat-amplitude [°]

% Fixed
FLY.time_fixed      = TRIG_fixed.time_sync; % video time
FLY.Fs_fixed        = round(1/mean(diff(FLY.time_fixed))); % video sampling rate
[b_fixed,a_fixed]   = butter(2,FLY.Fc/(FLY.Fs_fixed/2),'low'); % make lpf
FLY.lwing_fixed     = rad2deg(hampel(FLY.time_fixed,benifly_fixed.LWing)); % left wing angles [°]
FLY.rwing_fixed 	= rad2deg(hampel(FLY.time_fixed,benifly_fixed.RWing)); % right wing angles [°]
FLY.lwing_fixed 	= filtfilt(b_fixed,a_fixed,FLY.lwing_fixed); % left wing angles [°]
FLY.rwing_fixed 	= filtfilt(b_fixed,a_fixed,FLY.rwing_fixed); % right wing angles [°]
FLY.wba_fixed   	= FLY.lwing_fixed - FLY.rwing_fixed; % delta wing-beat-amplitude [°]
FLY.wba_fixed       = filtfilt(b_fixed,a_fixed,FLY.wba_fixed); % delta wing-beat-amplitude [°]
% FLY.wba_fixed       = FLY.wba_fixed - mean(FLY.wba_fixed); % delta wing-beat-amplitude [°]

% Normalize fly kinematics for experimental window
FLY.int_time_free       = TRIG_free.time_sync_exp;
FLY.int_head_free       = FLY.head_free(TRIG_free.range);
FLY.int_lwing_free      = FLY.lwing_free(TRIG_free.range);
FLY.int_rwing_free      = FLY.rwing_free(TRIG_free.range);
FLY.int_wba_free        = FLY.wba_free(TRIG_free.range);

FLY.int_time_fixed      = TRIG_fixed.time_sync_exp;
FLY.int_time_fixed      = interp1(1:length(FLY.int_time_fixed),FLY.int_time_fixed, ...
                                1:0.5:length(FLY.int_time_fixed)-0.5)';

FLY.int_lwing_fixed 	= FLY.lwing_fixed(TRIG_fixed.range);
FLY.int_rwing_fixed     = FLY.rwing_fixed(TRIG_fixed.range);
FLY.int_wba_fixed   	= FLY.wba_fixed(TRIG_fixed.range);

FLY.int_lwing_fixed 	= interp1(TRIG_fixed.time_sync_exp,FLY.int_lwing_fixed,FLY.int_time_fixed);
FLY.int_rwing_fixed 	= interp1(TRIG_fixed.time_sync_exp,FLY.int_rwing_fixed,FLY.int_time_fixed);
FLY.int_wba_fixed       = interp1(TRIG_fixed.time_sync_exp,FLY.int_wba_fixed,FLY.int_time_fixed);

%% Get video data
FLY.raw_free = squeeze(vid_free.vidData(:,:,TRIG_free.range)); % raw free video data
FLY.nframe_free = size(FLY.raw_free,3);
FLY.raw_fixed = squeeze(raw_fixed.vidData(:,:,TRIG_fixed.range)); % raw fixed video data
FLY.nframe_fixed = size(FLY.raw_fixed,3);

FLY.raw_fixed_norm = FLY.raw_free;
shft = 2;
pp = 0;
for n = 1:FLY.nframe_fixed-1
    FLY.raw_fixed_norm(:,:,pp+(n:n+shft-1)) = repmat(FLY.raw_fixed(:,:,n),1,1,shft);
    pp = pp + shft-1;
end

[FLY.raw_yP,FLY.raw_xP,~] = size(FLY.raw_free); % get size of raw video
FLY.raw_center = [round(FLY.raw_xP/2) , 1.25*round(FLY.raw_yP/2)]; % center point for pattern & fly
radius = floor(max([FLY.raw_yP FLY.raw_xP])/2); % radius of pattern
thickness = 6; % radius display width

%% Get benifly parameters/mask
% Heading
fid = fopen(fullfile(PATH.beninfly_free, FILE.mask_free));
raw = fread(fid,inf);
str = char(raw'); 
fclose(fid); 
params_free = jsondecode(str);

fid = fopen(fullfile(PATH.beninfly_fixed, FILE.mask_fixed));
raw = fread(fid,inf);
str = char(raw'); 
fclose(fid); 
params_fixed = jsondecode(str);

FLY.wing_length = 150;

FLY.body_free = 90 + rad2deg(atan2(params_free.gui.head.hinge.y -params_free.gui.abdomen.hinge.y, ...
                            params_free.gui.head.hinge.x - params_free.gui.abdomen.hinge.x));
FLY.body_fixed = 90 + rad2deg(atan2(params_fixed.gui.head.hinge.y -params_fixed.gui.abdomen.hinge.y, ...
                            params_fixed.gui.head.hinge.x - params_fixed.gui.abdomen.hinge.x));
                        
FLY.lwing_hinge_free = [params_free.gui.left.hinge.x  , params_free.gui.left.hinge.y];
FLY.rwing_hinge_free = [params_free.gui.right.hinge.x , params_free.gui.right.hinge.y];
FLY.lwing_tip_free = FLY.lwing_hinge_free - FLY.wing_length*[cosd(FLY.int_lwing_free + FLY.body_free),  ...
                                                             sind(FLY.int_lwing_free + FLY.body_free)];
FLY.rwing_tip_free = FLY.rwing_hinge_free + FLY.wing_length*[cosd(FLY.int_rwing_free + FLY.body_free), ...
                                                            -sind(FLY.int_rwing_free + FLY.body_free)];

FLY.head_length = 40;
% FLY.head_hinge = [head_data.cPoint.X  , head_data.cPoint.Y];
FLY.head_hinge = [params_free.gui.head.hinge.x  , params_free.gui.head.hinge.y];
FLY.head_tip = FLY.head_hinge + FLY.head_length*[sind(FLY.int_head_free) , -cosd(FLY.int_head_free)];

FLY.lwing_hinge_fixed = [params_fixed.gui.left.hinge.x  , params_fixed.gui.left.hinge.y];
FLY.rwing_hinge_fixed = [params_fixed.gui.right.hinge.x , params_fixed.gui.right.hinge.y];
FLY.lwing_tip_fixed = FLY.lwing_hinge_fixed - FLY.wing_length*[cosd(FLY.int_lwing_fixed +  + FLY.body_fixed), ...
                                                               sind(FLY.int_lwing_fixed +  + FLY.body_fixed)];
FLY.rwing_tip_fixed = FLY.rwing_hinge_fixed + FLY.wing_length*[cosd(FLY.int_rwing_fixed +  + FLY.body_fixed), ...
                                                              -sind(FLY.int_rwing_fixed +  + FLY.body_fixed)];

%% Make Movie
% Create structure to store frames
MOV(1:FLY.nframe_free) = struct('cdata', [], 'colormap',[]);

% Create video object
if export
    % VID = VideoWriter(fullfile(PATH.mov,FILE.montage),'Uncompressed AVI');
    VID = VideoWriter(fullfile(PATH.mov,FILE.montage),'MPEG-4');
    % VID.LosslessCompression = true;
    VID.FrameRate = vidFs;
    open(VID)
end

FIG = figure (1); clf % main figure window for display & export
set(FIG, 'Color', 'k', 'Renderer', 'OpenGL', 'Units', 'inches', ...
    'Position', [-16 2.5 14 9]);
% set(FIG, 'Visible','off');
movegui(FIG, 'center')
linewidth = 1.25; % displayed line width
fontsize = 12;

gs = pattern_data.pattern.gs_val + 1;
cmap = [zeros(gs,1), linspace(0,1,gs)', zeros(gs,1)];
colormap(cmap)

clear ax
ax(1) = subplot(2,4,[1:2]) ; cla; hold on; axis square % for raw fly & pattern vid free
ax(2) = subplot(2,4,[5:6]) ; cla; hold on; axis square % for raw fly & pattern vid fixed
ax(3) = subplot(2,4,[3:4]) ; cla; hold on
        ylabel('Free Head & \DeltaWBA (°)','Color','w','FontSize',fontsize)
        h.head = animatedline('Color','c','LineWidth',linewidth); % for head angle
        h.wba_free = animatedline('Color','r','LineWidth',linewidth); % for dWBA angle free
ax(4) = subplot(2,4,[7:8]) ; cla; hold on
        ylabel('Fixed \DeltaWBA (°)','Color','w','FontSize',fontsize)
    	xlabel('Time (s)','Color','w','FontSize',fontsize)
        h.wba_fixed = animatedline('Color','m','LineWidth',linewidth); % for dWBA angle fixed

set(ax(3:end), 'FontSize', 12, 'Color', 'k', 'YColor', 'w', 'XColor', 'w', 'FontWeight', 'bold',...
    'LineWidth', 1.5,'XLim', [0 round(FLY.int_time_free(end))])
set(ax(end),'XTick', 0:2:round(FLY.time_free(end)))
set(ax(3), 'XTickLabel', [], 'XColor', 'none')

max_y = 5*ceil(max(abs([FLY.int_wba_free;FLY.int_wba_fixed]))/5);

set(ax(3:4), 'YLim', max_y*[-1 1], 'YTick', 20*[-1 0 1])

linkaxes(ax(3:end),'xy')
align_Ylabels_ax(ax(3:end)')

iter = round(FLY.Fs_free/vidFs); % # of frames to skip per iteration to acheive desired frame rate
expframe = circshift(mod(1:FLY.nframe_free,iter)==1,0); % which frames to export
disp('Exporting Video...')
tic
pat_image = 255*pattern_data.pattern.Pats(1,:,1,pat_ypos);
for jj = 1:FLY.nframe_free % for each frame
    if expframe(jj) % if we want to display this frame`
        % Get frames
        disp(jj)
        if jj >= iter
            win = jj-(iter-1):jj;
        else
            win = jj;
        end
        Frame.free  = median(FLY.raw_free(:,:,win),iter-1); % current raw frame median across frames
        Frame.fixed = 1.1*median(FLY.raw_fixed_norm(:,:,win),iter-1); % current raw frame median across frames
        
        % Display free video
        subplot(2,4,[1:2]); cla; hold on; axis image
            imshow(Frame.free)
          	plot([FLY.head_hinge(1), mean(FLY.head_tip(win,1))], ... % update line drawn to head
                 [FLY.head_hinge(2), mean(FLY.head_tip(win,2))], 'Color', 'c', 'LineWidth', 2.5)
                         
          	plot([FLY.lwing_hinge_free(1) , mean(FLY.lwing_tip_free(win,1))], ... % update line drawn to left wing
                [FLY.lwing_hinge_free(2) , mean(FLY.lwing_tip_free(win,2))], 'Color', 'r', 'LineWidth', 2.5)
            
        	plot([FLY.rwing_hinge_free(1) , mean(FLY.rwing_tip_free(win,1))], ... % update line drawn to right wing
                [FLY.rwing_hinge_free(2) , mean(FLY.rwing_tip_free(win,2))], 'Color', 'r', 'LineWidth', 2.5)
            
            plot(FLY.lwing_hinge_free(1), FLY.lwing_hinge_free(2), 'r.', 'MarkerSize',20)
            plot(FLY.rwing_hinge_free(1), FLY.rwing_hinge_free(2), 'r.', 'MarkerSize',20)
            
            % Make pattern ring
            ax_pat = axes; axis image
            set(ax_pat, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', ...
                            'Position', ax(1).Position)
           	cla
            
            pat_pos = 3.75*round(mean(PAT_free.pos_exp(win)));
            theta = -deg2rad(pat_pos) + -(pi/2 + linspace(0, 2*pi, size(pat_image,2)));
            x = radius * cos(theta) + FLY.raw_center(1);
            y = radius * sin(theta) + FLY.raw_center(2);
            z = zeros(1,length(x));
            hs = surface(ax_pat,[x;x],[y;y],[z;z],[pat_image;pat_image], ...
                'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', thickness);
            ax(1).XLim = [-100 600];
            ax(1).YLim = [-100 350];
            
        % Display raw video
        subplot(2,4,[5:6]); cla; hold on; axis image
            imshow(Frame.fixed)                         
          	plot([FLY.lwing_hinge_fixed(1) , mean(FLY.lwing_tip_fixed(win,1))], ... % update line drawn to left wing
                [FLY.lwing_hinge_fixed(2) , mean(FLY.lwing_tip_fixed(win,2))], 'Color', 'm', 'LineWidth', 2.5)
            
        	plot([FLY.rwing_hinge_fixed(1) , mean(FLY.rwing_tip_fixed(win,1))], ... % update line drawn to right wing
                [FLY.rwing_hinge_fixed(2) , mean(FLY.rwing_tip_fixed(win,2))], 'Color', 'm', 'LineWidth', 2.5)
            
            plot(FLY.lwing_hinge_fixed(1), FLY.lwing_hinge_fixed(2), 'm.', 'MarkerSize',20)
            plot(FLY.rwing_hinge_fixed(1), FLY.rwing_hinge_fixed(2), 'm.', 'MarkerSize',20)
            
            % Make pattern ring
            ax_pat = axes; axis image
            set(ax_pat, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', ...
                            'Position', ax(2).Position)
           	cla
            
            pat_pos = 3.75*round(mean(PAT_free.pos_exp(win)));
            theta = -deg2rad(pat_pos) + -(pi/2 + linspace(0, 2*pi, size(pat_image,2)));
            x = radius * cos(theta) + FLY.raw_center(1);
            y = radius * sin(theta) + FLY.raw_center(2);
            z = zeros(1,length(x));
            hs = surface(ax_pat,[x;x],[y;y],[z;z],[pat_image;pat_image], ...
                'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', thickness);
            ax(2).XLim = [-100 600];
            ax(2).YLim = [-100 350];
    end
    
    % Head plot & WBA plot free
    subplot(2,4,[3:4]); hold on
        addpoints(h.head, FLY.int_time_free(jj), FLY.int_head_free(jj))
        addpoints(h.wba_free, FLY.int_time_free(jj), FLY.int_wba_free(jj))
        
    % WBA plot fixed
    subplot(2,4,[7:8]); hold on
        addpoints(h.wba_fixed, FLY.int_time_fixed(jj), FLY.int_wba_fixed(jj))
        
    drawnow
    
    if export
        if expframe(jj)
            fig_frame = getframe(FIG);
         	writeVideo(VID,fig_frame);
        end
    end
    pause(0.001)
end
toc

if export
 	disp('Saving...')
    pause(1)
    close(VID) % close video
end
disp('DONE')
end