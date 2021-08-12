function [MOV] = make_montage_rigid_head_pat(rootdir,rootpat,vidFs,export)
%% make_montage_rigid_head_pat: makes movie for fly in rigid tether
%
% 	Includes fly video, head tracking & pattern position
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

% Example Input
clear ; clc ; close all 
export = true;
vidFs = 50;
amp = 11.25;
rootdir = ['E:\EXPERIMENTS\RIGID\Experiment_Sinusoid\' num2str(amp)];
rootpat = 'Q:\OneDrive - PSU\OneDrive - The Pennsylvania State University\Git\Arena\Patterns';
pat_ypos = 5;

if ~isfolder(rootdir)
    dirflag = false;
    [rootdir,mainfile,mainext] = fileparts(rootdir);
else
    dirflag = true;
end

if ~isfolder(rootpat)
	[PATH.pat,FILE.pat,patext] = fileparts(rootpat);
    FILE.pat = [FILE.pat , patext];
else
	% Select pattern file
    [FILE.pat, PATH.pat] = uigetfile({'*.mat'}, ...
        'Select pattern file', rootpat, 'MultiSelect','off');
end

% Create data paths
PATH.raw            = rootdir;
PATH.vid            = fullfile(PATH.raw,'Vid');
PATH.head_track     = fullfile(PATH.raw,'Vid','tracked_head');

if dirflag
    % Select tracked angle file (use head tracked files to select)
    [FILE.raw, PATH.head_track] = uigetfile({'*.mat'}, ...
        'Select fly file', PATH.head_track, 'MultiSelect','off');
else
    FILE.raw = [mainfile , mainext];
end

% Create movie output directory
PATH.mov = fullfile(PATH.raw,'movie'); % movie directory
mkdir(PATH.mov) % create directory for export images

% Set file names
[~,FILE.basename,~] = fileparts(FILE.raw);
FILE.montage    = [FILE.basename '_Montage_Pat_new.mp4'];
FILE.mask       = [FILE.basename '.json'];

% Load data
disp('Loading Data ...')
pattern_data    = load(fullfile(PATH.pat,FILE.pat),'pattern'); % load pattern
raw_data        = load(fullfile(PATH.raw,FILE.raw),'data','t_p'); % load DAQ pattern positions
vid_data        = load(fullfile(PATH.vid,FILE.raw),'vidData'); % load raw video
head_data    	= load(fullfile(PATH.head_track,FILE.raw),'hAngles','cPoint'); % load head angles
disp('DONE')

%% Get pattern data & sync with start of visual stimulus
pattern_total_time = 10; % [s]
reg = true; % use interp times
start_idx = nan; % use first frame
add1 = true; % add 1st frame
debug = false;

[TRIG,PAT] = sync_pattern_trigger(raw_data.t_p, raw_data.data(:,2), pattern_total_time, ... 
                                raw_data.data(:,1), reg, start_idx, add1, debug);

%% Get kinematics data
FLY.time    = TRIG.time_sync; % video time
FLY.Fs      = round(1/mean(diff(FLY.time))); % video sampling rate
FLY.Fc      = 30; % cut off frequency for lpf
[b,a]       = butter(2,FLY.Fc/(FLY.Fs/2),'low'); % make lpf
FLY.head    = filtfilt(b,a,head_data.hAngles); % head angles [deg]
% FLY.head    = filtfilt(b,a,rad2deg(benifly_data.Head)); % head angles [deg]

% Normalize fly kinematics for experimental window
FLY.int_time    = TRIG.time_sync_exp; % video time
FLY.int_head    = FLY.head(TRIG.range);
PAT.norm        = 3.75*(PAT.pos_exp - mean(PAT.pos_exp));

PAT.Fs = 1 / mean(diff(PAT.time_sync));
PAT.Fc = 2*1.8;
% PAT.Fc = 15;
[PAT.b,PAT.a] = butter(3, PAT.Fc / (PAT.Fs/2),'low');
PAT.pos_filt = filtfilt(PAT.b, PAT.a, PAT.pos);
PAT.int_pos_filt = interp1(PAT.time_sync, PAT.pos_filt, TRIG.time_sync_exp, 'pchip');
PAT.int_norm_filt = 3.75*(PAT.int_pos_filt - mean(PAT.int_pos_filt));
rshift = 2*amp / range(PAT.int_norm_filt(500:1000));
PAT.int_norm_filt = rshift * PAT.int_norm_filt;
PAT.int_norm_filt = PAT.int_norm_filt - mean(PAT.int_norm_filt);

%% Get video data
FLY.raw = squeeze(vid_data.vidData(:,:,TRIG.range)); % raw video data
FLY.raw_crop = FLY.raw;
FLY.nframe = size(FLY.raw,3);

[FLY.raw_yP,FLY.raw_xP,~] = size(FLY.raw_crop); % get size of raw video
FLY.raw_center = [round(FLY.raw_xP/2) , 1.25*round(FLY.raw_yP/2)]; % center point for pattern & fly
radius = floor(max([FLY.raw_yP FLY.raw_xP])/1.0); % radius of pattern
thickness = 8; % radius display width

%% Heading
imshow(FLY.raw(:,:,1))
roi = drawpoint;

FLY.head_length = 40;
FLY.head_hinge = roi.Position;
% FLY.head_hinge = [params.gui.head.hinge.x  , params.gui.head.hinge.y];
FLY.head_tip = FLY.head_hinge + FLY.head_length*[sind(FLY.int_head) , -cosd(FLY.int_head)];

close

%% Make Movie
% Create structure to store frames
MOV(1:FLY.nframe) = struct('cdata', [], 'colormap',[]);

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
    'Position', [-16 2.5 14 5]);
% set(FIG, 'Visible','off');
linewidth = 1.25; % displayed line width
fontsize = 12;

gs = pattern_data.pattern.gs_val + 1;
cmap = [zeros(gs,1), linspace(0,1,gs)', zeros(gs,1)];
colormap(cmap)

clear ax
ax(1) = subplot(2,4,[1:2,5:6]) ; cla; hold on; axis square % for raw fly & pattern vid
ax(2) = subplot(2,4,[3:4,7:8]) ; cla; hold on
        ylabel('(°)','Color','w','FontSize',fontsize)
        h.pat = animatedline('Color','g','LineWidth',1); % for pattern angle
        h.pat.Color(4) = 0.5;
        h.head = animatedline('Color','c','LineWidth',linewidth); % for head angle

set(ax(2:end), 'FontSize', 12, 'Color', 'k', 'YColor', 'w', 'XColor', 'w', 'FontWeight', 'bold',...
    'LineWidth', 1.5,'XLim', [0 round(FLY.int_time(end))])
set(ax(end),'XTick', 0:2:round(FLY.time(end)))
set(ax(2), 'XTickLabel', [], 'XColor', 'none')
set(ax(2), 'YLim', 21*[-1 1], 'YTick', 15*[-1 0 1])

linkaxes(ax(2:end),'x')
align_Ylabels_ax(ax(2:end)')

iter = round(FLY.Fs/vidFs); % # of frames to skip per iteration to acheive desired frame rate
expframe = circshift(mod(1:FLY.nframe,iter)==1,0); % which frames to export
disp('Exporting Video...')
tic
pat_image = 255*pattern_data.pattern.Pats(1,:,1,pat_ypos);
for jj = 1:FLY.nframe % for each frame
    if expframe(jj) % if we want to display this frame
        % Get frames
        disp(jj)
        if jj >= iter
            win = jj-(iter-1):jj;
        else
            win = jj;
        end
        Frame.raw = 1.5*median(FLY.raw(:,:,win),3); % current raw frame median across frames
        
        % Display raw video
        subplot(2,4,[1:2,5:6]); cla; hold on; axis image
            imshow(Frame.raw)
            %surf(Frame.raw,  'linestyle', 'none')
          	plot([FLY.head_hinge(1), mean(FLY.head_tip(win,1))], ... % update line drawn to head
                 [FLY.head_hinge(2), mean(FLY.head_tip(win,2))], 'Color', 'c', 'LineWidth', 2.5)
                         
            
            % Make pattern ring
            ax_pat = axes; axis image
            set(ax_pat, 'Color', 'none', 'XColor', 'none', 'YColor', 'none', ...
                            'Position', ax(1) .Position)
           	cla
            
            pat_pos = 3.75*round(mean(PAT.pos_exp(win)));
            theta = -deg2rad(pat_pos) + -(pi/2 + linspace(0, 2*pi, size(pat_image,2)));
            x = radius * cos(theta) + FLY.raw_center(1);
            y = radius * sin(theta) + FLY.raw_center(2);
            z = zeros(1,length(x));
            hs = surface(ax_pat,[x;x],[y;y],[z;z],[pat_image;pat_image], ...
                'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', thickness);
    end
    
    % Head plot
    subplot(2,4,[3:4,7:8]); hold on
        addpoints(h.head, FLY.int_time(jj), FLY.int_head(jj))
        addpoints(h.pat, FLY.int_time(jj), PAT.int_norm_filt(jj))
        
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