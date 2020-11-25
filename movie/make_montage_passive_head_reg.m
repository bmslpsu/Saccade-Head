function [MOV] = make_montage_passive_head_reg(rootdir,vidFs,export)
%% make_montage_passive_head_reg: makes movie for fly in rigid tether
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

% Example Input
clear ; clc ; close all 
export = true;
vidFs = 50;
rootdir = 'H:\EXPERIMENTS\RIGID\Passive Head Displacement\Vid\registered';
if ~isfolder(rootdir)
    dirflag = false;
    [rootdir,mainfile,mainext] = fileparts(rootdir);
else
    dirflag = true;
end

% Create data paths (static)
PATH.vid = rootdir;
PATH.head_track	= fullfile(PATH.vid, 'tracked_head');

if dirflag
    % Select tracked angle file (use head tracked files to select)
    [FILE.vid, PATH.head_track] = uigetfile({'*.mat'}, ...
        'Select fly file', PATH.head_track, 'MultiSelect','off');
else
    FILE.vid = [mainfile , mainext];
end

% Create movie output directory
PATH.mov = fullfile(PATH.vid,'movie'); % movie directory
mkdir(PATH.mov) % create directory for export images

% Set file names
[~,FILE.basename,~] = fileparts(FILE.vid);
FILE.montage = [FILE.basename '_Montage.mp4'];

% Load data
disp('Loading Data ...')
vid_data = load(fullfile(PATH.vid,FILE.vid),'regvid'); % load vid video
head_data = load(fullfile(PATH.head_track,FILE.vid),'hAngles','cPoint','t_v'); % load head angles
disp('DONE')

%% Get kinematics data
FLY.time    = head_data.t_v; % video time
FLY.Fs      = round(1/mean(diff(FLY.time))); % video sampling rate
FLY.Fc      = 40; % cut off frequency for lpf
[b,a]       = butter(2,FLY.Fc/(FLY.Fs/2),'low'); % make lpf
FLY.head    = filtfilt(b,a,head_data.hAngles); % head angles [°]

% Find pulse window
head.showplot = false;
head.Fc_detect = [10 nan];
head.Fc_ss = [nan nan];
head.amp_cut = 4;
head.dur_cut = inf;
head.thresh = [60 , 0, 0, 0];
head.true_thresh = 100;
head.direction = 0;
head.sacd_length = nan;
head.pks = [];
head.min_pkdist = 0.2;
head.min_pkwidth = 0.03;
head.min_pkprom = 75;
head.min_pkthresh = 0;
head.boundThresh = [0 50];

% Extract head pulse
head_saccade = saccade_all(FLY.head, FLY.time, head.thresh, head.true_thresh, head.Fc_detect, ...
                            head.Fc_ss, head.amp_cut, head.dur_cut, head.direction, head.pks, ...
                            head.sacd_length, head.min_pkdist, head.min_pkwidth, head.min_pkprom, ...
                            head.min_pkthresh, head.boundThresh, head.showplot);
t_win = [0.2 0.5]; % [before after]
scdI = head_saccade.saccades{1}.Index;
winI = scdI(1) - round(t_win(1)*FLY.Fs) : scdI(end) + ( round(t_win(2)*FLY.Fs) - length(scdI) );
% FLY.time_win = (1/FLY.Fs)*(0:length(winI)-1)';
FLY.time_win = (-t_win(1):(1/FLY.Fs):t_win(2))';
FLY.head_win = FLY.head(winI);

norm_win = scdI(end):length(FLY.head);
norm_head = median(FLY.head(norm_win));
FLY.head_win_norm = FLY.head_win - FLY.head_win(end);

%% Get video data
FLY.vid = squeeze(vid_data.regvid(:,:,winI)); % video data
FLY.vid = medfilt3(FLY.vid, [3 3 3]);
% [~,rectout] = imcrop(FLY.vid(:,:,1));
crop_r = [0.1 0.1];
yy = nan(1,2);
yy(1) = ceil(crop_r(1)*size(FLY.vid,1));
yy(2) = size(FLY.vid,1);

xx = nan(1,2);
xx(1) = ceil(crop_r(2)*size(FLY.vid,2));
xx(2) = size(FLY.vid,2) - xx(1);

FLY.vid = FLY.vid(yy(1):yy(2),xx(1):xx(2),:);

FLY.nframe = size(FLY.vid,3);
[FLY.vid_yP,FLY.vid_xP,~] = size(FLY.vid); % get size of vid video
FLY.vid_center = [round(FLY.vid_xP/2) , 1.25*round(FLY.vid_yP/2)]; % center point for pattern & fly
FLY.head_length = 40;
FLY.head_hinge = -[xx(1) yy(1)] + [head_data.cPoint.X  , head_data.cPoint.Y];
FLY.head_tip = FLY.head_hinge + FLY.head_length*[sind(FLY.head_win) , -cosd(FLY.head_win)];

%% Make Movie
% Create structure to store frames
MOV(1:FLY.nframe) = struct('cdata', [], 'colormap',[]);

% Create video object
if export
    VID = VideoWriter(fullfile(PATH.mov,FILE.montage),'MPEG-4');
    VID.FrameRate = vidFs;
    open(VID)
end

FIG = figure (1); clf % main figure window for display & export
set(FIG, 'Color', 'k', 'Renderer', 'OpenGL', 'Units', 'inches', ...
    'Position', [-16 2.5 5 7]);
linewidth = 1.25; % displayed line width
fontsize = 12;
clear ax
ax(1) = subplot(3,4,[1:8]) ; cla; hold on; axis square % for vid fly 
ax(2) = subplot(3,4,[9:12]) ; cla; hold on
        ylabel('Head (°)','Color','w','FontSize',fontsize)
        

set(ax(2:end), 'FontSize', 12, 'Color', 'k', 'YColor', 'w', 'XColor', 'w', 'FontWeight', 'bold',...
    'LineWidth', 1.5,'XLim', [-0.02 - t_win(1) , t_win(2)])
set(ax(end),'XTick', -t_win(1):0.1:t_win(2))
% set(ax(2), 'XTickLabel', [], 'XColor', 'none')
ylim_start = 10*floor(median(FLY.head_win(1:20))/10);
set(ax(2), 'YLim', [ylim_start 5], 'YTick', [ylim_start:10:0])
xlabel('Time (s)')

% linkaxes(ax(2:end),'x')
% align_Ylabels_ax(ax(2:end)')
n_rep = 3;
iter = round(FLY.Fs/vidFs); % # of frames to skip per iteration to acheive desired frame rate
expframe = circshift(mod(1:FLY.nframe,iter)==1,0); % which frames to export
disp('Exporting Video...')
tic
for n = 1:n_rep
    subplot(3,4,[9:12]) ; cla
    h.head = animatedline('Color','c','LineWidth',linewidth); % for head angle
    for jj = 1:FLY.nframe % for each frame
        if expframe(jj) % if we want to display this frame
            % Get frames
            disp(jj)
            if jj >= iter
                win = jj-(iter-1):jj;
            else
                win = jj;
            end
            Frame.vid = 0.7*median(FLY.vid(:,:,win),3); % current vid frame median across frames

            % Display vid video
            subplot(3,4,[1:8]); cla; hold on; axis image
                title('Real time', 'Color', 'w', 'FontSize', 14)
                imshow(Frame.vid)
                plot([FLY.head_hinge(1), mean(FLY.head_tip(win,1))], ... % update line dvidn to head
                     [FLY.head_hinge(2), mean(FLY.head_tip(win,2))], 'Color', 'c', 'LineWidth', 2.5)  
        end

        % Head plot
        subplot(3,4,[9:12]); hold on
            addpoints(h.head, FLY.time_win(jj), FLY.head_win_norm(jj))

        drawnow

        if export
            if expframe(jj)
                fig_frame = getframe(FIG);
                writeVideo(VID,fig_frame);
            end
        end
        pause(0.001)
    end
end

% Add slow motion
subplot(3,4,[9:12]) ; cla
h.head = animatedline('Color','c','LineWidth',linewidth); % for head angle
for jj = 1:FLY.nframe % for each frame
	Frame.vid = 0.7*median(FLY.vid(:,:,jj),3); % current vid frame median across frames

	% Display vid video
    subplot(3,4,[1:8]); cla; hold on; axis image
        imshow(Frame.vid)
        plot([FLY.head_hinge(1), FLY.head_tip(jj,1)], ... % update line dvidn to head
             [FLY.head_hinge(2), FLY.head_tip(jj,2)], 'Color', 'c', 'LineWidth', 2.5)  
       	title('Slowed down 6.2x', 'Color', 'w', 'FontSize', 14)
        
    % Head plot
    subplot(3,4,[9:12]); hold on
        addpoints(h.head, FLY.time_win(jj), FLY.head_win_norm(jj))

    drawnow

    if export
        fig_frame = getframe(FIG);
        writeVideo(VID,fig_frame);
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