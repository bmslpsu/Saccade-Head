%% Space Time
clear ; close all ; clc

% Set image properties
wave = 30; % spatial wavelength [°]
sz = [1 1000]; % image size [pixels]
style = 2; % square wave
res = 1; % pixel scale

% Make image
[image_data, spat_res] = make_image(wave, sz, style, res, false);
init_image = image_data(:,:,1,1);
sz = size(init_image);

% Set eye spatial paramters
interommatidial_angle = 4.6; % interommatidial angle [deg]
k = 1.1;
% acceptance_angle = k*interommatidial_angle; % acceptance angle [°]
acceptance_angle = 5; % acceptance angle [°]
n_ommatidia = 40; % # of ommatidia (use default if empty)

% Construct EYE object
Eye = eye_model(interommatidial_angle, acceptance_angle, sz(2), n_ommatidia, false);

% Define image motion paramaters
Fs = 100; % temporal sampling frequency [Hz]
f = 2; % sinusoidal motion frequency [Hz]
A = 3.75; % sinusoidal motion amplitude [°]
phi = 0; % sinusoidal motion phase [rad]
T = 10; % simulation time [s]
t = (0:(1/Fs):T)'; % time vector [s]
x = A*sin(2*pi*f.*t + phi); % % 
n = length(x); % length of motion vector

vel = 60;
x = vel * t;

% Create image stream from motion vector, spatially filter each image with eye
motion_image = zeros(sz(1),sz(2),n); % raw images
motion_image_filt = zeros(sz(1),Eye.n_ommatidia,n); % images filtered by eye
for k = 1:n
    motion_image(:,:,k) = move_image(init_image, x(k));
    motion_image_filt(:,:,k) = EyeFilt(Eye, motion_image(:,:,k));
end

%% Space-Time
fig = figure;
set(fig, 'Color', 'w', 'Name','Space-Time')

ax(1) = subplot(4,2,[1,3,5]); hold on ; axis tight ; title('Raw Space-Time')
    imagesc((squeeze(motion_image)'))
    plot([Eye.eye_image(1) Eye.eye_image(1) Eye.eye_image(end) Eye.eye_image(end)], ...
        [ax(1).YLim(1) ax(1).YLim(2) ax(1).YLim(2) ax(1).YLim(1)], 'r', 'LineWidth', 2)
ax(2) = subplot(4,2,[2,4,6]); hold on ; axis tight ; title('Filtered Space-Time')
    imagesc((squeeze(motion_image_filt)'))
ax(3) = subplot(4,2,7); hold on ; axis tight ; title('Raw Frame')
    imagesc(motion_image(:,:,1))
    plot([Eye.eye_image(1) Eye.eye_image(1) Eye.eye_image(end) Eye.eye_image(end)], ...
        [ax(3).YLim(1) ax(3).YLim(2) ax(3).YLim(2) ax(3).YLim(1)], 'r', 'LineWidth', 2)
ax(4) = subplot(4,2,8); hold on ; axis tight ; title('Filtered Frame')
    imagesc(motion_image_filt(:,:,1))

linkaxes(ax([1,3]),'x')
linkaxes(ax([2,4]),'x')
set(ax([3,4]),'YTick',[])
set(ax([1,2]),'YTick', unique(sort([ax(1).YTick , 0 ])))
ax(1).YTickLabels = ax(1).YTick / Fs;
ax(2).YTickLabels = ax(2).YTick / Fs;

XLabelHC = get(ax(3), 'XLabel');
set(XLabelHC, 'String', 'Space (pixel)')
XLabelHC = get(ax(4), 'XLabel');
set(XLabelHC, 'String', 'Receptor (pixel)')

YLabelHC = get(ax(1:2), 'YLabel');
set([YLabelHC{:}], 'String', 'Time (s)')

nc = 100;
cmap = [zeros(nc,1), linspace(0,1,nc)', zeros(nc,1)];
colormap(cmap)

%% Motion
fig = figure;
set(fig, 'Color', 'w', 'Units', 'inches')
fig.Position(3:4) = [6 5];
movegui(fig, 'center')
colormap(parula)
xtickmarks = -180:90:180;
clear ax
ax(1) = subplot(2,1,1); hold on ; axis tight
    xaxes = (xtickmarks ./ spat_res);
    xaxes = abs(xaxes(1)) + xaxes;
    xaxes(1) = 1;
    xticks(xaxes)
    xticklabels(xtickmarks)
    xlabel('(°)')
ax(2) = subplot(2,1,2); hold on ; axis tight
    xaxes = unique([1,10:10:Eye.n_ommatidia, Eye.n_ommatidia]);
    xticks(xaxes)
    xlabel('Ommatidia')
set(ax, 'FontWeight', 'bold', 'YTick', [], 'box', 'on')
    
for k = 1:n
    subplot(2,1,1) ; cla
        imagesc(motion_image(:,:,k))
        plot([Eye.eye_image(1) Eye.eye_image(1) Eye.eye_image(end) Eye.eye_image(end)], ...
            [ax(1).YLim(1) ax(1).YLim(2) ax(1).YLim(2) ax(1).YLim(1)], 'r', 'LineWidth', 3)
    subplot(2,1,2) ; cla
        imagesc(motion_image_filt(:,:,k))
    pause(1/Fs)
end





%%
r = 1;
center = [0 0];
theta = pi/4 + linspace(0,2*pi,sz(2));
x = r * cos(theta) + center(1);
y = r * sin(theta) + center(2);
% z = 10*ones(size(x));
z = linspace(0,10,length(x));
col = init_image;  % This is the color, vary with x in this case.

fig = figure;
set(fig, 'Color', 'w', 'Units', 'inches')
fig.Position(3:4) = [6 5];
movegui(fig, 'center')
ax(1) = subplot(1,1,1); hold on
axis(1.2*[-r r -r r])
% surface([x;x],[y;y],[z;z],[col;col], 'FaceColor','none','EdgeColor','interp', 'LineWidth', 10);
surface([x;x],[y;y],[z;z],[col;col], 'FaceColor','none','EdgeColor','flat', 'LineWidth', 10);
axis off
rotate3d on

nc = 10;
cmap = [zeros(nc,1), linspace(0,1,nc)', zeros(nc,1)];
colormap(cmap)














