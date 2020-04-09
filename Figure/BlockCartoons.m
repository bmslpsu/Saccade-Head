function [] = BlockCartoons()
%% BlockCartoons:

Fs = 1000;
t = 0:(1/Fs):1;
x = sin(2*pi*1*t);
y = sin(2*pi*10*t);

fig = figure(1); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [2 2 1 2])

ax(1) = subplot(2,1,1); hold on
    plot(t, x, 'b', 'LineWidth', 1)
    plot([0 1], [0 0], 'k', 'LineWidth', 0.5)

ax(2) = subplot(2,1,2); hold on
    plot(t, y, 'r', 'LineWidth', 1)
    plot([0 1], -0.5*[1 1], 'k', 'LineWidth', 0.5)
    
set(ax,'Visible','off')

end