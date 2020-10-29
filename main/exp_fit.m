function [fitresult, gof] = exp_fit(tt, xx, showplot)
% exp_fit (tt, xx)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : tt
%      Y Output: xx
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%

if nargin < 3
    showplot = false;
end

%% Fit
[xData, yData] = prepareCurveData( tt, xx );

% Set up fittype and options.
ft = fittype( 'a*exp(-b*x)+c', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.StartPoint = [30 30 0];
opts.StartPoint = [-10 30 0];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if showplot
    % Plot fit with data.
    figure( 'Name', 'Fit' );
    h = plot( fitresult, xData, yData );
    legend( h, 'xx vs. tt', 'Fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
    xlabel( 'tt', 'Interpreter', 'none' );
    ylabel( 'xx', 'Interpreter', 'none' );
    grid on
end

end
