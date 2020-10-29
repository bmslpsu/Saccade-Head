function [fitresult, gof] = test_fit(tt, xx, showplot)
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
ft = fittype( 'a*exp(-b*x) + c*exp(-d*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Lower = [-inf 1 -inf 1];
opts.Upper = [inf inf inf inf];
opts.Display = 'Off';
opts.StartPoint = [1 -1 1 -1];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if showplot
    % Plot fit with data.
    %figure( 'Name', 'Fit' );
    cla ; hold on
    h = plot( fitresult, xData, yData );
    legend( h, 'xx vs. tt', 'Fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
    xlabel( 'tt', 'Interpreter', 'none' );
    ylabel( 'xx', 'Interpreter', 'none' );
    grid on
end

end
