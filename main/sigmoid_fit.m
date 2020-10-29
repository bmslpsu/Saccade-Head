function [fitresult, gof] = sigmoid_fit(xfit, yfit, showplot)
% sigmoid_fit(xfit,yfit)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : xfit
%      Y Output: yfit
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%

if nargin < 3
    showplot = false;
end

%% Fit
[xData, yData] = prepareCurveData( xfit, yfit );

% Set up fittype and options.
ft = fittype( '1 / (1 + a*exp(b*x))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.291272365370147 0.443346056774894];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if showplot
    % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData );
    legend( h, 'yfit vs. xfit', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'xfit', 'Interpreter', 'none' );
    ylabel( 'yfit', 'Interpreter', 'none' );
    grid on
end

end


