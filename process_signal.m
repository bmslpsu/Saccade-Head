function [Signal] = process_signal(tt,xx,Fc,tintp,IOFreq)
%% process_signal: 
%   INPUTS:
%       tt      : time
%       xx      : raw signal
%       Fc      : cutoff frequency (optional)
%       tintp 	: interpolated time (optional)
%       IOFreq 	: discrete frequency inputs (optional)
%
%   OUTPUTS:
%       Signal	:   structure with signal attributes
%

% Get signal attributes
Signal.Time = tt(:);
Signal.X(:,1) = xx(:);
tintp = tintp(:);

% If interpolation is specified, replace signal & time with interpolated results
if ~isempty(tintp)
    Signal.isinterpolated = true;
    if min(tintp)<min(Signal.Time) || max(tintp)>max(Signal.Time)
       warning('Interpolation time outside of range')
    end
    
    Signal.X(:,1) = interp1(Signal.Time, Signal.X(:,1), tintp, ...
                            'linear','extrap'); % interpolate signal to new time
    Signal.Time = tintp(:);
else
    Signal.isinterpolated = false;
end

% Sampling rate
Signal.Ts = mean(diff(Signal.Time));
Signal.Fs = 1 ./ Signal.Ts;

% Filter signal if cut-off frequency if specified
if ~isempty(Fc)
    Signal.Fc = Fc;
    % LP Filter signal
    [b,a] = butter(2,Signal.Fc/(Signal.Fs/2),'low'); % 2nd-order low-pass butterworth filter
    Signal.X(:,1) = filtfilt(b,a,Signal.X(:,1)); % filtered signal
else
    Signal.Fc = [];
end

% Get 2 time derivatives of data
for drv = 1:2
    Signal.X(:,drv+1) = [0 ; diff(Signal.X(:,drv)) ./ Signal.Ts];
end

% Compute basic signal statistics
Signal.Mean         = mean(Signal.X,1);  	% mean: signal & derivatives
Signal.Median       = median(Signal.X,1); 	% median: signal & derivatives
Signal.STD          = std(Signal.X,[],1); 	% std: signal & derivatives
Signal.AbsMean      = mean(Signal.X,1);    	% mean of absolute value: signal & derivatives
Signal.AbsMedian    = median(Signal.X,1);  	% median of absolute value: signal & derivatives
Signal.AbsSTD       = std(Signal.X,[],1); 	% std of absolute value: signal & derivatives

% Transform signal into frequency domain
for drv = 1:size(Signal.X,2)
    [Signal.Fv, Signal.Mag(:,drv), Signal.Phase(:,drv), Signal.FREQ(:,drv)] ...
        = FFT(Signal.Time,Signal.X(:,drv));
    % Get IO frequencies if specified
    if ~isempty(IOFreq)
        [~,Signal.IOMag(:,drv),Signal.IOPhase(:,drv)] = ...
            getfreqpeaks(Signal.Fv,Signal.Mag(:,drv),Signal.Phase(:,drv),IOFreq,[],false);
    end
end

end