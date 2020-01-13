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
Signal.t = tt(:);
Signal.x(:,1) = xx(:);
tintp = tintp(:);

% If interpolation is specified, replace signal & time with interpolated
% results
if ~isempty(tintp)
    Signal.isinterpolated = true;
    if min(tintp)<min(Signal.t) || max(tintp)>max(Signal.t)
       warning('Interpolation time outside of range')
    end
    
    Signal.x(:,1) = interp1(Signal.t, Signal.x(:,1), tintp, ...
                            'linear','extrap'); % interpolate signal to new time
    Signal.t = tintp(:);
else
    Signal.isinterpolated = false;
end

% Sampling rate
Signal.Ts = mean(diff(Signal.t));
Signal.Fs = 1 ./ Signal.Ts;

% Filter signal if cut-off frequency if specified
if ~isempty(Fc)
    Signal.Fc = Fc;
    % LP Filter signal
    [b,a] = butter(2,Signal.Fc/(Signal.Fs/2),'low'); % 2nd-order low-pass butterworth filter
    Signal.x(:,1) = filtfilt(b,a,Signal.x(:,1)); % filtered signal
else
    Signal.Fc = [];
end

% Get 2 time derivatives of data
for drv = 1:2
    Signal.x(:,drv+1) = [0 ; diff(Signal.x(:,1)) ./ Signal.Ts];
end

% Compute basic signal statistics
Signal.mean         = mean(Signal.x,2);        	% mean: signal & derivatives
Signal.median       = median(Signal.x,2);   	% median: signal & derivatives
Signal.std          = std(Signal.x,[],2);   	% std: signal & derivatives
Signal.absmean      = mean(Signal.x,2);        	% mean of absolute value: signal & derivatives
Signal.absmedian    = median(Signal.x,2);   	% median of absolute value: signal & derivatives
Signal.absstd       = std(Signal.x,[],2);   	% std of absolute value: signal & derivatives

% Transform signal into frequency domain
for drv = 1:size(Signal.x,2)
    [Signal.Fv, Signal.Mag(:,drv), Signal.Phase(:,drv), Signal.FREQ(:,drv)] ...
        = FFT(Signal.t,Signal.x(:,drv));
end

if ~isempty(IOFreq)


end