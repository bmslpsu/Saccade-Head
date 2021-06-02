%% Head free
clear ; close all ; clc

% WING saccade detection parameters
Fs = 200;
wing.showplot = false;
wing.Fc_detect = [5 nan];
wing.Fc_ss = [5 nan];
wing.amp_cut = 2;
wing.dur_cut = inf;
wing.thresh = [0, 2, 1.25, 0];
wing.true_thresh = 100;
wing.sacd_length = nan;
wing.pks = [];
wing.min_pkdist = 0.3;
wing.min_pkwidth = 0.03;
wing.min_pkprom = 20;
wing.min_pkthresh = 0;
wing.boundThresh = 0.35;
wing.Fc = 20;
[wing.b, wing.a] = butter(3, wing.Fc / (Fs/2) ,'low');

std_thresh = [1 1.25 1.5 2];
n_sweep = length(std_thresh);
for n = 1:n_sweep
    wing.thresh(3) = std_thresh(n);
    tag = ['stdThresh=' num2str(wing.thresh(3)) '_ampCut=' num2str(wing.amp_cut)];
    Make_Ramp_HeadFree_head_wing_leg_sweep(wing, tag);
end

%% Head fixed
clear ; close all ; clc

% WING saccade detection parameters
Fs = 200;
wing.showplot = false;
wing.Fc_detect = [5 nan];
wing.Fc_ss = [5 nan];
wing.amp_cut = 2;
wing.dur_cut = inf;
wing.thresh = [0, 2, 1.25, 0];
wing.true_thresh = 100;
wing.sacd_length = nan;
wing.pks = [];
wing.min_pkdist = 0.3;
wing.min_pkwidth = 0.03;
wing.min_pkprom = 20;
wing.min_pkthresh = 0;
wing.boundThresh = 0.35;
wing.Fc = 20;
[wing.b, wing.a] = butter(3, wing.Fc / (Fs/2) ,'low');

std_thresh = [1 1.25 1.5 2];
n_sweep = length(std_thresh);
for n = 1:n_sweep
    wing.thresh(3) = std_thresh(n);
    tag = ['stdThresh=' num2str(wing.thresh(3)) '_ampCut=' num2str(wing.amp_cut)];
    Make_Ramp_HeadFixed_wing_sweep(wing, tag);
end