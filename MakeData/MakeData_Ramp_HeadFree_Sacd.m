function [] = MakeData_Ramp_HeadFree_Sacd(wave,match,Fc)
%% MakeData_Ramp_HeadFree_Sacd:
%   INPUTS:
%       wave    :   spatial wavelength of data
%       match   :   what saccades direction to 
%       Fc      :   head data low-pass filter cutoff frequency [Hz]
%
%   OUTPUTS:
%       -
%
wave = 30;
match = 0;
Fc = 40;

% Data location
rootdir = ['H:\EXPERIMENTS\RIGID\Experiment_Asymmetry_Control_Verification\HighContrast\' num2str(wave)];

% What saccades we will get
if      match==0
    clss = 'All';
elseif  match==1
    clss = 'CO';
elseif  match==-1
    clss = 'Anti';
elseif  match==2
    clss = 'Positive';
elseif  match==-2
    clss = 'Negative';
else
    error('Invalid match condition')
end

% Output file name
filename = ['NewRamp_HeadFree_SACCD_' clss '_filt=' num2str(Fc) '_Wave=' num2str(wave)];

% Setup Directories 
PATH.daq = rootdir; % DAQ data location
PATH.vid = fullfile(PATH.daq,'\Vid'); % video data location
PATH.ang = fullfile(PATH.daq,'\Vid\tracked'); % tracked kinematic data location

% Select files
[D,I,N,U,T,FILES,~,basename] = GetFileData(PATH.ang,'*.csv',false,'fly','trial','vel','wave');

%% Get Data %%
% disp('Loading...')
clc
close all
showplot = false;
tintp = (0:(1/200):(10 - 1/200))';
Vel = U.vel{1};
Stim = (Vel*tintp')';
SACCADE = [I , table(num2cell(zeros(N.file,1)))]; % store saccade objects
SACCADE_STATS = []; % store saccade stats
for kk = 1:N.file
    disp(kk)   
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    % load(fullfile(PATH.vid, [basename{kk} '.mat']),'t_v'); % load head angles % time arrays
    benifly = ImportBenifly(fullfile(PATH.ang, FILES{kk}));
    
    % Sync video with trigger & pattern
    trig.raw_time   = t_p; % DAQ raw times for trigger
    trig.pos        = round(data(:,1)); % trigger values
    trig.diff       = diff(trig.pos); % trigger derivative (rising edge triggers frame)
    [~,trig.locs]   = findpeaks(trig.diff); % where each frame starts
    trig.time       = [0;trig.raw_time(trig.locs+1)]; % where each frame starts
   	
    % Get head data
    benifly.Head(1) = benifly.Head(2);
    Head = process_signal(trig.time, rad2deg(benifly.Head), Fc, tintp, [4 8 16 32 64]);
    
    % Get Saccade Stats
    peaks = [];
    direction = -sign(D.vel(kk)); % only get saccades in the opposite direction of visual motion
    head_saccade = saccade(Head.X(:,1), Head.Time, 3.5, direction, peaks,  showplot);
    head_saccade = stimSaccade(head_saccade, Stim(:,I.vel(kk)), showplot);
    
    SACCADE{kk,5} = {head_saccade};
    
    if isnan(head_saccade.count)
        rep = 1;
    else
        rep = head_saccade.count;
    end
    ITable = repmat(I(kk,:),rep,1);
    SACCADE_STATS = [SACCADE_STATS ; [ITable , head_saccade.SACD]];
    
end

%%

% test = cellfun(@(x) x.SACD, SACCADE.Var1, 'UniformOutput',false);
% test = cat(1,test{:});
% test = [SACCADE,test];

test = cellfun(@(x) x.normpeak_saccade, SACCADE.Var1, 'UniformOutput',false);

% test = cellfun(@(x) x.time, cellfun(@(x) x.normpeak_saccade, SACCADE.Var1, ...
%     'UniformOutput',false),'UniformOutput',false);

%% Normalize Head Saccades
%---------------------------------------------------------------------------------------------------------------------------------
varnames = {'Time','Position','Velocity','Position_Error','Velocity_Error',...
                'Position_IntError','Velocity_IntError','Stimulus_Position','Stimulus_Velocity'};

clear SACCADE
SACCADE.Head = cell(N{1,3},9);
SACCADE.cIdx = cell(N{1,3},1);
dR = cell(N{1,3},1);
center = 0;
dim = 1;
for jj = 1:N{1,3}
    [SACCADE.Head{jj,1},SACCADE.cIdx{jj},~,~,dR{jj}] = nancat_center(SACD.Saccade.Head{jj}(:,1), center, dim, [], []);
    for ww = 2:size(SACD.Saccade.Head{jj},2)
        for kk = 1:size(SACD.Saccade.Head{jj},1)
            for ii = 1:size(SACD.Interval.Head{jj}{kk,ww},2)
                SACCADE.Head{jj,ww}{kk,1}(:,ii) = cat_pad(SACD.Saccade.Head{jj}{kk,ww}(:,ii), dR{jj}{kk}(:,ii),nan);
            end
        end
    	SACCADE.Head{jj,ww} = cat(2,SACCADE.Head{jj,ww}{:});
    end
end
SACCADE.Head = cell2table(SACCADE.Head,'VariableNames',varnames);
SACCADE.HeadStats = cell2table(cellfun(@(x) MatStats(x,2), table2cell(SACCADE.Head),...
                            'UniformOutput',false),'VariableNames',varnames);
%%
clear INTERVAL
INTERVAL.Head = cell(N{1,3},9);
dR = cell(N{1,3},1);
center = 0;
dim = 1;
for jj = 1:N{1,3}
    [INTERVAL.Head{jj,1},~,~,~,dR{jj}] = nancat_center(SACD.Interval.Head{jj}(:,1), center, dim);
    for ww = 2:size(SACD.Interval.Head{jj},2)
        for kk = 1:size(SACD.Interval.Head{jj},1)
            for ii = 1:size(SACD.Interval.Head{jj}{kk,ww},2)
                INTERVAL.Head{jj,ww}{kk,1}(:,ii) = cat_pad(SACD.Interval.Head{jj}{kk,ww}(:,ii), dR{jj}{kk}(:,ii),nan);
            end
        end
    	INTERVAL.Head{jj,ww} = cat(2,INTERVAL.Head{jj,ww}{:});
    end
end
INTERVAL.Head = cell2table(INTERVAL.Head,'VariableNames',varnames);
INTERVAL.HeadStats = cell2table(cellfun(@(x) MatStats(x,2), table2cell(INTERVAL.Head),...
                            'UniformOutput',false),'VariableNames',varnames);
                        
%% SAVE %%
disp('Saving...')
% save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
%     'SACD','SACCADE','INTERVAL','Stim','TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')

save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
      'PATH','SACD','SACCADE','INTERVAL','Stim','D','I','U','N','T','-v7.3')

disp('SAVING DONE')
end