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
filename = ['Ramp_HeadFree_SACCD_' clss '_filt=' num2str(Fc) '_Wave=' num2str(wave)];

% Setup Directories 
PATH.daq = rootdir; % DAQ data location
PATH.vid = fullfile(PATH.daq,'\Vid'); % video data location
PATH.ang = fullfile(PATH.daq,'\Vid\tracked'); % tracked kinematic data location

% Select files
[D,I,N,U,T,FILES,~,basename] = GetFileData(PATH.ang,'*.csv',false,'fly','trial','vel','wave');

%% Get Data %%
% disp('Loading...')
% TRIAL = cell(N{1,1},N{1,3});

% clear SACD
% SACD.Head = [];
% SACD.Wing = [];
% SACD.Saccade.Head = cell(N{1,3},1);
% SACD.Interval.Head = cell(N{1,3},1);
% SACD.Stimulus.Saccade.Head = cell(N{1,3},1);
% SACD.Stimulus.Interval.Head = cell(N{1,3},1);
% badtrial = {};
clc
close all
Vel = U.vel{1};
tintp = (0:(1/200):(10 - 1/200))';
Stim = (Vel*tintp')';
% bad = 1;
for kk = 1:N.file
    % disp(kk)
    % disp(basename{kk})
    
    % Load HEAD & DAQ data
	load(fullfile(PATH.daq, [basename{kk} '.mat']),'data','t_p'); % load head angles % time arrays
    load(fullfile(PATH.vid, [basename{kk} '.mat']),'t_v'); % load head angles % time arrays
    benifly = ImportBenifly(fullfile(PATH.ang, FILES{kk}));
    
    % Sync video with trigger & pattern
    Trig.raw_time   = t_p; % DAQ raw times for trigger
    Trig.pos        = round(data(:,1)); % trigger values
    Trig.diff       = diff(Trig.pos); % trigger derivative (rising edge triggers frame)
    [~,Trig.locs] = findpeaks(Trig.diff); % where each frame starts
    Trig.time       = [0;Trig.raw_time(Trig.locs+1)]; % where each frame starts
   	
    % Get head data
    benifly.Head(1) = benifly.Head(2);
    Head = process_signal(Trig.time, rad2deg(benifly.Head), Fc, tintp, [4 8 16 32 64]);
 	test = saccade(Head.X(:,1),Head.Time,3.5,true);

    % Get Saccade Stats   
    [head.SACD,head.thresh,head.count,head.rate,head.SACDRmv] = sacddetect(Head.X(:,1),Head.Time,400,true);
    
    % HeadRmv = Fly(head.SACDRmv,Head.Time,[],[],tt);
    % WingRmv = Fly(wing.SACDRmv,Wing.Time,[],[],tt);
    
    head.match = table(head.SACD.Direction*sign(D.vel(kk)));
    mIdx = 1:head.count;
    if ~isnan(match)
        mIdx = mIdx([head.match{:,1}]==match);
    end
    head.match.Properties.VariableNames = {'Match'};
  	wing.match = table(wing.SACD.Direction*sign(D.vel(kk)));
    wing.match.Properties.VariableNames = {'Match'};
    
    if isnan(head.count)
        head.Rate = table(0);
        wing.Rate = table(0);
    else
        head.Rate = table(nan(head.count,1));
        head.Rate{1,1} = head.rate;
    	wing.Rate = table(nan(wing.count,1));
        wing.Rate{1,1} = wing.rate;
    end
    head.Rate.Properties.VariableNames = {'Rate'};
  	wing.Rate.Properties.VariableNames = {'Rate'};
    
    head.SACD = [head.SACD , head.match];
    wing.SACD = [wing.SACD , wing.match];
     
   	Dir = table(sign(D{kk,3}),'VariableNames',{'Dir'});
    I_table = [I(kk,1:3) , rowfun(@(x) abs(x), D(kk,3)), Dir];
    I_table.Properties.VariableNames{3} = 'velIdx';
    I_table.Properties.VariableNames{4} = 'speed';
    
    [Saccade,Interval,Stimulus,Error,IntError,matchFlag] = SaccdInter(Head.X(:,1),Head.Time,head.SACD, ...
                                                                    match, Stim(:,I{kk,3}), false);
    
    var1 = {Saccade.Time, Saccade.Pos,Saccade.Vel, Error.Saccade.Pos, Error.Saccade.Vel,...
                IntError.Saccade.Pos, IntError.Saccade.Vel, Stimulus.Saccade.Pos , Stimulus.Saccade.Vel};
            
    var2 = {Interval.Time, Interval.Pos, Interval.Vel, Error.Interval.Pos, Error.Interval.Vel,...
                IntError.Interval.Pos, IntError.Interval.Vel, Stimulus.Interval.Pos , Stimulus.Interval.Vel};
    
    if isnan(head.count)
        Err_table = nan(1,5);
        loop = [];
        emptyFlag = true;
    else
        Err_table = nan(head.count,5);
        loop = size(Error.Interval.Pos,2);
       	if matchFlag
            emptyFlag = true;
        else
            emptyFlag = false;
        end
    end
    
    for jj = 1:loop
        pos_err = Error.Interval.Pos(:,jj);
        pos_err = pos_err(~isnan(pos_err));
     	vel_err = Error.Interval.Vel(:,jj);
        vel_err = vel_err(~isnan(vel_err));
        
        pos_int_err = IntError.Interval.Pos(:,jj);
        pos_int_err = pos_int_err(~isnan(pos_int_err));
      	vel_int_err = IntError.Interval.Vel(:,jj);
        vel_int_err = vel_int_err(~isnan(vel_int_err));
        
        stim_pos = Stimulus.Interval.Pos(:,jj);
     	stim_pos = stim_pos(~isnan(stim_pos));
        if ~isempty(pos_err)
            Err_table(mIdx(jj),1) = nanmean(pos_err(end-5:end),1);
            Err_table(mIdx(jj),2) = nanmean(vel_err(end-5:end),1);
            Err_table(mIdx(jj),3) = pos_int_err(end);
            Err_table(mIdx(jj),4) = vel_int_err(end);
            Err_table(mIdx(jj),5) = stim_pos(end);
        end
    end
    Err_table = splitvars(table(Err_table));
    Err_table.Properties.VariableNames = {'Position_Error','Velocity_Error','Position_IntError',...
                                                'Velocity_IntError','Stimulus_Position'};
    if isnan(head.count)
        head.I_table = I_table;
    else
        head.I_table = repmat(I_table,head.count,1);
        head.Rate{1,1} = head.Rate{1,1}*(size(Error.Interval.Pos,2)/head.count);
  	end
    head.SACD = [head.SACD , head.Rate];
    
    if isnan(wing.count)
        wing.I_table = I_table;
    else
        wing.I_table = repmat(I_table,wing.count,1);
    end
    
	head.SACD = [head.I_table , [head.SACD , Err_table]];
  	wing.SACD = [wing.I_table , wing.SACD];
    
    SACD.Head = [SACD.Head ; head.SACD];
    SACD.Wing = [SACD.Wing ; wing.SACD];
    
    if ~emptyFlag
        SACD.Saccade.Head{I{kk,3},1}  = [SACD.Saccade.Head{I{kk,3},1}  ; var1];
        SACD.Interval.Head{I{kk,3},1} = [SACD.Interval.Head{I{kk,3},1} ; var2];
    end
    
    % pause()
    close all
    % clc
end

clear jj ii kk pp qq ww n a b  t_v hAngles data head wing pat tt I_table Dir loop Saccade Interval Stimulus Error IntError...
    Head Pat Wing  vars PATH t_p var1 var2 Err_table pos_err vel_err pos_int_err vel_int_err stim_pos matchFlag emptyFlag

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
%---------------------------------------------------------------------------------------------------------------------------------
disp('Saving...')
% save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
%     'SACD','SACCADE','INTERVAL','Stim','TRIAL','FLY','GRAND','D','I','U','N','T','-v7.3')

save(['H:\DATA\Rigid_Data\' filename '_' datestr(now,'mm-dd-yyyy') '.mat'],...
    'SACD','SACCADE','INTERVAL','Stim','D','I','U','N','T','-v7.3')

disp('SAVING DONE')
end