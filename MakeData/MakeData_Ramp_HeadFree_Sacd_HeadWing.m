function [] = MakeData_Ramp_HeadFree_Sacd_HeadWing(wave,Fc)
%% MakeData_Ramp_HeadFree_Sacd_HeadWing: 
%   INPUTS:
%       wave    :   spatial wavelength of data
%       Fc      :   head data low-pass filter cutoff frequency [Hz]
%
%   OUTPUTS:
%       -
%   

Fc = 30;
wave = 30;
rootdir = ['H:\EXPERIMENTS\Experiment_Asymmetry_Control_Verification\HighContrast\' num2str(wave)];

% Setup Directories
PATH.daq  = rootdir; % DAQ data location
PATH.vid  = fullfile(PATH.daq,'Vid'); % video data location
PATH.ang  = fullfile(PATH.vid,'tracked'); % tracked kinematic data location
PATH.sacd = fullfile(PATH.ang,'SACD','Head'); % extracted head saccades

% Select files
[D,I,N,U,T,FILES,~,basename] = GetFileData(PATH.ang,'*.csv',false,'fly','trial','vel','wave');

%% Get Data %%
disp('Loading...')
tintp = (0:(1/200):10 - (1/200))';
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
    load(fullfile(PATH.vid, [basename{kk} '.mat']),'t_v'); % load head angles % time arrays
    benifly = ImportBenifly(fullfile(PATH.ang, FILES{kk}));
    disp(basename{kk})
   	
    % Get head data   
    Head = process_signal(t_v, rad2deg(benifly.Head), [], Fc)
    
    Head = Fly(head.Pos,head.Time,head.Fc,[],tt); % head object
    if max(abs(Head.X(3:end,2)))>900
        badtrial{bad,1} = basename{kk};
        bad = bad + 1;
        % error('Head Trial Error')
    end
  	
    % Get wing data
    wing.Time       = t_v; % wing time [s]
    wing.Fs         = 1/mean(diff(wing.Time)); % sampling frequency [Hz]
    wing.Fc         = Fc; % cutoff frequency [Hz]
    [b,a]           = butter(2,wing.Fc/(wing.Fs/2)); % butterworth filter
	wing.Left       = rad2deg(hampel(wing.Time,benifly.LWing)); % left wing
    wing.Right      = rad2deg(hampel(wing.Time,benifly.RWing)); % right wing
 	wing.Left       = filtfilt(b,a,wing.Left); % left wing
    wing.Right      = filtfilt(b,a,wing.Right); % right wing
    
    wing.WBA        = wing.Left - wing.Right; % dWBA (L-R)
    Wing            = Fly(wing.WBA,wing.Time,wing.Fc,[],tt); % head object
	    
    [SACD.Head,~,~,~,~,~,~] = Sacd_Manual(Head.X(:,1),Head.Time,false);
    pause
 	[SACD.Wing,~,~,~,~,~,~] = Sacd_Manual(Wing.X(:,1),Wing.Time,false);
    pause
    close all
    
    % save(fullfile(PATH.sacd, [basename{kk} '.mat']), 'SACD')
end
disp('DONE')

%%

fig = figure (2) ; clf
set(fig,'Color','k','Units','Inches','Position',[2 2 8 3])
movegui(fig,'center')
ax(1) = subplot(1,1,1); hold on
set(ax,'Color','k','XColor','w','YColor','w','FontSize',12,'YLim',22*[-1 1])
h(1) = plot(Head.Time,Head.X(:,1) - mean(Head.X(:,1)),'c','LineWidth',1);
h(2) = plot(Head.Time,Wing.X(:,1) - mean(Wing.X(:,1)),'r','LineWidth',1);
xlabel('Time (s)')
ylabel('Angle (°)')
leg = legend(h,'Head','\Delta WBA');
leg.Box = 'off';
leg.TextColor = 'w';


end