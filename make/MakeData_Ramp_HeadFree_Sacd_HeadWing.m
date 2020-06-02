function [] = MakeData_Ramp_HeadFree_Sacd_HeadWing(wave,Fc)
%% MakeData_Ramp_HeadFree_Sacd_HeadWing: 
%   INPUTS:
%       wave    :   spatial wavelength of data
%       Fc      :   head data low-pass filter cutoff frequency [Hz]
%
%   OUTPUTS:
%       -
%   

% Fc = 30;
% wave = 30;
rootdir = ['H:\EXPERIMENTS\Experiment_Asymmetry_Control_Verification\HighContrast\' num2str(wave)];

% Setup Directories
PATH.daq  = rootdir; % DAQ data location
PATH.vid  = fullfile(PATH.daq,'Vid'); % video data location
PATH.ang  = fullfile(PATH.vid,'tracked'); % tracked kinematic data location
PATH.sacd = fullfile(PATH.ang,'SACD','Head'); % extracted head saccades

% Select files
[D,I,N,U,T,FILES,~,basename] = GetFileData(PATH.ang,'*.csv',false,'fly','trial','vel','wave');

% Get Data %
disp('Loading...')
% tintp = (0:(1/200):10 - (1/200))';
for kk = 1:N.file
    disp(kk)
    % Load HEAD & DAQ data
    load(fullfile(PATH.vid, [basename{kk} '.mat']),'t_v'); % load head angles % time arrays
    benifly = ImportBenifly(fullfile(PATH.ang, FILES{kk}));
    disp(basename{kk})
   	
    % Get head data
    Head = process_signal(t_v, rad2deg(benifly.Head), Fc, [], [8 16 32]);
	    
    [SACD.Head,~,~,~,~,~,~] = sacd_detect_manual(Head.X(:,1),Head.Time,false);
    pause
    close all
    save(fullfile(PATH.sacd, [basename{kk} '.mat']), 'SACD')
end
disp('DONE')

%%

fig = figure (2) ; clf
set(fig,'Color','k','Units','Inches','Position',[2 2 8 3])
movegui(fig,'center')
ax(1) = subplot(1,1,1); hold on
set(ax,'Color','k','XColor','w','YColor','w','FontSize',12,'YLim',22*[-1 1])
h(1) = plot(Head.Time,Head.X(:,1) - mean(Head.X(:,1)),'c','LineWidth',1);
% h(2) = plot(Head.Time,Wing.X(:,1) - mean(Wing.X(:,1)),'r','LineWidth',1);
xlabel('Time (s)')
ylabel('Angle (°)')
% leg = legend(h,'Head','\Delta WBA');
% leg.Box = 'off';
% leg.TextColor = 'w';


end