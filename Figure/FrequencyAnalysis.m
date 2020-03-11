function [] = FrequencyAnalysis()
%% FrequencyAnalysis:
root = 'H:\DATA\Rigid_Data\';

[FILE,PATH] = uigetfile({'*.mat', 'DAQ-files'}, ...
    'Select head angle trials', root, 'MultiSelect','off');

load(fullfile(PATH,FILE),'PATH','COUNT','SACCADE','SACCADE_STATS','FLY','GRAND','Stim','D','I','U','N')

clms = N.vel/2;
CC = repmat(hsv(clms),2,1);

% Vel = U.vel{1};
% Speed = Vel(1:N.vel/2);

%% Frequency Domain
time = SACCADE.saccade{1}.time;
[Fv] = FFT(time,SACCADE.saccade{1}.position);
MAG = cell(N.vel,1);
for jj = 1:N.vel
    MAG{jj} = cell(N.fly,1);
end
for kk = 1:N.file
    obj = SACCADE.saccade{kk};
    if ~isempty(obj)
        pos = obj.position;
        % pos = obj.shift.IntrpPosition;
        [~,mag] = FFT(time, pos);
        MAG{I.vel(kk)}{I.fly(kk),1}(:,end+1) = mag;
    end 
end

Mag_med = cell(N.vel,1);
Mag_med_all = cell(N.vel,1);
for jj = 1:N.vel
    Mag_med{jj} = cellfun(@(x) median(x,2), MAG{jj}, 'UniformOutput', false);
    Mag_med_all{jj} = cat(2,Mag_med{jj}{:});
end
Mag_grand_med = cellfun(@(x) median(x,2), Mag_med_all, 'UniformOutput', false);

FIG = figure (1) ; clf
FIG.Color = 'w';
FIG.Units = 'inches';
FIG.Position = [2 2 clms*(4/3) 3];
movegui(FIG,'center')
ax = gobjects(N.vel,1);
alpha = 0.2;
for jj = 1:N.vel
    ax(jj) = subplot(2,N.vel/2,jj); hold on
    title([num2str(U.vel{1}(jj)) ' (�/s)'])
        plot(Fv,Mag_med_all{jj}, 'Color', [CC(jj,:) alpha])
        plot(Fv,Mag_grand_med{jj}, 'Color', 'k', 'LineWidth', 1)
end

set(ax,'LineWidth',1,'FontWeight','bold','XLim',[0 45],'YLim',[0 1])
linkaxes(ax,'xy')
XLabelHC = get(ax(clms+1:N.vel), 'XLabel');
set([XLabelHC{:}], 'String', 'Frequency (Hz)')
YLabelHC = get(ax([1,clms+1]), 'YLabel');
set([YLabelHC{:}], 'String', 'Head Magnitude (�)')
end