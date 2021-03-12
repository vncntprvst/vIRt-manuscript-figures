%% vIRt paper figures %%
%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
global CHRONUXGPU
CHRONUXGPU = 1;
%% define figure colormap
cmap=lines;cmap=[cmap(1:7,:);(lines+flipud(copper))/2;autumn];

% data should be in DJ
baseDir='D:\Vincent\';
% load(fullfile(baseDir,'Analysis','vIRt_sessions_analysis.mat'));
% Data already saved with vIRt_SaveData

load(fullfile(baseDir,'Analysis','Cell_List.mat'))
allCells=1:size(cellList,1);
tunedCells=find(cellList.Tuning==1);
PTCells=find(cellList.PT==1);

doPlot={'PTPlots'}; %'TuningPlots'

%% Figure 2. Population phase tuning
if any(contains(doPlot,'TuningPlots'))
    if exist(fullfile(baseDir,'Analysis','Cell_Tuning.mat'),'file')
        load(fullfile(baseDir,'Analysis','Cell_List.mat'));
    else
        cellTuning=struct('global',struct('phaseStats',[],'phaseTuning',[],...
            'thetas',[],'phaseCoherence',[]),...
            'epochs',struct('phaseStats',[],'phaseTuning',[],...
            'thetas',[],'phaseCoherence',[]));
        wS=struct('angle',[],'phase',[]);
        
        for cellNum=1:numel(tunedCells)
            uIdx=tunedCells(cellNum);
            sessID=[char(cellList.Session(uIdx)) '_' num2str(cellList.RecordingID(uIdx))];
            dataDir=fullfile(baseDir,'Analysis','Data',sessID);
            load(fullfile(dataDir,[sessID '_behavior.mat']),'whiskers','wEpochMask','bWhisk');
            %             load(fullfile(dataDir,[sessID '_ephys.mat']),'ephys');
            ephys=load(fullfile(dataDir,[sessID '_Unit' num2str(cellList.unitIndex(uIdx)) '.mat']));
            ephys.selectedUnits=ephys.unitId;
            load(fullfile(dataDir,[sessID '_recInfo.mat']),'recInfo');
            ephys.recInfo=recInfo;
            
            %Angle and Phase spectrums
            wS(cellNum).angle=vIRt_WhiskingSpectrum(whiskers(bWhisk).angle_BP,wEpochMask);
            wS(cellNum).phase=vIRt_WhiskingSpectrum(whiskers(bWhisk).phase,wEpochMask);
            
            % Tuning, Coherence, Stats thetas,phaseStats,phaseTuning,phaseCoherence
            [cellTuning(cellNum).global.thetas,cellTuning(cellNum).global.phaseStats,...
                cellTuning(cellNum).global.phaseTuning,cellTuning(cellNum).global.phaseCoherence]=...
                vIRt_PhaseTuning(whiskers(bWhisk).phase,ephys,wEpochMask,false);
            [cellTuning(cellNum).epochs.thetas,cellTuning(cellNum).epochs.phaseStats,...
                cellTuning(cellNum).epochs.phaseTuning,cellTuning(cellNum).epochs.phaseCoherence]=...
                vIRt_PhaseTuning(whiskers(bWhisk).phase,ephys,wEpochMask,true);
            %  circ_rtest(cellTuning(cellNum).thetas{24})
            
            clearvars 'whiskers' 'wEpochMask' 'bWhisk' 'ephys' 'recInfo'
        end
    end
    
    %     [rhos,thetas]=deal(nan(numel(tunedCells),1));
    
    for cellNum=1:numel(tunedCells)
        if cellList.tuningEpochs(cellNum) == 'all'
            %         case 'all'
            coherVals(cellNum)=cellTuning(cellNum).global.phaseCoherence;
            peakIdx=coherVals(cellNum).coherMag==max(coherVals(cellNum).coherMag);
            cellTuning(cellNum).rhos=coherVals(cellNum).coherMag(peakIdx);
            cellTuning(cellNum).thetas=cellTuning(cellNum).global.phaseStats.mean;
            %             cellTuning(cellNum).thetas=coherVals(cellNum).coherPhase(peakIdx);
        end
        if isnan(rhos(cellNum)) || rhos(cellNum)<0.3
            %         otherwise
            % Tuning epoch noted in cellList:
            % str2double(char(cellList.tuningEpochs(cellNum)))
            % Too abitrary. Average over all epochs with Coherence Magnitude > 0.3.
            rhos_allEpochs=[cellTuning(cellNum).epochs.phaseCoherence.peakCoherMag];
            cellTuning(cellNum).rhos=mean(rhos_allEpochs(rhos_allEpochs>0.3));
            thetas_allEpochs=[cellTuning(cellNum).epochs.phaseCoherence.peakCoherPhase];
            cellTuning(cellNum).thetas=circ_mean([cellTuning(cellNum).epochs.phaseStats(rhos_allEpochs>0.3).mean]');
            %             cellTuning(cellNum).thetas=circ_mean(thetas_allEpochs(rhos_allEpochs>0.3)');
        end
    end
    thetas=[cellTuning.thetas];
    rhos=[cellTuning.rhos];
    
    figure;
    polarplot(thetas,rhos,'o','LineWidth',2);
    paH = gca;
    paH.ThetaZeroLocation='left';
    paH.ThetaTickLabel={'max Protraction','','','Retraction','','',...
        'max Retraction','','','Protraction','',''};
    paH.ThetaDir = 'counterclockwise';
    
    hold on
    polarplot(thetas(1:5),rhos(1:5),'ob','MarkerFaceColor','b','LineWidth',2);
    polarplot(thetas(32),rhos(32),'or','LineWidth',2);
    polarplot(thetas(16),rhos(16),'vr','LineWidth',2);
    
    Pcells=contains(string(cellList.Tuning_Comment),{'P';'mid-P';'mid P'});
    
    polarplot(thetas(Pcells(1:41)),rhos(Pcells(1:41)),'dk','LineWidth',2);
    % polarplot(mean(thetas),mean(rhos),'db','LineWidth',2);
    % polarplot(mean(thetas(rhos>0.3)),mean(rhos(rhos>0.3)),'dr','LineWidth',2);
    
    figure;
    plot(cellTuning(32).global.phaseCoherence.freqVals,cellTuning(32).global.phaseCoherence.coherMag)
    
    %     title(['Unit ' num2str(ephysData.selectedUnits(unitNum)) ' - ' strrep(recName,'_','') ...
    %         ' - Tuning to ' labels ' phase'],'interpreter','none');
end

%% PT plots

% waveforms
% traces
% rasters
% latency histogram + SALT
% tuning

if any(contains(doPlot,'PTPlots'))
    taggedCells=ones(numel(allCells),1);
    for cellNum=1:numel(allCells) %PTCells
        uIdx=allCells(cellNum);
        sessID=[char(cellList.Session(uIdx)) '_' num2str(cellList.RecordingID(uIdx))];
        dataDir=fullfile(baseDir,'Analysis','Data',sessID);
        load(fullfile(dataDir,[sessID '_ephys.mat']),'ephys');
        load(fullfile(dataDir,[sessID '_pulses.mat']),'pulses');
        load(fullfile(dataDir,[sessID '_recInfo.mat']),'recInfo');
        ephys.recInfo=recInfo;
        ephys.selectedUnits=cellList.unitIndex(uIdx);
                
        taggedCells(cellNum) = FindPhototagged(ephys,pulses);
    end
    for cellNum=1:numel(allCells) %PTCells
        if taggedCells(cellNum)<0.01 && cellList.unitFrequency(cellNum)>=0.04
            %% Check Phototagging summary
            
            uIdx=allCells(cellNum);
            sessID=[char(cellList.Session(uIdx)) '_' num2str(cellList.RecordingID(uIdx))];
            dataDir=fullfile(baseDir,'Analysis','Data',sessID);
            load(fullfile(dataDir,[sessID '_ephys.mat']),'ephys');
            load(fullfile(dataDir,[sessID '_pulses.mat']),'pulses');
            load(fullfile(dataDir,[sessID '_recInfo.mat']),'recInfo');
            ephys.recInfo=recInfo;
            ephys.selectedUnits=cellList.unitIndex(uIdx);

            %load spikes data
            ephys.spikes=load(fullfile(dataDir,[sessID '_Unit' num2str(cellList.unitIndex(uIdx)) '.mat']));
            ephys.spikes.selectedUnits=ephys.spikes.unitId;
            ephys.spikes.times=ephys.spikes.spikeTimes;
            ephys.spikes.unitID=ones(numel(ephys.spikes.times),1)*ephys.spikes.selectedUnits;
            ephys.spikes.preferredElectrode=ephys.spikes.preferredEl;
            ephys.spikes.waveforms=ephys.spikes.waveForms;
            ephys.spikes=rmfield(ephys.spikes,{'spikeTimes','unitId','preferredEl','waveForms'});
            if numel(ephys.spikes.rasters)==1
                ephys.spikes.rasters=ephys.rasters(ephys.spikes.selectedUnits,:);
            end

            if ~isfield(pulses,'duration')
                sessInfoFile=fullfile(ephys.recInfo.dirName,[ephys.recInfo.baseName '_info.json']);
                    if exist(sessInfoFile,'file')
                        sessInfo = fileread(sessInfoFile);
                        sessInfo = jsondecode(sessInfo);
                        pulses.duration=sessInfo.photoStim.pulseDur;
                        ephys.spikes.bitResolution=sessInfo.bitResolution;
                    else                
                        pulses.duration=0.010; 
                    end
            end
            if ~isfield(ephys,'traces')
                tracefileIdx=cellfun(@(x) contains(x,'traces'),{ephys.recInfo.sessFiles.name});
                if ~any(tracefileIdx)
                    tracefileIdx=cellfun(@(x) contains(x,'rec.bin'),{ephys.recInfo.sessFiles.name});
                end
                traceFile = fopen(fullfile(ephys.recInfo.sessFiles(tracefileIdx).folder,...
                    ephys.recInfo.sessFiles(tracefileIdx).name), 'r');
                ephys.traces = fread(traceFile,[ephys.recInfo.numRecChan,Inf],'single');
                fclose(traceFile);
            end

            PhotoTagPlots(ephys,pulses); %Implement SALT test
            ephys=rmfield(ephys,'traces');
        end
    end
end

if any(contains(doPlot,'OverviewPlot'))
    %% Overview plot
    NBC_Plots_Overview(whiskers(bWhisk),whiskingEpochs,breathing,ephys,pulses.TTLTimes,false,false);
end


%% transition rythm plot

%% Supplementary
% Slow rhythm tunings: setpoint / breathing
% Proba tuning plot: main figure for now
% retraction cell video
