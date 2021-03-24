%% vIRt paper figures %%
%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
clear global
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

doPlot={'Spectrum'}; %'TuningPlots' PTPlots

%% Population phase tuning
if any(contains(doPlot,'TuningPlots'))
    if exist(fullfile(baseDir,'Analysis','Cell_Tuning.mat'),'file')
        load(fullfile(baseDir,'Analysis','Cell_Tuning.mat'));
    else
        
        cellTuning=struct('global',struct('phaseStats',[],'phaseTuning',[],...
            'meanFR',[],'edges',[],'phaseCoherence',[]),...
            'epochs',struct('phaseStats',[],'phaseTuning',[],...
            'meanFR',[],'edges',[],'phaseCoherence',[]));
        wS=struct('angle',[],'phase',[]);
        
        for cellNum=1:numel(tunedCells)
            uIdx=tunedCells(cellNum);
            sessID=[char(cellList.Session(uIdx)) '_' num2str(cellList.RecordingID(uIdx))];
            dataDir=fullfile(baseDir,'Analysis','Data',sessID);
            load(fullfile(dataDir,[sessID '_behavior.mat']),'whiskers','wEpochMask','bWhisk');
            load(fullfile(dataDir,[sessID '_ephys.mat']),'ephys');
            spikes=load(fullfile(dataDir,[sessID '_Unit' num2str(cellList.unitIndex(uIdx)) '.mat']));
            ephys.selectedUnits=spikes.unitId;
            load(fullfile(dataDir,[sessID '_recInfo.mat']),'recInfo');
            ephys.recInfo=recInfo;
            
            %Angle and Phase spectrums
            wS(cellNum).angle=vIRt_WhiskingSpectrum(whiskers(bWhisk).angle_BP,wEpochMask);
            wS(cellNum).phase=vIRt_WhiskingSpectrum(whiskers(bWhisk).phase,wEpochMask);
            
            
            %% return global coherence (with stats, etc)
            % Tuning, Coherence, Stats thetas,phaseStats,phaseTuning,phaseCoherence
            [r.meanFR,r.edges,r.phaseStats,r.phaseTuning,r.phaseCoherence]=...
                vIRt_PhaseTuning(whiskers(bWhisk).phase,ephys,wEpochMask,false);
            cellTuning(cellNum).global=r; clearvars r;
            %% same but for each epoch
            [r.meanFR,r.edges,r.phaseStats,r.phaseTuning,r.phaseCoherence]=...
                vIRt_PhaseTuning(whiskers(bWhisk).phase,ephys,wEpochMask,true);
            cellTuning(cellNum).epochs=r; clearvars r;
            %  circ_rtest(cellTuning(cellNum).meanFR{24})
            
            clearvars 'whiskers' 'wEpochMask' 'bWhisk' 'ephys' 'recInfo'
        end
    end
    
    % for each cell, get which has significant different PDF
    phaseDiffTest=struct('globalPhaseDiff',[],'epochPhaseDiffIdx',[]);
    for cellNum=1:numel(tunedCells)
        phaseDiffTest(cellNum).globalPhaseDiff=cellTuning(cellNum).global.phaseStats.spikePhaseStats(1);
        phaseDiffTest(cellNum).epochPhaseDiffIdx=...
            cellfun(@(x) x(1)<=0.05, {cellTuning(cellNum).epochs.phaseStats.spikePhaseStats});
    end
    
    %% plot bad ones
    noDiffPDF=[13,18]; %find(phaseDiffTest>0.01);
    numBins=32;
    for cellNum=1:numel(noDiffPDF)
        spikePhasePDF=cellTuning(noDiffPDF(cellNum)).global.phaseStats.spikePhasePDF;
        phasePDF=cellTuning(noDiffPDF(cellNum)).global.phaseStats.phasePDF;
        figure('position',[1202         323         583         487]); hold on ;
        plot(linspace(-pi,pi, numBins+1),spikePhasePDF,'linewidth',1.2,'Color', [0 0 0]); %centers
        plot(linspace(-pi,pi, numBins+1),phasePDF,'linewidth',1.2,'Color', [0 0 0 0.5]); %centers
        set(gca,'ytick',0:0.05:1,...
            'xlim',[-pi pi],'xtick',[-pi 0 pi],'xticklabel',{'-\pi','0','\pi'},...
            'tickdir','out');
        axis tight
        legend('P(\phi_k|spike)','P(\phi_k)','location','southeast')
        legend('boxoff')
        title({'Probability density function'; 'of phase for spiking events'})
    end
    
    % for each cell get which epochs has significant coherence
    propEpochCoh=struct('coherEpochIdx',[],'fractionCoherEpoch',[],'manualClass',[]);
    
    for cellNum=1:numel(tunedCells)
        epochCoh=cellfun(@(x,y) any(x>=y), {cellTuning(cellNum).epochs.phaseCoherence.coherMag},{cellTuning(cellNum).epochs.phaseCoherence.confC});
        propEpochCoh(cellNum).fractionCoherEpoch=sum(epochCoh)/numel(epochCoh);
        if cellList.tuningEpochs(cellNum) == 'all' % for comparison with manual classification
            propEpochCoh(cellNum).manualClass=1;
        else
            propEpochCoh(cellNum).manualClass=0;
        end
        propEpochCoh(cellNum).coherEpochIdx=epochCoh;
    end
    
    %     [rhos,thetas]=deal(nan(numel(tunedCells),1));
    
    for cellNum=1:numel(tunedCells)
        if propEpochCoh(cellNum).fractionCoherEpoch==1 && phaseDiffTest(cellNum).globalPhaseDiff<=0.05
            %         case 'all'
            coherVals(cellNum)=cellTuning(cellNum).global.phaseCoherence;
            peakIdx=coherVals(cellNum).coherMag==max(coherVals(cellNum).coherMag);
            
            cellTuning(cellNum).peakCMag=coherVals(cellNum).coherMag(peakIdx);
            cellTuning(cellNum).cMag=coherVals(cellNum).coherMag;
            
            cellTuning(cellNum).peakCPhase=coherVals(cellNum).coherPhase(peakIdx);
            cellTuning(cellNum).cPhase=coherVals(cellNum).coherPhase;
            
            cellTuning(cellNum).meanPhase=cellTuning(cellNum).global.phaseStats.mean;
            %         end
            %         if ~isfield(cellTuning,'peakCMag') || isempty(cellTuning(cellNum).peakCMag) || cellTuning(cellNum).peakCMag<0.3
        else
            %         otherwise
            % Tuning epoch noted in cellList:
            % str2double(char(cellList.tuningEpochs(cellNum)))
            % Too abitrary. Average over all epochs with Coherence Magnitude > 0.3.
            % even better: take only epochs with significant coherence
            
            epochIdx=propEpochCoh(cellNum).coherEpochIdx & phaseDiffTest(cellNum).epochPhaseDiffIdx;
            
            rhos_allEpochs=[cellTuning(cellNum).epochs.phaseCoherence.peakCoherMag];
            cellTuning(cellNum).peakCMag=mean(rhos_allEpochs(epochIdx));%rhos_allEpochs>0.3
            
            thetas_allEpochs=[cellTuning(cellNum).epochs.phaseCoherence.peakCoherPhase];
            cellTuning(cellNum).peakCPhase=circ_mean(thetas_allEpochs(epochIdx)');%rhos_allEpochs>0.3
            
            cellTuning(cellNum).meanPhase=circ_mean([cellTuning(cellNum).epochs.phaseStats(epochIdx).mean]');%rhos_allEpochs>0.3
        end
    end
    rhos=[cellTuning.peakCMag]; % #13 and 18 are out based on phase diff and coherence
    thetas=[cellTuning.peakCPhase];
    thetas(isnan(rhos))=NaN;
    
    % plot coherence in polar coordinates
    figure;
    %     polarplot(thetas,rhos,'o','LineWidth',2);
    
    polarplot(thetas(1:5),rhos(1:5),'o','MarkerFaceColor','b','MarkerEdgeColor','None','LineWidth',2.5);
    
    paH = gca;
    paH.ThetaZeroLocation='left';
    paH.ThetaTickLabel={'max Protraction','','','Retraction','','',...
        'max Retraction','','','Protraction','',''};
    paH.ThetaDir = 'counterclockwise';
    
    hold on
    % other marker
    P_group=thetas>=deg2rad(150) | thetas<deg2rad(-115);
    R_group=thetas>=deg2rad(-30) & thetas<deg2rad(65);
    midP_group=thetas>=deg2rad(-115) & thetas<deg2rad(-30);
    midR_group=thetas>=deg2rad(65) & thetas<deg2rad(150);
    
    polarplot(thetas(P_group),rhos(P_group),'o',...
        'MarkerEdgeColor',cmap(2,:),'MarkerFaceColor','None','LineWidth',2);
    polarplot(thetas(R_group),rhos(R_group),'o',...
        'MarkerEdgeColor',cmap(5,:),'MarkerFaceColor','None','LineWidth',2);
    polarplot(thetas(midP_group),rhos(midP_group),'o',...
        'MarkerEdgeColor',cmap(4,:),'MarkerFaceColor','None','LineWidth',2);
    polarplot(thetas(midR_group),rhos(midR_group),'o',...
        'MarkerEdgeColor',cmap(3,:),'MarkerFaceColor','None','LineWidth',2);
    
    %     polarplot(thetas(32),rhos(32),'or','LineWidth',2);
    %     polarplot(thetas(16),rhos(16),'vr','LineWidth',2);
    %
    %     Pcells=contains(string(cellList.Tuning_Comment),{'P';'mid-P';'mid P'});
    %
    %     polarplot(thetas(Pcells(1:41)),rhos(Pcells(1:41)),'dk','LineWidth',2);
    % polarplot(mean(thetas),mean(rhos),'db','LineWidth',2);
    % polarplot(mean(thetas(rhos>0.3)),mean(rhos(rhos>0.3)),'dr','LineWidth',2);
    
    %     figure;
    %     plot(cellTuning(32).global.phaseCoherence.freqVals,cellTuning(32).global.phaseCoherence.coherMag)
    
    %     title(['Unit ' num2str(ephysData.selectedUnits(unitNum)) ' - ' strrep(recName,'_','') ...
    %         ' - Tuning to ' labels ' phase'],'interpreter','none');
    
    %% get mean coherence and FR? for each group
    %   mean ICoherenceI ± SD =
    mean(rhos(P_group)) % 0.6016
    std(rhos(P_group))% 0.1136
    mean(rhos(R_group)) % 0.5535
    std(rhos(R_group))% 0.1396
    mean(rhos(midP_group)) % 0.4499
    std(rhos(midP_group))% 0.1311
    mean(rhos(midR_group)) % 0.5088
    std(rhos(midR_group))% 0.1941
    
    mean(rhos(1:5)) % 0.7098
    std(rhos(1:5))% 0.1537
    
    mean(thetas(P_group)+pi) % -1.7615
    std(thetas(P_group)+pi)% 0.1136
    mean(thetas(R_group)+pi) % 0.5535
    std(thetas(R_group)+pi)% 0.1396
    mean(thetas(midP_group)+pi) % 0.4499
    std(thetas(midP_group)+pi)% 0.1311
    mean(thetas(midR_group)+pi) % 0.5088
    std(thetas(midR_group)+pi)% 0.1941
    
    mean(thetas(1:5)+pi) % 0.7098
    std(thetas(1:5)+pi)% 0.1537
    
    %% plot tuning curves on top of each other
    gLabels={'Protraction','Retraction','mid Protraction','mid Retraction'};
    gIdx={find(P_group),find(R_group), find(midP_group), find(midR_group)};
    gCMap=[cmap(2,:); cmap(5,:);cmap(3,:);cmap(4,:)];
    figure('name','Tuning to Whisking Phase','Color','white','position',...
        [1054 279 678 612]);
    for gNum=1:4
        [phEdgeColor,phFaceColor]=deal(gCMap(gNum,:));
        
        [cThetas,cEdges]=deal(cell(numel(gIdx{gNum}),1));
        for gCellNum=1:numel(gIdx{gNum})
            cellNum=gIdx{gNum}(gCellNum);
            if propEpochCoh(cellNum).fractionCoherEpoch==1 && phaseDiffTest(cellNum).globalPhaseDiff<=0.05
                cThetas{gCellNum}=cellTuning(cellNum).global.meanFR{:};
                cEdges{gCellNum}=cellTuning(cellNum).global.edges{:};
            else
                epochIdx=propEpochCoh(cellNum).coherEpochIdx & phaseDiffTest(cellNum).epochPhaseDiffIdx;
                cThetas{gCellNum}=vertcat(cellTuning(cellNum).epochs.meanFR{epochIdx});
                cEdges{gCellNum}=cellTuning(cellNum).epochs.edges{1};
                
                %                 [cThetas,cEdges] = histcounts(vertcat(cellTuning(1).epochs.meanFR{epochIdx}),cellTuning(1).epochs.edges{1});
                
            end
            polarhistogram(cThetas{gCellNum},cEdges{gCellNum},'Displaystyle','bar',...
                'Normalization','probability','LineWidth',2,...
                'EdgeColor',phEdgeColor,'FaceColor',phFaceColor,...
                'FaceAlpha',0.2,'EdgeAlpha',0);
            hold on
        end
        
        paH = gca;
        paH.ThetaZeroLocation='left';
        paH.ThetaTickLabel={'max Protraction','','','Retraction','','',...
            'max Retraction','','','Protraction','',''};
        paH.ThetaDir = 'counterclockwise';
        
        % plot all group outline
        polarhistogram(vertcat(cThetas{:}),cEdges{1},'Displaystyle','stairs',...
            'Normalization','probability','LineWidth',2,...
            'EdgeColor','k','FaceColor','none',...
            'FaceAlpha',0,'EdgeAlpha',0.8);
        
        %         title(['Tuning to Whisking phase - group ' gLabels{gNum}] ,'interpreter','none');
        
    end
end

%% Population power spectrum. Classify fast/slow oscillation

if any(contains(doPlot,'Spectrum'))
    if exist(fullfile(baseDir,'Analysis','Cell_PSD.mat'),'file')
        load(fullfile(baseDir,'Analysis','Cell_PSD.mat'));
    else    
        spS=struct('spectrumVals',[],'freqVals',[],'R',[],'Serr',[],...
            'spectrumValsPSD',[],'SerrPSD',[],'RPSD',[],'StatSigIdx',[]);
        for cellNum=1:size(cellList,1)
            sessID=[char(cellList.Session(cellNum)) '_' num2str(cellList.RecordingID(cellNum))];
            dataDir=fullfile(baseDir,'Analysis','Data',sessID);
            load(fullfile(dataDir,[sessID '_behavior.mat']),'wEpochMask');
            %     load(fullfile(dataDir,[sessID '_ephys.mat']),'ephys');
            spikes=load(fullfile(dataDir,[sessID '_Unit' num2str(cellList.unitIndex(cellNum)) '.mat']));
            %     ephys.selectedUnits=spikes.unitId;
            %     load(fullfile(dataDir,[sessID '_recInfo.mat']),'recInfo');
            %     ephys.recInfo=recInfo;
            %         spikeData.spikeTimes=spikes.spikeTimes;
            try
                spS(cellNum)=vIRt_SpikingSpectrum(spikes.spikeTimes,wEpochMask);%wEpochMask
            catch
                continue
            end
        end
        
        for cellNum=1:size(cellList,1)
            %     figure;
            %     plot(spS(2).freqVals(spS(2).StatSigIdx),spS(2).spectrumValsPDF(spS(2).StatSigIdx),'k')
            %     plot(ans(1).freqVals(ans(1).StatSigIdx),ans(1).spectrumValsPSD(ans(1).StatSigIdx),'k')
            
            sigFreq=spS(cellNum).freqVals(spS(cellNum).StatSigIdx);
            sigPSD=spS(cellNum).spectrumValsPSD(spS(cellNum).StatSigIdx);
            spS(cellNum).peakFreq=sigFreq(find(sigPSD==max(sigPSD),1));
            spS(cellNum).peakPSD=sigPSD(find(sigPSD==max(sigPSD),1));
            spS(cellNum).meanFreq=sigFreq(find(sigFreq>=mean(sigFreq),1));
            spS(cellNum).meanPSD=sigPSD(find(sigFreq>=mean(sigFreq),1));
        end
        
        emptyValIdx=cellfun(@isempty, {spS.meanPSD});
        [spS(emptyValIdx).peakFreq,spS(emptyValIdx).peakPSD,...
            spS(emptyValIdx).meanFreq,spS(emptyValIdx).meanPSD]=deal(NaN);
    end
    
    meanValIdx=~cellfun(@isnan, {spS.meanPSD})';
    manualTuningClass=cellList.Tuning;
    PTClass=cellList.PT;
    tClass=unique(manualTuningClass);
    
    figure; hold on
    for tClassNum=4:-1:1
        cellIdx=manualTuningClass==tClass(tClassNum) & meanValIdx;
        plot([spS(cellIdx).meanFreq],...
            [spS(cellIdx).meanPSD],'.')
        
        %      plot([spS(cellIdx).peakFreq],...
        %          [spS(cellIdx).peakPSD],'.')
    end
    
    cellIdx=1:5; %PTClass==1;
    plot([spS(cellIdx).meanFreq],...
        [spS(cellIdx).meanPSD],'bd')
    
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
    for cellNum=1:5%numel(allCells) %PTCells
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
