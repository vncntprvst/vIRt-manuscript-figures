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
PTCells=find(cellList.PT==1);

doPlot={'transition'}; %'PTPlots TuningPlots' CV2 CV2_Slow Spectrum transition transition_Slow

if contains(doPlot,'Slow')
    %Slow oscillation analysis:
    load(fullfile(baseDir,'Analysis','SlowOscillationIdx.mat'))
    tunedCells= find(lowOscillationIdx); %see "Spectrum" analysis below
else
    tunedCells=find(cellList.Tuning==1);
end

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
    %     noDiffPDF=[13,18]; %find(phaseDiffTest>0.01);
    %     numBins=32;
    %     for cellNum=1:numel(noDiffPDF)
    %         spikePhasePDF=cellTuning(noDiffPDF(cellNum)).global.phaseStats.spikePhasePDF;
    %         phasePDF=cellTuning(noDiffPDF(cellNum)).global.phaseStats.phasePDF;
    %         figure('position',[1202         323         583         487]); hold on ;
    %         plot(linspace(-pi,pi, numBins+1),spikePhasePDF,'linewidth',1.2,'Color', [0 0 0]); %centers
    %         plot(linspace(-pi,pi, numBins+1),phasePDF,'linewidth',1.2,'Color', [0 0 0 0.5]); %centers
    %         set(gca,'ytick',0:0.05:1,...
    %             'xlim',[-pi pi],'xtick',[-pi 0 pi],'xticklabel',{'-\pi','0','\pi'},...
    %             'tickdir','out');
    %         axis tight
    %         legend('P(\phi_k|spike)','P(\phi_k)','location','southeast')
    %         legend('boxoff')
    %         title({'Probability density function'; 'of phase for spiking events'})
    %     end
    %
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
            % Even better: take only epochs with significant coherence
            
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
    
    % define groups
    P_group=thetas>=deg2rad(150) | thetas<deg2rad(-115);% & ~lowMedFreq;
    R_group=thetas>=deg2rad(-30) & thetas<deg2rad(65);% & ~lowMedFreq;
    midP_group=thetas>=deg2rad(-115) & thetas<deg2rad(-30);% & ~lowMedFreq;
    midR_group=thetas>=deg2rad(65) & thetas<deg2rad(150);% & ~lowMedFreq;
    
    % plot coherence in polar coordinates
    figure;
    %     polarplot(thetas,rhos,'o','LineWidth',2);
    polarplot(thetas(R_group),rhos(R_group),'o',...
        'MarkerFaceColor','k','MarkerEdgeColor','None','LineWidth',2); %cmap(5,:)
    
    paH = gca;
    paH.ThetaZeroLocation='left';
    paH.ThetaTickLabel={'max Protraction','','','Retraction','','',...
        'max Retraction','','','Protraction','',''};
    paH.ThetaDir = 'counterclockwise';
    
    hold on
    
    
    polarplot(thetas(P_group),rhos(P_group),'o',...
        'MarkerFaceColor',cmap(2,:),'MarkerEdgeColor','None','LineWidth',2);
    polarplot(thetas(midR_group),rhos(midR_group),'o',...
        'MarkerFaceColor',cmap(3,:),'MarkerEdgeColor','None','LineWidth',2);
    polarplot(thetas(midP_group),rhos(midP_group),'o',...
        'MarkerFaceColor',cmap(4,:),'MarkerEdgeColor','None','LineWidth',2);
    polarplot(thetas(1:5),rhos(1:5),'o',...
        'MarkerFaceColor','b','MarkerEdgeColor','None','LineWidth',2.5);
    
    legend({'Retraction group','Protraction group','Mid-retraction group',...
        'Mid-protraction group','Photo-tagged cells'},'FontSize',8);
    legend('boxoff')
    
    %% other markers
    %     % median Frequency PSD groups (see PSD analysis)
    %     lowMedFreq=[spS(1:41).medFreq]<8;
    %     polarplot(thetas(lowMedFreq),rhos(lowMedFreq),'o',...
    %         'MarkerEdgeColor','None','MarkerFaceColor','r','LineWidth',2);
    
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
                %                 spS(cellNum)=vIRt_FRSpectrum(spikes.spikeTimes,wEpochMask);%wEpochMask
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
            
            breathFreqIdx=spS(cellNum).freqVals>=3 & spS(cellNum).freqVals<=8;
            whiskFreqIdx=spS(cellNum).freqVals>=12 & spS(cellNum).freqVals<=17;
            spS(cellNum).breathPSDInt=trapz(spS(cellNum).spectrumValsPSD(breathFreqIdx));
            spS(cellNum).whiskPSDInt=trapz(spS(cellNum).spectrumValsPSD(whiskFreqIdx));
            spS(cellNum).BWratio=spS(cellNum).whiskPSDInt/spS(cellNum).breathPSDInt;
            try
                freqIdx=spS(cellNum).freqVals>=3 & spS(cellNum).freqVals<=20;
                spS(cellNum).medFreq=medfreq(spS(cellNum).spectrumVals(freqIdx),spS(cellNum).freqVals(freqIdx));
            catch
                continue
            end
        end
        
        emptyValIdx=cellfun(@isempty, {spS.medFreq});
        [spS(emptyValIdx).medFreq]=deal(NaN);
        
        emptyValIdx=cellfun(@isempty, {spS.meanPSD});
        [spS(emptyValIdx).peakFreq,spS(emptyValIdx).peakPSD,...
            spS(emptyValIdx).meanFreq,spS(emptyValIdx).meanPSD]=deal(NaN);
        
    end
    
    meanValIdx=~cellfun(@isnan, {spS.meanPSD})';
    medFreqIdx=~cellfun(@isempty, {spS.medFreq})';
    manualTuningClass=cellList.Tuning;
    PTClass=cellList.PT;
    tClass=unique(manualTuningClass);
    
    figure; hold on
    for tClassNum=4:-1:1
        cellIdx=manualTuningClass==tClass(tClassNum) & medFreqIdx; %meanValIdx;
        
        scatter([spS(cellIdx).medFreq],...
            [spS(cellIdx).R],'.')
        
        %         plot([spS(cellIdx).breathPSDInt],...
        %             [spS(cellIdx).whiskPSDInt],'.')
        
        %         plot([spS(cellIdx).meanFreq],...
        %             [spS(cellIdx).meanPSD],'.')
        %
        %      plot([spS(cellIdx).peakFreq],...
        %          [spS(cellIdx).peakPSD],'.')
    end
    
    cellIdx=PTClass==1;
    plot([spS(cellIdx).medFreq],...
        [spS(cellIdx).R],'kd')
    
    cellIdx=1:5; %PTClass==1;
    plot([spS(cellIdx).medFreq],...
        [spS(cellIdx).R],'bd')
    
    cellIdx=104;
    plot([spS(cellIdx).medFreq],...
        [spS(cellIdx).R],'rd')
    
    %
    %     cellIdx=PTClass==1;
    %     plot([spS(cellIdx).meanFreq],...
    %         [spS(cellIdx).meanPSD],'kd')
    %
    %     cellIdx=1:5; %PTClass==1;
    %     plot([spS(cellIdx).meanFreq],...
    %         [spS(cellIdx).meanPSD],'bd')
    %
    %     cellIdx=PTClass==1;
    %     plot([spS(cellIdx).breathPSDInt],...
    %         [spS(cellIdx).whiskPSDInt],'kd')
    %
    %     cellIdx=1:5; %PTClass==1;
    %     plot([spS(cellIdx).breathPSDInt],...
    %         [spS(cellIdx).whiskPSDInt],'bd')
    
    %% low oscillation - PT
    %   lowOscillationIdx=cellfun(@(x) x<10 ,{spS.medFreq},'UniformOutput',true)';
    %   Started with that lowOscillationIdx index and run CellTuning and CV2
    %   analyses to get those:
    load(fullfile(baseDir,'Analysis','SlowOscillationIdx.mat'))
    load(fullfile(baseDir,'Analysis','Cell_CV2_SO.mat'))
    
    lowOscillationIdx=find(lowOscillationIdx);
    lowOscillationIdx=lowOscillationIdx(~isnan([CV2.CV2_withinBurst_ISI]));
    
    srCells=cellList(lowOscillationIdx,:);
    srSpS=spS(lowOscillationIdx);
    cellIdx=srCells.PT==1;
    plot([srSpS(cellIdx).medFreq],...
        [srSpS(cellIdx).R],'ko')
    
    CV2wbISI=[CV2(~isnan([CV2.CV2_withinBurst_ISI])).CV2_withinBurst_ISI];
    
    figure; hold on
    plot([srSpS.medFreq],...
        CV2wbISI,'k.')
    plot([srSpS(cellIdx).medFreq],...
        CV2wbISI(cellIdx),'bo')
    
    %     #104 isn't part of it because meanFreq is NaN - check that
    
    PT_LO=lowOscillationIdx(cellIdx);
    %     PT_LO([srSpS(cellIdx).R]>20) %sufficient rate? Not appropriate
    
end

%% PT plots
if any(contains(doPlot,'PTPlots'))
    if exist(fullfile(baseDir,'Analysis','Cell_PT.mat'),'file')
        load(fullfile(baseDir,'Analysis','Cell_PT.mat'));
    else
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
    end
    PTCells=find(taggedCells(:,3));
    %     cellIdx=find(cellList.PT==1 & cellfun(@(x) x<10 ,{spS.medFreq},'UniformOutput',true)');
    PTmeasures=struct('latency',[],'jitter',[],...
        'ISI',struct('onPulse',[],'offPulse',[]),'waveform',[]);
    for cellNum=28:numel(PTCells) %1:5 %104 %3 %1:4 %numel(allCells) %PTCells tunedCells
        %         if taggedCells(cellNum)<0.01 && cellList.unitFrequency(cellNum)>=0.04
        %% Check Phototagging summary
        
        uIdx=PTCells(cellNum); %allCells PTCells cellIdx tunedCells
        sessID=[char(cellList.Session(uIdx)) '_' num2str(cellList.RecordingID(uIdx))];
        dataDir=fullfile(baseDir,'Analysis','Data',sessID);
        load(fullfile(dataDir,[sessID '_ephys.mat']),'ephys');
        load(fullfile(dataDir,[sessID '_pulses.mat']),'pulses');vIRt_TTL;
        load(fullfile(dataDir,[sessID '_recInfo.mat']),'recInfo');
        if ~exist('Z:/','dir')
            recInfo.dirName=strrep(recInfo.dirName,'Z:\Vincent\Ephys\','D:\Vincent\');
        end
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
        
        % benchmark with S2
        %         ephys=S2file(ephys,'load');
        %
        %         fldNames=fieldnames(ephys.spikes);
        %         s2SortIdx=cellfun(@(x) contains(x,[ephys.recInfo.baseName '_Ch']),fldNames);
        %         if any(s2SortIdx)
        %             s2Sort=fldNames{s2SortIdx};
        %             spikefileIdx=cellfun(@(x) contains(x,'_ephys'),{ephys.recInfo.sessFiles.name});
        %             originalSF=load(fullfile(ephys.recInfo.sessFiles(spikefileIdx).folder,...
        %                 ephys.recInfo.sessFiles(spikefileIdx).name),'ephys');
        %             ephys=S2sort(ephys,originalSF.ephys.spikes.times,pulses,s2Sort);
        %         end
        
        %         vIRt_TracesSpikesPulses(ephys,pulses,uIdx,'average',false); %overlay excerpt average
        
        latency=OptoJitter(ephys.spikes,pulses.TTLTimes,ephys.selectedUnits,pulses.duration,NaN);
        PTmeasures(uIdx).latency=mean(latency);
        PTmeasures(uIdx).jitter=std(latency);
        
        %         if false
        vIRt_PhotoTagPlots(ephys,pulses,uIdx,true);
        %             vIRt_CollisionTest(ephys,pulses,uIdx);
        %         end
        
        ephys=rmfield(ephys,'traces');
    end
    
    % Latency / Jitter plot
    figure('Position',[214   108   747   754],'Color','w');
    hold on
    plot([PTmeasures.latency],[PTmeasures.jitter],'o','MarkerFaceColor',...
        [0.5 0.5 0.5], 'MarkerEdgeColor','k','MarkerSize',8); %[0.3 0.75 0.93]
    plot([PTmeasures(1:5).latency],[PTmeasures(1:5).jitter],'o','MarkerFaceColor',...
        [0.3 0.75 0.93], 'MarkerEdgeColor','k','MarkerSize',8);
    
    box off; grid('on');
    set(gca,'xlim',[1 ceil(max(get(gca,'xlim'))/10)*10],...
        'ylim',[0.1 ceil(max(get(gca,'ylim'))/10)*10]);
    set(gca,'xscale','log','yscale','log','GridAlpha',0.25,'MinorGridAlpha',1)
    set(gca,'XTickLabel',get(gca,'XTick'),'YTickLabel',get(gca,'YTick'));
    set(gca,'Color','white','FontSize',10,'FontName','calibri','TickDir','out');
    %     legend('Latency','ISI','location','southeast','box','off')
    xlabel('Latency (ms)')
    ylabel('Latency jitter, SD (ms)');
    hold off
    
    % Latency histogram plot
    figure('Position',[214   108   747   754],'Color','w');
    hold on
    Lathist=histogram([PTmeasures(1:5).latency],-1:10);
    Lathist.FaceColor = [1.0000    0.6784    0.0980]; %[1.0000    0.8941    0.7686];
    Lathist.EdgeColor = 'k';
    
    xlabel('Latency distribution (ms)')
    ylabel('Count')
    axis('tight'); box off;
    %     set(gca,'xlim',[0 10^3],'Color','white','FontSize',10,'FontName','calibri','TickDir','out');
    %     set(gca,'XTick',[1,10,100,1000],'XTickLabel',[1,10,100,1000]);
    hold off
end

%% Overview plot
if any(contains(doPlot,'OverviewPlot'))
    NBC_Plots_Overview(whiskers(bWhisk),whiskingEpochs,breathing,ephys,pulses.TTLTimes,false,false);
end

%% CV2
if any(contains(doPlot,'CV2'))
    switch doPlot{contains(doPlot,'CV2')}
        case 'CV2'
            load(fullfile(baseDir,'Analysis','Cell_Tuning.mat'));
            CV2file='Cell_CV2.mat';
        case 'CV2_Slow'
            load(fullfile(baseDir,'Analysis','Cell_Tuning_SO.mat'));
            CV2file='Cell_CV2_SO.mat';
    end
    if exist(fullfile(baseDir,'Analysis',CV2file),'file')
        load(fullfile(baseDir,'Analysis',CV2file));
    else
        
        %% first compute propEpochCoh and phaseDiffTest to get epoch index
        
        % for each cell, get which has significant different PDF
        phaseDiffTest=struct('globalPhaseDiff',[],'epochPhaseDiffIdx',[]);
        for cellNum=1:numel(tunedCells)
            phaseDiffTest(cellNum).globalPhaseDiff=cellTuning(cellNum).global.phaseStats.spikePhaseStats(1);
            phaseDiffTest(cellNum).epochPhaseDiffIdx=...
                cellfun(@(x) x(1)<=0.05, {cellTuning(cellNum).epochs.phaseStats.spikePhaseStats});
        end
        % for each cell get which epochs has significant coherence
        propEpochCoh=struct('coherEpochIdx',[],'fractionCoherEpoch',[],'manualClass',[]);
        
        for cellNum=1:numel(tunedCells)
            epochCoh=cellfun(@(x,y) any(x>=y), {cellTuning(cellNum).epochs.phaseCoherence.coherMag},...
                {cellTuning(cellNum).epochs.phaseCoherence.confC});
            propEpochCoh(cellNum).fractionCoherEpoch=sum(epochCoh)/numel(epochCoh);
            if cellList.tuningEpochs(cellNum) == 'all' % for comparison with manual classification
                propEpochCoh(cellNum).manualClass=1;
            else
                propEpochCoh(cellNum).manualClass=0;
            end
            propEpochCoh(cellNum).coherEpochIdx=epochCoh;
        end
        
        CV2=struct('mVals',[],'wb_mVals',[],'CV2_all_ISI',[],'CV2_withinBurst_ISI',[]);
        
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
            
            if ~contains(doPlot,'Slow')
                wEpochMask.epochIdx=(propEpochCoh(cellNum).coherEpochIdx & phaseDiffTest(cellNum).epochPhaseDiffIdx)';
            else
                wEpochs=bwconncomp(wEpochMask.behav);
                wEpochMask.epochIdx=true(1,sum(cellfun(@(x) length(x),wEpochs.PixelIdxList)>=3000));
            end
            CV2(cellNum)=vIRt_CV2(whiskers(bWhisk).angle_BP,ephys,wEpochMask);
        end
    end
    %% plot CV2s
    
    rhos=[CV2.CV2_withinBurst_ISI]; % #13 and 18 are out based on phase diff and coherence
    thetas=[cellTuning.peakCPhase];
    thetas(isnan(rhos))=NaN;
    
    % plot CV2 in polar coordinates
    figure;
    polarplot(thetas,rhos,'o',...
        'MarkerFaceColor','k','MarkerEdgeColor','None','LineWidth',2); %cmap(5,:)
    
    paH = gca;
    paH.ThetaZeroLocation='left';
    paH.ThetaTickLabel={'max Protraction','','','Retraction','','',...
        'max Retraction','','','Protraction','',''};
    paH.ThetaDir = 'counterclockwise';
    paH.RLim = [0 2];
    %     hold on
    %     legend('','FontSize',8);
    %     legend('boxoff')
    
    
end

%% transition rythm plot
% transition: ISI / spectrogram should go from unimodal to bimodal
% Figure 2 panel: vIRt47_0805_5744 U 8/13/19 (R) 1 (P) Cell#3 : 19
if any(contains(doPlot,'transition'))
    
    switch doPlot{contains(doPlot,'transition')}
        case 'transition'
            load(fullfile(baseDir,'Analysis','Cell_Tuning.mat'));
            transitionfile='Cell_transition.mat';
        case 'transition_Slow'
            load(fullfile(baseDir,'Analysis','Cell_Tuning_SO.mat'));
            transitionfile='Cell_transition_SO.mat';
    end
    %     if exist(fullfile(baseDir,'Analysis',...),'file')
    %         load(fullfile(baseDir,'Analysis',...));
    %     else
    %
    
    %% first compute propEpochCoh and phaseDiffTest to get epoch index
    
    % for each cell, get which has significant different PDF
    phaseDiffTest=struct('globalPhaseDiff',[],'epochPhaseDiffIdx',[]);
    for cellNum=1:numel(tunedCells)
        phaseDiffTest(cellNum).globalPhaseDiff=cellTuning(cellNum).global.phaseStats.spikePhaseStats(1);
        phaseDiffTest(cellNum).epochPhaseDiffIdx=...
            cellfun(@(x) x(1)<=0.05, {cellTuning(cellNum).epochs.phaseStats.spikePhaseStats});
    end
    % for each cell get which epochs has significant coherence
    propEpochCoh=struct('coherEpochIdx',[],'fractionCoherEpoch',[],'manualClass',[]);
    
    for cellNum=1:numel(tunedCells)
        epochCoh=cellfun(@(x,y) any(x>=y), {cellTuning(cellNum).epochs.phaseCoherence.coherMag},...
            {cellTuning(cellNum).epochs.phaseCoherence.confC});
        propEpochCoh(cellNum).fractionCoherEpoch=sum(epochCoh)/numel(epochCoh);
        if cellList.tuningEpochs(cellNum) == 'all' % for comparison with manual classification
            propEpochCoh(cellNum).manualClass=1;
        else
            propEpochCoh(cellNum).manualClass=0;
        end
        propEpochCoh(cellNum).coherEpochIdx=epochCoh;
    end
    
   transition=struct('ISIs',[],'spS',[]);

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
        
        % check whisking epochs first
        
        if ~isfield(wEpochMask,'ampThd')
            ampThd=8; %12; %18 %amplitude threshold
            freqThld=3; %1 %frequency threshold
            minBoutDur=3000; %500; % 1000 % minimum whisking bout duration: 1s
            whiskingEpochs=WhiskingFun.FindWhiskingEpochs(...
                whiskers(bWhisk).amplitude,whiskers(bWhisk).frequency,...
                ampThd, freqThld, minBoutDur);
            whiskingEpochs(isnan(whiskingEpochs))=false; %just in case
            whiskingEpochsList=bwconncomp(whiskingEpochs);
            [~,wBoutDurSort]=sort(cellfun(@length,whiskingEpochsList.PixelIdxList),'descend');
            whiskingEpochsList.PixelIdxListSorted=whiskingEpochsList.PixelIdxList(wBoutDurSort);
            
            
            figure; hold on;
            plot(whiskers(bWhisk).angle);
            plot(wEpochMask.behav*nanstd(whiskers(bWhisk).angle)+nanmean(whiskers(bWhisk).angle))
            plot(whiskingEpochs*nanstd(whiskers(bWhisk).angle)+nanmean(whiskers(bWhisk).angle))
            
            if false
                % vIRt51_1201_5300: ampThd=15; freqThld=3; minBoutDur=3000;  
                % vIRt47_0805_5744: ampThd=10; freqThld=3; minBoutDur=3000;
                % vIRt44_1210_5450: ampThd=18; freqThld=3; minBoutDur=1500;
                % vIRt51_1202_5280: ampThd=14; freqThld=3; minBoutDur=3000;
                % save to wEpochMask
                wEpochMask.ampThd=ampThd;
                wEpochMask.freqThld=freqThld;
                wEpochMask.minBoutDur=minBoutDur;
                diffBE=find(wEpochMask.ephys,1)-find(wEpochMask.behav,1);
                modEMask=false(1,size(wEpochMask.ephys,2));
                modEMask([1:size(wEpochMask.behav,2)]+diffBE)=whiskingEpochs;
                wEpochMask.behav=whiskingEpochs;
                wEpochMask.ephys=modEMask;
                save(fullfile(dataDir,[sessID '_behavior.mat']),'wEpochMask','-append')%,'bWhisk'
            end
        end
%         if ~contains(doPlot,'Slow')
%             wEpochMask.epochIdx=(propEpochCoh(cellNum).coherEpochIdx & phaseDiffTest(cellNum).epochPhaseDiffIdx)';
%         else
            wEpochs=bwconncomp(wEpochMask.behav);
            wEpochMask.epochIdx=true(1,sum(cellfun(@(x) length(x),wEpochs.PixelIdxList)>=3000));
%         end
        
          [transition(cellNum).ISIs,transition(cellNum).spS]=...
              vIRt_transition(whiskers(bWhisk).angle,ephys,wEpochMask,uIdx,true);
    end
    %     end
end


%% Supplementary
% Proba tuning plot: main figure for now - done
% retraction cell video - done
% Slow rhythm tunings: setpoint / breathing

% Additionally, 2 out of the 7 photo-tagged vIRtPV neurons and 10 untagged neurons
% were found to have slower rhythms, with bursting bouts occurring across multiple
% whisking cycles (Extended Data Fig. Xa), potentially linked to whisking midpoint
% and/or the slower breathing rhythm (Extended Data Fig. Xb-Xc).
% PT Slow
% See PT_LO

% Slow Oscillation / putative setpoint Tuning
% vIRt47_0803_5900 U27 srCell 33/40
% vIRt44_1210_5450 U12 srCell 34/40

% See
% Overview Excerpt - vIRt47_0803_5900 - T188 - 194s U27 slow oscillation: great, also PT
% Overview Excerpt - vIRt44_1210_5450 - T124-132s- U12 slow oscillation: also PT
% SetpointTuning - vIRt49_0915_5400 - Unit 25 - late retraction tuning-01
% Overview_Excerpt 57 77 vIRt49_0915_5400 Unit 25
% Overview_Excerpt1_vIRt49_0915_5048 - Unit 58

% Breathing candidates:
% vIRt44_1210_5500
% vIRt44_1210_5100
% vIRt47_0805_5200
% vIRt50_1022_5250
% vIRt50_1022_5450
% vIRt51_1203_5650

% See:
% Overview Excerpt - vIRt44_1210_5100 - T92-94s - U11 U14 opposite breathing rhythm
% Overview excerpt - vIRt50_1022_5250 - T151-167s - nice breathing units
% Overview_Excerpt_vIRt51_1203_5650_Breathing unit
