function [thetas,edges,phaseStats,phaseTuning,phaseCoherence,rsbinMeanSpikeRate]=vIRt_PhaseTuning(whiskerPhase,ephysData,dataMask,splitEpochs)

%% Data masking: look at each whisking epoch
wEpochs.behav=bwconncomp(dataMask.behav);
% mask epochs with short whisking bouts
durationThd=cellfun(@(x) length(x),wEpochs.behav.PixelIdxList)>=3000;
if ~any(durationThd)
    durationThd=cellfun(@(x) length(x),wEpochs.behav.PixelIdxList)>=1000;
end
dataMask.behav(vertcat(wEpochs.behav.PixelIdxList{~durationThd}))=false;
wEpochs.behav.PixelIdxList=wEpochs.behav.PixelIdxList(durationThd);
if splitEpochs
    wEpochs.behav.NumObjects=sum(durationThd);
else
    wEpochs.behav.PixelIdxList={vertcat(wEpochs.behav.PixelIdxList{:})};
    wEpochs.behav.NumObjects=1;
end
% do the same for ephys data
wEpochs.ephys=bwconncomp(dataMask.ephys);
dataMask.ephys(vertcat(wEpochs.ephys.PixelIdxList{~durationThd}))=false;
wEpochs.ephys.PixelIdxList=wEpochs.ephys.PixelIdxList(durationThd);
if splitEpochs
    wEpochs.ephys.NumObjects=sum(durationThd);
else
    wEpochs.ephys.PixelIdxList={vertcat(wEpochs.ephys.PixelIdxList{:})};
    wEpochs.ephys.NumObjects=1;
end

spikeRasters = ephysData.rasters(ephysData.selectedUnits,:);
spikeRate=ephysData.spikeRate(ephysData.selectedUnits,:);
phaseStats=struct('mean',[],'median',[],'var',[],'std',[],...
    'std0',[],'skewness',[],'skewness0',[],'kurtosis',[],...
    'kurtosis0',[]);
phaseStatsPDF=struct('spikePhaseIdx',[],'whiskerPhase',[],'spikePhaseBins',[],'phaseBins',[],'spikePhasePDF',[],...
    'phasePDF',[],'spikePhaseStats',[]);
phaseCoherence=struct('coherMag',[],'coherPhase',[],'freqVals',[],...
    'peakCoherMag',[],'peakCoherPhase',[],'confC',[],'phistd',[],'Cerr',[]);
[thetas,edges]=deal(cell(size(spikeRasters,1),1));
    
for unitNum=1:size(spikeRasters,1)
    
    if mean(spikeRate(unitNum,:)) < 0.5
        continue
    end
    
    numEpochs=wEpochs.behav.NumObjects;
    phaseTuning=nan(1,numEpochs);

    for wEpochNum=1:numEpochs
        clearvars eWhiskerPhase
        eWhiskerPhase=whiskerPhase(wEpochs.behav.PixelIdxList{wEpochNum});
        if isempty(eWhiskerPhase); continue; end
        %% probability density function of phase for spiking events
        numBins=32; % each bin = pi/16 radians
        edges{unitNum,wEpochNum} = linspace(min(eWhiskerPhase), max(eWhiskerPhase), numBins*2+1);
        %         edges{unitNum,wEpochNum} = linspace(-pi-pi/numBins,pi+pi/numBins, numBins+1);
        centers = mean([ edges{unitNum,wEpochNum}(1:end-1); edges{unitNum,wEpochNum}(2:end) ]);
        [~,~, phaseBins ] = histcounts(eWhiskerPhase, edges{unitNum,wEpochNum});
        samplingRate=1000; %change in case this isn't at 1kHz SR
        
        try
            unitSpikeEvent=logical(spikeRasters(unitNum,wEpochs.ephys.PixelIdxList{wEpochNum}));
            attribPhaseBin = phaseBins(unitSpikeEvent);
            %     number of spikes in each phase bin N(?k|spike)
            spikePhaseBinCount=histcounts(attribPhaseBin,[1 1+unique(phaseBins)]);%         [spikePhaseBinCount,uniqueSpikePhaseBins]=hist(phaseVals,unique(phaseVals));
            %     probability density function P(?k|spike)
            spikePhasePDF=spikePhaseBinCount/sum(spikePhaseBinCount);
            spikePhasePDF=sum(reshape([spikePhasePDF(2:end),spikePhasePDF(1)],2,numBins))/2;
            spikePhasePDF=[spikePhasePDF(end) spikePhasePDF];
            %         spikePhasePDF=movsum(spikePhasePDF,6);spikePhasePDF=spikePhasePDF(6:6:end);
            %     number of phase occurence for each phase bin N(?k)
            phaseBinCount=histcounts(phaseBins,[1 1+unique(phaseBins)]);
            %     probability density function P(?k)
            phasePDF=phaseBinCount/sum(phaseBinCount);
            phasePDF=sum(reshape([phasePDF(2:end),phasePDF(1)],2,numBins))/2;
            phasePDF=[phasePDF(end) phasePDF];
            
            try
                [kuiperstest_p, kuiperstest_k, kuiperstest_K] = circ_kuipertest(eWhiskerPhase, eWhiskerPhase(unitSpikeEvent));
            catch
                [kuiperstest_p, kuiperstest_k, kuiperstest_K]=deal(NaN);
            end
            %keep data
            phaseStatsPDF(wEpochNum).spikePhaseBins=attribPhaseBin;
            phaseStatsPDF(wEpochNum).phaseBins=phaseBins;
            phaseStatsPDF(wEpochNum).spikePhasePDF=spikePhasePDF;
            phaseStatsPDF(wEpochNum).phasePDF=phasePDF;
            phaseStatsPDF(wEpochNum).whiskerPhase=eWhiskerPhase;
            phaseStatsPDF(wEpochNum).spikePhaseIdx=unitSpikeEvent;
            phaseStatsPDF(wEpochNum).spikePhaseStats=[kuiperstest_p, kuiperstest_k, kuiperstest_K];
            
            %         phasePDF=movsum(phasePDF,2);phasePDF=phasePDF(2:2:end);
            % mean spike rate for each phase bin ?[?k] = SR*N(?k|spike)/N(?k)
            meanPhaseSpikeRate=samplingRate*spikePhaseBinCount./phaseBinCount;
            %     fit sine wave
            %     modulation depth of the averaged whisking response
            
            phaseCoherence(wEpochNum)=vIRt_PhaseCoherence(eWhiskerPhase(1:min([90000 end])),...
                unitSpikeEvent(1:min([90000 end]))); %limit to 30s
%             phaseCoherence(wEpochNum)=vIRt_PhaseCoherence(eWhiskerPhase,unitSpikeEvent);

            clearvars unitSpikeRate
            eSpikeRate= spikeRate(unitNum,wEpochs.ephys.PixelIdxList{wEpochNum});
            
            clearvars binMeanSpikeRate binSESpikeRate
            
            % Bining firing rate and spikes
            for binNum = 1:length(edges{unitNum,wEpochNum})-1 % : -1 : 1
                %                     chunkSpikeRate=unitSpikeRate(chunkIndex);
                ratesVect= eSpikeRate(phaseBins == binNum);%(chunkIndex)
                numSample = numel(ratesVect);
                if numSample == 0
                    meanSpikeRate = 0;
                    steSpikeRate = 0;
                else
                    meanSpikeRate = nanmean(ratesVect);
                    steSpikeRate = MMath.StandardError(ratesVect);
                end
                binMeanSpikeRate(binNum) = meanSpikeRate;
                binSESpikeRate(binNum) = steSpikeRate;
            end
            
            %           rescale
            rsbinMeanSpikeRate{unitNum,wEpochNum}=sum(reshape([binMeanSpikeRate(2:end),binMeanSpikeRate(1)],2,numBins))/2;
            rsbinMeanSpikeRate{unitNum,wEpochNum}=[rsbinMeanSpikeRate{unitNum,wEpochNum}(end) rsbinMeanSpikeRate{unitNum,wEpochNum}];
            
            %% convert to thetas: make as many phase # as FR for that phase #
            thetas{unitNum,wEpochNum}=cell(numel(centers),1);
            for binNum=1:numel(centers)
                thetas{unitNum,wEpochNum}{binNum}=ones(round(binMeanSpikeRate(binNum)),1)*centers(binNum);
            end
            thetas{unitNum,wEpochNum}=vertcat(thetas{unitNum,wEpochNum}{:});
            if isempty(thetas{unitNum,wEpochNum})
                disp('not enough spikes')
                continue
            end
            % stats
            phaseStats(wEpochNum)=circ_stats(thetas{unitNum,wEpochNum});
            if  circ_rtest(thetas{unitNum,wEpochNum})<0.05 %((phaseStats.kurtosis>0.04 || phaseStats.skewness<-0.02) || ...
                phaseTuning(wEpochNum)=rad2deg(phaseStats(wEpochNum).median);
            end            
        catch
            continue
        end
    end
    
    %% check figure
    if false 
        labels = 'whisking';
        cmap=lines;cmap=[cmap(1:7,:);(lines+flipud(copper))/2;autumn];
        if isfield(ephysData.recInfo,'sessionName'); recName=ephysData.recInfo.sessionName;
        else; recName='PhaseTuning_polarPlot'; end

        figure('name',['Unit ' num2str(ephysData.selectedUnits(unitNum)) ' - ' recName ...
            ' - Tuning to ' labels ' phase'],'Color','white','position',...
            [1059 106 355 935]);
        %         end
        
        %         %         ptWhisks=bwconncomp(eWhiskerPhase>0);
        %         %% probability density function of phase for spiking events
        %         numBins=32; % each bin = pi/16 radians
        %         edges{unitNum,wEpochNum} = linspace(min(eWhiskerPhase), max(eWhiskerPhase), numBins*2+1);
        %         %         edges{unitNum,wEpochNum} = linspace(-pi-pi/numBins,pi+pi/numBins, numBins+1);
        %         centers = mean([ edges{unitNum,wEpochNum}(1:end-1); edges{unitNum,wEpochNum}(2:end) ]);
        %         [ ~, ~, phaseBins ] = histcounts(eWhiskerPhase, edges{unitNum,wEpochNum});
        %         samplingRate=1000; %change in case this isn't at 1kHz SR
        %         phaseTuning=nan(size(spikeRate,1),numEpochs);
        
        %% polar plot
        if ~isnan(phaseTuning(wEpochNum))
            phEdgeColor=cmap(unitNum,:);phFaceColor=cmap(unitNum,:);
        else
            phEdgeColor='k';phFaceColor='k'; %EdgeAlpha=0.5;
        end
        
        subplot(3,1,1);
        polarhistogram(thetas{unitNum,wEpochNum},edges{unitNum,wEpochNum},'Displaystyle','bar',...
            'Normalization','count','LineWidth',2,...
            'EdgeColor',phEdgeColor,'FaceColor',phFaceColor,...
            'EdgeAlpha',0);
        paH = gca;
        paH.ThetaZeroLocation='left';
        paH.ThetaTickLabel={'max Protraction','','','Retraction','','',...
            'max Retraction','','','Protraction','',''};
        paH.ThetaDir = 'counterclockwise';
        
        title(['Unit ' num2str(ephysData.selectedUnits(unitNum)) ' - ' strrep(recName,'_','') ...
            ' - Tuning to ' labels ' phase'],'interpreter','none');
        
        %% plot PDF across phase
        subplot(3,1,2); hold on ;
        plot(linspace(-pi,pi, numBins+1),spikePhasePDF,'linewidth',1.2,'Color', [0 0 0]); %centers
        plot(linspace(-pi,pi, numBins+1),phasePDF,'linewidth',1.2,'Color', [0 0 0 0.5]); %centers
        set(gca,'ytick',0:0.05:1,...
            'xlim',[-pi pi],'xtick',[-pi 0 pi],'xticklabel',{'-\pi','0','\pi'},...
            'tickdir','out');
        axis tight
        legend('P(\phi_k|spike)','P(\phi_k)','location','southeast')
        legend('boxoff')
        title({'Probability density function'; 'of phase for spiking events'})
        
        %% average firing rate across phase
        subplot(3,1,3);hold on;
        plot(linspace(-pi,pi, numBins+1), rsbinMeanSpikeRate{unitNum,wEpochNum}, 'LineWidth',2) %centers %,'color',cmap(unitNum,:));%'k'
        
        axis tight
        yl = ylim;
        ylim([0 yl(2)]);
        text(0,  yl(2)/10, '0 = Max protraction')
        set(gca,'xlim',[-pi pi],'xtick',[-pi 0 pi],'xticklabel',{'-\pi','0','\pi'},...
            'tickdir','out');
        box off
        if wEpochs.behav.NumObjects==1
            xlabel('Phase (radians)')
            ylabel('Spike rate (Hz)')
        end
        title({'Average spike rate'; ['across ' labels ' phase']})
        % set(gca,'xdir', 'reverse'); %, 'ydir', 'reverse')
    end
end

phaseStats=CatStruct(phaseStats,phaseStatsPDF);
