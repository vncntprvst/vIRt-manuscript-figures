function [thetas,phaseStats,phaseTuning]=vIRt_PhaseTuning(whiskerPhase,ephysData,dataMask,splitEpochs)

% if isfield(ephysData.recInfo,'sessionName'); recName=ephysData.recInfo.sessionName;
% else; recName='PhaseTuning_polarPlot'; end

%% Data masking: look at each whisking epoch
wEpochs.behav=bwconncomp(dataMask.behav);
% mask epochs with short whisking bouts
durationThd=cellfun(@(x) length(x),wEpochs.behav.PixelIdxList)>=3000;
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

for unitNum=1:size(spikeRasters,1)
    
    if mean(spikeRate(unitNum,:)) < 0.5
        continue
    end
    
    numEpochs=wEpochs.behav.NumObjects;
    phaseTuning=nan(1,numEpochs);
    thetas=cell(1,numEpochs);
    for wEpochNum=1:numEpochs
        clearvars eWhiskerPhase
        eWhiskerPhase=whiskerPhase(wEpochs.behav.PixelIdxList{wEpochNum});
        
        %% probability density function of phase for spiking events
        numBins=32; % each bin = pi/16 radians
        edges = linspace(min(eWhiskerPhase), max(eWhiskerPhase), numBins*2+1);
        %         edges = linspace(-pi-pi/numBins,pi+pi/numBins, numBins+1);
        centers = mean([ edges(1:end-1); edges(2:end) ]);
        [ ~, ~, phaseBins ] = histcounts(eWhiskerPhase, edges);
        samplingRate=1000; %change in case this isn't at 1kHz SR

        try
            unitSpikeEvent=spikeRasters(unitNum,wEpochs.ephys.PixelIdxList{wEpochNum});
            attribPhaseBin = phaseBins(logical(unitSpikeEvent));
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
            %         phasePDF=movsum(phasePDF,2);phasePDF=phasePDF(2:2:end);
            % mean spike rate for each phase bin ?[?k] = SR*N(?k|spike)/N(?k)
            meanPhaseSpikeRate=samplingRate*spikePhaseBinCount./phaseBinCount;
            %     fit sine wave
            %     modulation depth of the averaged whisking response
            
            
            clearvars unitSpikeRate
            eSpikeRate= spikeRate(unitNum,wEpochs.ephys.PixelIdxList{wEpochNum});
            
            clearvars binMeanSpikeRate binSESpikeRate
            
            % Bining firing rate and spikes
            for binNum = 1:length(edges)-1 % : -1 : 1
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
            rsbinMeanSpikeRate=sum(reshape([binMeanSpikeRate(2:end),binMeanSpikeRate(1)],2,numBins))/2;
            rsbinMeanSpikeRate=[rsbinMeanSpikeRate(end) rsbinMeanSpikeRate];
            
            %% convert to thetas: make as many phase # as FR for that phase #
            thetas{wEpochNum}=cell(numel(centers),1);
            for binNum=1:numel(centers)
                thetas{wEpochNum}{binNum}=ones(round(binMeanSpikeRate(binNum)),1)*centers(binNum);
            end
            thetas{wEpochNum}=vertcat(thetas{wEpochNum}{:});
            if isempty(thetas{wEpochNum})
                disp('not enough spikes')
                continue
            end
            % stats
            phaseStats(wEpochNum)=circ_stats(thetas{wEpochNum});
            if  circ_rtest(thetas{wEpochNum})<0.05 %((phaseStats.kurtosis>0.04 || phaseStats.skewness<-0.02) || ...
                phaseTuning(wEpochNum)=rad2deg(phaseStats(wEpochNum).median);
            end
            
        catch
            continue
        end
    end
end
