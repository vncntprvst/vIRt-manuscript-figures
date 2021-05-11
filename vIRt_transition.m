function vIRt_transition(wAngle,ephysData,dataMask)
% burstiness

%% Data masking: look at each whisking epoch
wEpochs.behav=bwconncomp(dataMask.behav);
% mask epochs with short whisking bouts
durationThd=cellfun(@(x) length(x),wEpochs.behav.PixelIdxList)>=3000;

%% keep epochs with significant phase tuning
epochID=find(durationThd);
durationThd(epochID(~dataMask.epochIdx))=false;

% % mask epochs with too much preceding whisking
% keepEpoch=false(numel(epochID),1);
% for epochNum=1:numel(epochID)
%     epochWindow=wEpochs.behav.PixelIdxList{epochID(epochNum)};
%     if epochWindow(1)<501; continue; end
%     if mean(abs(wAngle(epochWindow(1)-501:epochWindow(1)-1)))<1
%         keepEpoch(epochNum)=true;
%     end
% %     mean(abs(wAngle(epochWindow(1):epochWindow(1)+499)))
% end
% durationThd(epochID(~keepEpoch))=false;

dataMask.behav(vertcat(wEpochs.behav.PixelIdxList{~durationThd}))=false;
wEpochs.behav.PixelIdxList=wEpochs.behav.PixelIdxList(durationThd);
wEpochs.behav.NumObjects=sum(durationThd);
% wEpochs.behav.PixelIdxList={vertcat(wEpochs.behav.PixelIdxList{:})};
% wEpochs.behav.NumObjects=1;

% % do the same for ephys data
% wEpochs.ephys=bwconncomp(dataMask.ephys);
% dataMask.ephys(vertcat(wEpochs.ephys.PixelIdxList{~durationThd}))=false;
% wEpochs.ephys.PixelIdxList=wEpochs.ephys.PixelIdxList(durationThd);
% wEpochs.ephys.NumObjects=sum(durationThd);
% % wEpochs.ephys.PixelIdxList={vertcat(wEpochs.ephys.PixelIdxList{:})};
% % wEpochs.ephys.NumObjects=1;
%
spikeRasters = ephysData.rasters(ephysData.selectedUnits,:);
numEpochs=wEpochs.behav.NumObjects;

for unitNum=1:size(spikeRasters,1)  
    %     [mVals,wb_mVals]=deal(cell(numEpochs,1));
      
    unitSpikeEvent=spikeRasters(unitNum,:) ;%wEpochs.ephys.PixelIdxList{wEpochNum});
    ISI=diff(find(unitSpikeEvent))/1000;
    %       figure;  histogram(ISI,30)
    %         wb_ISI=ISI(ISI<=0.015); %within bursts

    spikeTimes=struct('interW',{},'intraW',{});
    for wEpochNum=1:numEpochs
        
        epochIdx=wEpochs.behav.PixelIdxList{wEpochNum};
        if epochIdx(1)<1000 ; continue; end
        
        %compare histogram inter and intra-epoch
        spikeTimes(wEpochNum).interW=find(unitSpikeEvent(epochIdx(1)-3000:epochIdx(1)-1));
        interISI=diff(spikeTimes(wEpochNum).interW);
        spikeTimes(wEpochNum).intraW=find(unitSpikeEvent(epochIdx(1:3000)));
        intraISI=diff(spikeTimes(wEpochNum).intraW);
        %             if numel(interISI)<20 || numel(intraISI)<20; continue; end
        
        %             figure; hold on
        %             histogram(interISI,2:15:150)
        %             histogram(intraISI,2:15:150)
        
        %compute instantaneous ISI
        % check continuous spectrum
        
        
    end
        
    % isi with Chronux
    [interISI_hist,interBins]=isi(rmfield(spikeTimes,'intraW'));
    [intraISI_hist,intraBins]=isi(rmfield(spikeTimes,'interW'));
    
    figure; hold on
    bar(intraBins,intraISI_hist,2)
    bar(interBins,interISI_hist,2)
    legend('intra-whisking ISI','inter-whisking ISI')
end
end


