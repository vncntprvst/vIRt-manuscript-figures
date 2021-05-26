function activation=vIRt_activation(whisker,ephysData,dataMask,uIdx,savePlots)
% activation plots

if nargin < 5
    savePlots =false;
end
fileName=ephysData.recInfo.sessionName;

%% Data masking: look at each whisking epoch
wEpochs.behav=bwconncomp(dataMask.behav);
% mask epochs with short whisking bouts
durationThd=cellfun(@(x) length(x),wEpochs.behav.PixelIdxList)>=3000;

%% keep epochs with significant phase tuning (or just duration threshold if Slow)
epochID=find(durationThd);
durationThd(epochID(~dataMask.epochIdx))=false;

dataMask.behav(vertcat(wEpochs.behav.PixelIdxList{~durationThd}))=false;
wEpochs.behav.PixelIdxList=wEpochs.behav.PixelIdxList(durationThd);
wEpochs.behav.NumObjects=sum(durationThd);
% wEpochs.behav.PixelIdxList={vertcat(wEpochs.behav.PixelIdxList{:})};
% wEpochs.behav.NumObjects=1;

% % do the same for ephys data
wEpochs.ephys=bwconncomp(dataMask.ephys);
dataMask.ephys(vertcat(wEpochs.ephys.PixelIdxList{~durationThd}))=false;
wEpochs.ephys.PixelIdxList=wEpochs.ephys.PixelIdxList(durationThd);
wEpochs.ephys.NumObjects=sum(durationThd);
% wEpochs.ephys.PixelIdxList={vertcat(wEpochs.ephys.PixelIdxList{:})};
% wEpochs.ephys.NumObjects=1;

spikeRasters = ephysData.rasters(ephysData.selectedUnits,:);
numEpochs=wEpochs.behav.NumObjects;
[activation,activationEpochs]=deal(struct('FR',struct('pre',[],'post',[]),'WF',struct('pre',[],'post',[])));

for unitNum=1:size(spikeRasters,1)   
    unitSpikeEvent=spikeRasters(unitNum,:) ;
    
    for wEpochNum=1:numEpochs
        eEpochIdx=wEpochs.ephys.PixelIdxList{wEpochNum};
        wEpochIdx=wEpochs.behav.PixelIdxList{wEpochNum};
        if eEpochIdx(1)<1000 || wEpochIdx(1)<3000; continue; end     
        wPreEpoch=wEpochIdx(1)-1000:wEpochIdx(1)-1;
        wPostEpoch=wEpochIdx(1:1000);
        
        % get pre activation values, keeping only quite periods
        nonQuiteIdx=whisker.amplitude(wPreEpoch)>6;
        if all(nonQuiteIdx); continue; end
        activationEpochs(wEpochNum,unitNum).FR.pre=unitSpikeEvent(eEpochIdx(1)-1000:eEpochIdx(1)-1);
        activationEpochs(wEpochNum,unitNum).FR.pre=activationEpochs(wEpochNum,unitNum).FR.pre(~nonQuiteIdx);
        activationEpochs(wEpochNum,unitNum).WF.pre=whisker.setPoint(wPreEpoch);
        activationEpochs(wEpochNum,unitNum).WF.pre=activationEpochs(wEpochNum,unitNum).WF.pre(~nonQuiteIdx);
        % get post activation values
        activationEpochs(wEpochNum,unitNum).FR.post=unitSpikeEvent(eEpochIdx(1:1000));
        activationEpochs(wEpochNum,unitNum).WF.post=whisker.setPoint(wPostEpoch);
    end
    
    emptyEpochs=cellfun(@isempty, {activationEpochs.FR});  
    activationEpochs=activationEpochs(~emptyEpochs,:);
    
    activFR=cellfun(@(x) [sum(x.pre)/numel(x.pre)*1000 sum(x.post)],{activationEpochs.FR}','UniformOutput',false);
    activFR=mean(vertcat(activFR{:}));
    activWF=cellfun(@(x) [mean(x.pre) mean(x.post)],{activationEpochs.WF}','UniformOutput',false);
    activWF=mean(vertcat(activWF{:}));
       
    activation.FR.pre=activFR(1);
    activation.FR.post=activFR(2);
    activation.WF.pre=activWF(1);
    activation.WF.post=activWF(2);
          
    if savePlots
        figDir='D:\Vincent\Figures\vIRt';
        savefig(gcf,fullfile(figDir, 'activation', ['Spectrogram - Cell' num2str(uIdx) fileName ' - Unit' num2str(unitNum) ' - all whisking epochs average']));
        print(gcf,fullfile(figDir, 'activation', ['Spectrogram - Cell' num2str(uIdx) fileName ' - Unit' num2str(unitNum) ' - all whisking epochs average']),'-dpng');
        close(gcf)
    end
end
end

