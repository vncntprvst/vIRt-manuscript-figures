function spS=vIRt_SpikingSpectrum_Epochs(spikeTimes,dataMask) %

if nargin>1
    %% Data masking
    wEpochs.behav=bwconncomp(dataMask.behav);
    % mask epochs with short whisking bouts
    durationThd=cellfun(@(x) length(x),wEpochs.behav.PixelIdxList)>=2000;
    dataMask.behav(vertcat(wEpochs.behav.PixelIdxList{~durationThd}))=false;
    wEpochs.behav.PixelIdxList=wEpochs.behav.PixelIdxList(durationThd);
    cumulDur=cumsum(cellfun(@numel, wEpochs.behav.PixelIdxList)/1000);
    timeLimitIdx=find(cumulDur>=30,1); %keep only 30s or so of wisking
    if isempty(timeLimitIdx) %then keep all
        timeLimitIdx=numel(wEpochs.behav.PixelIdxList);
    end
    
    % apply mask
    %     wEpochs.behav.PixelIdxList=vertcat(wEpochs.behav.PixelIdxList{1:timeLimitIdx});
    %     %% conversion to structure
    %     if ~isstruct(spikeTimes)
    %         spikeTimes=struct('times',spikeTimes);
    %     end
    %     for epochNum=numel(wEpochs.behav.PixelIdxList):-1:1
    %         epochIdx=spikeTimes(1).times>=wEpochs.behav.PixelIdxList{epochNum}(1)/1000 & ...
    %             spikeTimes(1).times<=wEpochs.behav.PixelIdxList{epochNum}(end)/1000;
    %         spikeTimes(epochNum).times=spikeTimes(1).times(epochIdx);
    %     end
    
    if false % code below to apply mask
        wEpochs.behav.PixelIdxList=wEpochs.behav.PixelIdxList(1:timeLimitIdx);
        for epochNum=1:numel(wEpochs.behav.PixelIdxList)-1
            interEpochIdx=spikeTimes>wEpochs.behav.PixelIdxList{epochNum}(end)/1000 & ...
                spikeTimes<wEpochs.behav.PixelIdxList{epochNum+1}(1)/1000;
            spikeTimes(interEpochIdx)=NaN;
        end
        
        spikeTimes(spikeTimes>wEpochs.behav.PixelIdxList{end}(end)/1000)=NaN;
        spikeTimes=spikeTimes(~isnan(spikeTimes));
    else % keep epochs separate
        epochList=wEpochs.behav.PixelIdxList(1:timeLimitIdx);
    end
else % keep first 30s
    %     spikeTimes=spikeTimes(spikeTimes<=30);
    epochList={spikeTimes};
end

%% Power spectrum
% set parameters
params.Fs=1000; % sampling frequency
params.fpass=[3 30]; % % band of frequencies to be kept
params.NW=min([100 floor(max(spikeTimes))*2.5]); %3
params.tapers=[params.NW params.NW*2-1]; % taper parameters
params.pad=2; % pad factor for fft
params.err=[2 0.01];
params.trialave=0;
movingwin=[0.5 0.05];

spS=struct('spectrumVals',[],'freqVals',[],'R',[],'Serr',[],...
    'spectrumValsPSD',[],'SerrPSD',[],'RPSD',[],'StatSigIdx',[]);

for epochNum=1:numel(epochList)
    epochIdx=spikeTimes>=epochList{epochNum}(1)/1000 & ...
        spikeTimes<=epochList{epochNum}(end)/1000;
    spkTimes=spikeTimes(epochIdx)-epochList{epochNum}(1)/1000;
    
    % tic
    [spS(epochNum).spectrumVals,spS(epochNum).freqVals,spS(epochNum).R,...
        spS(epochNum).Serr]=mtspectrumpt(spkTimes,params);
    % toc
    % [spS(epochNum).spectrumVals,spS(epochNum).times,spS(epochNum).freqVals,spS(epochNum).R,spS(epochNum).Serr]=mtspecgrampt_optimized(spikeTimes,movingwin,params);
    
    % convert to dB (power spectral density)
    spS(epochNum).spectrumValsPSD= 10*log10(spS(epochNum).spectrumVals);
    spS(epochNum).SerrPSD(1,:)=10*log10(spS(epochNum).Serr(1,:));
    spS(epochNum).SerrPSD(2,:)=10*log10(spS(epochNum).Serr(2,:));
    spS(epochNum).RPSD=10*log10(spS(epochNum).R);
    
    spS(epochNum).StatSigIdx=spS(epochNum).SerrPSD(1,:)<=spS(epochNum).RPSD &...
        spS(epochNum).SerrPSD(2,:)>=spS(epochNum).RPSD;
end

%% figures
if false
    figure;
    plot(spS(epochNum).freqVals,spS(epochNum).spectrumVals)
    
    for epochNum=1:numel(epochList)
        figure; hold on
        plot(spS(epochNum).freqVals,10*log10(spS(epochNum).spectrumVals),spS(epochNum).freqVals,10*log10(spS(epochNum).Serr(1,:)),...
            spS(epochNum).freqVals,10*log10(spS(epochNum).Serr(2,:)));
        line(get(gca,'xlim'),[10*log10(spS(epochNum).R) 10*log10(spS(epochNum).R)]);
    end
    plot(spS(epochNum).freqVals(spS(epochNum).StatSigIdx),spS(epochNum).spectrumValsPSD(spS(epochNum).StatSigIdx),'k')
    
end

