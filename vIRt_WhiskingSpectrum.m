function wS=vIRt_WhiskingSpectrum(whiskingTrace,dataMask)

whiskingTrace=whiskingTrace-mean(whiskingTrace);

%% Data masking
if nargin>1
    wEpochs.behav=bwconncomp(dataMask.behav);
    % mask epochs with short whisking bouts
    durationThd=cellfun(@(x) length(x),wEpochs.behav.PixelIdxList)>=3000;
    dataMask.behav(vertcat(wEpochs.behav.PixelIdxList{~durationThd}))=false;
    wEpochs.behav.PixelIdxList=wEpochs.behav.PixelIdxList(durationThd);
    cumulDur=cumsum(cellfun(@numel, wEpochs.behav.PixelIdxList)/1000);
    timeLimitIdx=find(cumulDur>=30,1); %keep only 30s or so of wisking
    if isempty(timeLimitIdx) %then keep all
        timeLimitIdx=numel(wEpochs.behav.PixelIdxList);
    end
    wEpochs.behav.PixelIdxList=vertcat(wEpochs.behav.PixelIdxList{1:timeLimitIdx});
    % apply mask
    whiskingTrace=whiskingTrace(wEpochs.behav.PixelIdxList);
end

%% Power spectrum
% set parameters
params.Fs=1000; % sampling frequency
params.fpass=[0 30]; % band of frequencies to be kept
params.NW=min([floor(numel(whiskingTrace)/1000) 50])*3;
params.tapers=[params.NW params.NW*2-1]; % taper parameters
params.pad=0; % pad factor for fft
params.err=[2 0.05];
if size(whiskingTrace,1)>1
   params.trialave=1;
else
   params.trialave=0;
end

try
    [wS.spectrumVals,wS.freqVals]=mtspectrumc(whiskingTrace',params);
    wS.spectrumVals=wS.spectrumVals';
    wS.peakFreq=wS.freqVals(wS.spectrumVals==max(wS.spectrumVals));
catch
    [wS.spectrumVals,wS.freqVals,wS.peakFreq]=deal(NaN);
end


%  figure; hold on
%  plot(whiskingTrace')
%  plot(wS.spectrumVals')
% 
%  plot(wS.freqVals,wS.spectrumVals)
%  plot(wS.freqVals,smooth(wS.spectrumVals,100))

