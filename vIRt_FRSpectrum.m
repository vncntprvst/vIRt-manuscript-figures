function spS=vIRt_FRSpectrum(spikeTimes,dataMask) %

%% Make rasters and compute spike density function
rasters=EphysFun.MakeRasters(spikeTimes,ones(numel(spikeTimes),1),1);
spikeRate=EphysFun.MakeSDF(rasters,3);

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
    wEpochs.behav.PixelIdxList=vertcat(wEpochs.behav.PixelIdxList{1:timeLimitIdx});
    % apply mask
    spikeRate=spikeRate(wEpochs.behav.PixelIdxList);
else % keep first 30s
    %     spikeTimes=spikeTimes(spikeTimes<=30);
end

%% Power spectrum
% set parameters
params.Fs=1000; % sampling frequency
params.fpass=[3 25]; % % band of frequencies to be kept
params.NW=min([20 length(spikeRate)*2.5]); %3
params.tapers=[params.NW params.NW*2-1]; % taper parameters
params.pad=2; % pad factor for fft
params.err=[2 0.01];
params.trialave=0;
movingwin=[0.5 0.05];

% tic
[spS.spectrumVals,spS.freqVals,spS.Serr]=mtspectrumc(spikeRate',params);
% toc

% [S,t,f,Serr]=mtspecgramc(spikeRate',movingwin,params)

% convert to dB (power spectral density)
spS.spectrumValsPSD= 10*log10(spS.spectrumVals);
spS.SerrPSD(1,:)=10*log10(spS.Serr(1,:));
spS.SerrPSD(2,:)=10*log10(spS.Serr(2,:));
spS.RPSD=10*log10(spS.R);

spS.StatSigIdx=spS.SerrPSD(1,:)<=spS.RPSD & spS.SerrPSD(2,:)>=spS.RPSD;

%% figures
if false
    figure;
    plot(spS.freqVals,spS.spectrumVals)
    
    figure
    plot(spS.freqVals,10*log10(spS.spectrumVals),spS.freqVals,10*log10(spS.Serr(1,:)),...
        spS.freqVals,10*log10(spS.Serr(2,:)));
    line(get(gca,'xlim'),[10*log10(spS.R) 10*log10(spS.R)]);
    hold on
    plot(spS.freqVals(spS.StatSigIdx),spS.spectrumValsPSD(spS.StatSigIdx),'k')
    
%     figure;imagesc(t,f,10*log10(S)'); axis xy; colorbar

end

