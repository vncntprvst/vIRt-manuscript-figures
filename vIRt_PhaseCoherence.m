function pCoh=vIRt_PhaseCoherence(whiskerPhase,unitRasters,mvtType)

if nargin<3
    params.fpass=[3 30]; % % band of frequencies to be kept
else
    switch mvtType
        case 'whisking'
            params.fpass=[6 25];
        case 'fast_whisking'
            params.fpass=[15 20];
        case 'slow_whisking'
            params.fpass=[6 16];
    end
end            

whiskerPhase=whiskerPhase-mean(whiskerPhase);

% figure; 
% phasePH=subplot(2,1,1); plot(whiskerPhase); axis tight
% rasterPH=subplot(2,1,2); imagesc(unitRasters)
% linkaxes([phasePH rasterPH],'x')

params.Fs=1000; % sampling frequency
params.NW=floor(numel(whiskerPhase)/1000)*2.5;
params.tapers=[params.NW params.NW*2-1]; % taper parameters 
params.pad=0; % pad factor for fft
params.err=[2 0.01];
params.trialave=0;

[pCoh.coherMag,pCoh.coherPhase,~,~,~,pCoh.freqVals,~,...
    pCoh.confC,pCoh.phistd,pCoh.Cerr]=coherencycpb(...
    whiskerPhase',unitRasters',params); % pCoh.crossSpectrum pCoh.wPhaseSpectrum pCoh.spikeSpectrum


% plotsig(pCoh.coherMag,pCoh.confC,t,pCoh.freqVals);

% O = ebi([pCoh.wPhaseSpectrum;pCoh.spikeSpectrum],...
%     [zeros(numel(pCoh.wPhaseSpectrum),1);ones(numel(pCoh.wPhaseSpectrum),1)],...
%     'Test','CohZ');
pCoh.peakCoherMag=pCoh.coherMag(find(pCoh.coherMag==max(pCoh.coherMag),1));
pCoh.peakCoherPhase=pCoh.coherPhase(find(pCoh.coherMag==max(pCoh.coherMag),1));

% figure;
% plot(pCoh.freqVals,pCoh.coherMag)
% figure;
% plot(pCoh.freqVals,pCoh.coherPhase)

% figure;
% plot(pCoh.freqVals,pCoh.Cerr)

% tested trial concat: doesn't work
% numTrial=floor(numel(whiskerPhase)/3000);
% [whiskerPhaseTrials,unitRastersTrials]=deal(nan(3000,numTrial));
% for trialNum=1:numTrial
%     whiskerPhaseTrials(:,trialNum)=whiskerPhase((trialNum-1)*3000+1:trialNum*3000);
%     unitRastersTrials(:,trialNum)=unitRasters((trialNum-1)*3000+1:trialNum*3000);
% end
% params.trialave=0;
% 
% [pCoh.coherMag,pCoh.coherPhase,~,~,~,pCoh.freqVals,~,confC,phistd,Cerr]=coherencycpb(...
%     whiskerPhaseTrials,unitRastersTrials,params);
% 
% figure;
% plot(pCoh.freqVals,pCoh.coherMag)
% figure;
% plot(pCoh.freqVals,pCoh.coherPhase)
% 
% figure;
% plot(pCoh.freqVals,pCoh.Cerr)
