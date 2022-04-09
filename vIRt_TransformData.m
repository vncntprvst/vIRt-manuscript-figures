function [whiskers,bWhisk,breathing,ephys]=vIRt_TransformData(behav,ephys)

whiskers=behav.whiskers;
try 
    bWhisk=behav.whiskerTrackingData.keepWhiskerIDs==behav.whiskerTrackingData.bestWhisker; %best whisker
catch
    bWhisk=find([whiskers.bestWhisker]);
end

%% compute whisking frequency (different from instantaneous frequency
for wNum=1:numel(bWhisk)
    whisksIdx = bwconncomp(diff(whiskers(bWhisk(wNum)).phase)>0);
    peakIdx = zeros(1,length(whiskers(bWhisk(wNum)).velocity));
    peakIdx(cellfun(@(whisk) whisk(1), whisksIdx.PixelIdxList))=1;
    whiskers(bWhisk(wNum)).frequency=movsum(peakIdx,behav.whiskerTrackingData.samplingRate);
end

%% other behavior data
if ~isempty(behav.breathing)
    breathing.data=double(behav.breathing');
    breathing.data=breathing.data*(range(whiskers(bWhisk(1)).setPoint)/range(breathing.data));
    breathing.ts=linspace(0,ephys.recInfo.duration_sec,numel(behav.breathing));
else
    breathing=[];
end

%% compute rasters
% aim for same length for ephys traces and behavior data
[ephys.rasters,unitList]=EphysFun.MakeRasters(ephys.spikes.times,ephys.spikes.unitID,...
    1,size(whiskers(bWhisk(1)).angle,2)); %int32(size(ephys.traces,2)/ephys.spikes.samplingRate*1000));


%% compute spike density functions
ephys.spikeRate=EphysFun.MakeSDF(ephys.rasters);

%%create timeline
ephys.timestamps = 0:0.001:ephys.recInfo.duration_sec;
