function [whiskers,bWhisk,breathing,ephys]=vIRt_TransformData(behav,ephys)

whiskers=behav.whiskers;
bWhisk=behav.whiskerTrackingData.keepWhiskerIDs==behav.whiskerTrackingData.bestWhisker; %best whisker

%% compute whisking frequency (different from instantaneous frequency
whisksIdx = bwconncomp(diff(whiskers(bWhisk).phase)>0);
peakIdx = zeros(1,length(whiskers(bWhisk).velocity));
peakIdx(cellfun(@(whisk) whisk(1), whisksIdx.PixelIdxList))=1;
whiskers(bWhisk).frequency=movsum(peakIdx,behav.whiskerTrackingData.samplingRate);

%% other behavior data
if ~isempty(behav.breathing)
    breathing.data=double(behav.breathing');
    breathing.data=breathing.data*(range(whiskers(bWhisk).setPoint)/range(breathing.data));
    breathing.ts=linspace(0,ephys.recInfo.duration_sec,numel(behav.breathing));
else
    breathing=[];
end

%% compute rasters
% aim for same length for ephys traces and behavior data
[ephys.rasters,ephys.unitList]=EphysFun.MakeRasters(ephys.spikes.times,ephys.spikes.unitID,...
    1,size(whiskers(bWhisk).angle,2)); %int32(size(ephys.traces,2)/ephys.spikes.samplingRate*1000));
%% compute spike density functions
ephys.spikeRate=EphysFun.MakeSDF(ephys.rasters);

%%create timeline
ephys.timestamps = 0:0.001:ephys.recInfo.duration_sec;
