function [ephys,knownTS,waveForms]=S2sort(ephys,oTS,pulses,s2Sort)
% see benchmark_and_fix to export traces and save S2 sort

spikeTimes=single(ephys.spikes.times);
TTLtimes=pulses.TTLTimes;
duration=pulses.duration;
pulseIdx=false(size(spikeTimes,1),size(TTLtimes,1));
for TTLNum=1:length(TTLtimes)
    pulseIdx(:,TTLNum)=spikeTimes>=TTLtimes(TTLNum) & spikeTimes<=TTLtimes(TTLNum)+max([duration 0.010]); %keep min window of 10ms
end
onSpikes=logical(sum(pulseIdx,2));

uCodes=unique(ephys.spikes.(s2Sort).codes(:,1));
keepIdx=ismember(ephys.spikes.(s2Sort).codes(:,1),uCodes); %[0,1,2,3] %([1:6])
spikeTimes=single(ephys.spikes.(s2Sort).times(keepIdx));
pulseIdx=false(size(spikeTimes,1),size(TTLtimes,1));
for TTLNum=1:length(TTLtimes)
    pulseIdx(:,TTLNum)=spikeTimes>=TTLtimes(TTLNum) & spikeTimes<=TTLtimes(TTLNum)+max([duration 0.010]); %keep min window of 10ms
end
onSpikesS2=logical(sum(pulseIdx,2));

%find offset
spktOffset= ephys.spikes.times(find(onSpikes,1))-ephys.spikes.(s2Sort).times(find(onSpikesS2,1));

% check spike times correspond to spike times sorted with JRC
s2spikeTimes=ephys.spikes.(s2Sort).times(onSpikesS2)+spktOffset;
knownTS=s2spikeTimes; %oTS(ismember(oTS,s2spikeTimes)); %ephys.timestamps

% get waveforms
keepTrace = mode(ephys.spikes.preferredElectrode);
waveForms=ExtractChunks(ephys.traces(keepTrace,:),... %foo = PreProcData(foo,30000,{'bandpass',[300 3000]});
    knownTS*ephys.recInfo.samplingRate,50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
waveForms=waveForms*ephys.recInfo.bitResolution;

% add spike times if needed (should be in the double digits)
s2onlyIndex=~ismember(knownTS,ephys.spikes.times);
knownTS=knownTS(s2onlyIndex);
waveForms=waveForms(s2onlyIndex,:);

if false
    ephys.spikes.times=[ephys.spikes.times;knownTS];
    [ephys.spikes.times,sortIdx]=sort(ephys.spikes.times);
    ephys.spikes.waveforms=[ephys.spikes.waveforms;waveForms];
    ephys.spikes.waveforms=ephys.spikes.waveforms(sortIdx,:);
    ephys.spikes.preferredElectrode=[ephys.spikes.preferredElectrode;...
        ones(numel(knownTS),1,'uint32')*mode(ephys.spikes.preferredElectrode)];
    ephys.spikes.unitID=[ephys.spikes.unitID;...
        ones(numel(knownTS),1)*mode(ephys.spikes.unitID)];
else
    ephys.rasters=EphysFun.MakeRasters([ephys.spikes.times;knownTS],...
        [ephys.spikes.unitID;ones(numel(knownTS),1)*mode(ephys.spikes.unitID)],...
        1,int32(size(ephys.traces,2)/ephys.spikes.samplingRate*1000));
end

if false
    wfSEM=std(waveForms)/ sqrt(size(waveForms,2)); %standard error of the mean
    wfSEM = wfSEM * 1.96; % 95% of the data will fall within 1.96 standard deviations of a normal distribution
    
    figure; hold on
    plot(mean(ephys.spikes.waveforms(onSpikes,:))*ephys.spikes.bitResolution)
    plot(mean(ephys.spikes.(s2Sort).values(onSpikesS2,10:end))*75)
    plot(mean(waveForms))
    patch([1:length(wfSEM),fliplr(1:length(wfSEM))],...
        [mean(waveForms)-wfSEM,fliplr(mean(waveForms)+wfSEM)],...
        'k','EdgeColor','none','FaceAlpha',0.2); %cmap(cellNum,:)
    
    figure; hold on
    plot(ephys.spikes.times(onSpikes),ones(numel(find(onSpikes),1)),'bd')
    plot(ephys.spikes.(s2Sort).times(onSpikesS2)+spktOffset,ones(numel(find(onSpikesS2),1)),'ro')
end
