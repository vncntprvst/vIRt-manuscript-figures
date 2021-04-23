function ephys=benchmark_and_fix(ephys,option)

if contains(option,'save')
    prefElec=double(ephys.spikes.preferredElectrode(ismember(...
        ephys.spikes.unitID,ephys.selectedUnits(1))));
    keepTrace=mode(prefElec);
    
    % vIRt44_1210_5300 -> keepTrace=2;
    % vIRt51_1201_5300 -> keepTrace=29;
    % vIRt47_0805_5744 -> keepTrace=23;
    % vIRt44_1210_5450 -> keepTrace=14;
    
    trace=ephys.traces(keepTrace,:);
    fileID = fopen(fullfile('D:\Vincent\Analysis\Traces', [sessInfo.baseName...
        'Ch' num2str(keepTrace) '.bin']),'w');
    fwrite(fileID,trace,'int16');
    fclose(fileID);
end

%%
%
% ephys.spikes.times=vIRt44_1210_5300Ch2_Ch3.times;
% ephys.spikes.unitID=vIRt44_1210_5300Ch2_Ch3.codes(:,1);
% ephys.spikes.unitID=ones(numel(ephys.spikes.times),1)*2;
% ephys.spikes.waveforms=vIRt44_1210_5300Ch2_Ch3.values;
% ephys.spikes.preferredElectrode=ones(numel(ephys.spikes.unitID),1)*2;

% s2Sort=load(fullfile('D:\Vincent\Analysis\Traces', [sessInfo.baseName...
%     'Ch' num2str(keepTrace) '_S2sorted.mat']));
% ephys.spikes.times=s2Sort.([sessInfo.baseName...
%     'Ch' num2str(keepTrace) '_Ch2']).times;
% ephys.spikes.unitID=s2Sort.([sessInfo.baseName...
%     'Ch' num2str(keepTrace) '_Ch2']).codes(:,1);
% ephys.spikes.waveforms=s2Sort.([sessInfo.baseName...
%     'Ch' num2str(keepTrace) '_Ch2']).values;
%
% ephys.spikes.times=ephys.spikes.times(ephys.spikes.unitID==1);
% ephys.spikes.waveforms=ephys.spikes.waveforms(ephys.spikes.unitID==1,:);
% ephys.spikes.unitID=ephys.spikes.unitID(ephys.spikes.unitID==1);
% ephys.spikes.preferredElectrode=ones(numel(ephys.spikes.unitID),1)*2;

if contains(option,'load')
    prefElec=double(ephys.spikes.preferredElectrode(ismember(...
        ephys.spikes.unitID,ephys.selectedUnits(1))));
    keepTrace=mode(prefElec);
    
    % load(fullfile('D:\Vincent\Analysis\Traces', [sessInfo.baseName...
    %     'Ch' num2str(keepTrace) '_S2sorted.mat']));
    % varnames=who;
    % s2Sort=varnames{cellfun(@(x) contains(x,ephys.recInfo.baseName),varnames)};
    ephys.spikes.([ephys.recInfo.baseName '_Ch' num2str(keepTrace)]) =...
        load(fullfile('D:\Vincent\Analysis\Traces', [ephys.recInfo.baseName...
        'Ch' num2str(keepTrace) '_S2sorted.mat']));
    s2SortName=fieldnames(ephys.spikes.([ephys.recInfo.baseName '_Ch' num2str(keepTrace)]));
    ephys.spikes.([ephys.recInfo.baseName '_Ch' num2str(keepTrace)]) =...
        ephys.spikes.([ephys.recInfo.baseName '_Ch' num2str(keepTrace)]).(s2SortName{1});
end

%%
if contains(option,'plots')
    %     use S2sort instead
    fldNames=fieldnames(ephys.spikes);
    s2Sort=fldNames{cellfun(@(x) contains(x,[ephys.recInfo.baseName '_Ch']),fldNames)};
    
    if ~isempty(s2Sort)
        % find spikes that occur during pulses
        
        spikeTimes=single(ephys.spikes.times);
        TTLtimes=pulses.TTLTimes;
        duration=pulses.duration;
        pulseIdx=false(size(spikeTimes,1),size(TTLtimes,1));
        for TTLNum=1:length(TTLtimes)
            pulseIdx(:,TTLNum)=spikeTimes>=TTLtimes(TTLNum) & spikeTimes<=TTLtimes(TTLNum)+max([duration 0.010]); %keep min window of 10ms
        end
        onSpikes=logical(sum(pulseIdx,2));
        
        spikeTimes=single(ephys.spikes.(s2Sort).times);
        pulseIdx=false(size(spikeTimes,1),size(TTLtimes,1));
        for TTLNum=1:length(TTLtimes)
            pulseIdx(:,TTLNum)=spikeTimes>=TTLtimes(TTLNum) & spikeTimes<=TTLtimes(TTLNum)+max([duration 0.010]); %keep min window of 10ms
        end
        onSpikesS2=logical(sum(pulseIdx,2));
        
        spktOffset= ephys.spikes.times(find(onSpikes,1))-ephys.spikes.(s2Sort).times(find(onSpikesS2,1));
        s2spikeTimes=ephys.spikes.(s2Sort).times(onSpikesS2)+spktOffset;
        ismember(s2spikeTimes,ephys.spikes.times(onSpikes));
        knownTS=ismember(ephys.timestamps,s2spikeTimes);
        waveForms=ExtractChunks(ephys.traces(keepTrace,:),... %foo = PreProcData(foo,30000,{'bandpass',[300 3000]});
            ephys.timestamps(knownTS)*ephys.recInfo.samplingRate,50,'tshifted'); %'tzero' 'tmiddle' 'tshifted'
        waveForms=waveForms*ephys.recInfo.bitResolution;
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
    
    
    
    
end
