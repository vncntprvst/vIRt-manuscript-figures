clearvars
% list recordings
files = dir('**/*.abf') ;

sliceRecs=struct('date',[],'filename',[],'trace',[],'interval',[],'header',[]);

for fileNum = 1:length(files)
    fName=fullfile(files(fileNum).folder,files(fileNum).name);
    [~,sliceRecs(fileNum).date]=fileparts(files(fileNum).folder);
    sliceRecs(fileNum).filename=files(fileNum).name;
    try
        [data,sliceRecs(fileNum).interval,sliceRecs(fileNum).header] = abfload(fName) ;
        sliceRecs(fileNum).trace=reshape(squeeze(data(:,1,:)),[numel(data(:,1,:)), 1])';
        clearvars data
    catch
        continue
    end
    % end
    
    if false
        figure; hold on
        plot(sliceRecs(fileNum).trace(1:1000))
        %     plot(Primary__.values(1:1000)) %Spike export - for comparison
    end
    
    %% create timeline
    tStep=sliceRecs(fileNum).interval/10^6;% sampling interval in us
    timestamps = 0:tStep:numel(sliceRecs(fileNum).trace)*tStep;
    %% get notes
    notes=sliceRecs(fileNum).header.tags;
    
    for evtNum=1:size(notes,2)
        eventTime=notes(evtNum).timeSinceRecStart;
        preTime= 5*60; %
        postTime= 5*60;
        evWindow=int32((eventTime-preTime)/tStep:(eventTime+postTime)/tStep);
        evWindow=evWindow(evWindow>=1);evWindow=evWindow(evWindow<=numel(sliceRecs(fileNum).trace));
        filtTrace=PreProcData(sliceRecs(fileNum).trace(evWindow),1/tStep,{'bandpass',[5 6000]});
        %         sliceRecs(fileNum).trace=[];
        wTimestamps=timestamps(evWindow);
        threshold = -rms(filtTrace);
        
        %% get spike times
        spikes = logical(filtTrace < threshold);
        spikeProp = regionprops(spikes, 'Area', 'PixelIdxList');
        spikeIdx = cellfun(@(timeIndex) timeIndex(1), {spikeProp.PixelIdxList});
        spikeTimes = wTimestamps(spikeIdx);
        
        %% compute CV2
        ISI_pre=diff(spikeTimes(spikeTimes<eventTime));
        CV2_pre=2*[abs(diff(ISI_pre)) NaN]./movsum(ISI_pre,2);
        ISI_post=diff(spikeTimes(spikeTimes>=eventTime));
        CV2_post=2*[abs(diff(ISI_post)) NaN]./movsum(ISI_post,2);
        
        %% plot trace around event times
        figure('color','white','position',[21 478 1866 444])
        plot(wTimestamps,filtTrace,'k')
        xline(eventTime,'--r',notes(evtNum).comment,'LineWidth',1);
        yline(threshold,'--g','spike threshold')
        xlabel('time (s)')
        ylabel(['current amplitude (' (sliceRecs(fileNum).header.recChUnits{1}) ')'])
        title({['recording day ' sliceRecs(fileNum).date...
            ' file ' sliceRecs(fileNum).filename(1:end-4)];...
            ['Condition: ' deblank(notes(evtNum).comment) ...
            ' CV2_pre: ' num2str(nanmean(CV2_pre))...
            ' CV2_post: ' num2str(nanmean(CV2_post))]},'interpreter','none')
        hold on
        plot(spikeTimes,ones(numel(spikeTimes),1)*(min(filtTrace)-1),'bd')
        
        %% save figure
        if ~exist(fullfile(cd, 'Figures'),'dir')
            mkdir('Figures')
        end
        savefig(gcf,fullfile(cd, 'Figures', [sliceRecs(fileNum).date '_'...
            sliceRecs(fileNum).filename(1:end-4)...
            '_' deblank(notes(evtNum).comment) '.fig']));
        print(gcf,fullfile(cd, 'Figures', [sliceRecs(fileNum).date '_'...
            sliceRecs(fileNum).filename(1:end-4)...
            '_' deblank(notes(evtNum).comment)]),'-dpng');
        close(gcf)
    end
end
save('sliceRecs','sliceRecs','-v7.3');

%% Find bursts with ISI distribution
% check 10/26/21, 11/03/21, 1/18/2022, 02/17/22
% rec #89 20220118_22118003_K9 + 4AP + Mg0 ACSF
fileNum=89;
data=sliceRecs(fileNum);
tStep=data.interval/10^6;% sampling interval in us
timestamps = 0:tStep:numel(data.trace)*tStep;
% get notes
notes=data.header.tags;

%%choose event
evtNum=1;
eventTime=notes(evtNum).timeSinceRecStart;
preTime= 5*60; %
postTime= 5*60;
evWindow=int32((eventTime-preTime)/tStep:(eventTime+postTime)/tStep);
evWindow=evWindow(evWindow>=1);evWindow=evWindow(evWindow<=numel(data.trace));
filtTrace=PreProcData(data.trace(evWindow),1/tStep,{'bandpass',[5 6000]});
%         data.trace=[];
wTimestamps=timestamps(evWindow);
threshold = -rms(filtTrace)*5; %adjust threshold

figure('color','white','position',[21 478 1866 444])
hold on
plot(wTimestamps,filtTrace,'k')
xline(eventTime,'--r',notes(evtNum).comment,'LineWidth',1);
yline(threshold,'--g','spike threshold')
xlabel('time (s)')
ylabel(['current amplitude (' (data.header.recChUnits{1}) ')'])


%% get spike times
spikes = logical(filtTrace < threshold);
spikeProp = regionprops(spikes, 'Area', 'PixelIdxList');
spikeIdx = cellfun(@(timeIndex) timeIndex(1), {spikeProp.PixelIdxList});
spikeTimes = wTimestamps(spikeIdx);

plot(spikeTimes,ones(numel(spikeTimes),1)*(min(filtTrace)-1),'bd')

%% compute ISI and CV2
ISI_pre=diff(spikeTimes(spikeTimes<eventTime)*1000);
CV2_pre=2*[abs(diff(ISI_pre)) NaN]./movsum(ISI_pre,2);
ISI_post=diff(spikeTimes(spikeTimes>=eventTime)*1000);
CV2_post=2*[abs(diff(ISI_post)) NaN]./movsum(ISI_post,2);

figure('color','white','position',[871   460   560   420])
subplot(1,2,1)
histogram(ISI_pre,0:5:200)
xlabel('ISI (ms)')
title(['pre ' deblank(notes(evtNum).comment)])
subplot(1,2,2)
histogram(ISI_post,0:5:200)
xlabel('ISI (ms)')
title(['post ' deblank(notes(evtNum).comment)])

figure('color','white','position',[871   460   560   420])
subplot(1,2,1)
histogram(CV2_pre,0:0.05:2)
xlabel('CV2')
title(['pre ' deblank(notes(evtNum).comment)])
subplot(1,2,2)
histogram(CV2_post,0:0.05:2)
xlabel('CV2')
title(['post ' deblank(notes(evtNum).comment)])

%% Get Rhythmicity
% Find first spike of each burst, then do FFT
postST=spikeTimes(spikeTimes>=eventTime);
burstIdx=bwconncomp([0 diff(postST)<0.04]);
burstFirstSpikeIdx=cellfun(@(spkt) spkt(1), burstIdx.PixelIdxList);

figure('color','white','position',[871   460   560   420])
histogram(diff(postST(burstFirstSpikeIdx)),0:0.05:4)
xlabel('Inter-burst interval (s)')
title(['post ' deblank(notes(evtNum).comment)])

postST=spikeTimes(spikeTimes>=eventTime)*1000;
burstIdx=bwconncomp([0 diff(postST)<40]);
burstFirstSpikeIdx=cellfun(@(spkt) spkt(1), burstIdx.PixelIdxList);
postST=int32(postST-postST(1));
spikeRaster=zeros(1, postST(end));
for burstNum=1:size(burstFirstSpikeIdx,2)
    spikeRaster(postST(burstFirstSpikeIdx(burstNum):burstFirstSpikeIdx(burstNum)+4))=[1 1 1 1 1] ;
end

figure('color','white','position',[871   460   560   420])
spectrogram(spikeRaster(1:60000),128,120,128,1000,'yaxis')
title(['post ' deblank(notes(evtNum).comment)])
