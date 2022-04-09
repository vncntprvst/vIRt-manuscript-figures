clearvars
% list recordings
files = dir('**/*.abf') ;

% slice recordings data structure
sr=struct('date',[],'filename',[],'event',[],'trace',[],...
    'spikeTimes',[],'eventTime',[],'ISI',{},'CV2',{},...
    'interval',[],'header',[]);
rNum=0;

if exist(fullfile(cd,'sliceRecs.mat'),'file')
    load('sliceRecs.mat');
else
    for fileNum = 1:length(files)
        fName=fullfile(files(fileNum).folder,files(fileNum).name);
        [~,recDate]=fileparts(files(fileNum).folder);
        filename=files(fileNum).name;
        try
            [data,interval,header] = abfload(fName) ;
            trace=reshape(squeeze(data(:,1,:)),[numel(data(:,1,:)), 1])';
            clearvars data
        catch
            continue
        end
        % end
        
        if false
            figure; hold on
            plot(trace(1:1000))
            %     plot(Primary__.values(1:1000)) %Spike export - for comparison
        end
        
        %% create timeline
        tStep=interval/10^6;% sampling interval in us
        timestamps = 0:tStep:numel(trace)*tStep;
        %% get notes
        notes=header.tags;
        
        for evtNum=1:size(notes,2)
            
            eventTime=notes(evtNum).timeSinceRecStart;
            preTime= 5*60; %
            postTime= 5*60;
            evWindow=int32((eventTime-preTime)/tStep:(eventTime+postTime)/tStep);
            evWindow=evWindow(evWindow>=1);evWindow=evWindow(evWindow<=numel(trace));
            filtTrace=PreProcData(trace(evWindow),1/tStep,{'bandpass',[5 6000]});
            wTimestamps=timestamps(evWindow);
            threshold = -rms(filtTrace)*5;
            
            %% get spike times
            spikes = logical(filtTrace < threshold);
            spikeProp = regionprops(spikes, 'Area', 'PixelIdxList');
            spikeIdx = cellfun(@(timeIndex) timeIndex(1), {spikeProp.PixelIdxList});
            spikeTimes = wTimestamps(spikeIdx);
            
            %% compute ISI CV2
            ISI_pre=diff(spikeTimes(spikeTimes<eventTime));
            CV2_pre=2*[abs(diff(ISI_pre)) NaN]./movsum(ISI_pre,2);
            ISI_post=diff(spikeTimes(spikeTimes>=eventTime));
            CV2_post=2*[abs(diff(ISI_post)) NaN]./movsum(ISI_post,2);
            
            %% store data
            rNum=rNum+1;
            sr(rNum).filename=filename;
            sr(rNum).date=recDate;
            sr(rNum).interval=interval;
            sr(rNum).header=header;
            sr(rNum).event=deblank(notes(evtNum).comment);
            sr(rNum).eventTime=eventTime;
            sr(rNum).trace=filtTrace;
            sr(rNum).spikeTimes=spikeTimes;
            sr(rNum).ISI={ISI_pre;ISI_post};
            sr(rNum).CV2={CV2_pre;CV2_post};
        end
        clearvars trace interval header notes
    end
    save('sliceRecs','sr','-v7.3');
end

%%
doPlots = false;
if doPlots
    recList=find(contains({sr.event},'4AP'));
    recList=[33,34,35,38,41,42,39,40,44,49,48,64,62,63]; %1:64; 
    sr_p=sr(recList);
    
    % check 10/26/21, 11/03/21, 1/18/2022, 02/17/22
    % rec #89 20220118_22118003_K9 + 4AP + Mg0 ACSF
    
    for rNum=1:size(sr_p,2)
    % rNum=89;
    data=sr_p(rNum);
%     tStep=data.interval/10^6;% sampling interval in us
%     timestamps = 0:tStep:numel(data.trace)*tStep;
%     % get notes
%     % notes=data.header.tags;
%    
%     %%choose event
%     % evtNum=1;
%     eventTime=notes(evtNum).timeSinceRecStart;
%     preTime= 5*60; %
%     postTime= 5*60;
%     evWindow=int32((eventTime-preTime)/tStep:(eventTime+postTime)/tStep);
%     evWindow=evWindow(evWindow>=1);evWindow=evWindow(evWindow<=numel(data.trace));
%     filtTrace=PreProcData(data.trace(evWindow),1/tStep,{'bandpass',[5 6000]});
%     %         data.trace=[];
%     wTimestamps=timestamps(evWindow);
%     threshold = -rms(filtTrace)*5; %adjust threshold
    
    ISI_pre=data.ISI{1}; ISI_post=data.ISI{2};

    %% Plot ISI and CV2
    
    figure('color','white','position',[871   460   560   420])
    subplot(10,2,1:2)
    title({['recording day ' sr_p(rNum).date...
        ' file ' sr_p(rNum).filename(1:end-4)];...
        ['Condition: ' deblank(sr_p(rNum).event)]},...
        'interpreter','none')
    axis off
    subplot(10,2,3:2:19)
    histogram(ISI_pre,0:0.005:0.2)
    xlabel('ISI (s)')
    title(['pre ' deblank(sr_p(rNum).event)])
    subplot(10,2,4:2:20)
    histogram(ISI_post,0:0.005:0.2)
    xlabel('ISI (s)')
    title(['post ' deblank(sr_p(rNum).event)])
    
%     figure('color','white','position',[871   460   560   420])
%     subplot(1,2,1)
%     histogram(CV2_pre,0:0.05:2)
%     xlabel('CV2')
%     title(['pre ' deblank(notes(evtNum).comment)])
%     subplot(1,2,2)
%     histogram(CV2_post,0:0.05:2)
%     xlabel('CV2')
%     title(['post ' deblank(notes(evtNum).comment)])
    
    end
    rNum=7;
    %% plot trace around event times
    tStep=sr_p(rNum).interval/10^6;% sampling interval in us
    timestamps = 0:tStep:numel(sr_p(rNum).trace)*tStep;
    eventTime=sr_p(rNum).eventTime;
    preTime= 5*60; %
    %     postTime= 5*60;
    %     evWindow=int32((eventTime-preTime)/tStep:(eventTime+postTime)/tStep);
    %     evWindow=evWindow(evWindow>=1);evWindow=evWindow(evWindow<=numel(sr_p(rNum).trace));
    wTimestamps=timestamps(1:end-1);
    ISI_pre=sr_p(rNum).ISI{1};
    ISI_post=sr_p(rNum).ISI{2};
    
    figure('color','white','position',[21 478 1866 444])
    plot(wTimestamps,sr_p(rNum).trace,'k')
    xline(preTime,'--r',sr_p(rNum).event,'LineWidth',1);
    %     yline(threshold,'--g','spike threshold')
    xlabel('time (s)')
    ylabel(['current amplitude (' (sr_p(rNum).header.recChUnits{1}) ')'])
    title({['recording day ' sr_p(rNum).date...
        ' file ' sr_p(rNum).filename(1:end-4)];...
        ['Condition: ' deblank(sr_p(rNum).event) ...
        ' ISI_pre: ' num2str(nanmean(ISI_pre))...
        ' ISI_post: ' num2str(nanmean(ISI_post))]},'interpreter','none')
%     hold on
%     plot(sr_p(rNum).spikeTimes,ones(numel(spikeTimes),1)*(min(sr_p(rNum).trace)-1),'bd')
    
    %% save figure
    if ~exist(fullfile(cd, 'Figures'),'dir')
        mkdir('Figures')
    end
    savefig(gcf,fullfile(cd, 'Figures', [sr_p(rNum).date '_'...
        sr_p(rNum).filename(1:end-4)...
        '_' deblank(notes(evtNum).comment) '.fig']));
    print(gcf,fullfile(cd, 'Figures', [sr_p(rNum).date '_'...
        sr_p(rNum).filename(1:end-4)...
        '_' deblank(notes(evtNum).comment)]),'-dpng');
    close(gcf)

    %% Plot Rhythmicity
    % Find first spike of each burst, plot distribution
    postST=spikeTimes(spikeTimes>=eventTime);
    burstIdx=bwconncomp([0 diff(postST)<0.04]); % Find bursts from ISI distribution
    burstFirstSpikeIdx=cellfun(@(spkt) spkt(1), burstIdx.PixelIdxList);
     
    figure('color','white','position',[871   460   560   420])
    histogram(diff(postST(burstFirstSpikeIdx)),0:0.05:4)
    xlabel('Inter-burst interval (s)')
    title(['post ' deblank(notes(evtNum).comment)])
    
    % Find first spike of each burst, then do FFT
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
end