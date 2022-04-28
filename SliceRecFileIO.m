function SliceRecFileIO(files, opt)

% slice recordings data structure
sr=struct('date',[],'filename',[],'event',[],'trace',[],...
    'spikeTimes',[],'eventTime',[],'ISI',{},'CV2',{},...
    'interval',[],'header',[]);

rNum=0;
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
    if ~isempty(notes)
        eventIdx=contains({notes.comment},opt.eventType,'IgnoreCase',true) &...
            ~contains({notes.comment},opt.eventExclude,'IgnoreCase',true);
        notes=notes(eventIdx);
        for evtNum=1:size(notes,2)
            
            eventTime=notes(evtNum).timeSinceRecStart;
            preTime= opt.preTime*60; %
            postTime= opt.postTime*60;
            evWindow=int32((eventTime-preTime)/tStep:(eventTime+postTime)/tStep);
            evWindow=evWindow(evWindow>=1);evWindow=evWindow(evWindow<=numel(trace));
            filtTrace=PreProcData(trace(evWindow),1/tStep,{'bandpass',opt.filtThld});
            wTimestamps=timestamps(evWindow);
            threshold = -rms(filtTrace)*opt.rmsThld;
            
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
            if opt.keepTrace; sr(rNum).trace=filtTrace; end
            sr(rNum).spikeTimes=spikeTimes;
            sr(rNum).ISI={ISI_pre;ISI_post};
            sr(rNum).CV2={CV2_pre;CV2_post};
        end
    end
    clearvars trace interval header notes
end
save(opt.outfile,'sr','-v7.3');
end