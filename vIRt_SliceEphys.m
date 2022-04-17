clearvars
% list recordings
fileDir = 'D:\Vincent\vIRt - Leeyup';
files = dir(fullfile(fileDir,('**/*.abf')));

% slice recordings data structure
sr=struct('date',[],'filename',[],'event',[],'trace',[],...
    'spikeTimes',[],'eventTime',[],'ISI',{},'CV2',{},...
    'interval',[],'header',[]);

eventList='K_added'; %'all'
switch eventList
    case 'all'
        recFile='sliceRecs.mat';
        opt.eventType={'any'};
        opt.eventExclude={};
    case 'K_added'
        files_currated=readtable('Cells_Of_Interest');
        files=files(ismember({files.name},[files_currated.filename]));
        recFile='sliceRecs_K.mat';
        opt.eventType={'4AP','K'};
        opt.eventExclude={'STR'};
end

if exist(fullfile(cd,recFile),'file')
    load(recFile);
else
    opt.preTime=5;
    opt.postTime=10; %5
    opt.rmsThld=5;
    opt.filtThld=[5 6000];
    opt.keepTrace=false;
    opt.outfile=recFile;
    SliceRecFileIO(files, opt);
end

%%
doPlots = false;
if ~exist(fullfile(cd, 'Figures'),'dir'); mkdir('Figures'); end
        
recList=find(contains({sr.event},opt.eventType,'IgnoreCase',true) &...
    ~contains({sr.event},opt.eventExclude,'IgnoreCase',true));
%     recList=[33,34,35,38,41,42,39,40,44,49,48,64,62,63]; %1:64;
sr_p=sr(recList);
%     fldNm=fieldnames(sr_p);
%     writetable(struct2table(rmfield(sr_p,fldNm(4:end))),fullfile(cd,'Cells_Of_Interest.xls'));

% check 10/26/21, 11/03/21, 1/18/2022, 02/17/22
% rec #89 20220118_22118003_K9 + 4AP + Mg0 ACSF

numISI=cellfun(@(x) numel(x), [sr_p.ISI]);
% only keep those with enough spikes detected
sr_p=sr_p(all(numISI>100));

%4,5,(8 same cell as 5),10,14
%     data=sr_p([4,5,10,14]);

varTypes={'double','double','logical','logical','logical'};
varNames={'BC_pre','BC_post','BF_pre','BF_post','newBimod'};
bimodTest=table('Size',[size(sr_p,2),5],'VariableTypes',varTypes,'VariableNames',varNames);
varNames={'burstyProp_pre','burstyProp_post','isBurstFreq_pre','isBurstFreq_post','newBursty'};
distribShift=table('Size',[size(sr_p,2),5],'VariableTypes',varTypes,'VariableNames',varNames);

for rNum=1:size(sr_p,2)
    % rNum=6;
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
    CV2_pre=data.CV2{1}; CV2_post=data.CV2{2};
    
    % check for bimodality
    %     [xpdf, n, b] = compute_xpdf(ISI_pre);
    %     [dip, p_value, xlow, xup] = HartigansDipSignifTest(xpdf, nboot);
    
    [BF_pre, BC_pre] = bimodalitycoeff(ISI_pre);
    [BF_post, BC_post] = bimodalitycoeff(ISI_post);
    [ISIpc_pre,edges_pre]=histcounts(ISI_pre,0:0.005:0.5,'Normalization','probability');
    [ISIpc_post,edges_post]=histcounts(ISI_post,0:0.005:0.5,'Normalization','probability');
    burstyProp_pre=sum(ISIpc_pre(edges_pre(2:end)<0.035))*100;
    burstyProp_post=sum(ISIpc_post(edges_post(2:end)<0.035))*100;
    isBurstFreq_pre=burstyProp_pre>5;
    isBurstFreq_post=burstyProp_post>5;
    
    % save values
    bimodTest(rNum,:)={BC_pre,BC_post,BF_pre,BF_post,~BF_pre && BF_post};
    distribShift(rNum,:)={burstyProp_pre,burstyProp_post,isBurstFreq_pre,...
        isBurstFreq_post,~isBurstFreq_pre && isBurstFreq_post};
    if distribShift.newBursty(rNum) & doPlots  %~BF_pre && BF_post
        %% Plot ISI
        figure('color','white','position',[871   460   560   420])
        subplot(3,2,1:2); %axis off
        title({['recording day ' sr_p(rNum).date...
            ' file ' sr_p(rNum).filename(1:end-4)];...
            ['Condition: ' deblank(sr_p(rNum).event)]},...
            'interpreter','none')
        %         plot(data.trace)
        swarmchart([ones(1,numel(ISI_pre)), ones(1,numel(ISI_post))*2],[ISI_pre,ISI_post])
        ylim([0 0.03])
        subplot(3,2,[3,5])
        histogram(ISI_pre,0:0.005:0.5)
        xlabel('ISI (s)')
        title({['pre ' deblank(sr_p(rNum).event)];['Bimodality coeff. = ' num2str(BC_pre, 2)]})
        subplot(3,2,[4,6])
        histogram(ISI_post,0:0.005:0.5)
        xlabel('ISI (s)')
        title({['post ' deblank(sr_p(rNum).event)];['Bimodality coeff. = ' num2str(BC_post, 2)]})
        
        %% save figure
        epIdx=sr_p(rNum).header.tags(vertcat(sr_p(rNum).header.tags.timeSinceRecStart)==...
            sr_p(rNum).eventTime).episodeIndex;
        savefig(gcf,fullfile(cd, 'Figures', [sr_p(rNum).date '_'...
            sr_p(rNum).filename(1:end-4) '_' num2str(epIdx) ...
            '_' deblank(sr_p(rNum).event) '_ISI.fig']));
        print(gcf,fullfile(cd, 'Figures', [sr_p(rNum).date '_'...
            sr_p(rNum).filename(1:end-4) '_' num2str(epIdx) ...
            '_' deblank(sr_p(rNum).event) '_ISI']),'-dpng');
        close(gcf)
    
    end
    fldNm=fieldnames(sr_p);fldIdx=~contains(fldNm,{'filename','date','event'});
    summaryTable=[struct2table(rmfield(sr_p,fldNm(fldIdx))) bimodTest distribShift];
    writetable(summaryTable,fullfile(cd,'BurstyCells.xls'));
end

%% keep unique recordings
% summaryTable=readTable(fullfile(cd,'BurstyCells.xls'));
keepRecIdx=summaryTable.newBursty;
% sr_b=sr_p(keepRecIdx);
sr_b=summaryTable(keepRecIdx,:);
%choose which file to keep for duplicates
keepFile = true(size(sr_b,1),1);
for fileNum=1:size(sr_b,1)
    if sum(contains(sr_b.filename,sr_b.filename(fileNum)))>1
        duplFileIdx=find(contains(sr_b.filename,sr_b.filename(fileNum)));
        tFidx=duplFileIdx(duplFileIdx==fileNum);
        oFidx=duplFileIdx(duplFileIdx~=fileNum);
        keepFile(fileNum)=sr_b.burstyProp_post(tFidx)/sr_b.burstyProp_pre(tFidx)>...
            sr_b.burstyProp_post(oFidx)/sr_b.burstyProp_pre(oFidx);
    end
end

sr_b=sr_b(keepFile,:);

if doPlots
    figure('color','white','position',[1209 201 581 468]); hold on
    plot([ones(size(sr_b,1),1),ones(size(sr_b,1),1)*2]',[sr_b.burstyProp_pre,sr_b.burstyProp_post]',...
        'color',[0.5 0.5 0.5 0.5],'marker','o')
    plot([1 2],mean([sr_b.burstyProp_pre,sr_b.burstyProp_post]),'k','LineWidth',1.5,'marker','s')
    set(gca,'xlim',[0 3],'XTickLabel',{'','','pre','','post','',''})
    xlabel('Event chronoly')
    ylabel('Percentage of ISI <35ms')
end


if false
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