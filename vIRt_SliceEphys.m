clearvars
% list recordings
fileDir = 'D:\Vincent\vIRt - Leeyup';
files = dir(fullfile(fileDir,('**/*.abf')));

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
    opt.keepTrace=true;
    opt.outfile=recFile;
    SliceRecFileIO(files(11,:), opt);
end

%%
doPlots = false;
if ~exist(fullfile(cd, 'Figures'),'dir'); mkdir('Figures'); end

if ~exist('sr_p','var')
    recList=find(contains({sr.event},opt.eventType,'IgnoreCase',true) &...
        ~contains({sr.event},opt.eventExclude,'IgnoreCase',true));
    %     recList=[33,34,35,38,41,42,39,40,44,49,48,64,62,63]; %1:64;
    sr_p=sr(recList);
    
    % or load it
    %     load('sliceRecs_K.mat')
    
end

% check 10/26/21, 11/03/21, 1/18/2022, 02/17/22
% rec #89 20220118_22118003_K9 + 4AP + Mg0 ACSF

numISI=cellfun(@(x) numel(x), [sr_p.ISI]);
% only keep those with enough spikes detected
sr_p=sr_p(all(numISI>100));

% save file list in spreadsheet
%     fldNm=fieldnames(sr_p); fldIdx=~contains(fldNm,{'fileIndex','filename','date','event'});
%     writetable(struct2table(rmfield(sr_p,fldNm(fldIdx))),fullfile(cd,'Cells_Of_Interest.xls'));

indivPlot=false;
if indivPlot
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
        burstCutoff=0.025; %0.035
        [BF_pre, BC_pre] = bimodalitycoeff(ISI_pre);
        [BF_post, BC_post] = bimodalitycoeff(ISI_post);
        [ISIpc_pre,edges_pre]=histcounts(ISI_pre,0:0.005:0.5,'Normalization','probability');
        [ISIpc_post,edges_post]=histcounts(ISI_post,0:0.005:0.5,'Normalization','probability');
        burstyProp_pre=sum(ISIpc_pre(edges_pre(2:end)<burstCutoff))*100;
        burstyProp_post=sum(ISIpc_post(edges_post(2:end)<burstCutoff))*100;
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
            histogram(ISI_pre,0:0.005:0.5); % get errorbars with [N,B,E] = isi(cumsum(ISI_pre),[],1);
            xlabel('ISI (s)')
            title({['pre ' deblank(sr_p(rNum).event)];['Bimodality coeff. = ' num2str(BC_pre, 2)]})
            subplot(3,2,[4,6])
            histogram(ISI_post,0:0.005:0.5)
            xlabel('ISI (s)')
            title({['post ' deblank(sr_p(rNum).event)];['Bimodality coeff. = ' num2str(BC_post, 2)]})
            
            %% save figure
            if doPlots
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
            
        end
        fldNm=fieldnames(sr_p);fldIdx=~contains(fldNm,{'fileIndex','filename','date','event'});
        summaryTable=[struct2table(rmfield(sr_p,fldNm(fldIdx))) bimodTest distribShift];
        writetable(summaryTable,fullfile(cd,'BurstyCells.xls'));
    end
end

%% keep unique recordings
% summaryTable=readtable(fullfile(cd,'BurstyCells.xls'));
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
%load('sliceRecs_K.mat')
sr_p=sr_p(ismember([sr_p.fileIndex],[sr_b.fileIndex]));

%% population plots
if doPlots
    %% change in proportion of ISI<25ms ("Bursty prop")
    figure('color','white','position',[1209 201 581 468]); hold on
    plot([ones(size(sr_b,1),1),ones(size(sr_b,1),1)*2]',[sr_b.burstyProp_pre,sr_b.burstyProp_post]',...
        'marker','o'); %'color',[0.5 0.5 0.5 0.5]
    plot([1 2],mean([sr_b.burstyProp_pre,sr_b.burstyProp_post]),'k','LineWidth',1.5,'marker','s')
    set(gca,'xlim',[0 3],'XTickLabel',{'','','pre','','post','',''})
    xlabel('Event chronoly')
    ylabel('Percentage of ISI <25 ms')
    
    %%
    cmap=lines;
    allISIs=reshape([sr_p.ISI],1,size(sr_p,2)*2);allISIs=cellfun(@(x) x(x>0.001), allISIs,'UniformOutput',false);
    scIdx=cellfun(@(x,y) ones(1,numel(x))*y, allISIs,num2cell(1:numel(allISIs)),'UniformOutput',false);
    colIdx=cellfun(@(x,y) ones(numel(x),3).*cmap(y,:), allISIs,...
        num2cell(reshape([1:numel(allISIs)/2;1:numel(allISIs)/2],1,numel(allISIs))),'UniformOutput',false);
    
    figure('color','white','position',[871   460   560   420])
    hold on
    swarmchart([scIdx{1:2:end}],[allISIs{1:2:end}],36,vertcat(colIdx{1:2:end}),'o','XJitterWidth',0.8);
    swarmchart([scIdx{2:2:end}],[allISIs{2:2:end}],36,vertcat(colIdx{2:2:end}),'d','XJitterWidth',0.8);
    
    yline(0.025,'--k','burst cutoff = 40Hz')
    set(gca, 'YScale', 'log'); set(gca,'YTickLabel',get(gca,'Ytick'));
    box off; axis tight
    xlabel('Pre/Post event')
    ylabel('ISI (seconds, log scale)')
    
%     all_ISIs=cell2table(cellfun(@(x) x', allISIs, 'UniformOutput', false));
    
    %% Plot Rhythmicity
    
    %% Plot ISI spectrum
    % set parameters
    params.Fs=1000; %20000; % sampling frequency
    params.fpass=[0 100]; % band of frequencies to be kept
    params.NW = 3; %min([100 floor(max(spikeTimes))*2.5]); %3 %time-bandwidth product
    params.tapers=[params.NW params.NW*2-1]; % taper parameters
    params.pad=2; %0; % pad factor for fft
    params.err=[2 0.01];
    
    %Use multi-taper segmented spectrum analysis
    %     [spS.spectrumVals,spS.freqVals]=mtspectrumpt(spikeTimes,params);
    segmentWin=10; % ten second window
    %     spikeTimes=round(cumsum([0 allISIs{1}]),3);
    %     [spS.spectrumVals,spS.freqVals,spS.R,~,~,~,spS.Serr]=mtspectrumsegpt(spikeTimes,segementWin,params);
    
    if exist('Spectrum_ISI.mat','file')
        load('Spectrum_ISI.mat');
    else
        [spS.spectrumVals,spS.freqVals,spS.R,~,~,~,spS.Serr]=...
            cellfun(@(ISIs) mtspectrumsegpt(round(cumsum([0 ISIs]),3),segmentWin,params), allISIs,'UniformOutput',false);
    end
    
    % convert to dB (power spectral density)
    spS.spectrumValsPSD=cellfun(@(spectrumVals) 10*log10(spectrumVals), spS.spectrumVals, 'UniformOutput', false);
    %     spS.spectrumValsPSD= 10*log10(spS.spectrumVals);
    %     spS.SerrPSD(1,:)=10*log10(spS.Serr(1,:));
    %     spS.SerrPSD(2,:)=10*log10(spS.Serr(2,:));
    %     spS.RPSD=10*log10(spS.R);
    %
    %     spS.StatSigIdx=spS.SerrPSD(1,:)<=spS.RPSD & spS.SerrPSD(2,:)>=spS.RPSD;
    
    spS.diffPSD = cellfun(@(PSD_pre,PSD_post) PSD_post-PSD_pre ,...
        spS.spectrumValsPSD(1:2:end),spS.spectrumValsPSD(2:2:end), 'UniformOutput', false);
    
    %% figures
    if false
        figure; %     plot(spS.freqVals,spS.spectrumVals)
        subplot(2,1,1)
        plot(spS.freqVals{1}, horzcat(spS.spectrumVals{1:2:end}))
        ylim([0 50])
        legend
        subplot(2,1,2)
        plot(spS.freqVals{1}, horzcat(spS.spectrumVals{2:2:end}))
        ylim([0 50])
        legend
        
        figure; %     plot(spS.freqVals,spS.spectrumVals)
        subplot(2,1,1)
        plot(spS.freqVals{1}, horzcat(spS.spectrumValsPSD{1:2:end}))
        %     ylim([0 50])
        legend
        subplot(2,1,2)
        plot(spS.freqVals{1}, horzcat(spS.spectrumValsPSD{2:2:end}))
        %     ylim([0 50])
        legend
        
        figure('color','white'); %     plot(spS.freqVals,spS.spectrumVals)
        plot(spS.freqVals{1}, horzcat(spS.diffPSD{:}))
        title('Diff PSDs (all spikes)')
        box off; axis tight
        xlabel('Frequencies (Hz)')
        ylabel('dB')
        
        %     figure
        %     plot(spS.freqVals,10*log10(spS.spectrumVals),spS.freqVals,10*log10(spS.Serr(1,:)),...
        %         spS.freqVals,10*log10(spS.Serr(2,:)));
        %     line(get(gca,'xlim'),[10*log10(spS.R) 10*log10(spS.R)]);
        %     hold on
        %     plot(spS.freqVals(spS.StatSigIdx),spS.spectrumValsPSD(spS.StatSigIdx),'k')
        
    end
    %% plot IBI histogram
    cmap=lines;
    
    burstCutoff=0.025; %from ISI distribution histograms above
    allISIs=reshape([sr_p.ISI],1,size(sr_p,2)*2);
    allIBIs=cellfun(@(cumsumISI,burstIdx) diff(cumsumISI(burstIdx)), ...
        cellfun(@(x) cumsum(x), allISIs,'UniformOutput',false),...
        cellfun(@(x) find(x>0.001 & x<burstCutoff)-1, allISIs,'UniformOutput',false),...
        'UniformOutput',false);
    allIBIs=cellfun(@(x) 1./(x(x>burstCutoff)), allIBIs,'UniformOutput',false);
    
    
    allIBIs(cellfun(@isempty, allIBIs))={-1};
    
    scIdx=cellfun(@(x,y) ones(1,numel(x))*y, allIBIs,num2cell(1:numel(allIBIs)),'UniformOutput',false);
    colIdx=cellfun(@(x,y) ones(numel(x),3).*cmap(y,:), allIBIs,...
        num2cell(reshape([1:numel(allIBIs)/2;1:numel(allIBIs)/2],1,numel(allIBIs))),'UniformOutput',false);
    
    figure('color','white','position',[871   460   560   420])
    hold on
    swarmchart([scIdx{1:2:end}],[allIBIs{1:2:end}],36,vertcat(colIdx{1:2:end}),'o','XJitterWidth',0.8);
    swarmchart([scIdx{2:2:end}],[allIBIs{2:2:end}],36,vertcat(colIdx{2:2:end}),'d','XJitterWidth',0.8);
    box off; axis tight
    set(gca,'YLim',[0 1], 'XLim', [0 11]);
    set(gca,'YTickLabel',get(gca,'Ytick'));
    xlabel('Pre/Post event')
    ylabel('IBI frequency (Hz)')
    
%     all_IBIs=cell2table(cellfun(@(x) x', allIBIs, 'UniformOutput', false));
    
    %% Plot IBI spectrum
    
    
    % Same parameters as for spiking spectrum analysis above
    params.Fs=1000; %20000; % sampling frequency
    params.fpass=[0 100]; % band of frequencies to be kept
    params.NW = 3; %min([100 floor(max(spikeTimes))*2.5]); %3 %time-bandwidth product
    params.tapers=[params.NW params.NW*2-1]; % taper parameters
    params.pad=2; %0; % pad factor for fft
    params.err=[2 0.01];
    
    %Use multi-taper segmented spectrum analysis
    segmentWin=10; % ten second window
    
    if exist('Spectrum_IBI.mat','file')
        load('Spectrum_IBI.mat');
    else
        [burstS.spectrumVals,burstS.freqVals,burstS.R,~,~,~,burstS.Serr]=...
            cellfun(@(ISIs) mtspectrumsegpt(round(cumsum([0 ISIs]),3),segmentWin,params), allIBIs,'UniformOutput',false);
    end
    
    % convert to dB (power spectral density)
    burstS.spectrumValsPSD=cellfun(@(spectrumVals) 10*log10(spectrumVals), burstS.spectrumVals, 'UniformOutput', false);
    %     burstS.spectrumValsPSD= 10*log10(burstS.spectrumVals);
    %     burstS.SerrPSD(1,:)=10*log10(burstS.Serr(1,:));
    %     burstS.SerrPSD(2,:)=10*log10(burstS.Serr(2,:));
    %     burstS.RPSD=10*log10(burstS.R);
    %
    %     burstS.StatSigIdx=burstS.SerrPSD(1,:)<=burstS.RPSD & burstS.SerrPSD(2,:)>=burstS.RPSD;
    
    burstS.diffPSD = cellfun(@(PSD_pre,PSD_post) PSD_post-PSD_pre ,...
        burstS.spectrumValsPSD(1:2:end),burstS.spectrumValsPSD(2:2:end), 'UniformOutput', false);
    
    %% figures
    if false
        figure; %     plot(burstS.freqVals,burstS.spectrumVals)
        subplot(2,1,1)
        plot(burstS.freqVals{1}, horzcat(burstS.spectrumVals{1:2:end}))
        ylim([0 50])
        legend
        subplot(2,1,2)
        plot(burstS.freqVals{4}, horzcat(burstS.spectrumVals{4:2:end}))
        ylim([0 50])
        legend
        
        figure; %     plot(burstS.freqVals,burstS.spectrumVals)
        subplot(2,1,1)
        plot(burstS.freqVals{1}, horzcat(burstS.spectrumValsPSD{1:2:end}))
        %     ylim([0 50])
        legend
        subplot(2,1,2)
        plot(burstS.freqVals{1}, horzcat(burstS.spectrumValsPSD{2:2:end}))
        %     ylim([0 50])
        legend
        
        figure('color','white'); %     plot(burstS.freqVals,burstS.spectrumVals)
        plot(burstS.freqVals{1}, horzcat(burstS.diffPSD{:}))
        title('Diff PSDs (all spikes)')
        box off; axis tight
        xlabel('Frequencies (Hz)')
        ylabel('dB')
        
    end
    
    
end

%% Compare with Matlab SDF functions
% fs = 1000;
% spikeRasters=EphysFun.MakeRasters(sr_p(2).spikeTimes*1000,ones(1,numel(sr_p(2).spikeTimes)),fs);
% figure('position',[1187 347 579 438]);
% pmtm(spikeRasters',3,length(spikeRasters),fs);
%
% ntp = 13;
% nw = 7.5;
% w = linspace(0,1,10240);
% pmtm(spikeRasters',{nw,ntp},w*(fs/2),fs)
%
% figure('position',[1187 347 579 438])
% % pwelch(spikeRasters',fs)
% [pxx,f] = pwelch(spikeRasters',500,300,500,fs);
% plot(f,10*log10(pxx))
%
% xlabel('Frequency (Hz)')
% ylabel('PSD (dB/Hz)')

if indivPlot
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
    % spikeTimes defined in SliceRecFileIO
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