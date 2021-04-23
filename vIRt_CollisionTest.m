function vIRt_CollisionTest(ephysData,pulses,uIdx,savePlots)
if nargin < 3
    uIdx=NaN;
end
if nargin < 4
    savePlots =false;
end

%% variables
fileName=ephysData.recInfo.sessionName; %'vIRt44_1210_5450';
TTLtimes=pulses.TTLTimes; %(1,:); TTLs.end=pulses.TTLTimes(2,:);
pulseDur=round(pulses.duration,3); %  min(mode(TTLs.end-TTLtimes));
IPI=mode(diff(TTLtimes));
delay=0.005;
preAlignWindow=0.050;
postAlignWindow=0.05; %0.20;
SRR=ephysData.recInfo.SRratio;
traceExcerpt.excerptSize=SRR;

% spikeData.selectedUnits=[7,18,24]-1;
if islogical(ephysData.selectedUnits) %logical array
    ephysData.selectedUnits=find(ephysData.selectedUnits);
end

%% compute rasters
if isfield(ephysData,'rasters')
    spikeRasters=ephysData.rasters;
else
    spikeRasters=EphysFun.MakeRasters(ephysData.spikes.times,ephysData.spikes.unitID,...
        1,int32(size(ephysData.traces,2)/ephysData.spikes.samplingRate*1000)); %ephysData.spikes.samplingRate
end
if size(spikeRasters,1)>1
    spikeRasters=spikeRasters(ephysData.selectedUnits,:);
end
alignedRasters=EphysFun.AlignRasters(spikeRasters,TTLtimes,preAlignWindow,postAlignWindow,1000);

%% compute spike density functions
% spikeRate=EphysFun.MakeSDF(spikeRasters);

%% Figures
% %if need to load ephys data:
% spikeSortingDir=[ephysData.recInfo.dirName filesep 'SpikeSorting' filesep ephysData.recInfo.sessionName];
% LoadSpikeData(fullfile(spikeSortingDir, [ephysData.recInfo.sessionName '_export_res.mat'])) ;

for tShift=0:0.0001:0.001
    TTLtimes=pulses.TTLTimes-tShift;
    
    traceExcerpt=nan(size(ephysData.traces,1),SRR/100+(0.01*SRR),length(TTLtimes));
    excerptSpikeTimes=cell(numel(TTLtimes),3);
    
    prefElec=double(ephysData.spikes.preferredElectrode(ismember(...
        ephysData.spikes.unitID,ephysData.selectedUnits)));
    try
        keepTrace=mode(prefElec(ephysData.spikes.times>(TTLtimes(1))));
    catch
        keepTrace=mode(prefElec);
    end
    
    for pulseNum=1:length(TTLtimes)
        excerptWindow=int32((TTLtimes(pulseNum)*SRR)-SRR/100:...
            (TTLtimes(pulseNum)*SRR)+(0.01*SRR)-1);
        traceExcerpt(:,:,pulseNum)=ephysData.traces(1:size(ephysData.traces,1),excerptWindow);
        
        [sTimes,excerptSpikeTimes{pulseNum,1}]=deal(double(ephysData.spikes.times(ephysData.spikes.times>=...
            double(excerptWindow(1))/SRR & ephysData.spikes.times<=double(excerptWindow(end))/SRR))...
            -TTLtimes(pulseNum));
        excerptSpikeTimes{pulseNum,2}=sTimes(find(sTimes<=0,1,'last'));
        excerptSpikeTimes{pulseNum,3}=sTimes(find(sTimes>0,1))-sTimes(find(sTimes<=0,1,'last'));
    end
    scaleF = max(max(abs(traceExcerpt(keepTrace,:,:))));
    traceExcerpt=traceExcerpt./scaleF;
    
    figure('Position',[214   108   747   754],'name',...
        [fileName ' Unit' num2str(ephysData.selectedUnits)] );
    hold on
    
    for pulseNum=1:size(traceExcerpt,3)
        if isempty(excerptSpikeTimes{pulseNum,3})
            yShift = 0;
        else
            yShift=floor(excerptSpikeTimes{pulseNum,2}*1000);
        end
%         yShift=0;
            plot(squeeze(traceExcerpt(keepTrace,:,pulseNum))-yShift,'color',[0 0 0 0.5]);
            plot(excerptSpikeTimes{pulseNum,1}*SRR+SRR/100,...
                zeros(1,numel(excerptSpikeTimes{pulseNum,1}))-yShift-0.7,...
                'o','MarkerFaceColor','r','MarkerEdgeColor','None','LineWidth',0.5);
    end
    
    xline(SRR/100);
    xTicks=0:SRR/1000:0.02*SRR;
    set(gca,'xtick',xTicks,'xticklabel',xTicks/30-10)
    axis tight
end

title(['Collision test - Cell #' num2str(uIdx) ' '...
    ephysData.recInfo.baseName ' Unit '...
    num2str(ephysData.selectedUnits)],...
    'interpreter','none')

% excerptTTLtimes=double(TTLtimes(TTLtimes>(traceExcerpt.location-...
%     traceExcerpt.excerptSize)/SRR &...
%     TTLtimes<(traceExcerpt.location+10*traceExcerpt.excerptSize)/SRR)-...
%     (traceExcerpt.location-traceExcerpt.excerptSize)/...
%     SRR)*SRR;
% 
% excerptSpikeTimes={double(ephysData.spikes.times(ephysData.spikes.times>(traceExcerpt.location-...
%     traceExcerpt.excerptSize)/SRR &...
%     ephysData.spikes.times<(traceExcerpt.location+10*traceExcerpt.excerptSize)/SRR)-...
%     (traceExcerpt.location-traceExcerpt.excerptSize)/...
%     SRR)*SRR};
% 
% 
% prefElec=double(ephysData.spikes.preferredElectrode(ismember(...
%     ephysData.spikes.unitID,ephysData.selectedUnits(cellNum))));
% try
%     keepTrace=mode(prefElec(ephysData.spikes.times>(traceExcerpt.location-...
%         traceExcerpt.excerptSize)/SRR));
% catch
%     keepTrace=mode(prefElec);
% end
% traceExcerpt.data=ephysData.traces(keepTrace,excerptWindow);
% 
% figure('Position',[214   108   747   754],'name',...
%     [fileName ' Unit' num2str(ephysData.selectedUnits(cellNum))] ); %Ch' num2str(spikeData.selectedUnits(cellNum))
% hold on
% plot(traceExcerpt.data,'color','k');
% 
% if ~isempty(excerptTTLtimes)
%     for TTLNum=1:length(excerptTTLtimes)
%         patch([excerptTTLtimes(TTLNum), excerptTTLtimes(TTLNum),...
%             excerptTTLtimes(TTLNum)+pulseDur*SRR, excerptTTLtimes(TTLNum)+pulseDur*SRR], ...
%             [get(gca,'ylim') fliplr(get(gca,'ylim'))], ...
%             [0 0 0 0],[0.3 0.75 0.93],'EdgeColor','none','FaceAlpha',0.7);
%     end
%     
% end
% box off;
% xlabel('Time (ms)');
% set(gca,'Color','white','FontSize',12,'FontName','Helvetica');
% 
% %%
% for cellNum=1:length(excerptSpikeTimes)
%     if ~isempty(excerptSpikeTimes{cellNum}) && all(~isnan(excerptSpikeTimes{cellNum}))
%         %plot spike id labels
%         spkLabelYLoc=ones(1,size(excerptSpikeTimes{cellNum},2))*(min(get(gca,'ylim'))/4*3);
%         plot(double(excerptSpikeTimes{cellNum}), ...%-double(traceExcerpt.location-traceExcerpt.excerptSize),...
%             spkLabelYLoc,'Color','k',...%[cmap(cellNum,:),0.4],...
%             'linestyle','none','Marker','x'); % 0 1 2 3 s ^ o x
%     end
% end

if savePlots
    if ~exist(fullfile(cd, 'Figures'),'dir')
        mkdir('Figures')
    end
    savefig(gcf,fullfile(cd, 'Figures', [fileName '_Unit' num2str(cellNum) '_Collision.fig']));
    print(gcf,fullfile(cd, 'Figures', [fileName '_Unit' num2str(cellNum) '_Collision']),'-dpng');
    close(gcf)
end

end





