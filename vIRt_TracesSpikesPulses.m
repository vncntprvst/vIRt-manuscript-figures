function vIRt_TracesSpikesPulses(ephysData,pulses,uIdx,plotType,savePlots)
if nargin < 3
    uIdx=NaN;
end
if nargin < 4
    plotType='excerpt';
end
if nargin < 5
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
spikeRasters=spikeRasters(ephysData.selectedUnits,:);
alignedRasters=EphysFun.AlignRasters(spikeRasters,TTLtimes,preAlignWindow,postAlignWindow,1000);

%% compute spike density functions
% spikeRate=EphysFun.MakeSDF(spikeRasters);

%% Figures
% some issue with ttl times from npy -> see CH29 from 'vIRt22_2018-10-16_18-43-54_5100_50ms1Hz5mW_nopp' KS

%if need to load ephys data:
spikeSortingDir=[ephysData.recInfo.dirName filesep 'SpikeSorting' filesep ephysData.recInfo.sessionName];
LoadSpikeData(fullfile(spikeSortingDir, [ephysData.recInfo.sessionName '_export_res.mat'])) ;

for cellNum=1:size(ephysData.selectedUnits,1)
    
    switch plotType
        case 'excerpt'
            
            if ~isfield(ephysData.recInfo,'SRratio')
                SRR=double(ephysData.spikes.samplingRate/1000);
            end
            % excerptTTLtimes=double(TTLtimes(TTLtimes>(traceExcerpt.location-traceExcerpt.excerptSize)/spikeData.recInfo.SRratio &...
            %     TTLtimes<(traceExcerpt.location+traceExcerpt.excerptSize)/spikeData.recInfo.SRratio)-...
            %     (traceExcerpt.location-traceExcerpt.excerptSize)/spikeData.recInfo.SRratio)*spikeData.recInfo.SRratio;
            % if ~isempty(excerptTTLtimes)
            % %     excerptTTLtimes=excerptTTLtimes(end); %if wants to keep only one pulse
            % else % check further out in the trace
            traceExcerpt.location=TTLtimes(1)*SRR;
            %     mod(winIdxStart,traceData.traceInfo.numChan)
            if exist('traceData','var') && isa(traceData,'memmapfile')
                winIdxStart=(traceExcerpt.location-traceExcerpt.excerptSize)*traceData.traceInfo.numChan+1;
                winSize=2; %default 1 pulse
                winIdxEnd=winIdxStart+(winSize*2*traceExcerpt.excerptSize*traceData.traceInfo.numChan);
            else
                winIdxStart=(traceExcerpt.location-traceExcerpt.excerptSize); %*traceData.traceInfo.numChan+1;
                winIdxEnd=traceExcerpt.location+10*traceExcerpt.excerptSize;
            end
            excerptWindow=int32(winIdxStart:winIdxEnd-1);%-SRR;
            %     size(excerptWindow,2)>(2*traceExcerpt.excerptSize*traceData.traceInfo.numChan)
            if exist('traceData','var') && isa(traceData,'memmapfile')
                traceExcerpt.data=traceData.allTraces.Data(excerptWindow);
                traceExcerpt.data=reshape(traceExcerpt.data,[traceData.traceInfo.numChan traceExcerpt.excerptSize*2*winSize]);
                preprocOption={'CAR','all'};
                traceExcerpt.data=PreProcData(traceExcerpt.data,30000,preprocOption);
                traceExcerpt.data=traceExcerpt.data(channelNum,:);%     figure; plot(dataExcerpt(11,:))
            else
                %Sometimes not the best trace. Find a way plot most relevant trace
                %         prefElec=double(ephysData.spikes.preferredElectrode(ismember(...
                %             ephysData.spikes.unitID,ephysData.selectedUnits(cellNum))));
                %         try
                %             [traceFreq,uniqueTraces]=hist(prefElec,unique(prefElec));
                %             keepTrace=uniqueTraces(end);
                %             traceExcerpt.data=ephysData.traces(keepTrace,excerptWindow);
                %         catch
                %         try
                %             keepTrace=mode(prefElec(ephysData.spikes.times>(traceExcerpt.location-...
                %                 traceExcerpt.excerptSize)/SRR));
                %         catch
                %             keepTrace=mode(prefElec);
                %         end
                %         traceExcerpt.data=ephysData.traces(:,excerptWindow);
                %         end
                %         figure; plot(traceExcerpt.data)
                %         figure; plot(ephysData.traces(keepTrace,:))
            end
            
            excerptTTLtimes=double(TTLtimes(TTLtimes>(traceExcerpt.location-...
                traceExcerpt.excerptSize)/SRR &...
                TTLtimes<(traceExcerpt.location+10*traceExcerpt.excerptSize)/SRR)-...
                (traceExcerpt.location-traceExcerpt.excerptSize)/...
                SRR)*SRR;
            
            try
                excerptSpikeTimes={double(ephysData.spikes.times(ephysData.spikes.times>(traceExcerpt.location-...
                    traceExcerpt.excerptSize)/SRR &...
                    ephysData.spikes.times<(traceExcerpt.location+10*traceExcerpt.excerptSize)/SRR)-...
                    (traceExcerpt.location-traceExcerpt.excerptSize)/...
                    SRR)*SRR};
            catch
                excerptSpikeTimes={NaN};
            end
            
            figure('Position',[214   108   747   754],'name',...
                [fileName ' Unit' num2str(ephysData.selectedUnits(cellNum))] ); %Ch' num2str(spikeData.selectedUnits(cellNum))
            traceIDs=1:size(ephysData.traces,1); %max([keepTrace-5 1]):min([keepTrace+5 size(ephysData.traces,1)]);
            %     for traceNum=1:numel(traceIDs)
            %         subplot(numel(traceIDs),1,traceNum)
            hold on
            %         keepTrace=traceIDs(traceNum);
            traceExcerpt.data=ephysData.traces(traceIDs,excerptWindow);
            
            OptoRawTrace(traceExcerpt,excerptSpikeTimes,...
                SRR,excerptTTLtimes,pulseDur,'',gca)
            
            title(['Cell list # ' uIdx ' Unit ' num2str(ephysData.selectedUnits(cellNum))])
            %     end
            
            if savePlots
                savefig(gcf,[cd filesep ' Cell_' ... %'Z:\Vincent\Figures\vIRt\PT'
                    num2str(PTCells(cellNum)) '_Traces_Spikes_Pulses'],'compact');
                close(gcf)
            end
            
        case 'average'
            %             TTLtimes=TTLtimes-0.001;
            traceExcerpt=nan(size(ephysData.traces,1),SRR/200+1.2*(pulseDur*SRR),length(TTLtimes));
            
            for TTLNum=1:length(TTLtimes)
                excerptWindow=int32((TTLtimes(TTLNum)*SRR)-SRR/200:...
                    (TTLtimes(TTLNum)*SRR)+1.2*(pulseDur*SRR)-1);
                traceExcerpt(:,:,TTLNum)=ephysData.traces(1:size(ephysData.traces,1),excerptWindow);
            end
            figure('Position',[214   108   747   754],'name',...
                [fileName ' Unit' num2str(ephysData.selectedUnits(cellNum))] ); %Ch' num2str(spikeData.selectedUnits(cellNum))
            hold on
            
            scaleF = max(max(max(abs(traceExcerpt))));
            traceExcerpt=traceExcerpt./scaleF;
            
            %             for TTLNum=1:length(TTLtimes)
            for traceNum=1:size(traceExcerpt,1)
                plot(squeeze(traceExcerpt(traceNum,:,:))+(traceNum-1),'color',[0 0 0 0.1]);
            end
            %             end
            
            %             traceExcerpt=mean(traceExcerpt,3);%(:,:,5)
            %
            %             traceExcerpt = array2table(traceExcerpt','VariableNames',...
            %                 cellfun(@(x) ['Ch' num2str(x)],num2cell(1:size(ephysData.traces,1)),'UniformOutput',false));
            %             stackedplot(traceExcerpt)
            
            %             for traceNum=1:size(traceExcerpt,1)
            %                 plot(zscore(traceExcerpt(traceNum,:))+15*(traceNum-1)+... %+max(get(gca,'ylim'))
            %                     max(abs(zscore(traceExcerpt(traceNum,:)))),'color','k');
            %             end
            
%             set(gca,'xticklabels',cellfun(@(x) round(str2num(x)/30,1)-5, get(gca,'xticklabels')))
            xline(SRR/200); xline(SRR/200+(pulseDur*SRR))
            xTicks=0:SRR/400:SRR/200+1.2*(pulseDur*SRR);
            set(gca,'xtick',xTicks,'xticklabel',xTicks/30-5)
            axis tight
            
            title(['Cell #' num2str(uIdx) ' '...
                ephysData.recInfo.baseName ' Unit '...
                num2str(ephysData.selectedUnits(cellNum))],...
                'interpreter','none')
            
            if savePlots
%                 savefig(gcf,[cd filesep 'Cell_' ... %'Z:\Vincent\Figures\vIRt\PT'
%                     num2str(uIdx) '_Overlaped_Triggered_Traces'],'compact');
                print(gcf,[cd filesep 'Cell_' num2str(uIdx) '_Overlaped_Triggered_Traces'],'-dpng')
                close(gcf)
            end
    end
    
    
    
    
    
end

end



