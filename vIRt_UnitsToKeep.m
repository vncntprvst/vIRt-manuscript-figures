function selectedUnits=vIRt_UnitsToKeep(ephys,keepWhat)

switch keepWhat
    case 'mostFreq' %most frequent units
        % mostFrqUnits=EphysFun.FindBestUnits(ephys.spikes.unitID,1);%keep ones over x% spikes
        %most frequent units during whisking periods
        reconstrUnits=ephys.rasters(:,whiskingEpochs).*(1:size(ephys.rasters,1))';
        reconstrUnits=reshape(reconstrUnits,[1,size(reconstrUnits,1)*size(reconstrUnits,2)]);
        reconstrUnits=reconstrUnits(reconstrUnits>0);
        mostFrqUnits=EphysFun.FindBestUnits(reconstrUnits,1);
        keepUnits=ismember(unitList,mostFrqUnits);
        keepTraces=unique(ephys.spikes.preferredElectrode(ismember(ephys.spikes.unitID,unitList(keepUnits))));
        selectedUnits=find(keepUnits);
    case 'all' %all of them
        selectedUnits=ephys.unitList; %all units
    case 'SU' % only Single units
        [unitQuality,RPVIndex]=SSQualityMetrics(ephys.spikes);
        unitQuality=[unique(double(ephys.spikes.unitID)),unitQuality];
        unitIdx=ismember(ephys.spikes.unitID,unitQuality(unitQuality(:,2)>0.6,1));
        unitQuality(unitQuality(:,2)>0.6,3)=hist(double(ephys.spikes.unitID(unitIdx)),...
            unique(double(ephys.spikes.unitID(unitIdx))))/sum(unitIdx);
        qualityUnits=unitQuality(unitQuality(:,2)>0.6 & unitQuality(:,3)>0.01,:);
        selectedUnits=qualityUnits(:,1);
    case 'handPick' % set numbers
        selectedUnits= ephys.keepList; %[11;5;4;14;16;24;26;28][12; 26; 37]; [33; 30; 29; 17; 45; 40; 4; 36; 6; 11; 14];
end

if nargin > 2 && ~isempty(ephys.keepList)% add manually, e.g.,selectedUnits=[selectedUnits;54]; 1;2;19];
    selectedUnits=[selectedUnits;ephys.keepList];
    selectedUnits=unique(selectedUnits);
end