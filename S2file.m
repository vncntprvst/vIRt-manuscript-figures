function ephys=S2file(ephys,option)

if contains(option,'save')
    prefElec=double(ephys.spikes.preferredElectrode(ismember(...
        ephys.spikes.unitID,ephys.selectedUnits(1))));
    keepTrace=mode(prefElec);
      
    trace=ephys.traces(keepTrace,:);
    fileID = fopen(fullfile('D:\Vincent\Analysis\Traces', [sessInfo.baseName...
        'Ch' num2str(keepTrace) '.bin']),'w');
    fwrite(fileID,trace,'int16');
    fclose(fileID);
end

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
