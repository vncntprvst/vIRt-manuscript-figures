function [keepList,qualityReport]=vIRt_CheckUnitSelection(whiskers,bWhisk,ephys,wEpochMask,pulses)

unitQuality=SSQualityMetrics(ephys.spikes);
unitIndex=unique(double(ephys.spikes.unitID));
unitIdx=ismember(ephys.spikes.unitID,unitIndex);%(unitQuality>0.6,1) no quality threshold here
unitFrequency=(hist(double(ephys.spikes.unitID(unitIdx)),...
    unique(double(ephys.spikes.unitID(unitIdx))))/sum(unitIdx))';
keepIndex=ismember(unitIndex,ephys.selectedUnits);

% make tuning and PT report
[unitTuning,tuningEpochs,unitPT]=deal(repmat("None",numel(unitIndex),1));
qualityReport=table(unitIndex,keepIndex,unitTuning,tuningEpochs,unitPT,unitFrequency,unitQuality);

keepList=ephys.selectedUnits;

for unitNum=1:numel(keepList)
    
    ephys.selectedUnits=keepList(unitNum);
    
    % PT
    if ~isfield(pulses,'duration'); pulses.duration=0.010; end
    PhotoTagPlots(ephys,pulses);
    
    % Tuning
    NBC_Plots_PhaseTuning(whiskers(bWhisk).angle,whiskers(bWhisk).phase,...
        ephys,wEpochMask,'whisking',false,false);
    answer = questdlg('Check individual epochs?','Phase Tuning', ...
        'Yes', 'No','No');
    switch answer
        case 'Yes'
            NBC_Plots_PhaseTuning(whiskers(bWhisk).angle,whiskers(bWhisk).phase,...
                ephys,wEpochMask,'whisking',true,false);
    end
    
    prompt = {'Whisker Phase Tuning','Best Epoch','Phototagged'};
    dlgtitle = ['Unit ' num2str(keepList(unitNum)) 'Quality'];
    dims = 1;
    definput = {'None','None','None'};
    UQvalues = inputdlg(prompt,dlgtitle,dims,definput);
    rowIdx=(qualityReport.unitIndex==keepList(unitNum));
    qualityReport.unitTuning(rowIdx)= UQvalues{1};
    qualityReport.tuningEpochs(rowIdx)= UQvalues{2};
    qualityReport.unitPT(rowIdx)= UQvalues{3};
    
    close all;
    
end

[keepList,ephys.selectedUnits]=deal(unitIndex(...
    ~contains(qualityReport.unitTuning,"None") |...
    ~contains(qualityReport.unitPT,"None")));



