function sessData=vIRt_GetSessionData(sessDir,sessBaseName)

processedDataFiles = dir(fullfile(sessDir, '**','*processedData.mat'));
sessFileIdx=cellfun(@(recName) contains(recName,sessBaseName),...
    {processedDataFiles.name});
sessData=load(fullfile(processedDataFiles(sessFileIdx).folder,processedDataFiles(sessFileIdx).name));
sessData.sessFiles=dir(processedDataFiles(sessFileIdx).folder);
if ~isfield(sessData.ephys,'traces')
    traceFile = fopen(fullfile(processedDataFiles(sessFileIdx).folder,...
        sessData.sessFiles(cellfun(@(x) contains(x,'traces'),{sessData.sessFiles.name})).name), 'r');
    if traceFile<0
        traceFile = fopen(fullfile(processedDataFiles(sessFileIdx).folder,...
            sessData.sessFiles(cellfun(@(x) contains(x,'rec.bin'),{sessData.sessFiles.name})).name), 'r');
    end
    sessData.ephys.traces = fread(traceFile,[sessData.ephys.recInfo.numRecChan,Inf],'single');
    fclose(traceFile);
end