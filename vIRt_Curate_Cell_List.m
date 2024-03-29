%% Curate vIRt recordings %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
%% define figure colormap
cmap=lines;cmap=[cmap(1:7,:);(lines+flipud(copper))/2;autumn];

% data should be in DJ
baseDir='D:\Vincent\';
load(fullfile(baseDir,'Analysis','vIRt_sessions_analysis.mat'));

vIRt_sess.("XP type")=categorical(vIRt_sess.("XP type"));
vIRt_sess.("Manuscript version")=categorical(vIRt_sess.("Manuscript version"));

cellQR=table('Size',[0,12],'VariableTypes',...
    {'categorical','categorical','double','double','logical','string','string',...
    'string','double','double','categorical','categorical'},...
    'VariableNames',...
    {'Subject','Session','RecordingID','unitIndex','keepIndex','unitTuning',...
    'tuningEpochs','unitPT','unitFrequency','unitQuality','XP type','Manuscript version'});

%% For Revision 1
manuscript_version='rev1';
rev1_vIRt_sessIdx=contains(vIRt_sess.KeepList,"1") &...
    vIRt_sess.("Manuscript version")==manuscript_version;
vIRt_sess=vIRt_sess(rev1_vIRt_sessIdx,:);

%% load data across sessions
for sessionNum=1:size(vIRt_sess,1)
    if ~(vIRt_sess.WR(sessionNum)== '')
        sessDir=fullfile(baseDir,char(vIRt_sess.Subject(sessionNum)),char(vIRt_sess.Session(sessionNum)));
        sessBaseName=[char(vIRt_sess.Session(sessionNum)) '_' num2str(vIRt_sess.RecordingID(sessionNum))];
        if ~exist(fullfile(baseDir,'Analysis','Data',sessBaseName),'dir')
            %% One time load and save data
            %
            sessData=vIRt_GetSessionData(sessDir,sessBaseName);
            [whiskers,bWhisk,breathing,ephys]=vIRt_TransformData(sessData.behav,sessData.ephys);
            
            % whisking epochs (based on first trace, if multiple whisker tracked)
            ampThd=18;%amplitude threshold
            freqThld=1; %frequency threshold
            minBoutDur=1000; %500; % 1000 % minimum whisking bout duration: 1s
            [whiskingEpochs,whiskingEpochsList]=vIRt_Epochs(whiskers,bWhisk,ampThd,freqThld,minBoutDur);
            
            % epoch mask (all tuning computation must be done outside of pulse stimulation)
            wEpochMask=vIRt_EpochMask(whiskers,bWhisk,sessData.pulses,sessData.behav.vidTimes,...
                whiskingEpochs,size(ephys.rasters,2));
            
            % Decide which units to keep based on notes and plots
            ephys.selectedUnits=str2num(vIRt_sess.KeepList(sessionNum));
            [ephys.keepList,ephys.qualityReport]=vIRt_CheckUnitSelection(whiskers,bWhisk,...
                ephys,wEpochMask,sessData.pulses);
            
            % Apply selection list
            keepWhat = 'handPick'; %keepList=[];
            ephys.selectedUnits=vIRt_UnitsToKeep(ephys,keepWhat);
            
            % save data
            ephys=vIRt_SaveData(sessBaseName,baseDir,...
                whiskers,bWhisk,whiskingEpochs,whiskingEpochsList,wEpochMask,...
                breathing,ephys,sessData);
            
            %Keep info on saved cells
            %Add cell session ID / remove bad ones
            unitIdx=ismember(ephys.qualityReport.unitIndex,ephys.selectedUnits);
            cellQR=vertcat(cellQR,...
                horzcat(repmat(vIRt_sess(sessionNum,1:3),sum(unitIdx),1),...
                ephys.qualityReport(unitIdx,:),repmat(vIRt_sess(sessionNum,10:11),sum(unitIdx),1)));
        end
        save(fullfile(baseDir,'Analysis',['Cell_List_' manuscript_version '.mat']),'cellQR');
    end
    clearvars whiskers bWhisk breathing ephys whiskingEpochs whiskingEpochsList wEpochMask
end

writetable(cellQR,fullfile(baseDir,'Analysis',['Cell_List_' manuscript_version '.xls']));
