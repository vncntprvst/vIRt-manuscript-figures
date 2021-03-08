%% vIRt paper figures %%
%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
%% define figure colormap
cmap=lines;cmap=[cmap(1:7,:);(lines+flipud(copper))/2;autumn];

% data should be in DJ
baseDir='D:\Vincent\';
% load(fullfile(baseDir,'Analysis','vIRt_sessions_analysis.mat'));

load(fullfile(baseDir,'Analysis','Cell_List.mat'))
tunedCells=find(cellList.Tuning==1);
PTCells=find(cellList.PT==1);

%% Figure 2. Population phase tuning
cellTuning=struct('phaseStats_all',[],'phaseStats_epochs',[],...
    'phaseTuning_all',[],'phaseTuning_epochs',[],'thetas',[]);

for cellNum=1:numel(tunedCells)
    uIdx=tunedCells(cellNum);
    sessID=[char(cellList.Session(uIdx)) '_' num2str(cellList.RecordingID(uIdx))];
    dataDir=fullfile(baseDir,'Analysis','Data',sessID);
    load(fullfile(dataDir,[sessID '_behavior.mat']),'whiskers','wEpochMask','bWhisk');
    load(fullfile(dataDir,[sessID '_ephys.mat']),'ephys');
    load(fullfile(dataDir,[sessID '_recInfo.mat']),'recInfo');
    ephys.selectedUnits=cellList.unitIndex(uIdx);
    ephys.recInfo=recInfo;
    [cellTuning(cellNum).phaseStats_all,cellTuning(cellNum).phaseTuning_all]=vIRt_PhaseTuning(whiskers(bWhisk).phase,ephys,wEpochMask,false);
    [cellTuning(cellNum).thetas,cellTuning(cellNum).phaseStats_epochs,cellTuning(cellNum).phaseTuning_epochs]=vIRt_PhaseTuning(whiskers(bWhisk).phase,ephys,wEpochMask,true);
%  circ_rtest(cellTuning(cellNum).thetas{24})
end

%% Overview plot
NBC_Plots_Overview(whiskers(bWhisk),whiskingEpochs,breathing,ephys,pulses.TTLTimes,false,false);

%% Check Phototagging summary
% ephys.selectedUnits=[60 23]; 10; 2; 37; %12;
if ~isfield(pulses,'duration'); pulses.duration=0.010; end
PhotoTagPlots(ephys,pulses); %Implement SALT test
% PTunits=[12,26,37];




%% population plots

%% PT plots

% waveforms
% traces
% rasters
% latency histogram + SALT
% tuning


%% transition rythm plot

%% Supplementary
% Slow rhythm tunings: setpoint / breathing
% Proba tuning plot