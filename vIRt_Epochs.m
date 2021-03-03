function [whiskingEpochs,whiskingEpochsList]=vIRt_Epochs(whiskers,bWhisk,ampThd,freqThld,minBoutDur)

whiskingEpochs=WhiskingFun.FindWhiskingEpochs(...
    whiskers(bWhisk).amplitude,whiskers(bWhisk).frequency,...
    ampThd, freqThld, minBoutDur);
whiskingEpochs(isnan(whiskingEpochs))=false; %just in case
whiskingEpochsList=bwconncomp(whiskingEpochs);
[~,wBoutDurSort]=sort(cellfun(@length,whiskingEpochsList.PixelIdxList),'descend');
whiskingEpochsList.PixelIdxListSorted=whiskingEpochsList.PixelIdxList(wBoutDurSort);

% check epochs
if false
    figure; hold on;
    plot(whiskers(bWhisk).angle);
    plot(whiskingEpochs*nanstd(whiskers(bWhisk).angle)+nanmean(whiskers(bWhisk).angle))
    % plot(whisker.phase*nanstd(whisker.angle)/2+nanmean(whisker.angle));
end