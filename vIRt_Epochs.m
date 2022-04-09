function [whiskingEpochs,whiskingEpochsList]=vIRt_Epochs(whiskers,bWhisk,ampThd,freqThld,minBoutDur)

[whiskingEpochs,whiskingEpochsList]=deal(cell(numel(bWhisk),1));
for wNum=1:numel(bWhisk)
    whiskingEpochs{wNum}=WhiskingFun.FindWhiskingEpochs(...
        whiskers(bWhisk(wNum)).amplitude,whiskers(bWhisk(wNum)).frequency,...
        ampThd, freqThld, minBoutDur);
    whiskingEpochs{wNum}(isnan(whiskingEpochs{wNum}))=false; %just in case
    whiskingEpochsList{wNum}=bwconncomp(whiskingEpochs{wNum});
    [~,wBoutDurSort]=sort(cellfun(@length,whiskingEpochsList{wNum}.PixelIdxList),'descend');
    whiskingEpochsList{wNum}.PixelIdxListSorted=whiskingEpochsList{wNum}.PixelIdxList(wBoutDurSort);
end

if false
    figure; hold on;
    for wNum=1:numel(bWhisk)
    plot(whiskers(bWhisk(1)).angle);
    plot(whiskingEpochs{wNum}*nanstd(whiskers(bWhisk(1)).angle)+nanmean(whiskers(bWhisk(1)).angle))
    end
    % plot(whisker.phase*nanstd(whisker.angle)/2+nanmean(whisker.angle));
end