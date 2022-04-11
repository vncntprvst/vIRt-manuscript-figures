function wEpochMask=vIRt_EpochMask(whiskers,bWhisk,pulses,vTimes,wEpochs,ephysArrayLength)
% Any tuning computation has to be done outside of stimulation periods
wEpochMask=cell(numel(bWhisk),1);
for wNum=1:numel(bWhisk)
    pulseMask=false(1,size(whiskers(bWhisk(wNum)).angle,2));
    pulseMask((round(pulses.TTLTimes(1)*1000):round(pulses.TTLTimes(end)*1000))-round(vTimes(1)*1000))=true;
    wEpochMask{wNum}.behav=wEpochs{wNum};wEpochMask{wNum}.behav(pulseMask)=false;
    wEpochMask{wNum}.ephys=false(1,ephysArrayLength);
    ephysMaskIdx=(0:numel(wEpochMask{wNum}.behav)-1)+round(vTimes(1)*1000);
    wEpochMask{wNum}.ephys(ephysMaskIdx)=wEpochMask{wNum}.behav;
end