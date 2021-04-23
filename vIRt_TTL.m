if isfield(pulses,'lag')
    pulses.TTLTimes=pulses.TTLTimes-pulses.lag;
end