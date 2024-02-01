function SpikeTrains = getSpikeTrains(data,th)
    %SpikeTrains = diff(repmat(sign(mean(data,2)),1,size(data,2)).*data > th,1,2) == 1;
    SpikeTrains = diff(data > th,1,2) == 1;
    SpikeTrains = [SpikeTrains zeros(size(SpikeTrains,1),1)];

end

