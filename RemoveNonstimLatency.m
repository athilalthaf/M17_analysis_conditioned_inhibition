function tble = RemoveNonstimLatency(tble,USincluded)
    longer_non_cond_latency_idx = tble.latency >=3 & tble.stage ~= "Abs_cond" ; % removing all the detection greater than stimulus presentation
    tble.latency(longer_non_cond_latency_idx) = nan;
    tble.response(longer_non_cond_latency_idx) = nan;


    if USincluded == false
       longer_cond_latency_idx =  tble.latency >= 2 & tble.stage == "Abs_cond";
           tble.latency(longer_cond_latency_idx) = nan;
           tble.response(longer_cond_latency_idx) = nan;
    end
end