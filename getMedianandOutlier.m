function [act_hampel,out_idx] = getMedianandOutlier(act,time_window,nsigma)
    
[~,out_idx,act_med] = hampel(act',time_window,nsigma);


act_hampel = act_med';
act_hampel(out_idx') = act(out_idx');
act_hampel(isnan(act)) = nan;
out_idx = out_idx';

end