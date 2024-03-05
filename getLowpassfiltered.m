function filtered = getLowpassfiltered(tble,freq,sample_fs)
    act = cell2mat(tble.act');
    act_zeros = act;
    act_zeros(isnan(act_zeros)) = 0;
    filtered = lowpass(act_zeros',freq,sample_fs);
    
    filtered = filtered';
    filtered(isnan(act)) = nan;
    

end