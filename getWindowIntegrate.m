function window_sum_series = getWindowIntegrate(out_idx,time_window,fs)
    
    time_sample  = time2sample(time_window,fs);
    unit_kernel = ones(1,time_sample);
    window_sum_series = conv2(unit_kernel,1,out_idx,'same');
    
    


end