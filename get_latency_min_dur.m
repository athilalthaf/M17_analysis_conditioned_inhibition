function tble_edt = get_latency_min_dur(tble,min_dur,timerange,max_latency)
% this function gives the response and latency from the power of the signal
% by excluding shorter singal crossing events. inputs are a table (that has
% the power, threshold and cutoff info ) , minimum duration in seconds that
% the signal has to stay crossed, timerange of the signal and maximum permi
% tted latency. Output is the table itself with the new columns for response and the onset of response.
tble_edt = tble;
RESP_arr = nan(1,height(tble)); % initialising the response arr
LATENCY_arr = nan(1,height(tble)); % initialising the onset arr

for i = 1:height(tble)
    pw = cell2mat(tble.power(i))'; % get the power
    th = tble.cutoff(i); % get the cutoff
    log_id = getthatrialid(tble,tble.bee_id(i),tble.stage(i),tble.stim(i),tble.ss_norm_num(i));
    thresh_crossing = pw>th; % logical indexing the thresh crossing events

    if sum(thresh_crossing) ~= 0 % all the signals where the thres is crossed

        thresh_boundaries = [0, diff(pw>th)];  % get the diff to get the 
                                               % onset and offset of thresh crossing info. positive values
                                               % respresent the onset of crossing and negatives vise versa.

        group_start = timerange(thresh_boundaries>0); % corresponding onset times  
        group_end = timerange(thresh_boundaries<0); % corresponding offset times
          
        duration = group_end - group_start(1:length(group_end)); % difference gives the duration
      
        longer_onsets = group_start(duration > min_dur); % get all the onsets
                                                         % that stays longer than minimum duration
        longer_onsets = longer_onsets(longer_onsets > 0 & longer_onsets < max_latency); % get the latencies after the 
                                                                                        % stim onset and before the maximum permitted latency

        if ~isempty(longer_onsets) % if longer onsets are present 
            first_onset = longer_onsets(1); % get the first of the longer onsets
            RESP_arr(i) = 1;                 
            LATENCY_arr(i) = first_onset;
        else % if longer onsets are not there
            RESP_arr(i) = nan; % no response is registered 
            LATENCY_arr(i) = nan; 
        end
    else % if there are no thresh crossing events
        RESP_arr(i) = nan; % no response is registered 
        LATENCY_arr(i) = nan; 

    end

    








%  tble(i,:) 
end
tble_edt.response = RESP_arr';

tble_edt.latency = LATENCY_arr';

tble_edt.Properties.Description = sprintf("latency window : [%d,%d] ,min duration for latency : %.2f",0,max_latency,min_dur);

end