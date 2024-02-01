function [resp, onset, th] = getResponse_test(data,time,stdTh,trial_dur)
    intPre = time < 0;
    data = data' ; %  ## indexing correction
%     trial_dur = 7   ; %## total conditioning time in s
%     intStart = find(time == 0);
    [~,intStart] = min(abs(time - 0)); % ## the nearest index the corresponds to 0
%     intStop = find(time == 1.5);
    [~,intStop] = min(abs(time - trial_dur )); %#  the nearest the integer corresponding to stim end
    intStim = intStart:intStop;

    avgPre = mean(data(:,intPre),2);
    stdPre= std(data(:,intPre),0,2);
    
    th = avgPre+stdTh.*stdPre;
    onset = NaN(size(data,1),1);
    resp = zeros(size(data,1),1);
    
    
    for i = 1:size(data,1)
        idxs = find(data(i,intStim) > th(i));
        if ~isempty(idxs)
            resp(i) = 1;
            onset(i) = time(intStim(1)+idxs(1)); 
        end
    end
end