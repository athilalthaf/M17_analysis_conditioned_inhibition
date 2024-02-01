function [resp, onset, th] = getResponse(data,time,stdTh)
    intPre = time < 0;
    data = data' ; %  ## indexing correction
    time_cond = 7   ; %## total conditioning time in s
%     intStart = find(time == 0);
    [~,intStart] = min(abs(time - 0)); % ## the nearest index the corresponds to 0
%     intStop = find(time == 1.5);
    [~,intStop] = min(abs(time - time_cond -1 )); %#  the nearest the integer corresponding to stim end
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