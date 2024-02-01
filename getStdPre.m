function stdPre = getStdPre(data,time,fs)
%     preTime = time <-1;     
    preTime = time < 0 ;            %### pretime changed
    data = data' ;%              ###changed to match the indexing 
    pool = data(:,preTime);
    totLen = size(data,1)*size(pool,2);
    pool = reshape(pool',1,totLen);
    
    spikes = diff([0 pool] > 5*std(pool))==1;
    
%     kernel = ones(1,0.01 * fs)+1; % kernel of 10ms around the spike 
    kernel = ones(1,round(0.01 * fs)) + 1; %  ### rounding fraction   
    aroundSpikes = conv(spikes,kernel);
%     aroundSpikes = aroundSpikes(length(kernel)/2:end-length(kernel)/2);
    aroundSpikes = aroundSpikes(round(length(kernel)/2):round(end-length(kernel)/2)); % ### rounding fractions

    filtPool = pool(~aroundSpikes);
%%
%     figure();
%     tplot = (0:totLen-1)/fs;
%         
%     subplot(211);
%     hold on;
%     plot(tplot,pool);
%     plot(tplot(spikes),800,'r*');
%     title(sprintf('pre-phase pooled signal, std = %.3f',std(pool)));
%     xlabel('time [s]');
%     ylabel('V [\muV]');
%     
%     subplot(212);
%     plot(tplot(~aroundSpikes),filtPool);
%     title(sprintf('pre-phase pooled signal after spike removal, std = %.3f',std(filtPool)));
%     xlabel('time [s]');
%     ylabel('V [\muV]');
%%    
    stdPre = std(filtPool);
end