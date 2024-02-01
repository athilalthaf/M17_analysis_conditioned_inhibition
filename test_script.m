load smple_data.mat

smple_per = smple_data(~isnan(smple_data.latency),:)

fs = 1000;
dt = 1/fs;
t = -1:dt:7;
%%

% a = cell2mat(smple_per.activity(2));
a = Bee_m17(1).raw_dat;
a_noart = removeArtifact(a,t);
% [a_std_pre,a_mean_pre] = getStdPre(a_noart,t,fs);
% a_spikes = getSpikeTrains(a_noart,a_mean_pre + a_std_pre * 15);
% 
% 
% a_power = getPower(b_down,dt,100);
% plot(a_power)

function data = removeArtifact(data,time)
    artTime = time>=0 & time<0.01;
    data(:,artTime) = repmat(mean(data,2),1,sum(artTime));
end


function stdPre = getStdPre(data,time,fs)
    preTime = time <-1; 
    pool = data(:,preTime);
    totLen = size(data,1)*size(pool,2);
    pool = reshape(pool',1,totLen);
    
    spikes = diff([0 pool] > 5*std(pool))==1;
    
    kernel = ones(1,0.01 * fs)+1; % kernel of 10ms around the spike 
    aroundSpikes = conv(spikes,kernel);
    aroundSpikes = aroundSpikes(length(kernel)/2:end-length(kernel)/2);
    filtPool = pool(~aroundSpikes);
    
    stdPre = std(filtPool);
end

function SpikeTrains = getSpikeTrains(data,th)
    %SpikeTrains = diff(repmat(sign(mean(data,2)),1,size(data,2)).*data > th,1,2) == 1;
    SpikeTrains = diff(data > th,1,2) == 1;
    SpikeTrains = [SpikeTrains zeros(size(SpikeTrains,1),1)];
end


function power = getPower(data,dt,tau)
   kwidth = 500; % kernel width in ms
   t = 0:dt*1000:kwidth;
   k = exp(-t / tau);
   k = k ./ sum(k);
   power=filter(k,1,data.^2')';  
end