% function filt = causalLowPass(data,dt,tau)
function filt = causalLowPass(data,dt,fs,tau)  %### added fs argument
%     kwidth = 500; % kernel width (ms)
    kwidth = .5; % ###kernel width in s
    
    
%     t = 0:dt*1000:kwidth;
%     t = 0:dt*fs:kwidth; % ## using the sampling frequency
    t = 0:kwidth*fs; % ## uses 
%     tau = tau * 10^-3;   % ## change the ms into second

    k = exp(-t / (tau * 10^-3 * fs));
    k = k ./ (sum(k));
    %figure;
    %plot(t,k);
    filt=[filter(k,1,data')]';%; zeros(1,size(data,1))]'; 
%     filt = filt.*(1/dt); % convert in Hz
    filt = filt.*fs; % ###convert in Hz   # in fs terms

end