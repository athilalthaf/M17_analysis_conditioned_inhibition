% function power = getPower(data,dt,tau)
function power = getPower(data,dt,fs,tau) %## added fs for  taking account of downsampling
%    kwidth = 500; % kernel width in ms
   kwidth = .5; %## kernel width in s 
%    t = 0:dt*1000:kwidth;
%     fs = 30303 ; % ### introducing sampling frequency
   t = 0:kwidth * fs; % ## change based on sampling frequency
%    tau = ;   % ## change the ms into second
%    k = exp(-t / tau);
   k = exp(-t /(tau * 10^-3 * fs)); % ### chage the ms into second and the time window
   k = k ./ sum(k);
   power=filter(k,1,data.^2);  
   
end








