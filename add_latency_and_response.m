function [tble_edt,act_pw] = add_latency_and_response(tble, th,tau,US_latency_included)
   fs = 30303; % sampling frequency
% dt = 1/fs ; % timestep
% time_non_cond = 5; % in  seconds [1 3 1] [pre stim post ] timings
time_cond = 7 ; % in  seconds [1 5 1]

time_cond_range = linspace(-1,time_cond-1,fs* time_cond);
% time_non_cond_range = linspace(-1,time_non_cond-1,fs* time_non_cond);

downsample_num = 20;
time_cond_range = downsample(time_cond_range,downsample_num);
% time_non_cond_range = downsample(time_non_cond_range,downsample_num);
fs_ds_cond = round(size(time_cond_range,2) / time_cond ); 
dt_ds = 1/fs_ds_cond;
% act = cell2mat(tble.act');
high_freq = 600;
act = getLowpassfiltered(tble,high_freq,fs_ds_cond);
% act_nart = removeArtifact(act,time_cond_range);




act_pw = getPower(act,fs_ds_cond,tau);
[response,latency,cutoff] = getResponse(act_pw,time_cond_range,th);
response(response == 0) = nan;

tble_edt = tble;
% tble_edt.power = act_pw;


tble_edt.response = response;
tble_edt.latency = latency;
tble_edt.cutoff = cutoff;
longer_cond_latency_idx = tble_edt.latency >= 2 & tble_edt.stage == "Abs_cond";
if US_latency_included == false
    tble_edt.latency(longer_cond_latency_idx) = nan;
    tble_edt.response(longer_cond_latency_idx) = nan;
end
longer_non_cond_latency_idx = tble_edt.latency >=4 & tble_edt.stage ~= "Abs_cond" ; % removing all the detection greater than stimulus presentation
tble_edt.latency(longer_non_cond_latency_idx) = nan;
tble_edt.response(longer_non_cond_latency_idx) = nan;

tble_edt.resp_num = tble_edt.response .* tble_edt.bee_num;

tble_edt.cutoff = cutoff;

tble_edt.power = num2cell(act_pw,1)';




end