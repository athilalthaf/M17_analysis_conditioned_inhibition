%% Script to test the m17 latency on small pilot experiment
load("m17_test_data.mat") % load the test data ; m17_conditioning only experiments
dt_test =    1.0200e-04; % timeperiods of sampling 

pre_onset_stim = 1; % pre stimulus time 
post_onset_stim = 6; % stimulus time the some stimulus goes less than 7 sec 
% that is why 6 was chosen


fs_test = round(1/dt_test); %sampling frequency
pre_stim = pre_onset_stim * fs_test ;  % getting the number of entries for pre stim time
stim_interval = post_onset_stim * fs_test; % getting same for stim duration
time_test_range = -pre_onset_stim:dt_test : post_onset_stim + dt_test; % the whole time range 
tau = 150;  % time constant for convolution
stdThPw = 15; % std thresh
trial_time = 8; % conditioning duration
stim_thresh = 1;  % thresholding to find the stim onset

ACTIVITY_arr = []; %activity array to get the activity based on the stim onset
% STIM_sum = [];
for i = 1:height(m17_test_data) 
    stim_data = cell2mat(m17_test_data.stim(i)); %get the stim data
    act_rawdata = cell2mat(m17_test_data.activity_raw(i)); % get the raw activity data
    onset_region = find(stim_data>stim_thresh); % get the stimulus timing  
    onset_flag = onset_region(1); %first element to get the onset
    
%     activity_slice =  act_rawdata(onset_flag - pre_stim: end);
    
    stim_slice = stim_data(onset_flag - pre_stim: end); % region where the stim 
    stim_dur = sum(stim_data > stim_thresh)/ fs_test; % get the stim duration
    if stim_dur > post_onset_stim % if the stim duration is long ie post_onset_stim
        activity_slice =  act_rawdata(onset_flag - pre_stim: onset_flag + stim_interval); %get the activity for that duration 
    else 
        activity_slice = NaN(size(time_test_range))'; % replace with nan
    end
%     
ACTIVITY_arr = [ACTIVITY_arr; {activity_slice}]; % append the activity array
%     STIM_sum = [STIM_sum ; stim_dur];
end


m17_test_data.activity = ACTIVITY_arr; % replace the activity with the 

% cs_plus_only = m17_test_data(m17_test_data.modality == "compound",:);

%% power method testing
test_dat = cell2mat(m17_test_data.activity');  % converting to matrix
test_dat_nart = removeArtifact_test(test_dat,time_test_range,[0 0.23]); % remove the artifact in this case  
 
test_dat_pw = getPower(test_dat_nart,dt_test,fs_test,tau); % get the rectified and convolved signal
[resp_pw,onset_pw, th_pw] = getResponse_test(test_dat_pw,time_test_range,stdThPw,trial_time); % get the onset response and threshold data

m17_test_data.latency_onset = onset_pw; % assign the onset to the table
m17_test_data.latency_th = th_pw; % assign the thresh
%%

% indices =1:30;
% % ids = 21:35;
% tiledlayout(10,3,"TileSpacing","tight","Padding","tight")
% for idx= indices
%     nexttile;
%     yyaxis left
%     plot(time_test_range,test_dat(:,idx));
%     xline(onset_pw(idx),"-k",'LineWidth',1.5);
%     xline(m17_test_data.latency(idx),"-m",'LineWidth',1.5)
%     yyaxis right
%     h = plot(time_test_range,test_dat_pw(:,idx));
% %     set(h,'Color',[h.Color, 0.3]);
%     yline(th_pw(idx),"-k",'LineWidth',1)
%     xlim([-1 6]);
% %     title("trial no: " + num2str(idx) + "stim: "+  Bee_m17_cond(fl).stim(idx))
%     title(num2str(idx)+ "; " + m17_test_data.modality(idx))
%     sgtitle(m17_test_data.filename(idx) +  "; tau =" + num2str(tau) + ";  pw thresh =" + num2str(stdThPw) )
% end

%%
grouped_trails_m17 = grpstats(m17_test_data,["modality","trial_num"], ...
    ["mean","std","numel"],"DataVars",["latency_onset","latency"]); % group and get the descriptive stats out of the true amd predicted latency
color = lines(7); % get the default colors
COLORMAT = color([4 5 7],:);  

figure(2)
pred_scatter = gscatter(grouped_trails_m17.trial_num,grouped_trails_m17.mean_latency_onset,grouped_trails_m17.modality,COLORMAT,"...");hold on;
true_scatter = gscatter(grouped_trails_m17.trial_num,grouped_trails_m17.mean_latency,grouped_trails_m17.modality,COLORMAT,"+++");

% errorbar(grouped_trails_m17.mean_latency_onset,grouped_trails_m17.std_latency_onset)
title("M17 Onset During Conditioning")
xlabel("Trial number")
ylabel("Mean Latency (s)")
set(pred_scatter,"LineStyle","-");
 set(pred_scatter,"LineWidth",2);
 set(true_scatter,"LineStyle","-");
 set(true_scatter,"LineWidth",1.5);

true_leg = unique(grouped_trails_m17.modality);
true_leg = true_leg([3 2 1]) + " real";

pred_leg = unique(grouped_trails_m17.modality); 

pred_leg = pred_leg([3 2 1]) + " prediction";   %%%%%%%%% order needs to be checked  by simple legend command

legend([pred_leg;true_leg])

hold off

%% 
% CountNonNan = @(x) sum(~isnan(x));

% grouped_trails_count = grpstats(m17_test_data,["modality","trial_num"],@(x) sum(~isnan(x)),"DataVars",["latency_onset","latency"]);
% same but for the counts 
% gives us an idea about the undetected and falsely detected the onsets
figure(3)
pred_count_scatter = gscatter(grouped_trails_m17.trial_num,grouped_trails_m17.numel_latency_onset,grouped_trails_m17.modality,COLORMAT,"...");hold on;
true_count_scatter = gscatter(grouped_trails_m17.trial_num,grouped_trails_m17.numel_latency,grouped_trails_m17.modality,COLORMAT,"+++");

true_leg = unique(grouped_trails_m17.modality);
true_leg = true_leg([3 2 1]) + " real";

pred_leg = unique(grouped_trails_m17.modality);
pred_leg = pred_leg([3 2 1]) + " prediction";   %%%%%%%%%

legend([pred_leg;true_leg])
set(true_count_scatter,"LineWidth",1);
set(true_count_scatter,LineStyle="-")
set(pred_count_scatter,"LineWidth",2);
set(pred_count_scatter,LineStyle="-")
title("M17 detection During Conditioning")
xlabel("Trial number")
ylabel("Count (#)")

