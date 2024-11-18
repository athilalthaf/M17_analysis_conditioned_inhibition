% load("Beem17_trim.mat");
% load("m17_conditioning_data.mat")
% load("m17_post_conditioning_data.mat")
load("m17_sep_tab.mat")
%%
stdThSt = 15; 
stdThFr = 5;
stdThPw = 5;

tau = 100;

fs = 30303; % sampling frequency
dt = 1/fs ; % timestep
time_non_cond = 5; % in  seconds [1 3 1] [pre stim post ] timings
time_cond = 7 ; % in  seconds [1 5 1]

time_cond_range = linspace(-1,time_cond-1,fs* time_cond);
time_non_cond_range = linspace(-1,time_non_cond-1,fs* time_non_cond);

downsample_num = 20;
time_cond_range = downsample(time_cond_range,downsample_num);
time_non_cond_range = downsample(time_non_cond_range,downsample_num);
fs_ds_cond = round(size(time_cond_range,2) / time_cond ); 
dt_ds = 1/fs_ds_cond;
cond_range = 20;
fig_pos = [-1279 50 1280 907]; %% the postion of the figure second monitor fullscreen
fig_pos_half  = [-1281 462 1279 486]; % half of the fullscreen
fig_pos_minimal = [-1116 316 699 603]; 
fig_pos_max = [1          41        1920         963];

addpath('C:\Users\athil\OneDrive - uni-bielefeld.de\Desktop\Codes\gramm_matlab\gramm\gramm\');
%%

% for fl =  1:length(Bee_m17_cond(1:10)) 
%     cond_data = Bee_m17_cond(fl).act(:,1:20);
%     cond_data_ds  = downsample(cond_data,downsample_num);
% 
%     cond_data_noArt = removeArtifact(cond_data_ds,time_cond_range);
% thr_std_pre = stdThSt * getStdPre(cond_data_noArt,time_cond_range,fs_ds_cond);
% sp_train = getSpikeTrains(cond_data_noArt,thr_std_pre);
% 
% sp_train_fr = causalLowPass(sp_train,dt_ds,fs_ds_cond,tau);
% 
% [resp_thresh,onset_thresh, th_thresh] = getResponse(cond_data_noArt,time_cond_range,stdThFr); %%
% 
% cond_dat_pw  = getPower(cond_data_noArt,dt_ds,fs_ds_cond,tau);
% [resp_pw,onset_pw, th_pw] = getResponse(cond_dat_pw,time_cond_range,stdThPw);
% 
% figure(fl)
% set(gcf,'Position',get(0,"ScreenSize"));
% 
% tiledlayout(10,2,"TileSpacing","tight","Padding","tight")
% for idx= 1: cond_range
%     nexttile;
%     yyaxis left
%     plot(time_cond_range,cond_data_ds(:,idx));
%     xline(onset_pw(idx),"-k",'LineWidth',1.5);
%     xline(onset_thresh(idx),"-m",'LineWidth',1.5)
%     yyaxis right
%     h = plot(time_cond_range,cond_dat_pw(:,idx));
% %     set(h,'Color',[h.Color, 0.3]);
%     yline(th_pw(idx),"-k",'LineWidth',1)
%     xlim([-1 6]);
% %     title("trial no: " + num2str(idx) + "stim: "+  Bee_m17_cond(fl).stim(idx))
%     title(num2str(idx)+ "; " + Bee_m17_cond(fl).stim(idx))
%     sgtitle(Bee_m17_cond(fl).bee_id +  "; tau =" + num2str(tau) + ";  pw thresh =" + num2str(stdThPw) )
% end
% 
% legend(["data","power onset", "std onset", " power", "power thresh" ])
% saveas(gcf,"plots\threshbee\thresh_45\" + Bee_m17_cond(fl).bee_id +"tau_150.png")
% close all
% end

%% Single snippet
% bee_num = randi(40);
% trial_num = randi(20);
% % bee_num = 1;
% % trial_num = 1;
% 
% 
% tau = 100;
%  cond_data = Bee_m17_cond(bee_num).act(:,1:20);
%     cond_data_ds  = downsample(cond_data,downsample_num);
% 
%     cond_data_noArt = removeArtifact(cond_data_ds,time_cond_range);
% thr_std_pre = stdThSt * getStdPre(cond_data_noArt,time_cond_range,fs_ds_cond);
% sp_train = getSpikeTrains(cond_data_noArt,thr_std_pre);
% 
% sp_train_fr = causalLowPass(sp_train,dt_ds,fs_ds_cond,tau);
% 
% [resp_thresh,onset_thresh, th_thresh] = getResponse(cond_data_noArt,time_cond_range,stdThFr); %%
% 
% cond_dat_pw  = getPower(cond_data_noArt,dt_ds,fs_ds_cond,tau);
% [resp_pw,onset_pw, th_pw] = getResponse(cond_dat_pw,time_cond_range,stdThPw);
% 
% yyaxis left
% plot(time_cond_range,cond_data_ds(:,trial_num),LineWidth=2);
% xline(onset_thresh(trial_num),"-m",'LineWidth',1.5,Alpha=.4)
% xline(onset_pw(trial_num),"-k",'LineWidth',1.5,Alpha=.4);
% yyaxis right
% h = plot(time_cond_range,cond_dat_pw(:,trial_num),LineWidth=2);
% set(h,'Color',[h.Color, 0.5]);
% yline(th_pw(trial_num),"-k",'LineWidth',1)
% xlim([-1 6]);
% %     title("trial no: " + num2str(idx) + "stim: "+  Bee_m17_cond(fl).stim(idx))
% sgtitle(Bee_m17_cond(bee_num).bee_id +  "; tau =" + num2str(tau) + ";  pw thresh =" + num2str(stdThPw) + ...
%         " trial no: " + num2str(trial_num) + "stim: "+  Bee_m17_cond(bee_num).stim(trial_num))

% figure(2)
% s =stackedplot(time_cond_range,cond_data_noArt,"HandleVisibility","off");
% ax = findobj(s.NodeChildren,"Type","Axes");
% arrayfun(@(h,x) xline(h,x,'color','g','LineWidth',1.5),ax,onset_pw)
% arrayfun(@(h,x) xline(h,x,'color','m',LineWidth=1.5),ax,onset_thresh)

% cellfun(@(h,x) xline(h,x),ax,onset_thresh)

%% 20 random plots with same tau and different stdPw

% stdThPw = [5 15 45];
% tau = 50;
% 
% % bee_num = randi(40,1,5);
% % trial_num = randi(20,1,4);
% bee_num =     [3    32    23    31    19];
% trial_num = [8    16    11    19];
% 
% 
% 
% fig_num = 1;
% for bn = bee_num
%     for tr = trial_num
%     cond_data = Bee_m17_cond(bn).act(:,1:20);
%     cond_data_ds  = downsample(cond_data,downsample_num);
% 
%     cond_data_noArt = removeArtifact(cond_data_ds,time_cond_range);
% thr_std_pre = stdThSt * getStdPre(cond_data_noArt,time_cond_range,fs_ds_cond);
% sp_train = getSpikeTrains(cond_data_noArt,thr_std_pre);
% 
% sp_train_fr = causalLowPass(sp_train,dt_ds,fs_ds_cond,tau);
% 
% [resp_thresh,onset_thresh, th_thresh] = getResponse(cond_data_noArt,time_cond_range,stdThFr); %%
% 
% cond_dat_pw  = getPower(cond_data_noArt,dt_ds,fs_ds_cond,tau);
% ONSET_PW = [];
% TH_PW = [];
% for thresh = stdThPw
% [resp_pw,onset_pw, th_pw] = getResponse(cond_dat_pw,time_cond_range,thresh);
% ONSET_PW = [ONSET_PW, onset_pw(tr)];
% TH_PW = [TH_PW, th_pw(tr)];
% 
% end
% 
% figure(fig_num)
% yyaxis left
% plot(time_cond_range,cond_data_ds(:,tr),LineWidth=2,DisplayName="m17 data");
% xline(onset_thresh(tr),"-m",'LineWidth',1.5,Alpha=.4,DisplayName="std onset")
% if sum(isnan(ONSET_PW)) ~= 3
% 
%     xline(ONSET_PW(~isnan(ONSET_PW)),"-k",string(stdThPw(~isnan(ONSET_PW))),'LineWidth',1.5,Alpha=.4,DisplayName="pw onset");
% end
% % scatter(ONSET_PW,[0 0 0]) 
%  yyaxis right
%  h = plot(time_cond_range,cond_dat_pw(:,tr),LineWidth=2,DisplayName="filtered m17");
% % set(h,'Color',[h.Color, 0.5]);
% %  yline(th_pw(tr),"-k",'LineWidth',1)
% xlim([-1 6]);
% %     title("trial no: " + num2str(idx) + "stim: "+  Bee_m17_cond(fl).stim(idx))
% sgtitle(Bee_m17_cond(bn).bee_id +  "; tau =" + num2str(tau) + ";" + ...
%         " trial no: " + num2str(tr) + " stim: "+  Bee_m17_cond(bn).stim(tr))
%     fig_num = fig_num + 1;
%     legend(Location="best");
%     saveas(gcf,"plots\diff_thresh_t50\" + Bee_m17_cond(bn).bee_id+"_tr_"+num2str(tr)+"_th5_10_15_t"+num2str(tau) +".png")
%     end
%     
% end


%% testing previous m17 data


%  figure(1)
% sp = stackedplot(time_test_range,test_dat(:,ids))
% ax = findobj(sp.NodeChildren,"Type","Axes");
% 
% % arrayfun(@(h,x)xline(h,x,'color','k','LineWidth',1.5),ax,onset_pw(ids))
% % arrayfun(@(h,x)yline(h,x,'color','k','LineWidth',1.5),ax,th_pw(ids))
% 
%  arrayfun(@(h,x)xline(h,x,'color','r','LineWidth',1.5),ax,m17_test_data.latency(ids))
% 
% % figure(2)
% % sp2= stackedplot(time_test_range,test_dat_nart(:,ids))
% % ax2 = findobj(sp2.NodeChildren,"Type","Axes");
% arrayfun(@(h,x)xline(h,x,'color','k','LineWidth',1.5),ax2,onset_pw(ids))
% arrayfun(@(h,x)xline(h,x,'color','r','LineWidth',1.5),ax2,m17_test_data.latency(ids))
% 
% % arrayfun(@(h,x)yline(h,x,'color','k','LineWidth',1.5),ax,th_pw(ids))
%% condtioning data arizona
mixed_beeid = m17_conditioning_data.bee_id(m17_conditioning_data.stim == "M");

% clean_tb = m17_conditioning_data(m17_conditioning_data.bee_id ~=
% mixed_beeid',:);
stdThPw = 5;
tau = 100;
clean_cond_tb = m17_conditioning_data(~ismember(m17_conditioning_data.bee_id,mixed_beeid),:);
clean_postcond_tb = m17_post_conditioning_data(~ismember(m17_post_conditioning_data.bee_id,mixed_beeid),:);

clean_cond_act = cell2mat(clean_cond_tb.act');
clean_postcond_act = cell2mat(clean_postcond_tb.act');

clean_cond_act_nart = removeArtifact(clean_cond_act,time_cond_range);
cond_dat_pw = getPower(clean_cond_act_nart,dt_ds,fs_ds_cond,tau);
[resp_pw,onset_pw,th_pw] = getResponse(cond_dat_pw,time_cond_range,stdThPw);
% clean_cond_tb.latency = onset_pw;
context = clean_cond_tb.stim == clean_cond_tb.PreExp;
context_arr = strings(size(context));
context_arr(context) = "familiar";
context_arr(~context) = "novel";
clean_cond_tb.latency = onset_pw;
clean_cond_tb.context = context_arr;
%%
h = boxplot(onset_pw,context_arr);
boxes = findobj(h,"Tag","Box");
cols = lines(5);
cols= cols(end-1:end,:)
for i=1:length(boxes)
    patch(get(boxes(i),"XData"),get(boxes(i),"YData"),cols(i,:),"FaceAlpha",0.85);
end
title("Latency during conditioning")
ylabel("time (s)")
% clean_act = 
% cond_dat = 
grouped_cond = grpstats(clean_cond_tb,["stim","trial_num"],["mean","median","numel"],'DataVars','latency');
gsc = gscatter(grouped_cond.trial_num,grouped_cond.mean_latency,grouped_cond.stim);

set(gsc,LineStyle="-")
set(gsc,"LineWidth",2);

gsc2 = gscatter(grouped_cond.trial_num,grouped_cond.numel_latency./grouped_cond.GroupCount,grouped_cond.stim)


gs = gscatter(clean_cond_tb.trial_num - 100,clean_cond_tb.latency,clean_cond_tb.context,cols)
xlim([0 22]);
% clean_cond_tb.trial_num  = clean_cond_tb.trial_num - 100
w= linspace(min(clean_cond_tb.trial_num-100),max(clean_cond_tb.trial_num -100))
fit = fitlm(clean_cond_tb,"latency ~ trial_num + context",CategoricalVars="context");
line(w,feval(fit,w,"familiar"),'Color',cols(1,:),"LineWidth",2)
line(w,feval(fit,w,"novel"),'Color',cols(2,:),'LineWidth',2);

title("Latency during conditioning")
xlabel("trial num (#)")
ylabel("time (s)")

%% stacked plot for each beeid

% bees_id = unique(clean_cond_tb.bee_id);
% cond_trial_start = 100;
% % 
% % for id = 2%length(bees_id)
% %     figure(id)
% %     set(gcf, 'Position', get(0, 'Screensize'));
% % 
% %     sp = stackedplot(time_cond_range,clean_cond_act_nart(:,clean_cond_tb.bee_id==bees_id(id)), ...
% %         DisplayLabels=string(clean_cond_tb.trial_num(clean_cond_tb.bee_id==bees_id(id)) - cond_trial_start))
% %     ax = findobj(sp.NodeChildren,"Type","Axes")
% %     arrayfun(@(h,x) xline(h,x,'color','r','LineWidth',1.5),ax,clean_cond_tb.latency(clean_cond_tb.bee_id ==bees_id(id)))
% %     saveas(gcf,"plots\thresh_5_t150_cond\" + bees_id(id) +"_t5_tau_100.png")
% %     title(bees_id(id))
% % 
% %     close all
% % end
% 
% % tiledlayout(10,2,"TileSpacing","tight","Padding","tight")
% thresh_color = lines(7);
% thresh_color = thresh_color(7,:);
% for id = 1:length(bees_id)
%     figure(id);
% 
%     single_bee_condnum = length(clean_cond_tb.bee_id(clean_cond_tb.bee_id == bees_id(id)));
%     tiledlayout(single_bee_condnum,1,'TileSpacing','none','Padding','tight');
%     cond_trial_start = min(clean_cond_tb.trial_num(clean_cond_tb.bee_id == bees_id(id))) - 1;
%     for trl = 1:single_bee_condnum
%         nexttile;
%         single_trial = clean_cond_tb.bee_id==bees_id(id) & clean_cond_tb.trial_num - cond_trial_start==trl;
%         ylabel(string(clean_cond_tb.trial_num(single_trial)))
%         plot(time_cond_range,clean_cond_act_nart(:,single_trial));
%         xline(clean_cond_tb.latency(single_trial),Color=thresh_color,LineWidth=1.5,Alpha=0.9)
%     end
%             sgtitle(bees_id(id))
% %     saveas(gcf,"plots\thresh_5_t150_cond\" + bees_id(id) +"_t5_tau_100.png")
%     close all
% end



%% post conditioning

% clean_postcond_act_nart = removeArtifact(clean_postcond_act,time_cond_range);
% postcond_dat_pw = getPower(clean_postcond_act_nart,dt_ds,fs_ds_cond,tau);
% [resp_post_pw,onset_post_pw,th_post_pw] = getResponse(postcond_dat_pw,time_cond_range,stdThPw);
% % clean_cond_tb.latency = onset_pw;
% context = clean_postcond_tb.stim == clean_postcond_tb.PreExp;
% context_arr = strings(size(context));
% context_arr(context) = "familiar";
% context_arr(~context) = "novel";
% context_arr(clean_postcond_tb.stim == "M") = "mix";
% clean_postcond_tb.latency = onset_post_pw;
% clean_postcond_tb.context = context_arr;
% % clean_postcond_tb = clean_postcond_tb(~isnan(clean_postcond_tb.latency),:);
% h = boxplot(clean_postcond_tb.latency,clean_postcond_tb.context);
% boxes = findobj(h,"Tag","Box");
% cols = lines(5);
% cols= cols(end-2:end,:)
% for i=1:length(boxes)
% 
%     patch(get(boxes(i),"XData"),get(boxes(i),"YData"),cols(i,:),"FaceAlpha",0.85);
% end
% title("Latency during conditioning")
% ylabel("time (s)") 

%% All data

context = m17_sep_tab.stim == m17_sep_tab.PreExp;
context_arr = strings(size(context));
context_arr(context) = "familiar";
context_arr(~context) = "novel";
context_arr(m17_sep_tab.stim == "M") = "mix";



m17_sep_tab.context = context_arr;

m17_sep_tab.context = categorical(m17_sep_tab.context);
m17_sep_tab.stage = categorical(m17_sep_tab.stage);
m17_sep_tab.stim = categorical(m17_sep_tab.stim);
m17_sep_tab.bee_id = categorical(m17_sep_tab.bee_id);
m17_sep_tab.PreExp = categorical(m17_sep_tab.PreExp);
% clean_m17 = m17_sep_tab(~isnan(m17_sep_tab.latency),:);



m17_sep_tab.bee_num =  grp2idx(m17_sep_tab.bee_id);

mixed_stim_ids = m17_sep_tab.bee_id(m17_sep_tab.stage == "Abs_cond" & m17_sep_tab.stim == "M"); % The stims with M during conditions
m17_sep_tab_clean = m17_sep_tab(~ismember(m17_sep_tab.bee_id,mixed_stim_ids),:);

m17_sep_tab_clean = m17_sep_tab_clean(~(m17_sep_tab_clean.ss_norm_num >10 & m17_sep_tab_clean.stage ~= "Pre_exp"),:);

m17_sep_tab_raw = m17_sep_tab; %% just for all the m17data to be avaiable to cross check

m17_sep_tab = m17_sep_tab_clean;
addpath('C:\Users\athil\OneDrive - uni-bielefeld.de\Desktop\Codes\gramm_matlab\gramm')
whole_activity = cell2mat(m17_sep_tab.act');
whole_activity_nart = removeArtifact(whole_activity,time_cond_range);

whole_activity_pw = getPower(whole_activity_nart,fs_ds_cond,tau);
THRESH = 4:2:20;
% figure() %% geting the dimension and position in a desired format 
% fig_pos_half = get(gcf,'Position') %% after setting the command we will
% get the bounding box (i guess ) use it for the figure position


fig_pos = [-1279 50 1280 907]; %% the postion of the figure second monitor fullscreen
fig_pos_half  = [-1281 462 1279 486]; % half of the fullscreen
fig_pos_minimal = [-1116 316 699 603]; 
plot_path  = "plots\PER_thresh\post_troubleshoot\"; 
tau = 100;
%%
for th = 5%THRESH
[response,latency,cutoff] = getResponse(whole_activity_pw,time_cond_range,th);
m17_sep_tab.latency = latency;
m17_sep_tab.response = response;
longer_cond_latency_idx = m17_sep_tab.stage == 'Abs_cond' & m17_sep_tab.latency >=2;
m17_sep_tab.latency(longer_cond_latency_idx) = nan;

m17_sep_tab.response(longer_cond_latency_idx) = nan;

m17_sep_tab.resp_num = m17_sep_tab.response .* m17_sep_tab.bee_num; 
m17_sep_tab.resp_num(m17_sep_tab.resp_num == 0) = nan;
cond_tab  = m17_sep_tab(m17_sep_tab.stage == "Abs_cond",:);
postcond_tab  = m17_sep_tab(m17_sep_tab.stage == "Post_condtest",:);
precond_tab  = m17_sep_tab(m17_sep_tab.stage == "Post_exptest",:);

figure();
set(gcf,'Position',fig_pos);
g = gramm('x',m17_sep_tab.trial_num,'y',m17_sep_tab.resp_num,'color',m17_sep_tab.context);
g.facet_grid(m17_sep_tab.stim,[]);
g.geom_point();
g.set_title(sprintf("PER across trials and odour context (thresh= %d,tau= %d)",th,tau));
g.set_names('x','Trial num (#)','y',"Bee num (#)",'row','Stim','color','Context');
g.geom_vline('xintercept',[30.5,70.5,100.5,120.5],'style','--k')
g.draw(); %
export_path = "D:\ARIZONA BEES for Athil\learning\M17Analysis\M17_analysis_conditioned_inhibition\plots\PER_thresh";
saveas(gcf,plot_path + sprintf("overall_per_thresh_%d_tau_%d.png",th,tau))

figure();
set(gcf,'Position',fig_pos_minimal);
g = gramm('x',m17_sep_tab.ss_norm_num,'y',m17_sep_tab.resp_num,'color',m17_sep_tab.context,'subset',m17_sep_tab.stage~="Pre_exp" & m17_sep_tab.stage~="Pre_test");
g.facet_grid(m17_sep_tab.stim,m17_sep_tab.stage);
g.set_order_options("column",["Pre_test","Pre_exp","Post_exptest","Abs_cond","Post_condtest"])
g.geom_point();
g.set_title(sprintf("PER across trials and odour context (thresh= %d,tau= %d)",th,tau));
g.set_names('x','Trial num (#)','y',"Bee num (#)",'row','Stim','color','Context','column',"stage");
% g.geom_vline('xintercept',[30.5,70.5,100.5,120.5],'style','--k')
g.draw(); %
% export_path = "D:\ARIZONA BEES for Athil\learning\M17Analysis\M17_analysis_conditioned_inhibition\plots\PER_thresh";
% 
saveas(gcf,plot_path + sprintf("phase2_per_thresh_%d_tau_%d.png",th,tau))


figure()
set(gcf,'Position',fig_pos_minimal);

agg = grpstats(m17_sep_tab,["stage","stim","context","ss_norm_num"],["mean","numel"],DataVars="resp_num");
gr = gramm('x',agg.ss_norm_num,'y',agg.numel_resp_num,'color',agg.context,'subset',agg.stage~="Pre_exp" & agg.stage~="Pre_test");

gr.set_order_options('column',["Post_exptest","Abs_cond","Post_condtest"]);

gr.set_names('x','Trial num (#)','y',"Count (#)",'row','Stim','color','Context','column',"Stage");
gr.facet_grid(categorical(agg.stim),categorical(agg.stage));

gr.set_title(sprintf("PER Response  (thresh= %d, tau= %d)",th,tau));
% gr.set_text_options(facet_scaling=1)
gr.geom_point();
gr.geom_line();
gr.draw();
saveas(gcf,plot_path + sprintf("per_count_th_%d_tau_%d.png",th,tau));



figure()
set(gcf,'Position',fig_pos_half);

agg = grpstats(m17_sep_tab,["stage","context","ss_norm_num"],["mean","numel"],DataVars="resp_num");
gr = gramm('x',agg.ss_norm_num,'y',agg.numel_resp_num,'color',agg.context,'subset',agg.stage~="Pre_exp" & agg.stage~="Pre_test");

gr.set_order_options('column',["Post_exptest","Abs_cond","Post_condtest"]);

gr.set_names('x','Trial num (#)','y',"Count (#)",'color','Context','column',"Stage");
gr.facet_grid([],categorical(agg.stage));

gr.set_title(sprintf("PER Response  (thresh= %d,tau= %d)",th,tau));
% gr.set_text_options(facet_scaling=1)
gr.geom_point();
gr.geom_line();
gr.set_text_options('font','Helvetica');
gr.draw();
saveas(gcf,plot_path + sprintf("per_count_total_th_%d_tau_%d.png",th,tau));



end

%% troubleshooting 

% tab_check = m17_sep_tab(m17_sep_tab.latency <= 0.5 & m17_sep_tab.stage == "Abs_cond",:);
cutoff_trial_num = 5;
% tab_check = m17_sep_tab(m17_sep_tab.ss_norm_num <= cutoff_trial_num & m17_sep_tab.stage == "Abs_cond" ...
%                         & ~isnan(m17_sep_tab.resp_num) ,:); 
% tab_check = tab_check(1:40,:);
% th = 15;
% tau = 100;
% 
onset_set_path = "plots\onset_detection_abscond\";

% m17_sep_dat = cell2mat(m17_sep_tab.act');
% 
% check_data_nart = removeArtifact(m17_sep_dat,time_cond_range);
% % whole_activity_pw = getPower(whole_activity_nart,dt_ds,fs_ds_cond,tau);
% m17_sep_pw = getPower(check_data_nart,dt_ds,fs_ds_cond,tau);
% [response,latency,cutoff] = getResponse(m17_sep_pw,time_cond_range,th);
% m17_sep_tab.latency = latency;
% m17_sep_tab.response = response;
% m17_sep_tab.response(m17_sep_tab.response == 0) = nan;
% longer_cond_latency_idx = m17_sep_tab.latency >= 2 & m17_sep_tab.stage == "Abs_cond";
% 
% m17_sep_tab.latency(longer_cond_latency_idx) = nan;
% m17_sep_tab.response(longer_cond_latency_idx) = nan;
% m17_sep_tab.resp_num  = m17_sep_tab.bee_num .* m17_sep_tab.response ;
% 
% stage_checked = "Abs_cond";
% check_conditions = m17_sep_tab.stage == stage_checked & isnan(m17_sep_tab.latency);
% tab_check = m17_sep_tab(check_conditions,:);
% m17_sep_dat = m17_sep_dat(:,check_conditions);
% m17_sep_pw = m17_sep_pw(:,check_conditions);
% cutoff = cutoff(check_conditions);
% sample_plot_num = 20;
stage_needed = 'Abs_cond';
bee_id_needed = 'bee146';

stage_and_id = m17_sep_tab.bee_id == bee_id_needed & m17_sep_tab.stage == stage_needed;
ss_norm_num_max = 10;
subplot_num = 20;
fg1= figure(1);
set(fg1,'Position',fig_pos);
% tiledlayout(sample_plot_num,1,"TileSpacing",'none');
tiledlayout(ss_norm_num_max,2,'TileSpacing','compact');

bee_tab =sortrows(m17_sep_tab(stage_and_id,:),{'ss_norm_num','stim'});
current_bee_num = bee_tab.bee_num(1);
current_bee_exp = bee_tab.PreExp(1);
th = 3;

[bee_tab_processed,pw] = add_latency_and_response(bee_tab,th,tau,true);
% bee_tab_processed = sortrows(bee_tab_processed,"stim");
for idx= 1:subplot_num
    nexttile;

    sp = plot(time_cond_range,cell2mat(bee_tab_processed.act(idx))',LineWidth=1);
    xline(bee_tab_processed.latency(idx),"-r",'LineWidth',2);
%     xline(onset_thresh(idx),"-m",'LineWidth',1.5)
    ylabel(sprintf('stim %s, trial num %d',bee_tab_processed.stim(idx),bee_tab_processed.ss_norm_num(idx)),HorizontalAlignment="right",Rotation=0);
%     ylh = get(sp,'ylabel'); 
%     ylp = get(ylh,'Position');
%     set(ylh,Posi)
    ylim([-2 * 10^(-3) 2*10^(-3)]);
    yyaxis right
    h = plot(time_cond_range,pw(:,idx)');
    set(h,'Color',[h.Color, 1]);
    yline(bee_tab_processed.cutoff(idx),"-k",'LineWidth',1)
    xlim([-1 6]);

%     ylabel()
%     xticks([])
%     title("trial no: " + num2str(idx) + "stim: "+  Bee_m17_cond(fl).stim(idx))
%     title(num2str(idx)+ "; " + Bee_m17_cond(fl).stim(idx))
%     sgtitle(Bee_m17_cond(fl).bee_id +  "; tau =" + num2str(tau) + ";  pw thresh =" + num2str(th) )
    
end

% f1 =figure()
% set(gcf,'Position',fig_pos)
% sp = stackedplot(time_cond_range,act_data)
% sp.DisplayLabels = string(tab_check.bee_id) +" "+ string(tab_check.stim) +" " +string(tab_check.ss_norm_num) 
% 

sgtitle(sprintf("id = %s stage= %s,  thresh= %d  , tau= %d",bee_id_needed,replace(stage_checked,"_"," "),th,tau))
xlabel("time (s)")

onset_detected = length(bee_tab_processed.latency(~isnan(bee_tab_processed.latency)));
per_detected = length(bee_tab_processed.latency(bee_tab_processed.latency<2));
sprintf("onsets detected : %d \nper detected : %d",onset_detected,per_detected)
saveas(fg1,onset_set_path + sprintf("%s_stage_%s.png",stage_needed, bee_id_needed));


mini_fig_pos = [-771   448   717   148];
fg2 = figure(2);
set(fg2,'Position',mini_fig_pos);
single_bee_tab = m17_sep_tab(m17_sep_tab.bee_id == bee_id_needed,:);
single_bee_tab_process = add_latency_and_response(single_bee_tab,th,tau,false);
gscatter(single_bee_tab_process.trial_num,single_bee_tab_process.response,single_bee_tab_process.context);

title(sprintf("id : %s, th = %d, tau = %d, exp = %s",bee_id_needed,th,tau,current_bee_exp));
yticks([]);
ylabel(sprintf("# %d",current_bee_num));
xlim([0,151]);
xline([30.5,70.5,100.5,120.5])
legend("Location","eastoutside");
saveas(fg2,onset_set_path + sprintf("allstage_%s.png", bee_id_needed));

% ax = findobj(sp.NodeChildren,"Type","Axes")
% arrayfun(@(h,x) xline(h,x,'color','r','LineWidth',1.5,'Alpha',.6),ax,tab_check.latency)
%% threshold tuned
thresh_list = readtable("thresh_list.xlsx");
thresh_list.bee_id = categorical(thresh_list.bee_id);
vnames = ["bee_id"	"stim"	"trial_num"	"stage"	"stim_num"	"ss_norm_num"	"time_stamp"	"act"	"PreExp"	"context"	"bee_num"	"latency"	"response"	"resp_num"	"cutoff" "power"];
allbee_processed_tab  = array2table(nan(0,length(vnames)),'VariableNames',vnames);
allbee_processed_tab = convertvars(allbee_processed_tab,["bee_id","stim","stage","PreExp","context"],'categorical');

for entry = 1:height(thresh_list)
   single_bee_tab = m17_sep_tab(m17_sep_tab.bee_id == thresh_list.bee_id(entry),:);
   single_bee_tab_process  = add_latency_and_response(single_bee_tab,thresh_list.th(entry),tau,true);
   allbee_processed_tab  = [allbee_processed_tab ; single_bee_tab_process];
end 
OutlierCount = 5;
% allbee_processed_tab.latency = NaN;
high_freq = 600;
nsigma = 3;
act = getLowpassfiltered(allbee_processed_tab,high_freq,fs_ds_cond);

[medout,outidx] = getMedianandOutlier(act,500,nsigma);
stim_time = time_cond_range >= 0 & time_cond_range <= 2;
stim_time_cond = time_cond_range >= 0 & time_cond_range <= 2;
win_sum_sr = getWindowIntegrate(outidx,500,fs_ds_cond);
% % 
win_sum_sr(~stim_time,:) = nan; 
win_sum_sr(~stim_time_cond,allbee_processed_tab.stage=="Abs_cond") = nan;
% % 
artifact_filter = max(win_sum_sr,[],1) < OutlierCount;
% % % act_int  = act(:,artifact_filter);

artifact_tab = allbee_processed_tab(artifact_filter,:); % useful while checking whether the artifact filtering is working or not;
artifact_tab = artifact_tab(artifact_tab.response == 1,:)

allbee_processed_tab.latency(allbee_processed_tab.latency>2) = nan;
allbee_processed_tab.response(isnan(allbee_processed_tab.latency)) = nan;
allbee_processed_tab.resp_num(isnan(allbee_processed_tab.latency)) = nan;
allbee_processed_tab.response(artifact_filter) = nan;
allbee_processed_tab.resp_num(artifact_filter) = nan;
allbee_processed_tab.latency(artifact_filter) = nan;

save("allbee_processed_tab.mat","allbee_processed_tab")



% allbee_processed_tab.
% save(allbee_processed_tab,)


%% artifact filter 



%% sort by latent inhibition score Overall count
load("allbee_processed_tab.mat")
% LatInhib_agg = grpstats(allbee_processed_tab,["bee_id","context"],"numel",DataVars="latency");
LatInhib_agg = grpstats(allbee_processed_tab,"bee_id","numel",DataVars="latency");
% allbee_processed_tab.li_score = []; % use this if the table already has these variables
% allbee_processed_tab.li_sort = []; 
% allbee_processed_tab.sort_num = [];

% LI_SCORE = nan(size(unique(allbee_processed_tab.bee_id)))
LI_SCORE = [];
TOT_RESP = [];
for bee = unique(allbee_processed_tab.bee_id)'
%     f_count = LatInhib_agg.numel_latency(LatInhib_agg.bee_id == bee & LatInhib_agg.context=="familiar");
%     n_count = LatInhib_agg.numel_latency(LatInhib_agg.bee_id == bee & LatInhib_agg.context=="novel");
%     m_count = LatInhib_agg.numel_latency(LatInhib_agg.bee_id == bee & LatInhib_agg.context=="mix");
    counts = LatInhib_agg.numel_latency(LatInhib_agg.bee_id == bee);
%     li_score =  (f_count - n_count)/(f_count + n_count); 
    LI_SCORE(end + 1) = counts; %% count

%     LI_SCORE(end + 1) = li_score;
end

li_tab = table(categorical(unique(allbee_processed_tab.bee_id)),LI_SCORE',VariableNames=["bee_id","li_score"]);
li_tab = sortrows(li_tab,"li_score","ascend");
li_tab.li_sort = (1:height(li_tab))';
allbee_processed_tab = join(allbee_processed_tab,li_tab);

allbee_processed_tab.sort_num = allbee_processed_tab.response .* allbee_processed_tab.li_sort;


%% 




%% sort li score by stage


% latency
load("allbee_processed_tab.mat")
% allbee_processed_tab.latency(allbee_processed_tab.latency>2) = nan;
% allbee_processed_tab.response(isnan(allbee_processed_tab.latency)) = nan;
allbee_processed_tab.resp_num(isnan(allbee_processed_tab.latency)) = nan;
% allbee_processed_tab.li_score = [];
% allbee_processed_tab.li_sort = [];
% allbee_processed_tab.sort_num = [];
LatInhib_agg = grpstats(allbee_processed_tab,["bee_id","stage","context"],"numel",DataVars="latency");
% LI_SCORE = nan(size(unique(allbee_processed_tab.bee_id)))
LI_SCORE = [];
% li_score_sort_stage = "Post_exptest";
li_score_sort_stage = "Post_condtest";

for bee = unique(allbee_processed_tab.bee_id)'
    f_count = LatInhib_agg.numel_latency(LatInhib_agg.bee_id == bee & LatInhib_agg.context=="familiar" & LatInhib_agg.stage==li_score_sort_stage);
    n_count = LatInhib_agg.numel_latency(LatInhib_agg.bee_id == bee & LatInhib_agg.context=="novel" & LatInhib_agg.stage==li_score_sort_stage);
    li_score =  (f_count - n_count)/(f_count + n_count);
    LI_SCORE(end + 1) = li_score;
end

li_tab = table(categorical(unique(allbee_processed_tab.bee_id)),LI_SCORE',VariableNames=["bee_id","li_score"]);
li_tab = sortrows(li_tab,"li_score","descend");
li_tab.li_sort = (1:height(li_tab))';
allbee_processed_tab = join(allbee_processed_tab,li_tab);

allbee_processed_tab.sort_num = allbee_processed_tab.response .* allbee_processed_tab.li_sort;
% allbee_gramm = gramm('x',allbee_processed_tab.trial_num,'y',allbee_processed_tab.resp_num,'color',allbee_processed_tab.context);
%%
close all
fg1 =  figure(1);
set(fg1,'Position',fig_pos_max);

% allbee_gramm = gramm('x',allbee_processed_tab.trial_num,'y',allbee_processed_tab.sort_num,'color',allbee_processed_tab.context,marker=allbee_processed_tab.stim,subset=allbee_processed_tab.stim~="M");
allbee_gramm = gramm('x',allbee_processed_tab.trial_num,'y',allbee_processed_tab.sort_num,'color',allbee_processed_tab.context,marker=allbee_processed_tab.stim);

% allbee_gramm.facet_grid(allbee_processed_tab.stim,[]);
allbee_gramm.geom_point();
allbee_gramm.set_title(sprintf("PER where LI score is sorted by '%s' stage",strrep(li_score_sort_stage,"_"," ")));
% allbee_gramm.set_title(sprintf("PER sorted by no. of responses including US"));

allbee_gramm.set_names('x','Trial num (#)','y',"Bee ID (#)",'row','Stim','color','Context','marker',"stim");
allbee_gramm.set_point_options("base_size",8)
stage_breaks = [30.5,70.5,100.5,120.5];
allbee_gramm.geom_vline('xintercept',[30.5,70.5,100.5,120.5],'style','-k');
% allbee_gramm.geom_vline('xintercept',0.5:1:150.5,'style','-k');
% allbee_gramm.geom_hline('yintercept',0.5:1:max(allbee_processed_tab.bee_num)+ .5,'style','-k');
allbee_gramm.axe_property('XLim',[.5 150.5],YLim=[.5 max(allbee_processed_tab.li_sort)+ .5])
allbee_gramm.axe_property("YTick",1:max(allbee_processed_tab.li_sort),"TickDir","out" ,"Yticklabel",li_tab.bee_id')
allbee_gramm.draw(); %
xline(0.5:1:150.5,Color=[.5 .5 .5],Parent=allbee_gramm.facet_axes_handles(1))
yline(0.5:1:max(allbee_processed_tab.bee_num)+ .5,Color=[.3 .3 .3],Parent=allbee_gramm.facet_axes_handles(1))

stage_midpoints = filter([.5 .5],1,[0 stage_breaks 150])/ 150;
stage_midpoints = stage_midpoints(2:end);

fig = gcf;
% annotation(0,50,{"sd"},'Parent',allbee_gramm.facet_axes_handles(1))
% annotation(fig,"textbox",[0 .8 .1 .02], 'String',"a","EdgeColor",'none')
allstage_names = ["Pre_test","Pre_exp","Post_exptest","Abs_cond","Post_condtest"];
allstage_names = strrep(allstage_names,"_"," ");
pos_info = allbee_gramm.facet_axes_handles.Position;
x0 = pos_info(1);
w = pos_info(3);
tb_x = 0.0844166666666666;
for i  = 1:5
annotation(fig,'textbox',...
    [x0 + w * stage_midpoints(i)- tb_x/2 0.900311526479751 tb_x 0.0197300103842158],...
    'String',{allstage_names(i)},...
    'HorizontalAlignment','center','FitBoxToText','off',EdgeColor="none");
end
% subtitle("a")
set(fig,'Renderer','painters')
% exportgraphics(fig,"plots\overall_plots_seperate_thresholds\" + "allbeedotplots.png")
saveas(gcf,"plots\overall_plots_seperate_thresholds\mindurmethod\overal_per_"+li_score_sort_stage+"us_included.png")
% saveas(gcf,"plots\overall_plots_seperate_thresholds\mindurmethod\overal_per_count_us_included.png")
% agg_allbee =  



%%
fg2 = figure(2);
set(gcf,'Position',fig_pos_minimal);

agg = grpstats(allbee_processed_tab,["stage","stim","context","ss_norm_num"],["mean","numel"],DataVars="resp_num");
gr = gramm('x',agg.ss_norm_num,'y',(agg.numel_resp_num ./ agg.GroupCount) * 100,'color',agg.context,'subset',agg.stage~="Pre_exp" & agg.stage~="Pre_test");

gr.set_order_options('column',["Post_exptest","Abs_cond","Post_condtest"]);

gr.set_names('x','Trial num (#)','y',"% response",'row','Stim','color','Context','column',"Stage");
gr.facet_grid(categorical(agg.stim),categorical(agg.stage));

gr.set_title(sprintf("PER Response"));
% gr.set_text_options(facet_scaling=1)
gr.geom_point();
gr.geom_line();
gr.draw();
% saveas(gcf,"plots\overall_plots_seperate_thresholds\outlier removed\learning_curve_separate.png")

fg3 = figure(3);
set(fg3,'Position',fig_pos_half);

agg2 = grpstats(allbee_processed_tab,["stage","context","ss_norm_num"],["mean","numel"],DataVars="resp_num");

g3 = gramm('x',agg2.ss_norm_num,'y',(agg2.numel_resp_num ./ agg2.GroupCount) * 100,'color',agg2.context,'subset',agg2.stage~="Pre_exp" & agg2.stage~="Pre_test");

g3.set_order_options('column',["Post_exptest","Abs_cond","Post_condtest"]);

g3.set_names('x','Trial num (#)','y'," % response",'color','Context','column',"Stage");
g3.facet_grid([],categorical(agg2.stage));

g3.set_title(sprintf("PER Response"));
% gr.set_text_options(facet_scaling=1)
g3.geom_point();
g3.geom_line();
% gr.set_text_options('font','Helvetica');
g3.draw();
% saveas(gcf,"plots\overall_plots_seperate_thresholds\outlier removed\learning_curve_combined.png")

%% just the precond cond post cond dotplots for per 

figure();
set(gcf,'Position',fig_pos_minimal);
g = gramm('x',allbee_processed_tab.ss_norm_num,'y',allbee_processed_tab.resp_num,'color',allbee_processed_tab.context,'subset',allbee_processed_tab.stage~="Pre_exp" & allbee_processed_tab.stage~="Pre_test");
g.facet_grid(allbee_processed_tab.stim,allbee_processed_tab.stage);
g.set_order_options("column",["Pre_test","Pre_exp","Post_exptest","Abs_cond","Post_condtest"])
g.geom_point();
g.set_title(sprintf("PER across trials and odour context "));
g.set_names('x','Trial num (#)','y',"Bee num (#)",'row','Stim','color','Context','column',"stage");
% g.geom_vline('xintercept',[30.5,70.5,100.5,120.5],'style','--k')
g.draw(); %
% export_path = "D:\ARIZONA BEES for Athil\learning\M17Analysis\M17_analysis_conditioned_inhibition\plots\PER_thresh";
% 
saveas(gcf,"plots\overall_plots_seperate_thresholds\outlier removed\per_condstages.png")

%%
fg4 = figure(4);
set(fg4,'Position',fig_pos_minimal)
g4 = gramm('x',allbee_processed_tab.ss_norm_num,'y',allbee_processed_tab.latency,'color',allbee_processed_tab.context,'subset',allbee_processed_tab.stage~="Pre_exp" & allbee_processed_tab.stage~="Pre_test");
g4.set_order_options('column',["Post_exptest","Abs_cond","Post_condtest"]);

g4.facet_grid(allbee_processed_tab.stim,allbee_processed_tab.stage);
g4.set_names('x','Trial num (#)','y'," latency (s)",'color','Context','column',"Stage",'row',"stim");
% g4.geom_point();
% g4.stat_summary('type','ci','geom','line');
g4.geom_point('alpha',.5);
g4.stat_glm()
g4.draw()
saveas(gcf,"plots\overall_plots_seperate_thresholds\outlier removed\latency_separate.png")
g4.results.stat_glm.model
%% 
fg5 = figure(5);
set(fg5,'Position',fig_pos_half);
g5 = gramm('x',allbee_processed_tab.ss_norm_num,'y',allbee_processed_tab.latency,'color',allbee_processed_tab.context,'subset',allbee_processed_tab.stage~="Pre_exp" & allbee_processed_tab.stage~="Pre_test");
g5.set_order_options('column',["Post_exptest","Abs_cond","Post_condtest"]);

g5.facet_grid([],allbee_processed_tab.stage);
g5.set_names('x','Trial num (#)','y'," latency (s)",'color','Context','column',"Stage");
% g4.geom_point();
% g5.stat_summary('type','ci','geom','area');
g5.stat_glm()
g5.geom_point('alpha',.5);
% g5.stat_summary('geom','errorbar');
g5.draw()
saveas(gcf,"plots\overall_plots_seperate_thresholds\outlier removed\latency_combined.png")
g5.results.stat_glm.model
%%
fg5 = figure(6);
set(fg5,'Position',fig_pos_half);

first_last_dat = allbee_processed_tab(allbee_processed_tab.stage~= "Pre_exp" & allbee_processed_tab.stage~= "Pre_test",:);
first_last_dat = first_last_dat(first_last_dat.ss_norm_num ==1 | first_last_dat.ss_norm_num ==10,:);
first_last_dat = first_last_dat(~(first_last_dat.ss_norm_num == 1 & first_last_dat.stage== "Post_exptest"),:);


g6 = gramm('x', first_last_dat.ss_norm_num,'y',first_last_dat.latency,'color',first_last_dat.context);
g6.set_order_options('column',["Post_exptest","Abs_cond","Post_condtest"]);
g6.set_names('x','Trial num (#)','y'," latency (s)",'color','Context','column',"Stage");
g6.facet_grid([],first_last_dat.stage);
g6(1).axe_property('XTick',10)
g6(1).axe_property('XTick',[1 10])
% g6(1,3).axe_property('XTick',[1 10])
 g6.stat_boxplot()
g6.geom_jitter(dodge=.7)


latcounts = grpstats(first_last_dat,["stage","context","ss_norm_num"],["numel","max"],DataVars="latency")
stage_names = ["Post_exptest","Abs_cond","Post_condtest"];
cntxt_names = ["familiar","mix","novel"]

% g5.geom_point('alpha',.5);
g6.draw()
for i = 1:3
    if i ~=2
        for cntxt = cntxt_names
            
                box_entry = latcounts.stage == stage_names(i) & latcounts.context == cntxt;
                shift_pos = -2:2:2;
                x_pos = latcounts.ss_norm_num(box_entry) +shift_pos(i);
                y_pos = latcounts.max_latency(box_entry) 
                nums = latcounts.numel_latency(box_entry) ;
                text(x_pos,y_pos,"N="+num2str(nums),'Parent',g6.facet_axes_handles(i))
            
                
        end
    else
        for cntxt = ["familiar","novel"]
            
                box_entry = latcounts.stage == stage_names(i) & latcounts.context == cntxt;
                shift_pos = -1:1;
                x_pos = latcounts.ss_norm_num(box_entry)
                y_pos = latcounts.ss_norm_num(box_entry)
                nums = latcounts.numel_latency(box_entry) ;
                text(x_pos+shift_pos(i),y_pos,"N="+num2str(nums),'Parent',g6.facet_axes_handles(i))
            
                
        end
    end
end
latcounts.xpos_shift = latcounts.ss_norm_num - grp2idx(latcounts.stage)
% for i = 1:height(latcounts)
%     x_pos = latcounts.xpos_shift(i);
%     y_pos = latcounts.max_latency(i);
%     facet_info = grp2idx(latcounts.stage);
%    text(x_pos,y_pos,num2str(latcounts.numel_latency(i)),'Parent',g6.facet_axes_handles(facet_info(i)));
% end

