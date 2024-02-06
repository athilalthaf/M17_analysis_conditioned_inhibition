load("Beem17_trim.mat");
load("m17_conditioning_data.mat")
load("m17_post_conditioning_data.mat")
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

fs_ds_cond = round(size(time_cond_range,2) / time_cond ); 
dt_ds = 1/fs_ds_cond;
cond_range = 20;
%%

for fl =  1:length(Bee_m17_cond(1:10)) 
    cond_data = Bee_m17_cond(fl).act(:,1:20);
    cond_data_ds  = downsample(cond_data,downsample_num);

    cond_data_noArt = removeArtifact(cond_data_ds,time_cond_range);
thr_std_pre = stdThSt * getStdPre(cond_data_noArt,time_cond_range,fs_ds_cond);
sp_train = getSpikeTrains(cond_data_noArt,thr_std_pre);

sp_train_fr = causalLowPass(sp_train,dt_ds,fs_ds_cond,tau);

[resp_thresh,onset_thresh, th_thresh] = getResponse(cond_data_noArt,time_cond_range,stdThFr); %%

cond_dat_pw  = getPower(cond_data_noArt,dt_ds,fs_ds_cond,tau);
[resp_pw,onset_pw, th_pw] = getResponse(cond_dat_pw,time_cond_range,stdThPw);

figure(fl)
set(gcf,'Position',get(0,"ScreenSize"));

tiledlayout(10,2,"TileSpacing","tight","Padding","tight")
for idx= 1: cond_range
    nexttile;
    yyaxis left
    plot(time_cond_range,cond_data_ds(:,idx));
    xline(onset_pw(idx),"-k",'LineWidth',1.5);
    xline(onset_thresh(idx),"-m",'LineWidth',1.5)
    yyaxis right
    h = plot(time_cond_range,cond_dat_pw(:,idx));
%     set(h,'Color',[h.Color, 0.3]);
    yline(th_pw(idx),"-k",'LineWidth',1)
    xlim([-1 6]);
%     title("trial no: " + num2str(idx) + "stim: "+  Bee_m17_cond(fl).stim(idx))
    title(num2str(idx)+ "; " + Bee_m17_cond(fl).stim(idx))
    sgtitle(Bee_m17_cond(fl).bee_id +  "; tau =" + num2str(tau) + ";  pw thresh =" + num2str(stdThPw) )
end

legend(["data","power onset", "std onset", " power", "power thresh" ])
saveas(gcf,"plots\threshbee\thresh_45\" + Bee_m17_cond(fl).bee_id +"tau_150.png")
close all
end

%% Single snippet
bee_num = randi(40);
trial_num = randi(20);
% bee_num = 1;
% trial_num = 1;


tau = 100;
 cond_data = Bee_m17_cond(bee_num).act(:,1:20);
    cond_data_ds  = downsample(cond_data,downsample_num);

    cond_data_noArt = removeArtifact(cond_data_ds,time_cond_range);
thr_std_pre = stdThSt * getStdPre(cond_data_noArt,time_cond_range,fs_ds_cond);
sp_train = getSpikeTrains(cond_data_noArt,thr_std_pre);

sp_train_fr = causalLowPass(sp_train,dt_ds,fs_ds_cond,tau);

[resp_thresh,onset_thresh, th_thresh] = getResponse(cond_data_noArt,time_cond_range,stdThFr); %%

cond_dat_pw  = getPower(cond_data_noArt,dt_ds,fs_ds_cond,tau);
[resp_pw,onset_pw, th_pw] = getResponse(cond_dat_pw,time_cond_range,stdThPw);

yyaxis left
plot(time_cond_range,cond_data_ds(:,trial_num),LineWidth=2);
xline(onset_thresh(trial_num),"-m",'LineWidth',1.5,Alpha=.4)
xline(onset_pw(trial_num),"-k",'LineWidth',1.5,Alpha=.4);
yyaxis right
h = plot(time_cond_range,cond_dat_pw(:,trial_num),LineWidth=2);
set(h,'Color',[h.Color, 0.5]);
yline(th_pw(trial_num),"-k",'LineWidth',1)
xlim([-1 6]);
%     title("trial no: " + num2str(idx) + "stim: "+  Bee_m17_cond(fl).stim(idx))
sgtitle(Bee_m17_cond(bee_num).bee_id +  "; tau =" + num2str(tau) + ";  pw thresh =" + num2str(stdThPw) + ...
        " trial no: " + num2str(trial_num) + "stim: "+  Bee_m17_cond(bee_num).stim(trial_num))

% figure(2)
% s =stackedplot(time_cond_range,cond_data_noArt,"HandleVisibility","off");
% ax = findobj(s.NodeChildren,"Type","Axes");
% arrayfun(@(h,x) xline(h,x,'color','g','LineWidth',1.5),ax,onset_pw)
% arrayfun(@(h,x) xline(h,x,'color','m',LineWidth=1.5),ax,onset_thresh)

% cellfun(@(h,x) xline(h,x),ax,onset_thresh)

%% 20 random plots with same tau and different stdPw

stdThPw = [5 15 45];
tau = 50;

% bee_num = randi(40,1,5);
% trial_num = randi(20,1,4);
bee_num =     [3    32    23    31    19];
trial_num = [8    16    11    19];



fig_num = 1;
for bn = bee_num
    for tr = trial_num
    cond_data = Bee_m17_cond(bn).act(:,1:20);
    cond_data_ds  = downsample(cond_data,downsample_num);

    cond_data_noArt = removeArtifact(cond_data_ds,time_cond_range);
thr_std_pre = stdThSt * getStdPre(cond_data_noArt,time_cond_range,fs_ds_cond);
sp_train = getSpikeTrains(cond_data_noArt,thr_std_pre);

sp_train_fr = causalLowPass(sp_train,dt_ds,fs_ds_cond,tau);

[resp_thresh,onset_thresh, th_thresh] = getResponse(cond_data_noArt,time_cond_range,stdThFr); %%

cond_dat_pw  = getPower(cond_data_noArt,dt_ds,fs_ds_cond,tau);
ONSET_PW = [];
TH_PW = [];
for thresh = stdThPw
[resp_pw,onset_pw, th_pw] = getResponse(cond_dat_pw,time_cond_range,thresh);
ONSET_PW = [ONSET_PW, onset_pw(tr)];
TH_PW = [TH_PW, th_pw(tr)];

end

figure(fig_num)
yyaxis left
plot(time_cond_range,cond_data_ds(:,tr),LineWidth=2,DisplayName="m17 data");
xline(onset_thresh(tr),"-m",'LineWidth',1.5,Alpha=.4,DisplayName="std onset")
if sum(isnan(ONSET_PW)) ~= 3

    xline(ONSET_PW(~isnan(ONSET_PW)),"-k",string(stdThPw(~isnan(ONSET_PW))),'LineWidth',1.5,Alpha=.4,DisplayName="pw onset");
end
% scatter(ONSET_PW,[0 0 0]) 
 yyaxis right
 h = plot(time_cond_range,cond_dat_pw(:,tr),LineWidth=2,DisplayName="filtered m17");
% set(h,'Color',[h.Color, 0.5]);
%  yline(th_pw(tr),"-k",'LineWidth',1)
xlim([-1 6]);
%     title("trial no: " + num2str(idx) + "stim: "+  Bee_m17_cond(fl).stim(idx))
sgtitle(Bee_m17_cond(bn).bee_id +  "; tau =" + num2str(tau) + ";" + ...
        " trial no: " + num2str(tr) + " stim: "+  Bee_m17_cond(bn).stim(tr))
    fig_num = fig_num + 1;
    legend(Location="best");
    saveas(gcf,"plots\diff_thresh_t50\" + Bee_m17_cond(bn).bee_id+"_tr_"+num2str(tr)+"_th5_10_15_t"+num2str(tau) +".png")
    end
    
end


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

bees_id = unique(clean_cond_tb.bee_id);
cond_trial_start = 100;
% 
% for id = 2%length(bees_id)
%     figure(id)
%     set(gcf, 'Position', get(0, 'Screensize'));
% 
%     sp = stackedplot(time_cond_range,clean_cond_act_nart(:,clean_cond_tb.bee_id==bees_id(id)), ...
%         DisplayLabels=string(clean_cond_tb.trial_num(clean_cond_tb.bee_id==bees_id(id)) - cond_trial_start))
%     ax = findobj(sp.NodeChildren,"Type","Axes")
%     arrayfun(@(h,x) xline(h,x,'color','r','LineWidth',1.5),ax,clean_cond_tb.latency(clean_cond_tb.bee_id ==bees_id(id)))
%     saveas(gcf,"plots\thresh_5_t150_cond\" + bees_id(id) +"_t5_tau_100.png")
%     title(bees_id(id))
% 
%     close all
% end

% tiledlayout(10,2,"TileSpacing","tight","Padding","tight")
thresh_color = lines(7);
thresh_color = thresh_color(7,:);
for id = 1:length(bees_id)
    figure(id);

    single_bee_condnum = length(clean_cond_tb.bee_id(clean_cond_tb.bee_id == bees_id(id)));
    tiledlayout(single_bee_condnum,1,'TileSpacing','none','Padding','tight');
    cond_trial_start = min(clean_cond_tb.trial_num(clean_cond_tb.bee_id == bees_id(id))) - 1;
    for trl = 1:single_bee_condnum
        nexttile;
        single_trial = clean_cond_tb.bee_id==bees_id(id) & clean_cond_tb.trial_num - cond_trial_start==trl;
        ylabel(string(clean_cond_tb.trial_num(single_trial)))
        plot(time_cond_range,clean_cond_act_nart(:,single_trial));
        xline(clean_cond_tb.latency(single_trial),Color=thresh_color,LineWidth=1.5,Alpha=0.9)
    end
            sgtitle(bees_id(id))
%     saveas(gcf,"plots\thresh_5_t150_cond\" + bees_id(id) +"_t5_tau_100.png")
    close all
end



%% post conditioning

clean_postcond_act_nart = removeArtifact(clean_postcond_act,time_cond_range);
postcond_dat_pw = getPower(clean_postcond_act_nart,dt_ds,fs_ds_cond,tau);
[resp_post_pw,onset_post_pw,th_post_pw] = getResponse(postcond_dat_pw,time_cond_range,stdThPw);
% clean_cond_tb.latency = onset_pw;
context = clean_postcond_tb.stim == clean_postcond_tb.PreExp;
context_arr = strings(size(context));
context_arr(context) = "familiar";
context_arr(~context) = "novel";
context_arr(clean_postcond_tb.stim == "M") = "mix";
clean_postcond_tb.latency = onset_post_pw;
clean_postcond_tb.context = context_arr;
% clean_postcond_tb = clean_postcond_tb(~isnan(clean_postcond_tb.latency),:);
h = boxplot(clean_postcond_tb.latency,clean_postcond_tb.context);
boxes = findobj(h,"Tag","Box");
cols = lines(5);
cols= cols(end-2:end,:)
for i=1:length(boxes)

    patch(get(boxes(i),"XData"),get(boxes(i),"YData"),cols(i,:),"FaceAlpha",0.85);
end
title("Latency during conditioning")
ylabel("time (s)") 

%% All data
addpath('C:\Users\athil\OneDrive - uni-bielefeld.de\Desktop\Codes\gramm_matlab\gramm')
whole_activity = cell2mat(m17_sep_tab.act');
whole_activity_nart = removeArtifact(whole_activity,time_cond_range);
whole_activity_pw = getPower(whole_activity_nart,dt_ds,fs_ds_cond,tau);
[response,latency,~] = getResponse(whole_activity_pw,time_cond_range,stdThPw);
m17_sep_tab.latency = latency;
m17_sep_tab.response = response;

m17_sep_tab.latency(m17_sep_tab.stage == 'Abs_cond' & m17_sep_tab.latency >=2) = nan;

m17_sep_tab.response(m17_sep_tab.stage == 'Abs_cond' & m17_sep_tab.latency >=2) = 0;

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
clean_m17 = m17_sep_tab(~isnan(m17_sep_tab.latency),:);

cond_tab  = m17_sep_tab(m17_sep_tab.stage == "Abs_cond",:);
postcond_tab  = m17_sep_tab(m17_sep_tab.stage == "Post_condtest",:);


%%
figure(1)
g = gramm('x',m17_sep_tab.trial_num,'y',m17_sep_tab.latency,'color',m17_sep_tab.stage);
g.facet_grid(m17_sep_tab.stim,[]);
g.geom_point();
g.set_title("Latency across trials and phases");
g.set_names('x','Trial num (#)','y',"Time (s)",'row','Stim','color','Stages');

g.draw(); % 


cond_agg = grpstats(cond_tab,["context","trial_num"],["mean","numel"],DataVars="latency");
cond_agg = cond_agg(cond_agg.context~='mix',:);

postcond_agg = grpstats(postcond_tab,["context","trial_num"],["mean","numel"],DataVars="latency");
% postcond_agg = cond_agg(postcond_agg.context~='mix',:);

figure(2)
gr = gramm('x',cond_agg.trial_num,'y',cond_agg.numel_latency,'color',cond_agg.context);
gr.set_title('PER response during conditioning');
gr.set_names('x','Trial num (#)','y',"Count (#)",'color','Odour');
gr.geom_point();
gr.geom_line();
gr.draw()


figure(3)
gp = gramm('x',postcond_agg.trial_num,'y',postcond_agg.numel_latency ,'color',postcond_agg.context);
gp.set_title('PER response after conditioning');
gp.set_names('x','Trial num (#)','y',"Count (#)",'color','Odour');
gp.geom_point();
gp.geom_line();
gp.draw()
