%% subplots 

[resp,onset,th] = getResponse(b2_pw,time_cond_range,5);
IDX = 1:20;
fig = figure(1);
set(fig,'Position',fig_pos);
for i = IDX
    plot(time_cond_range,resp_act(:,i));
    xline(onset(th),'-r');
    

    yyaxis left
    plot(time_cond_range,b2_pw(:,i))
    yline(th(i),'-r')
    
end
%%

% resp_on = allbee_processed_tab(allbee_processed_tab.response == 1,:);
% filt = filt(:,~isnan(allbee_processed_tab.response));
% act_hampel = act_hampel(:,~isnan(allbee_processed_tab.response));

%%
for n = 41:50
figure(n)
idx = (n-1) * 10 + 1: n * 10;
tiledlayout(10,1,"TileSpacing","none");
for i = idx
    nexttile;
    plot(time_cond_range,filt(:,i)',DisplayName='act',Color=[.7 .7 .7]);
    ylabel(sprintf("stim %s, trial num %d",allbee_processed_tab.stim(i),allbee_processed_tab.ss_norm_num(i)), ...
                HorizontalAlignment='right',Rotation=0);
    hold on;
    if isnan(allbee_processed_tab.response(i))
    plot(time_cond_range,act_hampel(:,i),Color='k',DisplayName='median+outlier');
    else
        plot(time_cond_range,act_hampel(:,i),Color='g',DisplayName='median+outlier');
    end
    xline(allbee_processed_tab.latency(i),'-r','LineWidth',2);
    
    hold off
    box off
end
legend()
end
%%

resp_non_cond = ~isnan(allbee_processed_tab.response) & allbee_processed_tab.stage ~= "Abs_cond";

non_cond_resp_tab = allbee_processed_tab(resp_non_cond,:);

act = cell2mat(allbee_processed_tab.act');
filt = getLowpassfiltered(allbee_processed_tab,600,fs_ds_cond);

[act_hampel,outmat] = getMedianandOutlier(filt,round(fs_ds_cond/2),5);

outmat(isnan(filt)) = 0;

integ_out = getWindowIntegrate(outmat,500,fs_ds_cond);

integ_out(isnan(filt)) = nan;

integ_trim = integ_out(time_cond_range>=0 & time_cond_range<3,:); % reducing the time window to 0 - 3;
max_outliers = 10;
cond_idx = allbee_processed_tab.stage == "Abs_cond";
integ_out_cond  = integ_out(:,cond_idx);
integ_trim_cond = integ_out_cond(time_cond_range>=0 & time_cond_range<2,:);

allbee_processed_tab.response = double((max(integ_trim,[],1)>max_outliers))';
allbee_processed_tab.response(cond_idx) = double((max(integ_trim_cond,[],1)>max_outliers))';


allbee_processed_tab.response(allbee_processed_tab.response == 0) = nan;
allbee_processed_tab.resp_num = allbee_processed_tab.bee_num .* allbee_processed_tab.response;

% gscatter(allbee_processed_tab.trial_num,allbee_processed_tab.resp_num,allbee_processed_tab.context)

allbee_integ_response = allbee_processed_tab(~isnan(allbee_processed_tab.response),:);

[allbee_latencytab,pw] = add_latency_and_response(allbee_integ_response,5,100,false);

allbee_processed_tab(~isnan(allbee_processed_tab.response),:) = allbee_latencytab;

allbee_processed_tab.latency(isnan(allbee_processed_tab.response)) = nan;










