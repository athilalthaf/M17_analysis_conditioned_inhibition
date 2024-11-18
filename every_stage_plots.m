% plotting every stage and their respective thresholds 
load("allbee_processed_tab.mat");
thresh_list = readtable("thresh_list.xlsx");
nsigma = 3;

% plot_path  = "plots\PER_thresh\post_troubleshoot\"; 
%%
STAGES = ["Pre_test","Pre_exp","Post_exptest","Abs_cond","Post_condtest"];      %stage name
rownum = 10;
tau = 100;
% plot_path_allstages = "plots\allbee_allstages\with_hampel_filtering\";
plot_path_allstages = "plots\allbee_allstages\mindurmethod\us_included\";
min_dur  = .4;
max_latency = 2;
allbee_processed_tab_edt = get_latency_min_dur(allbee_processed_tab,min_dur,time_cond_range,max_latency);

max_latency_cond = 4;
allbee_processed_tab_cond = get_latency_min_dur(allbee_processed_tab_edt(allbee_processed_tab_edt.stage == "Abs_cond",:),min_dur,time_cond_range,max_latency_cond);
allbee_processed_tab_edt(allbee_processed_tab.stage == "Abs_cond",:) = allbee_processed_tab_cond;
%% 
for bee_idx = 3:height(thresh_list)
    for stg_idx = 1:length(STAGES)
        stage_and_id = allbee_processed_tab_edt.bee_id == thresh_list.bee_id(bee_idx) & allbee_processed_tab_edt.stage == STAGES(stg_idx);
        sub_tab = allbee_processed_tab_edt(stage_and_id,:);
%         act_hmp = 
        subplot_num = height(sub_tab);
        fig = figure();
        set(fig,'position',fig_pos);
        if stg_idx == 2             % prexposure phase
            col_num = 4;
        elseif stg_idx == 4         % cond phase
            col_num = 2;
        else
            col_num = 3;
        end
        

        tiledlayout(rownum,col_num,'TileSpacing','compact');
        
        sub_tab = sortrows(sub_tab,{'ss_norm_num','stim'});
        current_bee_num = sub_tab.bee_num(1); 
        current_bee_exp = sub_tab.PreExp(1);
        th = thresh_list.th(bee_idx);
%         [sub_tab_processed,pw] = add_latency_and_response(sub_tab,th,tau,true);
%         [act_hampel,out_idx_mat] = getMedianandOutlier(cell2mat(sub_tab_processed.act'),500,5);
        for sp_id = 1:subplot_num
            nexttile;
%             sp = plot(time_cond_range,cell2mat(sub_tab_processed.act(sp_id))');
            sp = plot(time_cond_range,cell2mat(sub_tab.act(sp_id))');% hold on
%             act_hmp = getMedianandOutlier(cell2mat(sub_tab.act(sp_id))',500,nsigma);
%             plot(time_cond_range, act_hmp,Color='k'); hold off;
            
            xline(sub_tab.latency(sp_id),'-r','LineWidth',2);
            ylabel(sprintf("stim %s, trial num %d",sub_tab.stim(sp_id),sub_tab.ss_norm_num(sp_id)), ...
                HorizontalAlignment='right',Rotation=0);
%             ylim([-2 * 10^(-3) 2*10^(-3)]);
            yyaxis right
            plot(time_cond_range,cell2mat(sub_tab.power(sp_id))');
%             set(h,'Color',[])
            yline(sub_tab.cutoff(sp_id),'-k');
%           
%             if stg_idx == 4
%                 xlim([-1 4]);
%             else
%                 xlim([-1 4]);
%             end

             xlim([-1 4])

        end

        stagename  = replace(STAGES(stg_idx),"_"," ");
        sgtitle(sprintf("id = %s, stage = %s, thresh = %d, tau = %d ",string(cell2mat(thresh_list.bee_id(bee_idx))),stagename,th,tau))
        xlabel("time (s)")
    
        saveas(fig,plot_path_allstages + sprintf("%s_stage_%d.png",string(cell2mat(thresh_list.bee_id(bee_idx))),stg_idx));
        close all

    end

end
%%  single_stage single bee plot



