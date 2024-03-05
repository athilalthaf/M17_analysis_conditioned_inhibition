% plotting every stage and their respective thresholds 
load("allbee_processed_tab.mat");
thresh_list = readtable("thresh_list.xlsx");
%%
STAGES = ["Pre_test","Pre_exp","Post_exptest","Abs_cond","Post_condtest"];
rownum = 10;
tau = 100;
plot_path_allstages = "plots\allbee_allstages\with_hampel_filtering";
for bee_idx = 1%:height(thresh_list)
    for stg_idx = 1:length(STAGES)
        stage_and_id = allbee_processed_tab.bee_id == thresh_list.bee_id(bee_idx) & allbee_processed_tab.stage == STAGES(stg_idx);
        sub_tab = allbee_processed_tab(stage_and_id,:);
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
        [sub_tab_processed,pw] = add_latency_and_response(sub_tab,th,tau,true);
        
        for sp_id = 1:subplot_num
            nexttile;
            sp = plot(time_cond_range,cell2mat(sub_tab_processed.act(sp_id))');
            xline(sub_tab_processed.latency(sp_id),'-r','LineWidth',2);
            ylabel(sprintf("stim %s, trial num %d",sub_tab_processed.stim(sp_id),sub_tab_processed.ss_norm_num(sp_id)), ...
                HorizontalAlignment='right',Rotation=0);
            ylim([-2 * 10^(-3) 2*10^(-3)]);
            yyaxis right
            plot(time_cond_range,pw(:,sp_id));
%             set(h,'Color',[])
            yline(sub_tab_processed.cutoff(sp_id),'-k');
            
            if stg_idx == 4
                xlim([-1 6]);
            else
                xlim([-1 4]);
            end



        end

        stagename  = replace(STAGES(stg_idx),"_"," ");
        sgtitle(sprintf("id = %s, stage = %s, thresh = %d, tau = %d",thresh_list.bee_id(bee_idx),stagename,th,tau))
        xlabel("time (s)")
        saveas(fig,plot_path_allstages + sprintf("%s_stage_%d.png",thresh_list.bee_id(bee_idx),stg_idx));
    end

end
