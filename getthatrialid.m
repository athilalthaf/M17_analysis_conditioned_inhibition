function id = getthatrialid(tab,beeid,stage,stim,ss_normnum) 
% function to return the data index from beeid , stage, stim and stim stage
% normalised trial number
% dat = allbee_processed_tab.act(allbee_processed_tab.bee_id == beeid & allbee_processed_tab.stage == stage & allbee_processed_tab.stim == stim &allbee_processed_tab.ss_norm_num ==ss_normnum);
log_id = tab.bee_id == beeid & tab.stage == stage & tab.stim == stim & tab.ss_norm_num ==ss_normnum;

id  = find(log_id);
end