param_tab = readtable("power_method_calib.xlsx",Sheet="Sheet3");
tab_stats = grpstats(param_tab,["tau","pw_thresh"],{'numel'});
%% overall onset detection
detect_tab = param_tab(:,["tau","pw_thresh","onset_correct"]);

detect_tab = grpstats(detect_tab,["tau","pw_thresh"],"numel");
detect_tab.detect = detect_tab.numel_onset_correct./ detect_tab.GroupCount * 100 ;
% figure(1)
% heatmap(detect_tab,"tau","pw_thresh",ColorVariable="detect",Title="percentage detection");
%% no signal ones omitted
signal_param = param_tab(param_tab.description~="no signal",:);
detect_tab = signal_param(:,["tau","pw_thresh","onset_correct"]);

detect_tab = grpstats(detect_tab,["tau","pw_thresh"],"numel");
detect_tab.detect = detect_tab.numel_onset_correct./ detect_tab.GroupCount * 100 ;
figure(2)
heatmap(detect_tab,"tau","pw_thresh",ColorVariable="detect",Title="percentage detection");

%% correctly identified onsets onset labeled ones
detect_tab = signal_param(:,["tau","pw_thresh","onset_correct"]);
detect_tab = detect_tab(~isnan(detect_tab.onset_correct),:);
breakdown_tab = groupcounts(detect_tab,["tau","pw_thresh","onset_correct"]);
breakdown_tab_condition = groupcounts(detect_tab,["tau","pw_thresh"]);

true_detect = breakdown_tab(breakdown_tab.onset_correct==1,:);
true_detect.detect = true_detect.GroupCount ./ breakdown_tab_condition.GroupCount * 100;

false_detect = breakdown_tab(breakdown_tab.onset_correct==0,:);
false_detect.detect = false_detect.GroupCount ./ breakdown_tab_condition.GroupCount * 100;

% undetect = breakdown_tab(isnan(breakdown_tab.onset_correct),:);
% undetect.detect = undetect.GroupCount ./ breakdown_tab_condition.GroupCount * 100;
% 

figure(3)
heatmap(true_detect,"tau","pw_thresh",ColorVariable="detect",Title=" Correctly detected onsets");

figure(4)
heatmap(false_detect,"tau","pw_thresh",ColorVariable="detect",Title=" Incorrectly detected onsets");

% figure(5)
% heatmap(undetect,"tau","pw_thresh",ColorVariable="detect",Title=" undetected onsets");
%% accuracy 
param_tab_red = param_tab(:,["tau","pw_thresh","onset_correct","description"]);
doubt_entries = param_tab_red.description == "doubt";
param_tab_red.description(doubt_entries) = {''};


param_tab_red.description  = string(param_tab_red.description);
param_condensed = groupcounts(param_tab_red,["tau","pw_thresh","onset_correct","description"])

efficiency_tab = groupcounts(param_tab_red,["tau","pw_thresh"]);
accuracy= [];
precision = [];
recall = [];
for i = 1:height(efficiency_tab)
    t = efficiency_tab.tau(i);
    pw = efficiency_tab.pw_thresh(i);
    
    tab_sub = param_condensed(param_condensed.tau == t & param_condensed.pw_thresh == pw,:)
    
    true_pos = tab_sub.GroupCount(tab_sub.onset_correct == 1);
    true_neg = tab_sub.GroupCount(isnan(tab_sub.onset_correct) & tab_sub.description == "no signal");
    false_pos = sum(tab_sub.GroupCount(tab_sub.onset_correct ==0));
    false_neg = tab_sub.GroupCount(isnan(tab_sub.onset_correct) & tab_sub.description == "");
    acc = (true_pos + true_neg) / efficiency_tab.GroupCount(i);
    prec = true_pos/(true_pos + false_pos);
    
    rec = true_pos/ (true_pos + false_neg);
    accuracy = [accuracy; acc];
    precision = [precision; prec];
    recall = [recall; rec ];
end

efficiency_tab.accuracy = accuracy;
efficiency_tab.precision = precision;
efficiency_tab.recall = recall;

figure(4)
heatmap(efficiency_tab,"tau","pw_thresh",ColorVariable="accuracy",Title=" Accuracy");

figure(5)
heatmap(efficiency_tab,"tau","pw_thresh",ColorVariable="precision",Title="Precision");

figure(6)
heatmap(efficiency_tab,"tau","pw_thresh",ColorVariable="recall",Title="Recall");

