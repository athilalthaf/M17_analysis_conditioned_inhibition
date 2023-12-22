% script to collect all the m17 data and annotate with the stim ids
% Written by Athil Althaf 



% cedpath = 'C:\Users\athil\CEDMATLAB\CEDS64ML';
% addpath(cedpath);
% CEDS64LoadLib(cedpath);
% rel_m17_path = fileparts("../rawdataAllBeesPreAnalysis/");

m17_dir = dir("..\ThiagoNeuralynx\m17_exported\"); % directory of the raw m17 mat file
stim_dir = dir("..\..\analysisData\"); % directory containing the stim id info

stim_path_info = stim_dir(1).folder; % path of stim id files 
stim_files = string({stim_dir.name}); % list containing all the stim file name
stim_files = stim_files(contains(stim_files,"Stim") & ~contains(stim_files,"Mix")); % only taking the files named stim and excluding the mix
 
fs = 30303; % sampling frequency
dt = 1/fs ; % timestep
time_non_cond = 5; % in  seconds [1 3 1] [pre stim post ] timings
time_cond = 7 ; % in  seconds [1 5 1]


m17_path_info = m17_dir(1).folder ; % path info for the m17 mat files
m17_files = string({m17_dir.name}); % getting the list of mat file names 
m17_files = m17_files(contains(m17_files,".mat")); % isolating mat files only



bee_id  = []; % initialise bee id list for the stim table
stim_id  = []; % initialise stimtype list for the stim table
timestamp = []; %initialise timestamp list for the stim table
for st_fl = stim_files % for each stimulus files
    
    beeid_stimid = split(st_fl,"_"); % splitting the file string by _ ; it has both bee id and stim type
    
    stim_data = load(stim_path_info + "\" + st_fl); % load the stim file 
    bee_id = [bee_id; repmat(beeid_stimid(1),size(stim_data,1),1)]; % creating the list with beeid as the size of the file  and appending to the big list 
    stim_type = char(beeid_stimid(2)); % getting the stim type
    stim_type = string(stim_type(5));  % isolating the stimulus char
    stim_id = [stim_id;repmat(stim_type,size(stim_data,1),1) ];  % creating the list with stim type as the size of the file and appending to the big list 
    timestamp = [timestamp; stim_data];  % appending the timestamp to the bigger list
end
stim_info_tab = table(bee_id,stim_id,timestamp); % making a table with the id,stim, actual timestamp
stim_info_tab = sortrows(stim_info_tab,[1,3]); % sorting the list based on bee id and then with the timestamp order

% 
Bee_m17 = struct(); % initialising the struct for incorporating the m17 data
 
for fl = 1:length(m17_files) % for each m17 file   %% check whole list is called or not
    file_name = m17_files(fl); % get the file name
    raw_m17 = load(m17_path_info + "\" + file_name); % load the m17 mat file
    m17_field = fieldnames(raw_m17); % get the fieldname of the nested struct
    raw_m17_act = raw_m17.(m17_field{1}).values;  %access the m17 activities within that struct
    raw_m17_times = raw_m17.(m17_field{1}).times; %access the m17 times within that struct
    
    bee_name = char(file_name); % isolating the beeid
    bee_name = bee_name(1:end -7); 
    Bee_m17(fl).bee_id = string(bee_name) ; % attributing the bee ids
    Bee_m17(fl).stim  = stim_info_tab.stim_id(stim_info_tab.bee_id == bee_name,:)'  ; % getting the stim info the the previous stim table created
    Bee_m17(fl).timestamp  = stim_info_tab.timestamp(stim_info_tab.bee_id == bee_name,:)' ; % gettomg the timestamp info
    time_st = Bee_m17(fl).timestamp; % getting the timestamp list for indexing 
    % note that the activity has the length of longest recordings that is
    % the conditioning phase 
    Bee_m17(fl).act = NaN(time_cond*fs,length(time_st)); % pre initialising the activity by nan
    

    for ts = 1:length(time_st) % for each timestamp
        [~,stim_flag] = min(abs(time_st(ts) - raw_m17_times)); % finding the nearest index of the timestamp
        
        if (ts > 100 ) && (ts <= 120) % trials during conditioning 
            % use the flag to get 1 second before stim and 1 second after
            % stim , here the total stim time is taken as 7s
            Bee_m17(fl).act(1:fs*time_cond,ts) = raw_m17_act(stim_flag-(1*fs):stim_flag + (time_cond - 1)*fs - 1);
        else  % all the other trials
            % here the total stim time is taken as 5s
            Bee_m17(fl).act(1:fs*time_non_cond,ts) = raw_m17_act(stim_flag-(1*fs):stim_flag + (time_non_cond - 1)*fs - 1);
        end

    end
          
    sprintf(bee_name + " done , percentage : " + num2str(fl/length(m17_files) * 100)) % status of the progress
end


% [~,rel_stim_path] = fileparts("..\..\..\");
% 

% stim_dir = fullfile(pwd,rel_stim_path)



% unloadlibrary ceds64int

 
stdThSt = 15; 
stdThFr = 5;
stdThPw = 45;

tau = 100; %time constant for convolve function
% fignum = 1;

time_cond_range = linspace(-1,time_cond-1,fs* time_cond); % time range for conidtioning data
time_non_cond_range = linspace(-1,time_non_cond-1,fs* time_non_cond); %time range of all other cases
%%%%%%%%%%%
% cond_data = Bee_m17(1).act(:,101:120);
% post_cond_data = Bee_m17(1).act(:,121:150);
%%
downsample_num = 20; 
time_cond_range = downsample(time_cond_range,downsample_num); % downsampling the time range
 
fs_ds_cond = round(fs / downsample_num );  % new sampling frequency for downsampled data 
dt_ds = 1/fs_ds_cond; % corresponding time period
% cond_range = 20; 
%%

% cond_data_noArt = removeArtifact(cond_data,time_cond_range);
% thr_std_pre = stdThSt * getStdPre(cond_data_noArt,time_cond_range,fs);
% sp_train = getSpikeTrains(cond_data_noArt,thr_std_pre);
% 
% sp_train_fr = causalLowPass(sp_train,dt_ds,tau);
% 
% [resp_thresh, onset_thresh, th_thresh] = getResponse(cond_data_noArt,time_cond_range,stdThFr); %%
% 
% cond_dat_pw  = getPower(cond_data_noArt,dt_ds,tau);
% [resp_pw,onset_pw, th_pw] = getResponse(cond_dat_pw,time_cond_range,stdThPw);
% figure(fignum)
% 
% tiledlayout(10,2,"TileSpacing","tight","Padding","tight")
% for idx= 1: 20
%     nexttile;
%     yyaxis left
%     plot(time_cond_range,cond_data(:,idx));
%     xline(onset_pw(idx),"-k",'LineWidth',1.5);
%     xline(onset_thresh(idx),"-m",'LineWidth',1.5)
%     yyaxis right
%     h = plot(time_cond_range,cond_dat_pw(:,idx));
%     set(h,'Color',[h.Color, 0.3]);
%     yline(th_pw(idx),"-k",'LineWidth',1)
%     xlim([-1 6]);
%     sgtitle("tau =" + num2str(tau) + ",  pw thresh =" + num2str(stdThPw) )
% end
% 
% legend(["data","power onset", "std onset", " power", "power thresh" ])
% % fignum = fignum + 1;
%% 
Bee_m17_tab = struct2table(Bee_m17); 
trial_len = cellfun(@length, Bee_m17_tab.stim); % get the total trial numbers of each experiment
Bee_m17_tab.trial_num = trial_len; % get a new column with trial length
 
Bee_m17_tab = Bee_m17_tab(Bee_m17_tab.trial_num>=150,:); % exclude every thing that has not made till the end
cond_num = 20;
time_diff_cutoff = 200; % time threshold to get the different phases
long_interval = 30*60; % 
POST_COND_START_IDX = zeros(height(Bee_m17_tab),1); % initialise post condition start index array
COND_START_IDX = zeros(height(Bee_m17_tab),1); % initialise condition phase  starting index array

COND_ITEMS= strings(height(Bee_m17_tab),3); % initialse array for conditioning stims 
COND_ITEM_COUNTS = zeros(height(Bee_m17_tab),3); % respective counts
COND_TRIAL_NUM = zeros(height(Bee_m17_tab),1); % check the number of elements in conditioning phase
for bee = 1:height(Bee_m17_tab) % 
    bee_timestamps = cell2mat(Bee_m17_tab.timestamp(bee)); %convert every timestamp to matrix
    time_stamp_diff = [0 diff(bee_timestamps)]; % get the first order difference 
    time_stamp_stage_idx = find(time_stamp_diff >time_diff_cutoff); % find all the difference that goes above a threshold gives us the starting indices of each respective stages
    COND_TRIAL_NUM(bee) = time_stamp_stage_idx(4) - time_stamp_stage_idx(3); % number of conditioning trials
    COND_START_IDX(bee) = time_stamp_stage_idx(3); % starting index for conditioning phase
%     POST_COND_START_IDX(bee) = find(time_stamp_diff == time_stamp_diff(time_stamp_diff> long_interval)); % 
    POST_COND_START_IDX(bee) = time_stamp_stage_idx(4);  % post conditioning phase starting index
    bee_stim_types = Bee_m17_tab.stim(bee); % get the stim types
    bee_stim_types = bee_stim_types{1}; % access the cell contents 
    bee_stim_cond = bee_stim_types(POST_COND_START_IDX(bee) - cond_num: POST_COND_START_IDX(bee) - 1); % not really relevant
    [counts, items] = groupcounts(bee_stim_cond'); % not relevant
    COND_ITEMS(bee,1:length(items)) = items';  
    COND_ITEM_COUNTS(bee,1:length(counts)) = counts'; 
    
%     plot(time_stamp_diff);     %%%%%% to the timestamp difference if
%     needed
    
%     hold on

end
% yline(200);
% hold off

Bee_m17_tab.post_cond_start_idx = POST_COND_START_IDX; % append a column for post condition start index 
Bee_m17_tab.cond_start_idx = COND_START_IDX;  % append a column for conditioning start index
Bee_m17_tab.cond_num = COND_TRIAL_NUM;  % append a column for 
 
Bee_m17_tab.cond_items  = COND_ITEMS; % not really relevant but can see the the conditioning phase diagnosis
Bee_m17_tab.cond_count  = COND_ITEM_COUNTS;


%% seperate the timestamps and stim type for the ease of access also adding the pre exposure info


bee_data_time_stamps = Bee_m17_tab(:,["bee_id","stim","timestamp","trial_num","cond_items","cond_count"]); % get the table without activity for easier manipulation 

% total_entry_length = sum(bee_data_time_stamps.trial_num); 
time_stamp_arr = []; % initialising a set of arrays for the restructured table not elegant way but it works
stim_type_arr = []; 
bee_id_arr = [];
stim_num_arr = [];
trial_num_arr = [];
act_arr = [];

for entry = 1:height(bee_data_time_stamps)
    bee_id_arr = [bee_id_arr ; repmat(bee_data_time_stamps.bee_id(entry), bee_data_time_stamps.trial_num(entry),1)]; %repeating the beeid to match length
    stim_cell = bee_data_time_stamps.stim(entry); % get the stim types
    stim_arr_sub = stim_cell{:}; % unpack the cell contents
    stimnum_arr_sub = zeros(size(stim_arr_sub));  % trial number of each stim array initialising
    trial_num_arr_sub = 1:size(stim_arr_sub,2); % trial number of overall experiment of each bee
    for typ = ["A","B","M"] % for each stim types
        stim_arr_sub_id = stim_arr_sub == typ; % get the elements where the stim type matches
        stimnum_arr_sub(stim_arr_sub_id) = 1:sum(stim_arr_sub_id); %  fill the corresponding indices with the range 1:maximum num
    end
    timestamp_cell = bee_data_time_stamps.timestamp(entry); % get the timestamp cell
    stim_type_arr = [stim_type_arr ; stim_arr_sub']; % appending the arrays
    time_stamp_arr = [time_stamp_arr ; timestamp_cell{:}'];
    stim_num_arr = [stim_num_arr ; stimnum_arr_sub'];
    trial_num_arr = [trial_num_arr; trial_num_arr_sub'];
    act_sub = Bee_m17_tab.act(entry); % get the activity cell
%     act_sub = act_sub{:};
    act_sub_ds = downsample(act_sub{:},downsample_num); % downsample it 
    act_sub = num2cell(act_sub_ds,1); % convert back to cell
    act_arr = [act_arr; act_sub']; % append the activity array

end

m17_sep_tab = table(bee_id_arr,stim_type_arr,trial_num_arr,stim_num_arr,time_stamp_arr,act_arr, ...
    VariableNames=["bee_id","stim","trial_num","stim_num","time_stamp","act"]); % get a new table with all the new column in the new format
% save("m17_sep_tab.mat","m17_sep_tab","-v7.3")

bee_ids = unique(m17_sep_tab.bee_id); % get the unique bee ids
pre_exp_arr = []; % initialise the preexposure array 
for bee = bee_ids' % for all the ids
  single_bee_tab  = m17_sep_tab(m17_sep_tab.bee_id == bee,:); % get the single bee talbe
  [counts, items] = groupcounts(single_bee_tab.stim); % count and factor the stim types
%   items == ["A","B","M"]' %checking the order is correct or not  
   if counts(1) - counts(2) > 0  % if number of A > B
       PreExp = "A"; % assign the preExp as A
   else 
       PreExp = "B"; % else otherwise
   end
   pre_exp_arr = [pre_exp_arr ; repmat(PreExp, height(single_bee_tab),1)]; % append the preExp array
%   pre_exp_tab = single_bee_tab(single_bee_tab.trial_num >30 & single_bee_tab.trial_num <=70,:);
%   unique(pre_exp_tab.stim);
end
m17_sep_tab.PreExp = pre_exp_arr; % add the info to the table 

%% get the 


% Bee_m17_cond = struct();
% 
% post_trial_end = 150;
% cond_start_trial = 101;
% 
% for cond_fl = 1:length(Bee_m17)
%     if length(Bee_m17(cond_fl).stim) >= post_trial_end
%         Bee_m17_cond(cond_fl).bee_id = Bee_m17(cond_fl).bee_id;
%         
%         stim_trim = Bee_m17(cond_fl).stim;
%         stim_trim = stim_trim(cond_start_trial:post_trial_end);
%         Bee_m17_cond(cond_fl).stim = stim_trim; 
% 
%         timest_trim = Bee_m17(cond_fl).timestamp;
%         timest_trim = timest_trim(cond_start_trial:post_trial_end);
%         Bee_m17_cond(cond_fl).timestamp = timest_trim;
%         
%         act_trim = Bee_m17(cond_fl).act;
%         act_trim = act_trim(:,cond_start_trial:post_trial_end);
%         Bee_m17_cond(cond_fl).act = act_trim;
%     else
%         continue;
%     end
% 
% end
% 
% non_empty_rows = ~cellfun(@isempty,{a.bee_id});
% Bee_m17_cond = Bee_m17_cond(non_empty_rows);
% 

% figure(2)
% tiledlayout(10,2)
% for idx= 1: 20
%     nexttile;
%     plot(time_cond_range,cond_dat_pw(:,idx));
% %     xline(onset_pw(idx),"-g",'LineWidth',2);
%     yline(th_pw(idx),"-g",'LineWidth',2,Alpha=.5)
%     xlim([-1 6]);
% end
m17_conditioning_data = []; % conditioning data array
m17_post_conditioning_data = []; % post conditioning data
for i = 1:height(Bee_m17_tab) % for each bee
   bee_tab = m17_sep_tab(m17_sep_tab.bee_id == Bee_m17_tab.bee_id(i),:); % get single bee tab
   single_cond_bee_tab = bee_tab(bee_tab.trial_num >= Bee_m17_tab.cond_start_idx(i) & ...
       bee_tab.trial_num < Bee_m17_tab.post_cond_start_idx(i),:); % get the conditioning trials in that bee
    m17_conditioning_data = [m17_conditioning_data ; single_cond_bee_tab];  % append to conditioning data

    single_post = bee_tab(bee_tab.trial_num>=Bee_m17_tab.post_cond_start_idx(i),:); % get all the post condiitoning trials in that bee
    m17_post_conditioning_data = [m17_post_conditioning_data;single_post]; % append the post conditioning 
end 
save("m17_post_conditioning_data.mat","m17_post_conditioning_data","-v7.3");  % save the data for conditioning and post conditioning
save("m17_conditioning_data.mat","m17_conditioning_data","-v7.3");
