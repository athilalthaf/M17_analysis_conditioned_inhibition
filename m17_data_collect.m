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
 
for fl = 1:length(m17_files(1)) % for each m17 file   %% check whole list is called or not
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
          

end


% [~,rel_stim_path] = fileparts("..\..\..\");
% 

% stim_dir = fullfile(pwd,rel_stim_path)



% unloadlibrary ceds64int