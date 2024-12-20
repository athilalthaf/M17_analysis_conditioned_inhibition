load smple_data.mat

smple_per = smple_data(~isnan(smple_data.latency),:)

fs = 1000;
dt = 1/fs;
th = 5;

%% check a particular trial onset works or not;

% act = cell2mat(allbee_processed_tab.act');
% log_id = getthatrialid(allbee_processed_tab,"bee102","Post_exptest","B",10);
log_id = 101
pw = cell2mat(allbee_processed_tab.power(log_id))';
act = cell2mat(allbee_processed_tab.act(log_id))';
th = allbee_processed_tab.cutoff(log_id);

min_dur = .4;
potential_onsets = find(pw>th);
thresh_crossing = pw>th;

if sum(thresh_crossing) ~= 0
thresh_boundaries = [0,diff(pw>th)];
group_start = time_cond_range(thresh_boundaries>0);
group_end = time_cond_range(thresh_boundaries<0);
duration = group_end - group_start(1:length(group_end));
longer_onsets = group_start(duration > min_dur);
else
    group_start = [];
    group_end = [];
    duration = group_end(1:length(group_start)) - group_start;
end
% if ~isempty(potential_onsets)
%    onset_diff = diff(potential_onsets);
%    group_start = [potential_onsets(1), potential_onsets(find(onset_diff > 1)+1)];
%    group_end = [potential_onsets(find(onset_diff > 1)),potential_onsets(end)];
%    
%    valid_onset = [];
%    
%    start_time = time_cond_range(group_start);
%    end_time = time_cond_range(group_end);
%    duration = end_time - start_time;
%    longer_onsets = start_time(duration > min_dur);
  if ~isempty(longer_onsets)
   first_onset = longer_onsets(1);
  end
%    for i = 1%:length(group_start)
%        start_time = time_cond_range(group_start(i));
%        end_time = time_cond_range(group_end(i));
%        
%        duration = (end_time - start_time);
% 
%        if duration >=min_dur
%            valid_onset(end+1) = start_time;
% 
%        end

%    end


% end

plot(time_cond_range,act)

yyaxis right
plot(time_cond_range,pw);

% title()
yline(th)
if ~isempty(longer_onsets) && sum(thresh_crossing) ~= 0
% xline(start_time,"g")
xline(first_onset,"r")

end
%% testing interactive function

fig =figure();
ax= axes(fig);
y_data =sin(-pi:0.001:pi) ;
plot(ax,y_data)
hold(ax,"on");


slider_wf = uicontrol(fig,"Style",'slider',"Min",-2,"Max",2,"Value",0,'Position',[75 3 440 20], ...
    "Sliderstep",[.01 .01],"Callback",@slidercallback);


function slidercallback(hobject,eventdata,handles,y_data)

    persistent lineHandle
    if isempty(lineHandle) || ~isvalid(lineHandle)
        lineHandle = yline(gca,0);
    end

    sliderval = get(hobject,'Value');
    lineHandle.Value = sliderval;
    if lineHandle.Value < y_data
        onset = find(lineHandle.Value < y_data);
        xline(onset(1));
       
    end
end
% 
% h = images.roi.Line(ax,'Position',[xlim(ax);[0 0]]', ...
%     'Label',"cutoff", ...
%     'InteractionsAllowed',"translate")
% addlistener(h,'MovingROI',@linemovingfcn)


% function linemovingfcn(obj,event)
% ax = ancestor(obj,'axes');
% ax.Title.String = sprintf("y = %.3f",event.CurrentPosition(1));
% end
% function updateLines(fig,signals,x)
% for i = 1:
% 
% end


%%
function data = removeArtifact(data,time)
    artTime = time>=0 & time<0.01;
    data(:,artTime) = repmat(mean(data,2),1,sum(artTime));
end


function stdPre = getStdPre(data,time,fs)
    preTime = time <-1; 
    pool = data(:,preTime);
    totLen = size(data,1)*size(pool,2);
    pool = reshape(pool',1,totLen);
    
    spikes = diff([0 pool] > 5*std(pool))==1;
    
    kernel = ones(1,0.01 * fs)+1; % kernel of 10ms around the spike 
    aroundSpikes = conv(spikes,kernel);
    aroundSpikes = aroundSpikes(length(kernel)/2:end-length(kernel)/2);
    filtPool = pool(~aroundSpikes);
    
    stdPre = std(filtPool);
end

function SpikeTrains = getSpikeTrains(data,th)
    %SpikeTrains = diff(repmat(sign(mean(data,2)),1,size(data,2)).*data > th,1,2) == 1;
    SpikeTrains = diff(data > th,1,2) == 1;
    SpikeTrains = [SpikeTrains zeros(size(SpikeTrains,1),1)];
end


function power = getPower(data,dt,tau)
   kwidth = 500; % kernel width in ms
   t = 0:dt*1000:kwidth;
   k = exp(-t / tau);
   k = k ./ sum(k);
   power=filter(k,1,data.^2')';  
end