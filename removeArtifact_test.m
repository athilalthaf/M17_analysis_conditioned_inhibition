function data = removeArtifact(data,time,time_range)
     data  = data'; %### transpose to get the indexing correct
    artTime = time>=time_range(1) & time<time_range(2);
    data(:,artTime) = repmat(mean(data,2),1,sum(artTime));
    data = data' ;  %## transponse to get back the original dimensions
    
end