function dataout = reformatICGCdata(donor_ids, sample_ids, feature_ids, data)
%
% dataout = reformatICGCdata(donor_ids, sample_ids, feature_ids, data)
%returns the ICGC data reformated into a cellarray. The ICGC 
%data comes in the format where the data for samples is stored in a vector
%format, all features for one sample are followed by features for another
%sample.  
%The function rearranges the data into a cellarray  where the first
%row stores donor_ids, second row stores samples ids and the rest of the
%rows store the features. This way, the data for each sample is stored in
%its own column.
%
%Input parameters:
% donor_ids = a vector of donor ids. Pass the entire column 'icgc_donor_id'
% as is
% sample_ids = a vector of sample ids. Pass the entire column 'icgc_sample_id'
% as is
% feature_ids = a vector of feature ids. Pass the entire column of
% features as is, these could be: mirna_id or ind_ensid, etc.
% data - a numeric matrix that stores the data that you wish to use. For
% eamples the 'normalized_read_count' column or the 'raw_read_count'
% column.
%
%output parameters:
%dataout - the reformated cellarray
%
%example:
% celldata_raw = reformatICGCdata(mirnaseq.icgc_donor_id, mirnaseq.icgc_sample_id, mirnaseq.mirna_id, mirnaseq.raw_read_count);
% celldata_norm = reformatICGCdata(mirnaseq.icgc_donor_id, mirnaseq.icgc_sample_id, mirnaseq.mirna_id, mirnaseq.normalized_read_count);
%
%author: Kathrin Tyryshkin
%date: January 2018

%initialize
dataout = {}; %empty cell array
 
%check the input
if length(donor_ids) ~= length(sample_ids) || ...
   length(donor_ids) ~= length(feature_ids)|| ...
   length(donor_ids) ~= length(data)
    disp('Input vectors must be of the same length');
    return
elseif ~isnumeric(data)
    disp('The input data vector must be double (numeric)');
    return
end
%convert strings to character arrays if necessary
if isstring(donor_ids)
    donor_ids = cellstr(donor_ids);
end
if isstring(feature_ids)
    feature_ids = cellstr(feature_ids);
end
if isstring(sample_ids)
    sample_ids = cellstr(sample_ids');
end


uniqueIDs = unique(donor_ids); %get unique donor ids
uniqueFeatures = unique(feature_ids); %get unique feature ids

%initialize the output cell array - first column is the feature names, the
%rest of the columns are samples. all values are initialized to 0.
dataout = ['donour IDs' uniqueIDs']; 
dataout(3:length(uniqueFeatures)+2, 1) = uniqueFeatures;
dataout(2,1) = {'sample IDs'};
dataout(3:end, 2:end) = num2cell(zeros(length(uniqueFeatures),length(uniqueIDs)));


%assemble the data for all samples (unique ids)
for i=1:length(uniqueIDs)
    curr_donor_id = uniqueIDs(i);
    %find all indices to the current donor id
    indA = find(strcmp(donor_ids, curr_donor_id));
    curr_sample_ids = sample_ids(indA);
    curr_feature_ids = feature_ids(indA);
    curr_data = data(indA);
    
    %match the feature ids in the output cellarray and current features
    %use 'stable' to keep the order of the features (otherwise they will be
    %sorted)
    [~, indA, indB] = intersect(dataout(:, 1), curr_feature_ids,'stable');
    %fill into the array
    dataout(2, i+1) = curr_sample_ids(1);
    dataout(indA, i+1) = num2cell(curr_data(indB));    
end

