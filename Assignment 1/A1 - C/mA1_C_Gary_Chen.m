% Gary Chen - 20359001

% Loading the data
load("ForMiniAssignments.mat");

% Cell array case

cell_datain = cell2mat(celldata_ICGCdata_norm(3:end, 2:end)); % Rows 3:end skip the first two rows, Columns 2:end skip the first column
% Convert a subset of the cell array celldata_ICGCdata_norm into a numeric matrix

cellResult = MarkLowCounts(cell_datain, 0.9);

% Table case

tbl_datain = tbl_ICGCdata_norm{:, 2:end};
% Columns 2:end skip the first column, and use the curly brackets for only
% the raw data and not a subtable

tableResult = MarkLowCounts(cell2mat(tbl_datain), 0.9);
% input the cleaned and processed data to the SampleCorrelation
% function, after normalizing the data into an array of doubles

% Pretests
try 
    datain = (celldata_ICGCdata_norm(3:end,2:end));

result = MarkLowCounts(datain, -1);
    
catch e
   if strcmp(e.identifier,"MATLAB:quantile:BadProbs")
       error("If this message is seen, then the function is not filtering out the incorrect values");
   end
end


try 
    datain = (celldata_ICGCdata_norm(3:end,2:end));

result = MarkLowCounts(datain, 1.4);
    
catch e
   if strcmp(e.identifier,"MATLAB:quantile:BadProbs")
       error("If this message is seen, then the function is not filtering out the incorrect values");
   end
end

try 
    cell_datain = (celldata_ICGCdata_norm(3:end,2:end));

    result = MarkLowCounts(cell_datain,'a');
    
catch e
   if strcmp(e.identifier,"MATLAB:quantile:BadProbs")
       error("If this message is seen, then the function is not filtering out the incorrect values");
   end
end

try 
    datain = (celldata_ICGCdata_norm(3:end,2:end));

result = MarkLowCounts(datain, 0.9);
    
catch e
   if strcmp(e.identifier,"MATLAB:prctile:InvalidData")
       error("If this message is seen, then the function is not filtering out the incorrect values");
   end
end

try 
    datain = (celldata_ICGCdata_norm(3:end,2:end));

result = MarkLowCounts(datain);
    
catch e
   if strcmp(e.identifier,"MATLAB:minrhs")
       error("If this message is seen, then the function is not filtering out the incorrect values");
   end
end