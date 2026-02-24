% Gary Chen - 20359001

% Loading the data
load("ForMiniAssignments.mat");

% Cell array case

cell_datain = cell2mat(celldata_ICGCdata_norm(3:end, 2:end)); % Rows 3:end skip the first two rows, Columns 2:end skip the first column
% Convert a subset of the cell array celldata_ICGCdata_norm into a numeric matrix

cellSpear = SampleCorrelation(cell_datain, "Spearman");
cellPearson = SampleCorrelation(cell_datain, "Pearson");

% Table case

tbl_datain = tbl_ICGCdata_norm{:, 2:end};
% Columns 2:end skip the first column, and use the curly brackets for only
% the raw data and not a subtable

cellRes = SampleCorrelation(cell2mat(tbl_datain), "Spearman");
tablePearson = SampleCorrelation(cell2mat(tbl_datain), "Pearson");
% input the cleaned and processed data to the SampleCorrelation
% function, after normalizing the data into an array of doubles

% Test Cases
% Testing for non numeric data
try 
    datain = (celldata_ICGCdata_norm(3:end,2:end));

result = SampleCorrelation(datain, "Pearson");
    
catch e
   if strcmp(e.identifier,"MATLAB:UndefinedFunction")
       error("Non Numeric Data Inputed")
   end
end

% Testing if a string was inputted
try 
    cell_datain = (celldata_ICGCdata_norm(3:end,2:end));

    result = SampleCorrelation(cell_datain,'');
    
catch e
   if strcmp(e.identifier,"stats:corr:UnknownType")
       error("Strings were passed to the function")
   end
end

% Testing for empty inputs
try 
    datain = {};

result = SampleCorrelation(datain, "Pearson");
    
catch e
   if strcmp(e.identifier,"stats:corr:TooFewInputs")
       error("There were no values passed to the function")
   end
end

% input missing
try 
    datain = cell2mat(celldata_ICGCdata_norm(3:end,2:end));

result = SampleCorrelation(datain);
    
catch e
   if strcmp(e.identifier,"MATLAB:minrhs")
       error("There was nothing passed to the function")
   end
end