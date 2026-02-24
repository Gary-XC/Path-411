function lowExpression = MarkLowCounts(inputMatrix, quantileLevel)
% MarkLowCounts.m
%
% This function identifies lowly expressed features from the quantile
% threshold. For each feature we check if it has an expression above 
% quantileLevel and apply inverse boolean to find true values of low
% expression from high expression in genes.
%
% Input:
%   inputMatrix - matrix of type double
%   quantileLevel - double for where 0 <= quantileLevel <= 1
%
% Output:
%   lowExpression - boolean vector where T has no samples above threshold,
%   and F has at least one sample above threshold

% Checks for double type and quantileLevel in range (assign .5 default)
arguments
    inputMatrix double {mustBeNonempty}
    quantileLevel double {mustBeInRange(quantileLevel, 0, 1)} = 0.5
end

% Extra matrix check
if ~ismatrix(inputMatrix)
    error("Input must be a 2D matrix.");
end

% Values that we want from our matrix (for size and shape)
rowValue = 1;
colValue = 1;

% Grab the number of rows in the data from 1st dimension (rows)
featureNum = size(inputMatrix, rowValue);

% Create a column vector of #rows x 1 col filled with false boolean
lowExpression = false(featureNum, colValue);

% Standard dimension values of a 2D matrix
rowDimension = 1;
colDimension = 2;

% Reshape the original matrix so it is flattened to a 1D array as the
% quantile functions requires a vector and not a matrix
matrixVector = reshape(inputMatrix, [rowDimension, size(inputMatrix, rowDimension)*size(inputMatrix, colDimension)]);

% Calculate threshold number based on our matrixVector for low expression
% cutoff
threshold = quantile(matrixVector, quantileLevel);

% For each feature, check if any sample > threshold across COLUMNS, and return True for
% rows with at least one sample > threshold (returns column vector with 1
% result per row)
highExpression = any(inputMatrix > threshold, colDimension);

% Return lowExpression genes as True if NOT(FALSE) from highExpression
lowExpression = ~highExpression;

end




