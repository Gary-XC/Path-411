% SampleNormalizationRF.m 
% 
% Normalizes a matrix by computing the sum
% of all values in a column and dividing each value in the column of the
% matrix by its sum.
%
% Input: inputMatrix - A matrix of type double
% Output: normalizedMatrix - A matrix normalized from column sums

function normalizedMatrix = SampleNormalizationRF(inputMatrix)

% Check input is of type double and not empty
arguments
    inputMatrix double {mustBeNonempty}
end

% Check if input is a matrix if not throw error
if ~ismatrix(inputMatrix)
    error("The input must be a 2D matrix.");
end

% Compute the sum for each column in the matrix
columnSums = sum(inputMatrix);

% Divide each element in the matrix by the column sum
normalizedMatrix = inputMatrix ./ columnSums;

end


