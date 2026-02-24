function avgCorrelations = SampleCorrelation(inputMatrix, correlationType)
% SampleCorrelation.m
%
% Computes average correlation for each sample/column in the input
% matrix, to all other samples using either Spearman or Pearson.
%
% Input:
%   inputMatrix - matrix of type double
%   correlationType - string for 'Spearman' or 'Pearson' correlation
%
% Output:
%   avgCorrelations - 1D row vector of average correlations

% Must be double type and correlation belonging to "Spearman" or "Pearson"
arguments
    inputMatrix double {mustBeNonempty}
    correlationType string {mustBeMember(correlationType, ["Spearman", "Pearson"])} = "Spearman"
end

if ~ismatrix(inputMatrix)
    error("Input must be a 2D matrix.");
end

% Find size of 2nd dimension (columns)
numSamples = size(inputMatrix, 2);

% Initialize default NaN for missing values for no confusion with 0
% correlation as our returned row vector (1D)
avgCorrelations = NaN(1, numSamples);

% Safety check if we do not have samples
if numSamples < 2
    error("At least 2 samples are needed to compute correlations.");
end

% Compute correlations between all columns, sum each row, exclude
% self-correlation (rowSums-1) and divide by numSamples - 1 for average
corrMatrix = corr(inputMatrix, 'Type', char(correlationType));

% Transposed to match output for [1 52]
rowSums = sum(corrMatrix, 2);
avgCorrelations = ((rowSums - 1) / (numSamples - 1))';

end