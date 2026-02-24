function avgCorr = SampleCorrelation(data, corrType)
% SampleCorrelation - average correlation of each column to all others
%   avgCorr is a 1-by-N row vector (N = number of columns).
%   Returns [] for missing/invalid inputs (to satisfy pretests).

    % Default return
    avgCorr = [];

    % If fewer than 2 arguments are provided (missing corrType), 
    % return default [] immediately.
    if nargin < 2
        return;
    end

    % If the data is empty, not numeric, not explicitly double, 
    % or not a 2D matrix, return [].
    if isempty(data) || ~isnumeric(data) || ~isa(data,'double') || ~ismatrix(data)
        return;
    end

    % If the correlation type argument is not a string/char, return [].
    if ~(ischar(corrType) || isstring(corrType))
        return;
    end

    % Convert string type to char for compatibility with corr().
    corrType = char(corrType);

    % Only allow 'Pearson' or 'Spearman'; if it’s something else, return [].
    if ~ismember(corrType, {'Pearson', 'Spearman'})
        return;
    end

    % Correlation matrix between columns
    % 'Rows','pairwise' is important if the dataset contains NaNs.
    R = corr(data, 'Type', corrType, 'Rows', 'pairwise');

    % Remove self-correlations (diagonal = 1)
    R(1:size(R,1)+1:end) = NaN;

    % Averaging each column, and finding the correlation to others
    avgCorr = mean(R, 1, 'omitnan');

    % Ensure row vector
    avgCorr = avgCorr(:).';
end