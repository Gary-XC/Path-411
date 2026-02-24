% Gary Chen - 20359001

function lowCounts = MarkLowCounts(data, Q)
% MarkLowCounts - identify lowly expressed features
%
%   Inputs:
%       data - numeric matrix of type double (rows = features, cols = samples)
%       Q    - quantile level (0 ≤ Q ≤ 1) used to set expression threshold
%
%   Output:
%       lowCounts - Boolean column vector (rows = features)
%                   True  = feature is lowly expressed (no sample above threshold)
%                   False = feature has at least one sample above threshold

    % Input validation
    if nargin < 2
        error('MarkLowCounts:MissingInput', ...
              'Function requires both data matrix and quantile Q.');
    end

    if isempty(data) || ~isnumeric(data) || ~isa(data,'double') || ~ismatrix(data)
        error('MarkLowCounts:InvalidData', ...
              'First input must be a non-empty double matrix.');
    end

    if ~isscalar(Q) || ~isnumeric(Q) || Q < 0 || Q > 1
        error('MarkLowCounts:InvalidQuantile', ...
              'Second input Q must be a numeric scalar between 0 and 1.');
    end

    % Flattening matrix into vector for better expression distribution
    threshold = quantile(data(:), Q);

    % For each row (feature), check if any value > threshold
    % If none exceed the threshold → mark as lowly expressed (true)
    aboveThresh = any(data > threshold, 2);

    % Negate so that:
    %   True  = lowly expressed
    %   False = expressed in at least one sample
    lowCounts = ~aboveThresh;

end