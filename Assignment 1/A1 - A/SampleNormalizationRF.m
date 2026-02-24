% Gary Chen - 20359001

function normData = SampleNormalizationRF(data)
    % SampleNormalizationRF: Performs relative frequency (RF) normalization
    % Input:  data - numeric matrix (double)
    % Output: normData - RF-normalized matrix (each column sums to 1)

    % making sure that the data is the proper type, checking first that it
    % is numeric, and then checking the type is double
    % otherwise it throws an error and stops the function
    if ~isnumeric(data) || isa(data, 'double')
        error('SampleNormalizationRF: Invalid Type');
    end

    % checking if the data is empty and if it is, then the function stops
    % and it throws an error and stops the function
    if isempty(data)
        error('SampleNormalizationRF: EmptyInput');
    end

    % Calculating column sums
    colSums = sum(data, 1);

    % Avoiding divide-by-zero by setting values that equal to 0 to 1 
    colSums(colSums == 0) = 1;

    % Normalizing each column
    normData = data ./ colSums;
end


