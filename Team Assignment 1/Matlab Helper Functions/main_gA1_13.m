% PATH 411: Group Assignment 1
% main_gA1_13.m
% Author(s): Sebastian, Lenard, Gary
%
% The program processes the MALY-DE (ICGC) miRNA-seq for DLBCL only
% end to end for quality control. We import donor + mirna_seq, subset DLBCL by ICD-10 with C83.3 
% diagnosis code, reformat the matrices to wide format, do an initiali visualization
% for expression distribution box plots, IQR, and total counts for raw and normalized data,
% log2-transform the data for a more normal distribution, visualize
% totals/IQR/boxplots, compare provied normalized data to RF normalized data, 
% compute the average Spearman correlation, detect
% outliers where IQR = 0 or where samples are 0.2 or more below the alpha
% threshold for correlation analysis, and we repeat the QC until the
% dataset is stable and no outliers are present. Finally, we apply a
% threshold filtering for highly-expressed genes.
%
% In general we:
% 1. Load data into the correct format for preprocessing.
% 2. Visualize initial data for raw/norm counts with scatterplots/boxplots for total samples,
% IQR, expression distribution.
% 3. Normalize and log 2 transform our data and check normal distribution.
% 4. Compare normalization techniques from norm counts vs RF normalization.
% 5. Perform QC, checking for IQR, total reads, Spearman correlation, and
% expression distribution while removing any outliers (repeat).
% 6. Filter out lowly-expressed miRNA genes under a quantile threshold.
% 7. Final visualization of filtered data with IQR, expression
% distribution, total reads, Spearman correlation, and a
% clustergam/heatmap.

% Identifying Diffuse Large B-Cell Lymphoma (DLBCL) with the ICD-10 code: C83.3
dlbcl_index = donor.donor_diagnosis_icd10 == 'C83.3';
dlbcl_tbl = donor(dlbcl_index, 'icgc_donor_id');

% Filter miRNA-seq data to keep just DLBCL samples only
mirnaseq_dlbcl = innerjoin(mirna_seq, dlbcl_tbl, 'Keys', 'icgc_donor_id');

% Reformat data from long to wide format with miRNAs as rows, samples as
% colums

celldata_raw = reformatICGCdata(mirnaseq_dlbcl.icgc_donor_id, mirnaseq_dlbcl.icgc_sample_id, mirnaseq_dlbcl.mirna_id, mirnaseq_dlbcl.raw_read_count);
celldata_norm = reformatICGCdata(mirnaseq_dlbcl.icgc_donor_id, mirnaseq_dlbcl.icgc_sample_id, mirnaseq_dlbcl.mirna_id, mirnaseq_dlbcl.normalized_read_count);

% Grab donors, samples, data as rows from cell arrays and feature names and
% sample data as columns

% Sample IDs and miRNA IDs
sample_ids = celldata_raw(2, 2:end);
mirna_ids = celldata_raw(3:end, 1);

% Expression matrices for raw and normlized counts
expression_raw = cell2mat(celldata_raw(3:end, 2:end));
expression_norm = cell2mat(celldata_norm(3:end, 2:end));

% Check for int data type for raw
is_integer_data = all(mod(expression_raw(:), 1) == 0);

% Check small values for normalized data
max(expression_norm(:))

% Check for matching dimensions
size(expression_raw)
size(expression_norm)

% Initial visualization for raw data QC
total_counts_raw = sum(expression_raw, 1);
iqr_raw = iqr(expression_raw, 1);

% Replace any zeroes with low values to perform the log transforms for
% skewed data
expression_raw_replace_zeros = replaceZeros(expression_raw, 'lowval');
expression_raw_log = log2(expression_raw_replace_zeros);

n_samples_raw = length(sample_ids);
alpha_raw = computeAlphaOutliers(n_samples_raw);

% Plot Raw Data QC

% Testing plot for raw total counts (not used in report)
figure('Position', [100 100 1000 500]);
plot(1:length(total_counts_raw), total_counts_raw, 'bo-', ...
     'MarkerSize', 8, 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
xlabel('Sample Index', 'FontSize', 12);
ylabel('Total Read Counts', 'FontSize', 12);
title('miRNA-Seq DLBCL Raw Data: Total Counts per Sample', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Testing plot for IQR per sample where IQR = 0 is no variation (not used
% in report)
figure('Position', [100 100 1000 500]);
plot(1:length(iqr_raw), iqr_raw, 'ro-', ...
     'MarkerSize', 8, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
xlabel('Sample Index', 'FontSize', 12);
ylabel('IQR', 'FontSize', 12);
title('miRNA-Seq DLBCL Raw Data: IQR per Sample', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Boxplot for raw expression distribution 
figure('Position', [100 100 1200 500]);
boxplot(expression_raw, 'PlotStyle', 'traditional', ...
        'Colors', 'k', 'Whisker', 1.5);
ylabel('Linear Raw Counts', 'FontSize', 12);
xlabel('Samples', 'FontSize', 12);
title('miRNA-Seq DLBCL Raw Data: Expression Distribution', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Boxplot for log transformed raw expression distribution
figure('Position', [100 100 1200 500]);
boxplot(expression_raw_log, 'PlotStyle', 'traditional', ...
        'Colors', 'k', 'Whisker', 1.5);
xlabel('Samples');
ylabel('Log2(Normalized Expression)');
title('Provided miRNA-Seq DLBCL Raw Data: Expression Distribution'); 
grid on;

% Using scatterplot functions to visualize total raw counts and raw IQR
[outliers_total_raw, outlier_ids_total_raw, ~, ~] = ...
    scatterplotMarkOutliers(total_counts_raw, ...
                           'ColumnLabels', sample_ids, ...
                           'PlotTitle', 'miRNA-Seq DLBCL Raw Data: Total Counts per Sample', ...
                           'xlabel', 'Sample Index', ...
                           'ylabel', 'Total Read Counts', ...
                           'ShowXTickLabel', false);

% Raw data IQR
[outliers_iqr_raw, outlier_ids_iqr_raw, ~, ~] = ...
    scatterplotMarkOutliers(iqr_raw, ...
                           'ColumnLabels', sample_ids, ...
                           'PlotTitle', 'miRNA-Seq DLBCL Raw Data: IQR per Sample', ...
                           'xlabel', 'Sample Index', ...
                           'ylabel', 'IQR', ...
                           'ShowXTickLabel', false);

% Plot Normalized Data QC
total_norm = sum(expression_norm, 1);
iqr_norm = iqr(expression_norm, 1);

% Log transforms to handle skewed data distribution
expression_norm_replace_zeros = replaceZeros(expression_norm, 'lowval');
expression_norm_log = log2(expression_norm_replace_zeros);

% Scatterplot for total norm counts
[outliers_total_norm, outlier_ids_total_norm, ~, ~] = ...
    scatterplotMarkOutliers(total_norm, ...
        'ColumnLabels', sample_ids, ...
        'PlotTitle', 'Provided miRNA-Seq DLBCL Norm Counts: Total per Sample', ...
        'xlabel', 'Sample Index', ...
        'ylabel', 'Total', ...
        'ShowXTickLabel', false);

% Scatterplot for IQR norm counts
[outliers_iqr_norm, outlier_ids_iqr_norm, ~, ~] = ...
    scatterplotMarkOutliers(iqr_norm, ...
        'ColumnLabels', sample_ids, ...
        'PlotTitle', 'Provided miRNA-Seq DLBCL Norm Counts: IQR per Sample', ...
        'xlabel', 'Sample Index', ...
        'ylabel', 'IQR', ...
        'ShowXTickLabel', false);

% Boxplot for norm expression distribution (not log transformed)
figure('Position', [100 100 1200 500]);
boxplot(expression_norm, 'PlotStyle', 'traditional', ...
        'Colors', 'k', 'Whisker', 1.5);
xlabel('Samples');
ylabel('Linear Normalized Expression');
title('Provided miRNA-Seq DLBCL Norm Counts: Expression Distribution'); 
grid on;

% Log transformed expression distribution boxplot
figure('Position', [100 100 1200 500]);
boxplot(expression_norm_log, 'PlotStyle', 'traditional', ...
        'Colors', 'k', 'Whisker', 1.5);
xlabel('Samples');
ylabel('Log2(Normalized Expression)');
title('Provided miRNA-Seq DLBCL Norm Counts: Expression Distribution'); 
grid on;

% Peform RF normalization on the raw samples with log2 transformation and
% choose best normalization
expression_rf = SampleNormalizationRF(expression_raw);

% Check if normalization is summed to 1.0
column_sums = sum(expression_rf, 1);

% Log transformation
expression_rf_replace_zeros = replaceZeros(expression_rf, 'lowval');
expression_rf_log = log2(expression_rf_replace_zeros);

% Comparison of before and after data skew with log transformation for more
% normal distribution

% Raw histograms
figure('Position', [100 100 1200 500]);
subplot(1, 2, 1);
histogram(expression_raw(:), 50, 'FaceColor', 'b', 'EdgeColor', 'k');
xlabel('Raw Counts (linear scale)', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('Before Log: Highly Skewed', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

subplot(1, 2, 2);
histogram(expression_rf_log(:), 50, 'FaceColor', 'g', 'EdgeColor', 'k');
xlabel('Log2(RF Normalized)', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('After Log: More Normal', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Norm histograms
figure('Position', [100 100 1200 500]);
subplot(1, 2, 1);
histogram(expression_norm(:), 50, 'FaceColor', 'b', 'EdgeColor', 'k');
xlabel('Norm Counts (linear scale)', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('Before Log: Highly Skewed', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

subplot(1, 2, 2);
histogram(expression_norm_log(:), 50, 'FaceColor', 'g', 'EdgeColor', 'k');
xlabel('Log2(RF Normalized)', 'FontSize', 12);
ylabel('Frequency', 'FontSize', 12);
title('After Log: More Normal', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Comparing both provided normalized counts and RF normalized counts to
% check better alignment
figure('Name', 'Normalization Comparison', 'Position', [100 100 1200 500]);

% RF Normalization boxplot
subplot(1, 2, 1);
boxplot(expression_rf_log, 'PlotStyle', 'traditional', ...
        'Colors', 'k', 'Whisker', 1.5);
ylabel('Log2(RF Normalized Expression)', 'FontSize', 12);
xlabel('Samples', 'FontSize', 12);
title('RF Normalization', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
% Store y-limits for comparison scale
ylim_rf = ylim;  

% Provided Normalization boxplot
subplot(1, 2, 2);
boxplot(expression_norm_log, 'PlotStyle', 'traditional', ...
        'Colors', 'k', 'Whisker', 1.5);
ylabel('Log2(Normalized Expression)', 'FontSize', 12);
xlabel('Samples', 'FontSize', 12);
title('Original Normalization', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
ylim_prov = ylim;

% Match y-axis limits for fair comparison so scales are consistent
ylim_common = [min(ylim_rf(1), ylim_prov(1)), max(ylim_rf(2), ylim_prov(2))];
subplot(1, 2, 1); ylim(ylim_common);
subplot(1, 2, 2); ylim(ylim_common);

% expression_rf_log was deemed the better normalization method
normalization_data = expression_rf_log;

% Perform quality control with RF normalization with IQR and Spearman
% Correlation
sample_iqr = iqr(normalization_data);
average_correlation = SampleCorrelation(normalization_data, "Spearman");

% Calculating sample size and alpha threshold check
n_samples = size(normalization_data, 2);
alpha = computeAlphaOutliers(n_samples);

% Look for outliers in IQR
[outliers_iqr_qc1, outlier_ids_iqr_qc1, low_outliers_iqr1, high_outliers_iqr1] = ...
    scatterplotMarkOutliers(sample_iqr, ...
                           'ColumnLabels', sample_ids, ...
                           'PlotTitle', 'RF Normalized Quality Control: Sample IQR', ...
                           'xlabel', 'Sample Index', ...
                           'ylabel', 'IQR', ...
                           'ShowXTickLabel', false);

% Look for outliers from correlation
[outliers_corr1, outlier_ids_corr1, low_outliers_corr1, high_outliers_corr1] = ...
    scatterplotMarkOutliers(average_correlation, ...
                           'ColumnLabels', sample_ids, ...
                           'PlotTitle', 'RF Normalized Quality Control: Average Spearman Correlation', ...
                           'xlabel', 'Sample Index', ...
                           'ylabel', 'Average Correlation', ...
                           'ShowXTickLabel', false);

% Plotting expression distribution of RF normalized data to see alignment
figure('Position', [100 100 1200 500]);
boxplot(normalization_data, 'PlotStyle', 'traditional', ...
        'Colors', 'k', 'Whisker', 1.5);
xlabel('Samples');
ylabel('Log2(RF Normalized Expression)');
title('miRNA-Seq DLBCL RF Norm Counts: Expression Distribution'); 
grid on;

% Removing samples of technical failure where IQR = 0
donor_ids_header = celldata_raw(1, 2:end);
sample_ids_header = sample_ids;

% Samples where IQR = 0 indicate technical failure and should be removed
% A boolean mask is used and we remove sample IDs/donor IDs corresponding
% to the True boolean
failed_samples = (sample_iqr == 0);
remove_sample_ids = sample_ids_header(failed_samples);
remove_donor_ids = donor_ids_header(failed_samples);

% Check if samples are not in the remove_sample_ids, and keep them in our
% new matrix
keep_by_id = ~ismember(sample_ids_header, remove_sample_ids);

% Update the new matrix and new ids for labelling and QC later
updated_normalized_matrix = normalization_data(:, keep_by_id);
updated_sample_ids = sample_ids_header(keep_by_id);

sample_iqr1 = iqr(updated_normalized_matrix, 1);
average_correlation1 = SampleCorrelation(updated_normalized_matrix, "Spearman");

% Samples and alpha threshold check
n_samples1 = size(updated_normalized_matrix, 2);
alpha1 = computeAlphaOutliers(n_samples1);

% Look for IQR outliers
[outliers_iqr_qc2, outlier_ids_iqr_qc2, low_outliers_iqr2, high_outliers_iqr2] = ...
    scatterplotMarkOutliers(sample_iqr1, ...
                           'ColumnLabels', updated_sample_ids, ...
                           'PlotTitle', 'Quality Control after IQR Outlier Removal: Sample IQR', ...
                           'xlabel', 'Sample Index', ...
                           'ylabel', 'IQR', ...
                           'ShowXTickLabel', false);

% Look for outliers from correlation
[outliers_corr2, outlier_ids_corr2, low_outliers_corr2, high_outliers_corr2] = ...
    scatterplotMarkOutliers(average_correlation1, ...
                           'ColumnLabels', updated_sample_ids, ...
                           'PlotTitle', 'Quality Control after IQR Outlier Removal: Average Spearman Correlation', ...
                           'xlabel', 'Sample Index', ...
                           'ylabel', 'Average Correlation', ...
                           'ShowXTickLabel', false);

% Expression distribution for new matrix after outlier removal
figure('Position', [100 100 1200 500]);
boxplot(updated_normalized_matrix, 'PlotStyle', 'traditional', ...
        'Colors', 'k', 'Whisker', 1.5);
xlabel('Samples');
ylabel('Log2(RF Normalized Expression)');
title('miRNA-Seq DLBCL RF Norm Counts after IQR Outlier Removal: Expression Distribution'); 
grid on;

% After QC, repeat process and filter outliers 0.2 below correlation
% threshold with IQR calculations
Q1_corr = quantile(average_correlation1, 0.25);
Q3_corr = quantile(average_correlation1, 0.75);
IQR_corr = Q3_corr - Q1_corr;
lower_boundary_corr = Q3_corr - alpha * IQR_corr;
outlier_threshold = lower_boundary_corr - 0.2;

% Remove outliers more than 0.2 less than the alpha threshold with boolean
% mask
severe_outliers = (average_correlation1 < outlier_threshold);
to_remove = updated_sample_ids(severe_outliers);

keep = ~ismember((updated_sample_ids), to_remove);

final_data = updated_normalized_matrix(:, keep);  % chosen log2 matrix after removal
sample_ids_final = updated_sample_ids(keep);

final_data_iqr = iqr(final_data, 1);
final_data_corr = SampleCorrelation(final_data, 'Spearman');

% Quick post-removal QC 
[outliers_iqr_qc3, outlier_ids_iqr_qc3, low_outliers_iqr3, high_outliers_iqr3] = ...
    scatterplotMarkOutliers(final_data_iqr, ...
                           'ColumnLabels', sample_ids_final, ...
                           'PlotTitle', 'Quality Control after Correlation Outlier Removal: Sample IQR', ...
                           'xlabel', 'Sample Index', ...
                           'ylabel', 'IQR', ...
                           'ShowXTickLabel', false);

% Look for outliers from correlation
[outliers_corr3, outlier_ids_corr3, low_outliers_corr3, high_outliers_corr3] = ...
    scatterplotMarkOutliers(final_data_corr, ...
                           'ColumnLabels', sample_ids_final, ...
                           'PlotTitle', 'Quality Control after Correlation Outlier Removal: Average Spearman Correlation', ...
                           'xlabel', 'Sample Index', ...
                           'ylabel', 'Average Correlation', ...
                           'ShowXTickLabel', false);

% Boxplot expression distribution for final data
figure('Position', [100 100 1200 500]);
boxplot(final_data, 'PlotStyle', 'traditional', ...
        'Colors', 'k', 'Whisker', 1.5);
xlabel('Samples');
ylabel('Log2(RF Normalized Expression)');
title('miRNA-Seq DLBCL RF Norm Counts after Correlation Outlier Removal: Expression Distribution'); 
grid on;

% Filter lowly expressed genes with 0.90 quantile threshold and keep miRNAs
% that are not below the threshold (0.90 for strictness)
quantile_threshold = 0.90;

low_expression_mask = MarkLowCounts(final_data, quantile_threshold);

% Boolean mask to keep miRNAs not lowly expressed
keep_mirnas = ~low_expression_mask;

expression_final = final_data(keep_mirnas, :);

final_total = sum(expression_final, 1);
final_iqr = iqr(expression_final, 1);
final_corr = SampleCorrelation(expression_final, 'Spearman');

% Total expression per sample
figure('Name', 'Final Total Expression', 'Position', [100 100 800 500]);
plot(1:length(final_total), final_total, 'o-', ...
     'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
xlabel('Sample Index', 'FontSize', 12);
ylabel('Total Expression (Log2)', 'FontSize', 12);
title('Final Data: Total Expression per Sample', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% IQR per sample
figure('Name', 'Final IQR', 'Position', [100 100 800 500]);
plot(1:length(final_iqr), final_iqr, 'o-', ...
     'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('Sample Index', 'FontSize', 12);
ylabel('IQR', 'FontSize', 12);
title('Final Data: IQR per Sample', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Final distribution boxplot
figure('Name', 'Final Distribution', 'Position', [100 100 1200 500]);
boxplot(expression_final, 'PlotStyle', 'traditional', ...
        'Colors', 'k', 'Whisker', 1.5);
ylabel('Log2(Expression)', 'FontSize', 12);
xlabel('Samples', 'FontSize', 12);
title('Final Data: Expression Distribution', 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Final IQR Scatterplot
[outliers_iqr_qc4, outlier_ids_iqr_qc4, low_outliers_iqr4, high_outliers_iqr4] = ...
    scatterplotMarkOutliers(final_iqr, ...
                           'ColumnLabels', sample_ids_final, ...
                           'PlotTitle', 'Final Quality Control Sample IQR', ...
                           'xlabel', 'Sample Index', ...
                           'ylabel', 'IQR', ...
                           'ShowXTickLabel', false);

% Final Spearman Correlation Scatterplot
[outliers_corr4, outlier_ids_corr4, low_outliers_corr4, high_outliers_corr4] = ...
    scatterplotMarkOutliers(final_corr, ...
                           'ColumnLabels', sample_ids_final, ...
                           'PlotTitle', 'Final Quality Control: Average Spearman Correlation', ...
                           'xlabel', 'Sample Index', ...
                           'ylabel', 'Average Correlation', ...
                           'ShowXTickLabel', false);

% Final clustergram/heatmap for viz, where we calculate variance across
% each miRNA to see differential expresssion and sort with highest variance
% first, selecting the top 50 and grabbing our sample/expression indices
% for display with euclidean distance
mirna_variance = var(expression_final, 0, 2);
[~, top_var_idx] = sort(mirna_variance, 'descend');
top_n = min(50, size(expression_final, 1)); 
top_mirnas_idx = top_var_idx(1:top_n);

expression_top = expression_final(top_mirnas_idx, :);
mirna_ids_top = mirna_ids(keep_mirnas);
mirna_ids_top = mirna_ids_top(top_mirnas_idx);

figure('Position', [100 100 1400 800]);

% Create clustergram with proper labeling
cg = clustergram(expression_top, ...
                 'Rowlabels', mirna_ids_top, ...
                 'ColumnLabels', sample_ids_final, ...
                 'Colormap', redbluecmap, ...
                 'DisplayRange', 10, ...
                 'RowPDist', 'euclidean', ...
                 'ColumnPDist', 'correlation', ...
                 'Linkage', 'average', ...
                 'OptimalLeafOrder', true, ...
                 'Standardize', 'row');

addTitle(cg, 'Hierarchical Clustering: Top 50 Variable miRNAs in DLBCL Samples');


 

















