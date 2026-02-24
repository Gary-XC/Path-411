% PATH 411: Group Assignment 2
% main_gA2_13.m
% Author(s): Sebastian, Lenard, Gary
%
% This assignment identifies discriminative miRNA features to distinguish
% between gastrointestinatl neuroendocrine tumours (GI-NET) locations and
% develop classification models for foregut and midgut tumours.
%
% The dataset is GI-NET miRNA data with 69 patients consisting of
% 24 pancreatic, 26 ileal, 12 appendiceal, 7 rectal with locations where
% 24 samples are from the foregut, 38 are from the midgut, and 7 are from
% hindgut.
%
% We assume the data has already been preprocessed with outliers removed, 
% batch effects corrected, low-expressed miRNAs filtered at 90th threshold.
%
% Process: 
%
% 1. EDA with t-SNE dimensionality reduction to visualize clustering by
% cancer type.
% 2. Hierarchical clustering with all miRNAs first, comparing by type
% (pancreatic, ileal, appendiceal, rectal), and by grade (low-grade,
% intermediate-grade) with multiple distance metrics tested (correlation,
% euclidean, cosine) and multiple linkages (single, average, complete) to
% view baseline clustering. Data was first log transformed and median
% centered.
% 3. Comparing feature selection algorithms with relieff, fscchi2, and
% fscmrmr with only forgut and midgut location samples, where data has been
% log transformed and median centered.
% 4. Choosing an appropriate number of features from the selection by 
% looking at feature importance and
% score distribution to balance discrimination and prevent overfitting.
% 5. Hierarchical clustering performed with only the new selected feature
% to compare clustering and validate that these features improved class
% seperation. 
% 6. Creating the classification models through the Classification App with
% various machine learning techniques (SVM, KNN, Decision Trees, Ensembles,
% etc.) with/without PCA and evaluating validation accuracy, confusion
% matrices, ROC-AUC, PR-AUC, etc. to find best models.

% DATASET LOADED AS A WORKSPACE ADDED TO PATH
load('GINETdata_forA2.mat')

% Grab miRNA expression data (miRNA (rows) x samples (columns))
miRNA_names = celldata_miRNAs_4analysis(2:end, 1);
miRNA_data = cell2mat(celldata_miRNAs_4analysis(2:end, 2:end));

% Setting seed for reproducibility for t-SNE
rng(42);

% Calculate t-SNE embeddings with barneshut algorithm and spearman distance
% perplexity set to 20 to determine effective neighbours and balance
% local/global data
tSNE_miRNA = tsne(miRNA_data', ...
    'Algorithm', 'barneshut', ...
    'Distance', 'spearman', ...
    'Perplexity', 20); 

% Extract cancer type from clinical data
cancer_type = tbl_clinical.Type;

% Plot t-SNE
figure('Position', [100, 100, 800, 600]);
gscatter(tSNE_miRNA(:,1), tSNE_miRNA(:,2), cancer_type);
xlabel('t-SNE Dimension 1');
ylabel('t-SNE Dimension 2');
title('t-SNE Visualization of GI-NET miRNA Data by Type');
legend('Location', 'best');

% Hierarchical clustering with log transform and median centering
miRNA_replace_zeros = replaceZeros(miRNA_data, 'lowval');
miRNA_log = log2(miRNA_replace_zeros);

% Calculate median of each miRNA across all samples and center
median_value = median(miRNA_log, 2);
miRNA_median_centered = miRNA_log - median_value;

% Extracting labels for clustering
sample_labels_type = tbl_clinical.Type;
sample_labels_grade = tbl_clinical.Grade;

% HIERARCHICAL CLUSERTING FOR TYPE AND GRADE: Testing distances 
% for correlation (captures similarity in expression regardless or
% magnitude), euclidean (standard distance sensitive to magnitude
% differences), cosine (measuring angles between vectors; good for
% high-dimensional data)
% Linkages methods tested are single (minimum distance between clusters
% which may create long chains), average (balanced for average distance
% between all pairs), complete (max distance between clusters for compact
% clustering)

% HIERARCHICAL CLUSTERING BY CANCER TYPE

% Correlation Distance + Single Linkage
figure('Position', [100 100 1400 800]);
cg_type_corr_single = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_type, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'correlation', ...
    'Linkage', 'single');
addXLabel(cg_type_corr_single, 'Cancer Type');
addYLabel(cg_type_corr_single, 'miRNA');
addTitle(cg_type_corr_single, 'Hierarchical Clustering by Type: Correlation Distance + Single Linkage');

% Correlation Distance + Average Linkage
figure('Position', [100 100 1400 800]);
cg_type_corr_avg = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_type, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'correlation', ...
    'Linkage', 'average');
addXLabel(cg_type_corr_avg, 'Cancer Type');
addYLabel(cg_type_corr_avg, 'miRNA');
addTitle(cg_type_corr_avg, 'Hierarchical Clustering by Type: Correlation Distance + Average Linkage');

% Correlation Distance + Complete Linkage
figure('Position', [100 100 1400 800]);
cg_type_corr_comp = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_type, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'correlation', ...
    'Linkage', 'complete');
addXLabel(cg_type_corr_comp, 'Cancer Type');
addYLabel(cg_type_corr_comp, 'miRNA');
addTitle(cg_type_corr_comp, 'Hierarchical Clustering by Type: Correlation Distance + Complete Linkage');

% Euclidean Distance + Single Linkage
figure('Position', [100 100 1400 800]);
cg_type_euc_single = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_type, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'euclidean', ...
    'Linkage', 'single');
addXLabel(cg_type_euc_single, 'Cancer Type');
addYLabel(cg_type_euc_single, 'miRNA');
addTitle(cg_type_euc_single, 'Hierarchical Clustering by Type: Euclidean Distance + Single Linkage');

% Euclidean Distance + Average Linkage
figure('Position', [100 100 1400 800]);
cg_type_euc_avg = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_type, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'euclidean', ...
    'Linkage', 'average');
addXLabel(cg_type_euc_avg, 'Cancer Type');
addYLabel(cg_type_euc_avg, 'miRNA');
addTitle(cg_type_euc_avg, 'Hierarchical Clustering by Type: Euclidean Distance + Average Linkage');

% Euclidean Distance + Complete Linkage
figure('Position', [100 100 1400 800]);
cg_type_euc_comp = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_type, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'euclidean', ...
    'Linkage', 'complete');
addXLabel(cg_type_euc_comp, 'Cancer Type');
addYLabel(cg_type_euc_comp, 'miRNA');
addTitle(cg_type_euc_comp, 'Hierarchical Clustering by Type: Euclidean Distance + Complete Linkage');

% Cosine Distance + Single Linkage
figure('Position', [100 100 1400 800]);
cg_type_cos_single = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_type, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'cosine', ...
    'Linkage', 'single');
addXLabel(cg_type_cos_single, 'Cancer Type');
addYLabel(cg_type_cos_single, 'miRNA');
addTitle(cg_type_cos_single, 'Hierarchical Clustering by Type: Cosine Distance + Single Linkage');

% Cosine Distance + Average Linkage
figure('Position', [100 100 1400 800]);
cg_type_cos_avg = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_type, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'cosine', ...
    'Linkage', 'average');
addXLabel(cg_type_cos_avg, 'Cancer Type');
addYLabel(cg_type_cos_avg, 'miRNA');
addTitle(cg_type_cos_avg, 'Hierarchical Clustering by Type: Cosine Distance + Average Linkage');

% Cosine Distance + Complete Linkage
figure('Position', [100 100 1400 800]);
cg_type_cos_comp = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_type, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'cosine', ...
    'Linkage', 'complete');
addXLabel(cg_type_cos_comp, 'Cancer Type');
addYLabel(cg_type_cos_comp, 'miRNA');
addTitle(cg_type_cos_comp, 'Hierarchical Clustering by Type: Cosine Distance + Complete Linkage');

% HIERARCHICAL CLUSTERING BY CANCER GRADE

% Correlation Distance + Single Linkage
figure('Position', [100 100 1400 800]);
cg_grade_corr_single = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_grade, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'correlation', ...
    'Linkage', 'single');
addXLabel(cg_grade_corr_single, 'Tumor Grade');
addYLabel(cg_grade_corr_single, 'miRNA');
addTitle(cg_grade_corr_single, 'Hierarchical Clustering by Grade: Correlation Distance + Single Linkage');

% Correlation Distance + Average Linkage
figure('Position', [100 100 1400 800]);
cg_grade_corr_avg = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_grade, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'correlation', ...
    'Linkage', 'average');
addXLabel(cg_grade_corr_avg, 'Tumor Grade');
addYLabel(cg_grade_corr_avg, 'miRNA');
addTitle(cg_grade_corr_avg, 'Hierarchical Clustering by Grade: Correlation Distance + Average Linkage');

% Correlation Distance + Complete Linkage
figure('Position', [100 100 1400 800]);
cg_grade_corr_comp = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_grade, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'correlation', ...
    'Linkage', 'complete');
addXLabel(cg_grade_corr_comp, 'Tumor Grade');
addYLabel(cg_grade_corr_comp, 'miRNA');
addTitle(cg_grade_corr_comp, 'Hierarchical Clustering by Grade: Correlation Distance + Complete Linkage');

% Euclidean Distance + Single Linkage
figure('Position', [100 100 1400 800]);
cg_grade_euc_single = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_grade, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'euclidean', ...
    'Linkage', 'single');
addXLabel(cg_grade_euc_single, 'Tumor Grade');
addYLabel(cg_grade_euc_single, 'miRNA');
addTitle(cg_grade_euc_single, 'Hierarchical Clustering by Grade: Euclidean Distance + Single Linkage');

% Euclidean Distance + Average Linkage
figure('Position', [100 100 1400 800]);
cg_grade_euc_avg = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_grade, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'euclidean', ...
    'Linkage', 'average');
addXLabel(cg_grade_euc_avg, 'Tumor Grade');
addYLabel(cg_grade_euc_avg, 'miRNA');
addTitle(cg_grade_euc_avg, 'Hierarchical Clustering by Grade: Euclidean Distance + Average Linkage');

% Euclidean Distance + Complete Linkage
figure('Position', [100 100 1400 800]);
cg_grade_euc_comp = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_grade, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'euclidean', ...
    'Linkage', 'complete');
addXLabel(cg_grade_euc_comp, 'Tumor Grade');
addYLabel(cg_grade_euc_comp, 'miRNA');
addTitle(cg_grade_euc_comp, 'Hierarchical Clustering by Grade: Euclidean Distance + Complete Linkage');

% Cosine Distance + Single Linkage
figure('Position', [100 100 1400 800]);
cg_grade_cos_single = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_grade, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'cosine', ...
    'Linkage', 'single');
addXLabel(cg_grade_cos_single, 'Tumor Grade');
addYLabel(cg_grade_cos_single, 'miRNA');
addTitle(cg_grade_cos_single, 'Hierarchical Clustering by Grade: Cosine Distance + Single Linkage');

% Cosine Distance + Average Linkage
figure('Position', [100 100 1400 800]);
cg_grade_cos_avg = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_grade, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'cosine', ...
    'Linkage', 'average');
addXLabel(cg_grade_cos_avg, 'Tumor Grade');
addYLabel(cg_grade_cos_avg, 'miRNA');
addTitle(cg_grade_cos_avg, 'Hierarchical Clustering by Grade: Cosine Distance + Average Linkage');

% Cosine Distance + Complete Linkage
figure('Position', [100 100 1400 800]);
cg_grade_cos_comp = clustergram(miRNA_median_centered, ...
    'RowLabels', miRNA_names, ...
    'ColumnLabels', sample_labels_grade, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 3, ...
    'RowPDist', 'euclidean', ...
    'ColumnPDist', 'cosine', ...
    'Linkage', 'complete');
addXLabel(cg_grade_cos_comp, 'Tumor Grade');
addYLabel(cg_grade_cos_comp, 'miRNA');
addTitle(cg_grade_cos_comp, 'Hierarchical Clustering by Grade: Cosine Distance + Complete Linkage');

%% PERFORMING FEATURE SELECTION WITH RELIEFF to identify discriminative
% miRNAs from foregut and midgut

% Relieff is good for high-D data and ranks features by ability to seperate
% classes while also looking out for nearest neighours

% Grabbing location labels
location = tbl_clinical.Location;

% Filter and search for strings to match that type and apply the boolean to
% our data subset so it only uses foregut and midgut locations
filter_location_names = strcmp(location, 'foregut') | strcmp(location, 'midgut');
location_subset = miRNA_data(:, filter_location_names);
labels = location(filter_location_names);

% Log-transform and median center
relieff_replace_zeros = replaceZeros(location_subset, 'lowval');
relieff_log = log2(relieff_replace_zeros);

% Calculate median of each miRNA across all samples
relieff_median_value = median(relieff_log, 2);
relieff_median_centered = relieff_log - relieff_median_value;

% Transpose the data for the relieff method with k = 10
X = relieff_median_centered';
[~, weights] = relieff(X, labels, 10);

% Sort by feature importance in descending order
[sorted_weights, sorted_index] = sort(weights, 'descend');

% Plotting first 100 to see changes in importance to grab top n features
figure;
plot(sorted_weights(1:100), 'b-', 'LineWidth', 2);
xlabel('Feature Rank', 'FontSize', 12);
ylabel('ReliefF Score', 'FontSize', 12);
title('Top 100 Features with ReliefF for Foregut and Midgut miRNA GI-NET Samples', 'FontSize', 14);
yticks(0:0.05:0.3);
grid on;

% Selecting features
n_features = 10;

% Extract the indices and selected features for clustering
top_index = sorted_index(1:n_features);

% Subset to top
select_data = relieff_median_centered(top_index, :);

% Grabbing names
selected_names = miRNA_names(top_index);

% Seeing scores of selected features
figure;
bar(sorted_weights(1:n_features));
xlabel('Selected Features');
ylabel('ReliefF Score');
title('Top 10 Selected Features with ReliefF for Foregut and Midgut miRNA GI-NET Samples');

% Hierarchical Clustering with selected features to validate better
% grouping
figure('Position', [100 100 1400 800]);
cg = clustergram(select_data, ...
    'RowLabels', selected_names, ...
    'ColumnLabels', labels, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 4, ...
    'ColumnPDist', 'correlation', ...
    'Linkage', 'average');
addXLabel(cg, 'Location (Foregut vs Midgut)');
addYLabel(cg, 'Selected miRNAs');
addTitle(cg, 'Hierarchical Clustering: Top 10 Features (ReliefF)');

% Creating table for the Classification App to run the machine learning
% with our selected feature data

% Transposing the data for samples x features
class_data = select_data';

% Table with miRNA names as the variables
class_table = array2table(class_data, 'VariableNames', selected_names);

% Location labels for the response
class_table.Location = labels;


%% PERFORMING FEATURE SELECTION WITH fscchi2 to identify discriminative
% miRNAs from foregut and midgut

raw = celldata_miRNAs_4analysis;
miRNA_names = raw(2:end, 1); % (264 × 1) list of miRNA identifiers
exprMatrix = cell2mat(raw(2:end, 2:end)); % (264 × 69) numeric matrix (miRNAs × samples)

% Extracting the clinical label for machine learning: Location (foregut / midgut)

locAll = tbl_clinical.Location;

isForegut = strcmp(locAll, 'foregut');
isMidgut = strcmp(locAll, 'midgut');
keepIdx = isForegut | isMidgut;

expr_sub = exprMatrix(:, keepIdx); % subset of samples
loc_sub = locAll(keepIdx); % labels for selected samples
Y = categorical(loc_sub); % convert to categorical that is required by fscchi2

% Feature selection using fscchi2 (Filter Method)

X = expr_sub.';     % Converting to (nSub × 264) for fscchi2

[idxFeat, scores] = fscchi2(X, Y);
% idxFeat to sort feature indices sorted best → worst
% scores are chi-square values the listing of feature importance

%  Displaying the top 20 features
topK = 20; 
topIdx = idxFeat(1:topK);
topNames = miRNA_names(topIdx);
topScores = scores(topIdx);

T = table(topNames(:), topScores(:));  

figure;
bar(topScores); % bar heights
set(gca, 'XTick', 1:topK); % set x-axis ticks
set(gca, 'XTickLabel', topNames); % label with miRNA names
xtickangle(45); % rotating the labels
ylabel('Chi-square Score');
title('Top 20 miRNAs Ranked by Chi-square Score');
grid on;

% Ploting the feature importance curve for the first 150 features(to avoid
% the trailing tail for features later on
figure;
n = 150;  % number of points to plot

line(1:n, scores(idxFeat(1:n)), 'LineWidth', 2);
xlabel('Feature Rank');
ylabel('Chi-square Score');
title('Feature Importance using fscchi2 For Location (Foregut vs Midgut)');
grid on;

% Creating a dataset of selected miRNAs for the top 20 features

K = 20; % Selecting the top 20 features
selectedFeatureIdx   = idxFeat(1:K); % indexing the features from 1 to K (20)
selectedFeatureNames = miRNA_names(selectedFeatureIdx); % grabbing the miRNAs names that were selected

expr_sel = expr_sub(selectedFeatureIdx, :);   % K × nSub

% Log transform to reduce the skewness for gene expression preprocessing
log_sel = log2(expr_sel + 1e-6);

% Median-centering: subtracting per-miRNA median across samples
med_sel = log_sel - median(log_sel, 2);

% Hierarchical clustering using selected miRNAs
% Input for clustergram: (features × samples)

selectedNames = miRNA_names(selectedFeatureIdx); % RowLabels

cg_selected = clustergram(med_sel, ...
    'RowLabels',   selectedNames, ...  % Selected miRNA names
    'ColumnLabels', cellstr(Y_loc_final), ... % foregut / midgut
    'Standardize', 'none', ... % data has already centered
    'RowPDist',    'correlation', ... % correlation distance
    'ColumnPDist', 'correlation', ...
    'Colormap',    redbluecmap); % cluster samples by similarity

addTitle(cg_selected,'Hierarchical Clustergram Based on Top 20 Chi-Square Ranked miRNAs (Location)'); 
% Naming the graph

% Classification Model Data
X_sel_final = med_sel.'; % Features
Y_loc_final = Y; % Label that the classification model is guessing

%% PERFORMING FEATURE SELECTION WITH fscmrmr to identify discriminative

cell_clinical = table2cell(tbl_clinical);

% Extract names of miRNAs
miRNA_names = celldata_miRNAs_4analysis (2:end,1);

% Expression matrix for miRNA seq data
miRNA_data = cell2mat(celldata_miRNAs_4analysis(2:end, 2:end));

% Extracting samples by location (foregut or midgut)

% Extract samples with cancer located in either the foregut or midgut
cancer_location = cell_clinical(:, 14);
index_cancer_foregut_midgut = ismember (cancer_location, {'foregut', 'midgut'});

% Filter GI-NET miRNA data for samples with cancer in the foregut or midgut
fore_mid_miRNA_data = miRNA_data (:, index_cancer_foregut_midgut');

% Filter clinical information for samples with cancer in the foregut or
% midgut
fore_mid_clinical_info = cell_clinical(index_cancer_foregut_midgut, :);
fore_mid_data = fore_mid_clinical_info (:, 14);

% Transformation and median centering of the data

% Replace zeros and log2 transform the foregut and midgut data
transform_fore_mid_miRNA_data = log2(replaceZeros(fore_mid_miRNA_data, 'lowval'));

% Median center the foregut and midgut data
fore_mid_miRNA_median_centered = transform_fore_mid_miRNA_data - median (transform_fore_mid_miRNA_data, 2);

% Feature selection using fscmrmr

% Perform feature selection
[fore_mid_feature_select, feature_scores] = fscmrmr (fore_mid_miRNA_median_centered', fore_mid_data); 

% Plot the predictor importance scores for miRNA features
bar(feature_scores(fore_mid_feature_select))
xlabel ('Predictor Rank')
ylabel ('Predictor Importance Score')
title ('Feature Selection using fscmrmr')
grid on;

% Sorting predictor importance scores in descending order
[sorted_feature_scores, importance_index] = sort (feature_scores, 'descend');

% Plotting first 50 features to identify those of greatest importance
figure;
bar(sorted_feature_scores(1:50));
xlabel ('Predictor Rank')
ylabel ('Predictor Importance Score')
title ('Top 50 Features using fscmrmr')
grid on;

% Selecting the top 11 features
top_features = 11;

% Plotting 11 features to identify those of greatest importance
figure;
bar(sorted_feature_scores(1:top_features));
xlabel ('Predictor Rank')
ylabel ('Predictor Importance Score')
title ('Top 11 Features using fscmrmr')
grid on;

% Creating an index to extract top features
top_11_scores_index = importance_index (1:top_features);

% Extracting miRNA data from top 11 features
top_11_scores_data = fore_mid_miRNA_median_centered (top_11_scores_index, :);

% Extracting miRNA name of top 11 features
top_11_scores_names = miRNA_names (top_11_scores_index);

% Hierarchical Clustering of top 11 features using fscmrmr 
cg_fscmrmr = clustergram(top_11_scores_data, ...
    'RowLabels', top_11_scores_names, ...
    'ColumnLabels', fore_mid_data, ...
    'Colormap', redbluecmap, ...
    'DisplayRange', 4, ...
    'ColumnPDist', 'correlation', ...
    'Linkage', 'average');
addXLabel(cg_fscmrmr, 'Location');
addYLabel(cg_fscmrmr, 'miRNAs');
addTitle(cg_fscmrmr, 'Hierarchical Clustering of Top 11 Features using fscmrmr');

% Prepare data for classification model
miRNA_data_classifier = top_11_scores_data;
miRNA_data_classifier_labels = fore_mid_data;