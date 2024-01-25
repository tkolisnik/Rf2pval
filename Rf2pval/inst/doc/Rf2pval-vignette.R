## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----eval=FALSE---------------------------------------------------------------
#  Installation:
#  1. Install devtools if not already installed:
#  install.packages("devtools")
#  library(devtools)
#  
#  2. Use devtools to install Rf2pval
#  devtools::install_github("tkolisnik/Rf2pval", build_vignettes = TRUE)
#  
#  3. Load Rf2pval
#  library(Rf2pval)

## ----eval=FALSE---------------------------------------------------------------
#  if (!requireNamespace("reticulate")) {
#      install.packages("reticulate")
#  }
#  library(reticulate)
#  
#  reticulate::use_python("/usr/local/bin/python3")
#  reticulate::conda_create(envname = "rf2pval-conda")
#  reticulate::conda_install(envname = "rf2pval-conda", packages = c("scikit-learn", "numpy","shap"))
#  reticulate::use_condaenv("rf2pval-conda", required = TRUE)
#  
#  # You can check if this has worked by running
#  reticulate::py_config()
#  
#  # If conda installation has worked, the previous command will result in an output message similar to:
#  python:         /opt/homebrew/Caskroom/miniconda/base/envs/rf2pval-conda/bin/python
#  libpython:      /opt/homebrew/Caskroom/miniconda/base/envs/rf2pval-conda/lib/libpython3.9.dylib
#  pythonhome:     /opt/homebrew/Caskroom/miniconda/base/envs/rf2pval-conda:/opt/homebrew/Caskroom/miniconda/base/envs/rf2pval-conda
#  version:        3.9.18 | packaged by conda-forge | (main, Aug 30 2023, 03:53:08)  [Clang 15.0.7 ]
#  numpy:          /opt/homebrew/Caskroom/miniconda/base/envs/rf2pval-conda/lib/python3.9/site-packages/numpy
#  numpy_version:  1.26.0
#  sklearn:        /opt/homebrew/Caskroom/miniconda/base/envs/rf2pval-conda/lib/python3.9/site-packages/sklearn
#  
#  NOTE: Python version was forced by use_python() function
#  
#  # You may now proceed to using the package

## ----eval=FALSE---------------------------------------------------------------
#  load("/path/to/Rf2pval/data/demo_rnaseq_data.RData")

## ----dataframe-creation, echo=FALSE-------------------------------------------
library(knitr)
load("../data/demo_rnaseq_data.RData")

kable(head(demo_rnaseq_data$training_data[1:10,1:6], 10), 
      caption = "training_data", 
      align = 'l')

kable(head(demo_rnaseq_data$validation_data[1:10,1:6], 10), 
      caption = "testing_data", 
      align = 'l')

kable(head(demo_rnaseq_data$testing_data[1:10,1:6], 10), 
      caption = "validation_data", 
      align = 'l')

kable(head(demo_rnaseq_data$target_categories[1:10,1:2], 10), 
      caption = "target_categories", 
      align = 'l')

## ----eval=FALSE---------------------------------------------------------------
#  # Load target data:
#  training_data <- read_csv("/path/to/file/training_data.csv")
#  validation_data <- read_csv("/path/to/file/validation_data.csv")
#  testing_data <- read_csv("/path/to/file/testing_data.csv")
#  target_categories <- read_csv("/path/to/file/target_categories.csv")

## ----eval=FALSE---------------------------------------------------------------
#  # Load the tibble package
#  library(tibble)
#  
#  # Ensure each dataset is converted to a tibble and combine them into a list
#  demo_rnaseq_data <- list(
#    training_data = as_tibble(training_data),
#    validation_data = as_tibble(validation_data),
#    testing_data = as_tibble(testing_data),
#    target_categories = as_tibble(target_categories)
#  )

## ----eval=FALSE---------------------------------------------------------------
#  # This function prepares a dataset for scikit-learn random forest machine learning by creating a matrix from the dataset excluding certain columns, and extracting a target vector. It is flexible and can handle different types of data sets such as training, testing, or validation. Data must be in the correct input format, see vignette or example dataset for details.
#  processed_training_data <- Rf2pval::create_feature_matrix(demo_rnaseq_data$training_data,"training")
#  
#  processed_validation_data <- Rf2pval::create_feature_matrix(demo_rnaseq_data$validation_data,"validation")
#  
#  processed_testing_data <- Rf2pval::create_feature_matrix(demo_rnaseq_data$testing_data,"testing")

## ----eval=FALSE---------------------------------------------------------------
#  # This function uses scikit-learn's python based GridSearchCV to perform hyperparameter tuning and training of a RandomForestClassifier. It allows for customizable parameter grids and includes preprocessing steps of one-hot encoding and scaling. The function is designed to find the best hyperparameters based on accuracy. Please reference the scikit-learn GridSearchCV documentation for the full description of options, however our defaults are comprehensive.
#  
#  # If running defaults:
#  tuning_results <- Rf2pval::tune_and_train_rf_model(processed_training_data$X_training_mat, processed_training_data$y_training_vector)
#  print(tuning_results$grid_search$best_params_)
#  print(tuning_results$grid_search$best_score_)
#  
#  # If not running defaults and customization of parameter is desired:
#  custom_parameter_grid <- list(
#    bootstrap = list(TRUE),
#    class_weight = list(NULL),
#    max_depth = list(5L, 20L, NULL),
#    n_estimators = as.integer(seq(10, 90, 10)),
#    max_features = list("sqrt", 0.2),
#    criterion = list("gini"), # It is not recommended to change this, as it is used for Rf2pval's rank-based feature reduction approach.
#    warm_start = list(FALSE),
#    min_samples_leaf = list(1L, 50L),
#    min_samples_split = list(2L, 200L)
#  )
#  
#  # Train your model and perform hyperparameter tuning using Rf2pval. This can be extremely time consuming and scales rapidly with dataset size.
#  # A small dataset such as the demo data may take >5 minutes, however a larger 300 x 50,000 dataset may take days depending on the number of cores used.
#  # It is strongly recommended to increase the number of cores to reduce the time this takes.
#  
#  tuning_results <- Rf2pval::tune_and_train_rf_model(processed_training_data$X_training_mat, processed_training_data$y_training_vector,
#                                  scoring='f1',
#                                  seed = 123,
#                                  param_grid = custom_parameter_grid,
#                                  n_jobs = 1,
#                                  n_cores = 7)
#  
#  # After the model is tuned, print the best parameters and the best score. See scikit-learn documentation for interpreting scores. It is recommended to re-tune on different parameters if score less than 70.
#  
#  print(tuning_results$grid_search$best_params_)
#  print(tuning_results$grid_search$best_score_)

## ----eval=FALSE---------------------------------------------------------------
#  # This function fits a Random Forest model using the provided hyperparameters and training data, then evaluates its performance on a validation set.
#  
#  # Fit the results to the validation set to check performance
#  fitting_results <- Rf2pval::fit_and_evaluate_rf(tuning_results$best_params,processed_training_data$X_training_mat,processed_training_data$y_training_vector,processed_validation_data$X_validation_mat,processed_validation_data$y_validation_vector)
#  # Print the fitting results, provides accuracy, f1 score, precision, recall and roc_auc scores on the model as fitted to the validation set
#  print(fitting_results)
#  
#  # If performance is not adequate, you may want to reconstruct your model, tweak your data, and rerun step 7 before proceeding.
#  # This function is also used to evaluate the testing set when the model is finalized.
#  fitting_results <- Rf2pval::fit_and_evaluate_rf(tuning_results$best_params,processed_training_data$X_training_mat,processed_training_data$y_training_vector,processed_testing_data$X_testing_mat,processed_testing_data$y_testing_vector)
#  # Print the fitting results, provides accuracy, f1 score, precision, recall and roc_auc scores on the model as fitted to the validation set
#  print(fitting_results)

## ----eval=FALSE---------------------------------------------------------------
#  # This function fits a Random Forest model and calculates the true and permuted feature importances. It performs permutations on the target variable to generate permuted importances for comparison.
#  # This step is also time consuming and scales with the number of permutations used.
#  feat_importances <- Rf2pval::calculate_feature_importances(fitting_results$model,processed_training_data$X_training_mat,processed_training_data$y_training_vector,n_permutations=1000)
#  
#  # Print top features from model pre-Rank-based feature reduction
#  print(feat_importances$top_features)

## ----eval = FALSE-------------------------------------------------------------
#  # Calculate quantiles
#  # This function computes the mean, lower, and upper quantiles of permuted feature importance scores and compares them with the observed importance scores from the true values. It is used to assess the significance of feature importances from random forest models by comparing them against a distribution of importances obtained through permutation.
#  quantile_data <- Rf2pval::calculate_quantiles(feat_importances$true_importances, feat_importances$permuted_importances)
#  
#  
#  # Calculate p-value for the entire feature set
#  # This function calculates the proportion-based p-value for the entire feature set by comparing the observed sum of absolute deviations from the mean feature importance to those from permuted data. It is particularly useful in permutation tests to assess the statistical significance of the observed data.
#  pvalue_set <- Rf2pval::calculate_full_set_pvalue(feat_importances$permuted_importances, quantile_data)
#  
#  # Calculate p-value for each rank
#  # This function calculates p-values for each feature rank in the dataset by comparing the observed feature importances against the distribution of importances in the permuted data. It returns ranks and their respective p-values and proportions up to a specified alpha threshold.
#  pvalues_ranks <- Rf2pval::calculate_ranked_based_pvalues(feat_importances$true_importances, feat_importances$permuted_importances, alpha=0.05)
#  
#  # Print top features from model post-Rank based feature reduction
#  print(pvalues_ranks)

## ----eval = FALSE-------------------------------------------------------------
#  # Generate Full feature importance plot
#  # This function generates a figure which shows the true observed data plotted against the permuted data, by rank. The intersection of the true data with the upper quartile is shown, which we recommend as a significance cutoff. Note: There are ample parameters for controlling the axes scale, label location, and zoom, because of data variability, you will almost certainly have to adjust these to fit your plot.
#  fiplot_full <- Rf2pval::generate_fi_rank_plot(feat_importances$permuted_importances, quantile_data, xlimitmin=1, xlimitmax = 10, ylimitmin= -5, ylimitmax = 0, labelhorizontaladjust = -0.05,labelverticaladjust = 1.5,focusedView = FALSE,logOn = TRUE)
#  
#  print(fiplot_full)
#  
#  # Same plot as above, but zoomed to show only the significant results
#  fiplot_focused <- Rf2pval::generate_fi_rank_plot(feat_importances$permuted_importances, quantile_data, xlimitmin=1, xlimitmax = 10, ylimitmin= -5, ylimitmax = 0, labelhorizontaladjust = 1.05,labelverticaladjust = 1.5,focusedView = TRUE,logOn = TRUE)
#  
#  print(fiplot_focused)
#  # Generate pi Histogram plot
#  # Creates a histogram showing the sum of absolute deviations by count, the measures used to calculate the pi statistics used in the p-value calculation for the set.
#  pihist_plot <- Rf2pval::generate_pi_histogram(feat_importances$permuted_importances, quantile_data)
#  
#  print(pihist_plot)
#  # Generate ECDF plot
#  # Creates a plot of the empirical cumulative density function showing the sum of absolute deviations by fraction of data.
#  ecdf_plot <- Rf2pval::generate_pi_ECDF_plot(feat_importances$permuted_importances, quantile_data)
#  
#  print(ecdf_plot)

## ----eval=FALSE---------------------------------------------------------------
#  # Produce a list of SHAP values
#  # This function calculates SHAP values for a given dataset using a provided model. It then identifies significant features based on the SHAP values for the specified class. Additionally, it prepares a long-format data frame of individual SHAP values suitable for visualization.
#  shapvals <- Rf2pval::calculate_SHAP_values(fitting_results$model, processed_training_data$X_training_mat, class_index = 1, shap_std_dev_factor = 0.5)
#  
#  # Produce the SHAP Plots
#  # Generates a combined bar and beeswarm plot to show global and local feature importance
#  shapplot <- Rf2pval::generate_shap_plots(mean_shap_values = shapvals$significant_features,long_shap_data = shapvals$long_shap_data)
#  
#  print(shapplot)

## ----eval=FALSE---------------------------------------------------------------
#  library(clusterProfiler)
#  library(org.Hs.eg.db) # Assuming your data is from human, if not see bioconductor for other organism dbs.
#  library(enrichplot)
#  
#  # We have 3 lists of significant features from our data that can be explored further
#  
#  # The 'raw' list of siginificant features from the RF model (pre-Rank based feature reduction):
#  print(feat_importances$top_features)
#  
#  # The list of significant features after rank-based permutation feature reduction (most conservative list):
#  print(pvalues_ranks)
#  
#  # The list of significant features from the SHAP analysis (can be further divided by directionality +/-):
#  print(shapvals$significant_features)
#  
#  # Subsequent steps can be ran on any of these feature lists, although if the list is too small these steps may not produce results.
#  features_list <- shapvals$significant_features$feature # Features from SHAP
#  
#  # Alternatively:
#  #features_list <- pvalues_ranks$feature # Post-feature reduction (most conservative but least amount of features so may not have desired results)
#  #features_list <- feat_importances$top_features # Top features from model pre-feature reduction (least robust)
#  
#  
#  # Run the GO Enrichment Analysis
#  # From clusterProfiler: GO Enrichment Analysis of a gene set. Given a vector of genes, this function will return the enrichment GO categories after FDR control.
#  ego <- enrichGO(gene         = features_list,
#                  OrgDb        = org.Hs.eg.db,
#                  keyType      = "ENSEMBL",
#                  ont          = "ALL", # consider BP, CC, and MF ontologies
#                  pAdjustMethod = "none", # adjust p-values using the BH method
#                  readable     = TRUE)  # convert Entrez IDs back to gene symbols
#  
#  # Plot the gene ontology
#  enrichplot::dotplot(ego, showCategory=20) + ggtitle("Gene Ontology for Top Features")
#  

## ----eval=FALSE---------------------------------------------------------------
#  library(gprofiler2)
#  
#  # Perform the enrichment analysis
#  # From gprofiler2: Interface to the g:Profiler tool g:GOSt (https://biit.cs.ut.ee/gprofiler/gost) for functional enrichments analysis of gene lists. In case the input 'query' is a list of gene vectors, results for multiple queries will be returned in the same data frame with column 'query' indicating the corresponding query name. If 'multi_query' is selected, the result is a data frame for comparing multiple input lists, just as in the web tool.
#  gprofiler_results <- gost(query = features_list, organism = "hsapiens", ordered_query = FALSE, user_threshold = 0.5)
#  
#  # Create gProfiler Manhattan plot
#  # From gprofiler2: This function creates a Manhattan plot out of the results from gprofiler2::gost().
#  gostplot_mh_plot<-gostplot(gprofiler_results, capped = FALSE, interactive = F)
#  print(gostplot_mh_plot)
#  
#  # Create gProfiler Manhattan plot with table of terms beneath it
#  # From gprofiler2: This function allows to highlight a list of selected terms on the Manhattan plot created with the gprofiler2::gostplot() function. The resulting plot is saved to a publication ready image if 'filename' is specified. The plot is very similar to the one shown in the g:GOSt web tool after clicking on circles.
#  gostplot_mh_plot_with_table<-publish_gostplot(p, highlight_terms = gprofiler_results$result$term_id,
#                         width = NA, height = NA, filename = NULL )
#  print(gostplot_mh_plot_with_table)
#  
#  # Assign gprofiler results
#  result_df <- gprofiler_results$result
#  
#  # Filter for significant results
#  significant_results <- result_df[result_df$significant, ]
#  
#  # Create a new column for -log10(p-value)
#  significant_results$log_p_value <- -log10(significant_results$p_value)
#  
#  # Print significant gprofiler results:
#  print(significant_results)
#  
#  # Plot significant gProfiler results
#  library(ggplot2)
#  ggplot(significant_results, aes(x = reorder(term_name, log_p_value), y = log_p_value)) +
#    geom_bar(stat = "identity") +
#    theme_minimal() +
#    coord_flip() +
#    labs(x = "Enriched Terms", y = "-log10(p-value)", title = "Enrichment Analysis Results from gProfiler") +
#    theme(axis.text.x = element_text(angle = 45, hjust = 1))

