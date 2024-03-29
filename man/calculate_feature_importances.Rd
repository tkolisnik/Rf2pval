% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelTrainingTuningFittingTesting.R
\name{calculate_feature_importances}
\alias{calculate_feature_importances}
\title{Calculate True and Permuted Feature Importances}
\usage{
calculate_feature_importances(model, X_train, y_train, n_permutations = 1000)
}
\arguments{
\item{model}{The Random Forest model to be used.}

\item{X_train}{Training feature matrix.}

\item{y_train}{Training target vector.}

\item{n_permutations}{The number of permutations to perform (default: 1000).}
}
\value{
A list containing three data frames: one for true feature importances, one for permuted importances, and one containing the top features (filtered non-zero true importances).
}
\description{
This function fits a Random Forest model and calculates the true and permuted feature importances.
It performs permutations on the target variable to generate permuted importances for comparison.
}
\examples{
feat_importances <- Rf2pval::calculate_feature_importances(fitting_results$model,processed_training_data$X_training_mat,processed_training_data$y_training_vector,n_permutations=1000)

Print the fitting results, provides accuracy, f1 score, precision, recall and roc_auc scores on the model as fitted to the validation set
print(top_features)
}
