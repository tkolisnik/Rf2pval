% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadDatasets.R
\name{load_true_fi}
\alias{load_true_fi}
\title{Load feature importance scores for true data
This function loads the first random forest results file}
\usage{
load_true_fi(pathtofile, featureImportanceColumnName, featureRankColumnName)
}
\arguments{
\item{pathtofile}{Path to the input file. Input file must be in .csv format as a table of feature importance scores or Gini Importance scores for each named feature. Must have a minimum of two columns, one named featurename and another named featureImportance.}

\item{featureImportanceColumnName}{The name of the column in the input .csv file which contains the feature importance scores.}

\item{featureRankColumnName}{The name of the column in the input .csv file which contains the ordered feature ranks. Option: enter "rownames" if the rank is stored as the rownames. Feature ranks should start at 1.}
}
\value{
A tibble containing all of the feature importance score data for the true (observed) dataset.
}
\description{
Load feature importance scores for true data
This function loads the first random forest results file
}
