# Rf2pval

Welcome to Rf2pval a comprehensive tool designed to revolutionize your approach to genomic data analysis using Random Forest Models in R. Tailored for expression data, such as RNA-seq or Microarray, Rf2pval is built for bioinformaticians and researchers looking to explore the relationship between biological features and a matched binary outcome variable using Random Forest models. 

Please see our [vignette](https://github.com/tkolisnik/Rf2pval/blob/main/vignettes/Rf2pval-vignette.pdf) for instructions that will guide you through Rf2pval's seamless integration of [scikit-learn's](https://scikit-learn.org/stable/common_pitfalls.html) Random Forest methodologies (imported to R via [reticulate](https://rstudio.github.io/reticulate/)) for model development, evaluation, and our custom feature reduction approach by way of rank-based permutation. You will also be directed through our integration with [Enrichr Enrichment Analysis & Gene Ontology](https://maayanlab.cloud/Enrichr/), [SHAP](https://shap.readthedocs.io/en/latest/) and [gProfiler](https://biit.cs.ut.ee/gprofiler/gost).

## Features

- Integration of Python's scikit-learn Random Forest models in R.
- Custom rank-based feature reduction methodology.
- Compatibility with RNA-seq and Microarray data.
- Integration with Gene Ontology, Enrichment Analysis, SHAP, and gProfiler.

## Installation

```r
# Install devtools if not already available

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install Rf2pval from GitHub

devtools::install_github("tkolisnik/Rf2pval")

# Load Rf2pval

library(Rf2pval)

Note: Rf2pval is developed on R version 4.3.1 for Apple Mac M1 arm64 architecture.

Python and Conda Setup

Before using Rf2pval, ensure Python 3.9.18 and Conda are installed and
set up in a Conda environment with the necessary packages:

# Install and load reticulate

if (!requireNamespace("reticulate", quietly = TRUE))
    install.packages("reticulate")
library(reticulate)

# Setup Conda environment

reticulate::use_python("/usr/local/bin/python3")
reticulate::conda_create(envname = "rf2pval-conda")
reticulate::conda_install(envname = "rf2pval-conda", packages = c("scikit-learn", "numpy", "shap"))
reticulate::use_condaenv("rf2pval-conda", required = TRUE)

# Verify Python configuration

reticulate::py_config()
```

## Data Preparation

Rf2pval works with structured genomic data. Your data should be formatted as a list containing four tibbles (training, validation, testing, and target categories). See the [vignette](https://github.com/tkolisnik/Rf2pval/blob/main/vignettes/Rf2pval-vignette.pdf) for detailed structuring instructions.

[Demo data](https://github.com/tkolisnik/Rf2pval/blob/main/data/demo_rnaseq_data.RData) is available using the command <code>data("demo_rnaseq_data", package = "Rf2pval")</code>

## Usage

Rf2pval includes functionalities for data preprocessing, model tuning, fitting, evaluation, and obtaining feature importances. Detailed usage examples and function documentation are available in the package [vignette](https://github.com/tkolisnik/Rf2pval/blob/main/vignettes/Rf2pval-vignette.pdf).

## Advanced Features

Explore advanced functionalities including SHAP value analysis and integration with Enrichment, Gene Ontology and gProfiler modules.

## Support and Contribution

For more information, detailed documentation, and how to contribute, visit the Rf2pval GitHub repository.

## Acknowledgments

Created by Tyler Kolisnik with support from Dr. Olin Silander, Dr. Adam Smith & Faeze Keshavarz.

## Contact

For queries, contact Tyler Kolisnik at tkolisnik@gmail.com.
