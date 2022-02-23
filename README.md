# Rf2pval
# Created by: Tyler Kolisnik under advisement of Dr. Adam Smith and Dr. Olin Silander of Massey University, Auckland, New Zealand
# Contact: tkolisnik@gmail.com

Description:
An R package for obtaining p-values for features and feature sets as well as suggested feature cutoff points from the resultant feature importance scores of random forest models. 
This package is compatible with any random forest model that produces gini impurity based feature importance scores.
The method used to calculate p-values in this package is rank and permutation based. Before running this package you must permute your random forest model a recommended minimum of 100 times by randomizing the outcome variables, and obtain feature importance scores for each permutation. This provides us with data that can be used to generate a null-distribution in from which to calculate rank-based p-values with. 
This package contains a series of 7 functions that are meant to be run in sequential order.

Notes:
The input files must be in a particular format and must be .csv

Input File Format:
There must be a column for feature importance and a column for feature rank, these can be named whatever you like as column name is specified in the data import functions.
Feature Importances must be ordered highest to lowest. DO NOT ORDER THEM BY FEATURE NAME, this is a rank based test.
Feature Ranks must start at 1.
Instead of using a column for feature rank you can alternatively use rownames by indicating "rownames" as the column name in the function input.
Additional (extra) column names in the input files will be kept and remain mapped throughout, such as feature name, gene symbol, ensemble id but are not necessary for calculations. 

Sample Input for load_true_fi():
  my-ml-model-true-results.csv
      rownames  GENEID     ENSEMBLEID      featureImportance        featureRank 
          <chr>        <chr>                       <dbl>              <dbl>                
      1   PDE3A        ENSG00000172572            0.0210                 1                
      2   MCMDC2       ENSG00000178460            0.0189                 2                

Sample Input for load_permutated_fi():
  ~desktop/my-ml-permutations-folder/
      /my-ml-model-permutated-results-1.csv
      rownames  GENEID     ENSEMBLEID      featureImportance        featureRank 
          <chr>        <chr>                       <dbl>              <dbl>                
      1   KLRG2        ENSG00000188883            0.0015                 1                
      2   FRS2         ENSG00000166225            0.0089                 2  
      /my-ml-model-permutated-results-2.csv
      rownames  GENEID     ENSEMBLEID      featureImportance        featureRank 
          <chr>        <chr>                       <dbl>              <dbl>                
      1   CD2AP       ENSG00000198087            0.0005                  1                
      2   MZB1        ENSG00000170476            0.0009                  2  
      



This package contains a series of 7 functions that are meant to be run in the following order:

##Load the true data for your random forest model. This will be assigned as "permutation 0" (although it is not a permuted dataset).
## This is a path to a file.
t_val<-load_true_fi("path/to/file","FeatureImportancesColname","FeatureRankColName") #path to a single file for the true feature importance scores

##Load Permuted Data:
## This is a path to a folder containing all the 100+ (recommended) permuted files. **Do not include the true data for your model in this folder!**
perm_val<-load_permuted_fi("path/to/folder","FeatureImportancesColname","FeatureRankColName") #path to a folder of files with feature importance scores

##Creates a tibble of quantiles.
q_val<-calculateQuantiles(t_val,perm_val)

##Calculate overall p-value for entire feature set.
p_val_set<-calculatePvalueforSet(perm_val,q_val)

##Generate a histogram showing the null (permuted) distribution of sum of absolute deviations compared to the sum of absolute deviations for the true set. Illustrates the pi statistics used to calculate p-values for the entire dataset used in the calculatePvalueforSet() function.
piHistogram(perm_val,q_val)

##Calculate p-value for each individual rank. You will likely want to export this output. Recommended cutoff is the first "x" features with a p-value less than 0.05. 
p_val_rank<-calculatePvalueforRank(t_val,perm_val,q_val)

##Generate a feature importance score-rank plot showing the  mean feature importance scores of each rank in the observed set compared to the mean feature importance score of each rank in the permuted set.

#with Log disabled
fiRankPlot(permutedvalues=perm_val,quantiledata=q_val,xlimitmin=0,xlimitmax=241,ylimitmin=0.0001,ylimitmax=0.0225,labelverticaladjust=1.3,labelhorizontaladjust=-0.02,indvPermScoresOn=FALSE,logOn=FALSE)

#with Log enabled (recommended)
fiRankPlot(permutedvalues=perm_val,quantiledata=q_val,xlimitmin=0,xlimitmax=250,ylimitmin=-8,ylimitmax=-3.5,labelverticaladjust=1.3,labelhorizontaladjust=-0.02,indvPermScoresOn=FALSE, logOn=TRUE)
