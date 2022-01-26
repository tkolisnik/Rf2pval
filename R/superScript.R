#This script shows an example of a possible workflow using this package
## double hashed denotes an alternate run

#Load Data:
#t_val<-load_true_fws("/Users/tyler/Desktop/2020 PAPER/14-04-2021 runs/plots/meanfeatureweights-resfile-side-permutated-0.csv","MeanFeatureWeight")
##t_val<-load_true_fws("/Users/tyler/Desktop/bestsinglecase-58-means-resfile-side-permutated-0.csv") #path to a single file for the true feature weights

#perm_val<-load_permutated_fws("/Users/tyler/Desktop/permutated_files","MeanFeatureWeight")
##perm_val<-load_permutated_fws("/Users/tyler/Desktop/SingleCase-Results-2","MeanFeatureWeight") #path to a folder of files with feature weights

#Create a tibble of quantiles.
#q_val<-calculateQuantiles(t_val,perm_val)

#Calculate overall p-value for entire feature set.
#p_val_set<-calculatePvalueforSet(perm_val,q_val)

#Calculate p-value for each individual rank.
#p_val_rank<-calculatePvalueforRank(t_val,perm_val,q_val)

#Generate a weight-rank plot showing the mean feature weights of each rank in the observed set compared to the mean feature weight of each rank in the permutated set.
#wrplot<-weightRankPlot(t_val,perm_val,q_val)

#Generate a log weight-rank plot showing the log mean feature weights of each rank in the observed set compared to the log mean feature weight of each rank in the permutated set.
#lwrplot<-logWeightRankPlot(t_val,perm_val,q_val)
##lwrplot<-logWeightRankPlot(t_val,perm_val,q_val,1,60,-7,-3,5,5)

#Generate a histogram showing the null (permutated) distribution of sum of absolute deviations compared to the sum of absolute deviations for the true set. Illustrates the pi statistics used to calculate p-values for the entire dataset.
#pi_hist<-piHistogram(perm_val,q_val)

