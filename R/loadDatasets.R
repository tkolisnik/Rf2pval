
#' Load feature importance scores for true data
#'
#' This function loads the first random forest results file
#'
#' @param pathtofile Path to the input file. Input file must be in .csv format as a table of feature importance scores or Gini Importance scores for each named feature. Must have a minimum of two columns, one named featurename and another named featureImportance.
#' @param featureImportanceColumnName The name of the column in the input .csv file which contains the feature importance scores.
#' @param featureRankColumnName The name of the column in the input .csv file which contains the ordered feature ranks. Option: enter "rownames" if the rank is stored as the rownames. Feature ranks should start at 1.
#' @return  A tibble containing all of the feature importance score data for the true (observed) dataset.
#' @export
load_true_fi <- function(pathtofile,featureImportanceColumnName,featureRankColumnName){
  if(grepl("\\.csv$", pathtofile)){
  d0 <- readr::read_csv(pathtofile,show_col_types=FALSE)
  if(featureRankColumnName=="rownames"){
  d2 <- d0 %>% dplyr::mutate(permutation = 0) %>% mutate(featureRank=as.numeric(rownames(d0))) %>% rename(featureImportance= all_of(featureImportanceColumnName))
  }else{
  d2 <- d0 %>% dplyr::mutate(permutation = 0) %>% mutate(featureRank=as.numeric(featureRankColumnName)) %>% rename(featureImportance= all_of(featureImportanceColumnName))
  }
  d3 <- d2 %>% mutate("LogfeatureImportance" = log(featureImportance))
  d3$featureRank<-d3$featureRank-1 #necessary to start feature ranking at 0 for plots.
  if(any(is.na(d3$featureRank))){
    stop("NA's detected in feature rank column ; Please check you used the right column name and try again")
  } else {
  return(d3) }
  }
  else{
    stop("Input file type is not .CSV ; Please check that it is .csv and try again")
  }
}



#' Load feature importance scores for permuted data
#'
#' This function loads the first random forest results files for the 100 randomized outcome variable null permutations
#'
#' @param pathtofile Path to a folder containing minimum 100 csv files (one for each null permutation). Input files must be in .csv format as a table of feature importance scores or Gini Importance scores for each named feature. Must have a minimum of two columns, one named featurename and another named featureImportance. The true non-randomized data .csv file must NOT be in this folder.
#' @param featureImportanceColumnName The name of the column in the input .csv file which contains the feature importance scores.
#' @param featureRankColumnName The name of the column in the input .csv file which contains the ordered feature ranks. Option: enter "rownames" if the rank is stored as the rownames. Feature ranks should start at 1.
#' @return A tibble containing all of the feature importance score data for the permuted datasets.
#' @export
load_permuted_fi <- function(pathtofile,featureImportanceColumnName,featureRankColumnName){
  if(length(list.files(path=pathtofile,pattern=".csv"))>0){
  allfiles<-list.files(path=pathtofile,pattern=".csv",full.names=TRUE)
  dd0 = readr::read_csv(allfiles[1],show_col_types = FALSE)
  if(featureRankColumnName=="rownames"){
  dd1 = dd0 %>% mutate(permutation = 1) %>% mutate(featureRank=as.numeric(rownames(dd0)))
  }else{
    dd1 = dd0 %>% mutate(permutation = 1) %>% mutate(featureRank=as.numeric(featureRankColumnName))
  }
  z<-1
  for (i in allfiles[2:length(allfiles)]){
    dd1 = dd1
    if(featureRankColumnName=="rownames"){
    dd2 = readr::read_csv(i,show_col_types = FALSE) %>% mutate(featureRank=as.numeric(rownames(dd0))) %>% dplyr::mutate(permutation = z)
    } else {
    dd2 = readr::read_csv(i,show_col_types = FALSE) %>% mutate(featureRank=as.numeric(featureRankColumnName)) %>% dplyr::mutate(permutation = z)
    }
    dd1 = dd1 %>% bind_rows(dd2)
    z<-z+1
  }
  dd1 = dd1 %>% rename(featureImportance= all_of(featureImportanceColumnName))
  dd1 = dd1 %>% dplyr::mutate(LogfeatureImportance = log(featureImportance))
  dd1$featureRank<-dd1$featureRank-1 #necessary to start feature ranking at 0 for plots.
  if(any(is.na(dd1$featureRank))){
    stop("NA's detected in feature rank column ; Please check you used the right column name and try again")
  } else {
  return(dd1) }
  }
  else{
    stop("Input file list is either empty or is not .CSV ; Please check and try again")
  }
}

