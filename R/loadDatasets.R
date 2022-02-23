
#' Load feature weight for true data
#'
#' This function loads the first random forest results file
#'
#' @param pathtofile Path to the input file. Input file must be in .csv format as a table of feature weights or Gini Importance scores for each named feature. Must have a minimum of two columns, one named featurename and another named featureweight.
#' @return  A tibble containing all of the feature weight data for the true (observed) dataset.
#' @export
load_true_fws <- function(pathtofile,featureWeightColumnName){
  if(grepl("\\.csv$", pathtofile)){
  d0 <- readr::read_csv(pathtofile,show_col_types=FALSE)
  d2 <- d0 %>% dplyr::mutate(permutation = 0) %>% mutate(featureRank=as.numeric(rownames(d0))) %>% rename(featureWeight= all_of(featureWeightColumnName))
  d3 <- d2 %>% mutate("LogFeatureWeight" = log(featureWeight))
  return(d3)
  }
  else{
    stop("Input file type is not .CSV ; Please check that it is .csv and try again")
  }
}



#' Load feature weights for permuted data
#'
#' This function loads the first random forest results files for the 100 randomized outcome variable null permutations
#'
#' @param pathtofile Path to a folder containing minimum 100 csv files (one for each null permutation). Input files must be in .csv format as a table of feature weights or Gini Importance scores for each named feature. Must have a minimum of two columns, one named featurename and another named featureweight. The true non-randomized data .csv file must NOT be in this folder.
#' @return A tibble containing all of the feature weight data for the permuted datasets.
#' @export
load_permuted_fws <- function(pathtofile,featureWeightColumnName){
  if(length(list.files(path=pathtofile,pattern=".csv"))>0){
  allfiles<-list.files(path=pathtofile,pattern=".csv",full.names=TRUE)
  dd0 = readr::read_csv(allfiles[1],show_col_types = FALSE)
  dd1 = dd0 %>% mutate(permutation = 1) %>% mutate(featureRank=as.numeric(rownames(dd0)))
  z<-1
  for (i in allfiles[2:length(allfiles)]){
    dd1 = dd1
    dd2 = readr::read_csv(i,show_col_types = FALSE) %>% mutate(featureRank=as.numeric(rownames(dd0))) %>% dplyr::mutate(permutation = z)
    dd1 = dd1 %>% bind_rows(dd2)
    z<-z+1
  }
  dd1 = dd1 %>% rename(featureWeight= all_of(featureWeightColumnName))
  dd1 = dd1 %>% dplyr::mutate(LogfeatureWeight = log(featureWeight))
  return(dd1)
  }
  else{
    stop("Input file list is either empty or is not .CSV ; Please check and try again")
  }
}

