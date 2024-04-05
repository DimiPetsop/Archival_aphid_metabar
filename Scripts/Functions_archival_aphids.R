# This function finds congruence with both definitions as percent and
# It makes certain assumptions for variable names and input dataframe sizes I do not have the time to fix these things I might at some point.
congruence_find <- function(df1,df2){
  df_common <-matrix(nrow=183, ncol=5)
  counter = 1
  
  for (i in unique(df1$date)){
    df_common[counter,] <- i
    df_common[counter,][2] <- length(intersect(df1[df1$date == i,]$genus,
                                               df2[df2$date == i,]$genus))
    df_common[counter,][3] <- length(setdiff(df1[df1$date == i,]$genus, 
                                             df2[df2$date == i,]$genus))
    df_common[counter,][4] <- length(setdiff(df2[df2$date == i,]$genus,
                                             df1[df1$date == i,]$genus))
    df_common[counter,][5] <- nrow(df2[df2$date == i,][2])
    counter = counter +1}

  df_common <- as.data.frame(df_common)
  colnames(df_common) <-c("Sample","Common","only_met","only_morph","total_morph")
  df_common[,-1] <- lapply(df_common[,-1],as.numeric)
  df_common <-transform(df_common,percent = 
                          round((100* df_common$Common)/df_common$total_morph))
  df_common$jaccard <- df_common$Common/((df_common$only_met+df_common$Common) +
                                           df_common$total_morph-df_common$Common)

return(df_common)
}


congruence_find_species <- function(df1,df2){
  df_common <-matrix(nrow=183, ncol=5)
  counter = 1
  
  for (i in unique(df1$date)){
    df_common[counter,] <- i
    df_common[counter,][2] <- length(intersect(df1[df1$date == i,]$species,
                                               df2[df2$date == i,]$species))
    df_common[counter,][3] <- length(setdiff(df1[df1$date == i,]$species, 
                                             df2[df2$date == i,]$species))
    df_common[counter,][4] <- length(setdiff(df2[df2$date == i,]$species,
                                             df1[df1$date == i,]$species))
    df_common[counter,][5] <- nrow(df2[df2$date == i,][2])
    counter = counter +1}
  
  df_common <- as.data.frame(df_common)
  colnames(df_common) <-c("Sample","Common","only_met","only_morph","total_morph")
  df_common[,-1] <- lapply(df_common[,-1],as.numeric)
  df_common <-transform(df_common,percent = 
                          round((100* df_common$Common)/df_common$total_morph))
  df_common$jaccard <- df_common$Common/((df_common$only_met+df_common$Common) +
                                           df_common$total_morph-df_common$Common)
  
  return(df_common)
}

calculate_sensitivity <- function(df){
  df <- df %>% mutate(sensitivity= df$Common/(df$Common + df$only_morph),
                      precision = df$Common/(df$Common + df$only_met))
  return(df)
}
