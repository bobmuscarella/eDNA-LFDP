########################################################
### COMPUTE CONFUSION MATRIX INDICES WITH THE INTENSIVE SUB-SAMPLE PLOTS
########################################################

### THE GOTU COLLAPSED DATASET OF THE G-PLOTS ARE PREPARED IN SCRIPT 2b
library(phyloseq)

### Get intensive sample points phyloseq object (list of 12, for each bioinfo)
glist <- readRDS("Processed_data/PR_eDNA-for-analysis_intenplot_2025-04-01.RData")

for(gdat in seq_along(glist)){
  
  message(paste("working on bioinformatic filtering", names(glist)[gdat]), ": ", gdat)
  
  g <- readRDS(paste0("Processed_data/stem-soil-intensive-data-", 
                      names(glist)[gdat], 
                      "-20250406.RDA"))
   
  # Function to compute standardized effect size
  ses <- function(obs, rand){
    return((obs - colMeans(rand))/apply(rand, 2, sd))
  }
  
  # Initialize holders for results
  conf_stats_obs <- matrix(nrow = 20, ncol = 13)
  conf_stats_ses <- matrix(nrow = 20, ncol = 13)
  
  for(r in 1:20){
    message(paste("working on radius", r))
    confus_obs <- caret::confusionMatrix(table(+(!g$dnamat.pa.G),
                                               +(!(1 * (g$stem.otu$abund[[r]][1,] > 0)))))
   
    conf_stats_obs[r, 1:12] <- c(seq(5,100,5)[r], confus_obs$byClass)
    
    conf_stats_obs[r, 13] <- mltools::mcc(TP=confus_obs$table[2,2],
                                          FP=confus_obs$table[2,1],
                                          TN=confus_obs$table[1,1],
                                          FN=confus_obs$table[1,2])
    
    # Randomization loop
    for(rand in 1:999){
      
      if(rand==1) {confus_rand <- list()}
      
      confus_rand[[rand]] <- caret::confusionMatrix(table(+(!g$dnamat.pa.G),
                                                 +(!(1 * (g$stem.otu$abund[[r]][rand,] > 0)))))
    }
      # Calculate SES for each class
      conf_stats_ses[r, 2:12] <- ses(conf_stats_obs[r, 2:12], 
                                    do.call(rbind, lapply(confus_rand, function(x) x$byClass)))
      
      
      mcc_all <- unlist(lapply(confus_rand, function(x) {
        mltools::mcc(TP=x$table[2,2],
                     FP=x$table[2,1],
                     TN=x$table[1,1],
                     FN=x$table[1,2])}))
      
      conf_stats_ses[r, 13] <- (conf_stats_obs[r, 13] - mean(mcc_all))/sd(mcc_all)

    }

  conf_stats_ses[, 1] <- seq(5,100,5)
  
  colnames(conf_stats_obs) <- c("r", names(confus_obs$byClass), "MCC")
  colnames(conf_stats_ses) <- c("r", names(confus_obs$byClass), "MCC")
  
  conf_stats_obs <- data.frame(conf_stats_obs)
  conf_stats_ses <- data.frame(conf_stats_ses)
  
  ### Save the data
  saveRDS(list(conf_stats_obs=conf_stats_obs, 
               conf_stats_ses=conf_stats_ses), 
          file=paste0("Processed_data/Conf_matrix_output-Gplots-", 
                 names(glist)[gdat], 
                 "-20250406.RDA"))
}
