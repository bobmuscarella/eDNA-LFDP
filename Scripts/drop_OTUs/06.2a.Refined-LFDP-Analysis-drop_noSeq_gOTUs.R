########################################################
### MAIN 40-POINT ANALYSES
########################################################

### Load libraries
library(vegan)
library(caret)
library(compositions)
library(viridis)
library(ecodist)
library(sf)
library(spdep)
library(parallel)

### Make labels for saving plots
label <- c("lenient_10k",
           "lenient_200k",
           "lenient_40k",
           "lenient_70k",
           "lenientf1",
           "rare_lenient_10k",
           "rare_lenient_200k",
           "rare_lenient_40k",
           "rare_lenient_70k",
           "rare_repfiltered",
           "repfiltered",
           "repfilteredf1")

### Name data files
outfiles <- c( "stem-soil-40pt-data-lenient_10k-drop_noSeq_gOTUs-20250406.RDA",
               "stem-soil-40pt-data-lenient_200k-drop_noSeq_gOTUs-20250406.RDA",
               "stem-soil-40pt-data-lenient_40k-drop_noSeq_gOTUs-20250406.RDA",
               "stem-soil-40pt-data-lenient_70k-drop_noSeq_gOTUs-20250406.RDA",
               "stem-soil-40pt-data-lenientf1-drop_noSeq_gOTUs-20250406.RDA",
               "stem-soil-40pt-data-rare_lenient_10k-drop_noSeq_gOTUs-20250406.RDA",
               "stem-soil-40pt-data-rare_lenient_200k-drop_noSeq_gOTUs-20250406.RDA",
               "stem-soil-40pt-data-rare_lenient_40k-drop_noSeq_gOTUs-20250406.RDA",
               "stem-soil-40pt-data-rare_lenient_70k-drop_noSeq_gOTUs-20250406.RDA",
               "stem-soil-40pt-data-rare_repfiltered-drop_noSeq_gOTUs-20250406.RDA",
               "stem-soil-40pt-data-repfiltered-drop_noSeq_gOTUs-20250406.RDA",
               "stem-soil-40pt-data-repfilteredf1-drop_noSeq_gOTUs-20250406.RDA")

### Select data file to load (for testing)
# data_selector <- 1

### Select data file to load (for full loop)
for(data_selector in seq_along(outfiles)){
  
  message(paste("Working on dataset", data_selector, ":", label[data_selector]))
  
  ### Load selected data file
  datafile <- paste0("Processed_data/", outfiles[data_selector])
  
  ### Load data
  dnamat <- readRDS(datafile)[[1]]
  stem.23 <- readRDS(datafile)[[2]]
  
  ### Count number of valid eDNA sample points (out of 40)
  (npts <- nrow(stem.23$abund[[1]]) - length(grep('random', rownames(stem.23$abund[[1]]))))
  
  ### Species code full plot summaries
  lfdp <- readRDS("Raw_data/LFDP2023-extract-v2-20240427.RDA")
  
  ### gOTU full plot summaries
  lfdp23 <- readRDS("Raw_data/LFDP2023-extract-v2-drop_noSeq_gOTUs-20240427-gOTUs.RDA")
  
  ### Make presence absence matrices
  dnamat.pa <- 1*(dnamat>0)
  stem.pa.list <- lapply(stem.23$abund, function(x) 1*(x>0))

  ##########################################################
  ######################################################
  ### 40 Point Analysis
  ########################################################
  ########################################################
  
  ########################################################
  ### 1. Correlation of species richness between DNA and stem data
  ########################################################
  corrs <- corrs_rare <- list()
  
  for(r in seq_along(stem.23$abund)){
    corrs_rare[[r]] <- cor.test(rarefy(stem.23$abund[[r]][1:npts,],
                                       min(rowSums(stem.23$abund[[r]][1:npts,]))),
                                rarefy(dnamat, min(rowSums(dnamat))))
    
    corrs[[r]] <- cor.test(rowSums(stem.23$abund[[r]][1:npts,]>0), rowSums(dnamat.pa))
  }

  ########################################################
  ### Figure 2
  ########################################################
  cp <- rev(viridis::viridis(20))
  cpsp <- viridis::viridis_pal(option = "A")(20)[cut(log10(lfdp23$total_abund), 20)]
  
  # pdf(paste0("Figures/Fig2.rank-abundance-plot_", 
  #            label[data_selector],
  #            "-drop_noSeq_gOTUs.pdf"), width = 6, height = 7.5)
  
  par(mfrow=c(3,2), mar=c(4,4,2,1))
  
  # TAXON ACCUMULATION IN STEM DATA WITH INCREASING RADII AROUND SAMPLE POINTS
  b <- boxplot(sapply(stem.23$abund, function(x) rowSums(x>0)[1:npts]), 
               col=cp,
               xlab="Radius (m)", 
               ylab="Taxon richness", 
               xlim=c(-0.5,20), ylim=c(0,no_gotu_stem))
  b2 <- boxplot(rowSums(dnamat.pa), add=T, at=-0.5, width=2, col=2)
  axis(1, at=-0.5, labels="DNA")
  abline(v=0.25)
  abline(h=no_gotu_stem, lty=2)
  mtext("A", adj=0, line=0.1)
  
  # Full plot abundance
  plot(rank(lfdp23$total_abund), rank(colSums(dnamat.pa)), 
       ylab="Rank abundance of DNA reads",
       xlab="Rank abundance of stems", 
       pch=21, cex=1.5, bg=cpsp)
  mtext("B", adj=0, line=0.1)
  
  # Rank stem abundance is significantly positively correlated with rank DNA read abundance
  (rankcor <- cor.test(rank(lfdp23$total_abund), rank(colSums(dnamat.pa))))
  rank_corr_est <- rankcor$estimate
  rank_corr_p <- rankcor$p.value
  
  plot(lfdp23$total_ba * 100,
       (colSums(dnamat.pa)/npts), log='x',
       ylab="Prop. of sites detected in DNA",
       xlab=bquote(log[10] ~ "Total basal area" ~ (m^2)),
       pch=21, bg=cpsp, ylim=c(0,1), cex=1.5)
  mtext("C", adj=0, line=0.1)
  
  legend_image <- as.raster(matrix(rev(viridis_pal(option = "A")(50)), ncol=1))
  rasterImage(legend_image, 0.0015, 0.5, 0.0025, 0.95)                       ## the gradient
  polygon(x=c(0.0015,0.0025,0.0025,0.0015), y=c(0.5,0.5,0.95,0.95))
  text(0.0025, 0.925, "High", adj=-0.1, cex=0.75)
  text(0.0025, 0.525, "Low", adj=-0.1, cex=0.75)
  text(0.0025, 0.725, "Abundance", adj=-0.1)
  
  # Total BA is significantly positively correlated with prop. of DNA samples where detected
  (bacor <- cor.test(lfdp23$total_ba * 100, colSums(dnamat.pa)/npts))
  ba_corr_est <- bacor$estimate
  ba_corr_p <- bacor$p.value

  ### Logistic regression with distance to nearest tree and detection
  y <- as.vector(dnamat.pa[1:npts,])
  x <- log10(as.vector(unlist(stem.23$nn[1:npts,])))
  
  plot(x, jitter(y, 0.1), pch=16, cex=0.75, col=rgb(0,0,0,0.05),
       ylab="Presence in DNA",
       xlab="Dist. to nearest tree (m)", #[log10]", 
       xlim=c(0,3), cex.lab=1.25, 
       axes=F)
  axis(1, c(0,1,2,3), labels=c(1, 10, 100,1000))
  axis(2)
  box()
  
  mod <- glm(as.vector(dnamat.pa) ~ x, family="binomial")
  log_p <- (summary(mod)$coefficients["x","Pr(>|z|)"])
  nd <- data.frame(x=seq(0, max(x), length.out=100))
  ypred <- predict(mod, nd, type='response')
  
  # Add overall fit
  lines(nd$x, ypred, col="blue", lwd=3)
  
  # Add taxon-specific fits
  coeffs <- vector()
  for(sp in 1:ncol(stem.23$nn)){
    xx <- log(unlist(stem.23$nn[1:npts,sp]))
    yy <- unlist(as.data.frame(1*(dnamat[,sp]>1)))
    mod <- glm(yy ~ xx, family=binomial)
    coeffs[sp] <- coef(mod)[2]
    nd <- data.frame(xx=seq(0, max(x), length.out=100))
    ypred <- predict(mod, nd, type='response')
    lines(nd$xx, ypred, 
          col=rgb(0,0,0,0.2),
          lwd=0.5)
    }
  mtext("D", adj=0, line=0.1)
  
  plot(seq_along(stem.23$abund), sapply(corrs_rare, function(x) x$estimate), 
       ylim=c(-0.5,0.75), axes=F,
       pch=21, bg='grey', xlab="Radius (m)",
       ylab="Pearson correlation",
  )
  segments(seq_along(stem.23$abund), y0=sapply(corrs_rare, function(x) x$conf.int[1]),
           y1=sapply(corrs_rare, function(x) x$conf.int[2]))
  abline(h=0, lty=2)
  points(seq_along(stem.23$abund), sapply(corrs_rare, function(x) x$estimate), 
         pch=21, bg=cp, cex=2)
  axis(1, labels=seq(5, 100, 5), at=1:20)
  axis(2)
  box()
  mtext("E", adj=0, line=0.1)
  
  dev.off()
  
  ########################################################
  ### SUMMARY DATA AT WHOLE PLOT LEVEL
  ########################################################
  ### Initiate a table to hold plot-level summary stats
  if(data_selector==1){
    plot_level_summary_table <- data.frame(dataset=NA,
                                           npts=NA,
                                           no_gotu_dna=NA,
                                           no_gotu_stem=NA,
                                           mdna=NA, sddna=NA, mindna=NA, maxdna=NA,
                                           m5=NA, sd5=NA, min5=NA, max5=NA,
                                           m100=NA, sd100=NA, min100=NA, max100=NA,
                                           rank_corr_est=NA, rank_corr_p=NA,
                                           ba_corr_est=NA, ba_corr_p=NA,
                                           log_p=NA)
  }
  
  # Total number of gOTUs in the DNA data
  (no_gotu_dna <- sum(colSums(dnamat.pa)>0))
  
  # Total gOTUs in the stem data:
  (no_gotu_stem <- sum(lfdp23$total_abund>0))
  
  ### Percent of gOTUs from stem data that are detected in DNA data
  round(100 * (no_gotu_dna/no_gotu_stem), 1)
  
  ### Mean and SD of taxon richness in eDNA data
  (mdna <- mean(rowSums(dnamat.pa)))
  (sddna <- sd(rowSums(dnamat.pa)))
  (mindna <- min(rowSums(dnamat.pa)))
  (maxdna <- max(rowSums(dnamat.pa)))
  
  ### Mean and SD of taxon richness in stem data at 5 m
  (m5 <- mean(rowSums(stem.23$abund[[1]]>0)[1:npts]))
  (sd5 <- sd(rowSums(stem.23$abund[[1]]>0)[1:npts]))
  (min5 <- min(rowSums(stem.23$abund[[1]]>0)[1:npts]))
  (max5 <- max(rowSums(stem.23$abund[[1]]>0)[1:npts]))
  
  ### Mean and SD of taxon richness in stem data at 100 m
  (m100 <- mean(rowSums(stem.23$abund[[20]]>0)[1:npts]))
  (sd100 <- sd(rowSums(stem.23$abund[[20]]>0)[1:npts]))
  (min100 <- min(rowSums(stem.23$abund[[20]]>0)[1:npts]))
  (max100 <- max(rowSums(stem.23$abund[[20]]>0)[1:npts]))
  
  plot_level_summary_table[data_selector,] <- c(label[data_selector],
                                                npts,
                                                no_gotu_dna,
                                                no_gotu_stem,
                                                mdna, sddna, mindna, maxdna,
                                                m5, sd5, min5, max5,
                                                m100, sd100, min100, max100,
                                                rank_corr_est, rank_corr_p,
                                                ba_corr_est, ba_corr_p,
                                                log_p)
  
  plot_level_summary_table[,-1] <- round(apply(plot_level_summary_table[,-1], 2, as.numeric), 3)
  
  if(data_selector==length(outfiles)){
    write.csv(plot_level_summary_table, 
              "Processed_data/plot_level_summary_table-drop_noSeq_gOTUs.csv", row.names = F)
  }
  
  
  ########################################################
  ### 2. Confusion matrix analyses
  ########################################################
  
  ### **Note that this takes while to run**
  message(paste("working on confusion matrix for", datafile, "; dataset", data_selector))
  
  # Set the number of cores to use for parallel processing
  num_cores <- parallel::detectCores()
  
  # Create clusters for parallel processing
  cl <- makeCluster(num_cores)
  
  # Initialize lists to store results
  conf_stats_obs_list <- vector("list", length = 20)
  conf_stats_ses_list <- vector("list", length = 20)
  
  nruns <- 999
  
  # Function to compute standardized effect size
  ses <- function(obs, rand){
    return((obs - colMeans(rand))/apply(rand, 2, sd))
  }
  
  # Export the caret package to all nodes
  clusterEvalQ(cl, library(caret))
  clusterEvalQ(cl, library(mltools))
  
  # Export necessary objects to all nodes
  clusterExport(cl, c("dnamat.pa", "stem.23", "nruns", "ses", "nruns", "npts"))
  
  # Parallel loop over radius
  results <- parLapply(cl, 1:20, function(r) {
    message(paste("radius =", r))
    
    # Initialize matrices to store observation and SES results
    conf_stats_obs <- matrix(nrow = npts, ncol = 12)
    conf_stats_ses <- matrix(nrow = npts, ncol = 12)
    
    # Perform computations for each site
    for (site in 1:npts) {
      confus_obs <- caret::confusionMatrix(table(+(!dnamat.pa[site,]), 
                                                 +(!(1*(stem.23$abund[[r]][site,]>0)))
      ))
      
      conf_stats_obs[site,1:11] <- confus_obs$byClass
      
      # Compute MCC
      conf_stats_obs[site,12] <- mltools::mcc(TP=confus_obs$table[2,2],
                                              FP=confus_obs$table[2,1],
                                              TN=confus_obs$table[1,1],
                                              FN=confus_obs$table[1,2])
      
      # Generate randomizations
      randomizations <- 1:nruns + npts
      
      confus_rands <- lapply(randomizations, function(s) {
        confusionMatrix(table(+(!dnamat.pa[site,]),
                              +(!(1 * (stem.23$abund[[r]][s,] > 0)))))
      })
      
      # Calculate SES for each class
      conf_stats_ses[site,1:11] <- ses(conf_stats_obs[site,1:11], 
                                       do.call(rbind, lapply(confus_rands, 
                                                             function(x) x$byClass)))
      
      # Compute MCC of randomized points
      rand_mcc <- unlist(lapply(confus_rands, function(tab) {
        mltools::mcc(confusionM = matrix(tab$table, nrow = 2))
      }))
      
      conf_stats_ses[site,12] <- (conf_stats_obs[site,12] - mean(rand_mcc))/sd(rand_mcc)
      
    }
    
    # Convert matrices to data frames
    conf_stats_obs <- as.data.frame(conf_stats_obs)
    conf_stats_ses <- as.data.frame(conf_stats_ses)
    names(conf_stats_ses) <- names(conf_stats_obs) <- c(names(confus_obs$byClass), "MCC")
    
    # Return results as a list
    list(conf_stats_obs = conf_stats_obs, conf_stats_ses = conf_stats_ses)
  })
  
  # Close clusters
  stopCluster(cl)
  
  # Unpack results
  for (i in 1:20) {
    conf_stats_obs_list[[i]] <- results[[i]]$conf_stats_obs
    conf_stats_ses_list[[i]] <- results[[i]]$conf_stats_ses
  }
  
  ### Median Sensitivity, Specificity & Balanced Accuracy at 5 m
  round(median(lapply(conf_stats_obs_list, function(x) x$Sensitivity)[[1]]), 2)
  round(median(lapply(conf_stats_obs_list, function(x) x$Specificity)[[1]]), 2)
  round(median(lapply(conf_stats_obs_list, function(x) x$MCC)[[1]]), 2)
  
  ### Median Sensitivity, Specificity & Balanced Accuracy at 100 m
  round(median(lapply(conf_stats_obs_list, function(x) x$Sensitivity)[[20]]), 2)
  round(median(lapply(conf_stats_obs_list, function(x) x$Specificity)[[20]]), 2)
  round(median(lapply(conf_stats_obs_list, function(x) x$`MCC`)[[20]]), 2)
  
  # Save confusion matrix
  saveRDS(list(conf_stats_obs_list=conf_stats_obs_list, 
               conf_stats_ses_list=conf_stats_ses_list), 
          file=paste0("Processed_data/Conf_matrix_output-", 
                      label[data_selector], 
                      "-drop_noSeq_gOTUs-20250406.RDA"))
  
}