##########################################
### SIMULATION TO VALIDATE NULL MODEL APPROACH
##########################################


#####################
### Load data
###### For some simulations, bioinformatic method doesn't matter since we simulate eDNA
###### For simulations based on empirical eDNA data, choose bioinformatic method below
#####################

# Load stem data
stem.23 <- readRDS("Processed_data/stem-soil-40pt-data-lenient_10k-20250405.RDA")[[2]]

# Get number of sample points to use as "observed", leaving 999 "random" points
nrand <- 999
npts <- sum(!grepl('random', rownames(stem.23$abund$`5`)))

# Convert stem data to presence absence
samp.abun.pa <- lapply(stem.23$abund, function(x){ 1 * (x[1:npts,] > 0)})

### gOTU full plot summaries
lfdp23 <- readRDS("Raw_data/LFDP2023-extract-v2-20240427-gOTUs.RDA")
lfdp23 <- lfdp23[order(lfdp23$total_abund, decreasing=F),]
lfdp23$cumprop <- round(cumsum(lfdp23$total_abund)/sum(lfdp23$total_abund), 3)
lfdp23$cumprop_ba <- round(cumsum(lfdp23$total_ba)/sum(lfdp23$total_ba), 3)

#####################
### Helper function
#####################
# Function to compute standardized effect size
ses <- function(obs, rand){
  return((obs - colMeans(rand))/apply(rand, 2, sd))
}

#####################
### Find OTUS to drop
#####################

### Don't drop anything (Sim 1a and 1b)
# drops <- NULL

### Drop the 5 most common gOTUs (Sim 2 - Dominant taxa removal)
# drops <- rownames(lfdp23)[order(lfdp23$total_abund, decreasing=T)][1:5]

### Drop OTUs that make up lower 10% of total basal area (Sim 3 - Rare taxa removal)
lfdp23 <- lfdp23[order(lfdp23$total_ba, decreasing=F),]
lfdp23$cumprop_ba <- round(cumsum(lfdp23$total_ba)/sum(lfdp23$total_ba), 3)
drops <- rownames(lfdp23)[lfdp23$cumprop_ba <= 0.10]

### What is % total basal area of dropped taxa?
round(sum(lfdp23$total_ba[rownames(lfdp23) %in% drops]) /  sum(lfdp23$total_ba), 2)
length(drops)


#####################
### Setup eDNA community
#####################

### Choose the starting DNA community (e.g., 100% match with the 5 m census)
samp.dna.pa <- samp.abun.pa$`25`

### Impose undetected species (drops from above)
samp.dna.pa[,drops] <- samp.dna.pa[,drops] * 0

#####################
### Run the simulation
#####################

# Initialize results holders
conf_stats_obs_list <- list()
conf_stats_ses_list <- list()

# Loop through radii
for (r in 1:20) {
  message(paste('radius', r))
  # Subset the stem and dna communities to the focal radius
  tmp.samp.abun.pa <- samp.abun.pa[[r]]
  
  # Initialize matrices to store observation and SES results
  conf_stats_obs <- matrix(nrow = npts, ncol = 12)
  conf_stats_ses <- matrix(nrow = npts, ncol = 12)
  
  # Loop through sites
  for (site in 1:npts) {
    message(paste('-site', site))
    # Initialize a results holder
    if(site == 1) { confus_rand <- confus_obs <- list() }

    confus_obs[[site]] <- caret::confusionMatrix(table(+(!samp.dna.pa[site,]),
                                                 +(!tmp.samp.abun.pa[site,])))

    # Loop through randomizations
    for (rand in 1:nrand){
      confus_rand[[rand]] <- caret::confusionMatrix(table(+(!samp.dna.pa[site,]),
                                                          +(!stem.23$abund[[r]][npts + rand,])))
      }
    
    # Calculate SES
    conf_stats_ses[site,1:11] <- ses(confus_obs[[site]]$byClass,
                                     do.call(rbind, lapply(confus_rand, function(x) x$byClass)))
    
    obs_mcc <- mltools::mcc(confusionM = matrix(confus_obs[[site]]$table, nrow=2))
    
    rand_mmc <- do.call(c, lapply(confus_rand, function(x)
      mltools::mcc(confusionM = matrix(x$table, nrow=2)
      )))
    
    conf_stats_ses[site,12] <- (obs_mcc - mean(rand_mmc))/(sd(rand_mmc))

  }

  conf_stats_obs[,1:11] <- do.call(rbind, lapply(confus_obs, function(x) x$byClass))
  
  conf_stats_obs[,12] <- do.call(rbind, lapply(confus_obs, function(x){
    mltools::mcc(confusionM = matrix(x$table, nrow=2))
  }))
  
  
  colnames(conf_stats_ses) <- c(names(confus_obs[[1]]$byClass), "MCC")
  conf_stats_ses <- as.data.frame(conf_stats_ses)

  colnames(conf_stats_obs) <- c(names(confus_obs[[1]]$byClass), "MCC")
  conf_stats_obs <- as.data.frame(conf_stats_obs)
  
  conf_stats_obs_list[[r]] <- conf_stats_obs
  conf_stats_ses_list[[r]] <- conf_stats_ses
  
}


#####################
### Save output
#####################

saveRDS(list(conf_stats_obs_list=conf_stats_obs_list, 
             conf_stats_ses_list=conf_stats_ses_list), 
        file="Processed_data/Simulation/Sim3-RareTaxaRemoval-20250609.RDA")



#####################
### Plot it!
#####################

files <- c("Processed_data/Simulation/Sim1a-5mRef-20250609.RDA",
           "Processed_data/Simulation/Sim1b-25mRef-20250609.RDA",
           "Processed_data/Simulation/Sim2-DomTaxaRemoval-20250609.RDA",
           "Processed_data/Simulation/Sim3-RareTaxaRemoval-20250609.RDA")

for(f in seq_along(files)){
  
  conf_stats_obs_list <- readRDS(files[f])[[1]]
  conf_stats_ses_list <- readRDS(files[f])[[2]]
  
  pdf(paste0("Figures/FigureS", 
             c(10:13)[f], ".",
             tools::file_path_sans_ext(basename(files[f])),
             ".pdf"), width = 8, height = 10)
  
  par(mfcol=c(3,2), mar=c(4,4,1,1), oma=c(1,1,1,1))
  
  cols <- rev(viridis::viridis(length(conf_stats_obs_list)))
  
  focdist <- c(1,5,5,5)[f]
  
  # "Observed" (raw simulated) values
  boxplot(lapply(conf_stats_obs_list, function(x) x$Sensitivity),
          ylab="Sensitivity", 
          xlab="Radius (m)", axes=F, col=cols)
  axis(1, labels=names(stem.23$abund), at=1:20); axis(2); graphics::box()
  mtext("A", adj=0, line=0.5)
  abline(v=focdist, col='red', lwd=2, lty=2)
  
  boxplot(lapply(conf_stats_obs_list, function(x) x$Specificity), 
          ylab="Specificity", 
          xlab="Radius (m)", axes=F, col=cols)
  axis(1, labels=names(stem.23$abund), at=1:20); axis(2); graphics::box()
  mtext("B", adj=0, line=0.5)
  abline(v=focdist, col='red', lwd=2, lty=2)
  
  boxplot(lapply(conf_stats_obs_list, function(x) x$MCC),
          ylab="MCC", 
          xlab="Radius (m)", axes=F, col=cols)
  axis(1, labels=names(stem.23$abund), at=1:20); axis(2); graphics::box()
  mtext("C", adj=0, line=0.5)
  abline(v=focdist, col='red', lwd=2, lty=2)
  
  # SES values
  boxplot(lapply(conf_stats_ses_list, function(x) x$Sensitivity), 
          ylab="Sensitivity (SES)", 
          xlab="Radius (m)", axes=F, col=cols)
  polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
  boxplot(lapply(conf_stats_ses_list, function(x) x$Sensitivity), 
          axes=F, col=cols, add=T)
  axis(1, labels=names(stem.23$abund), at=1:20); axis(2); graphics::box()
  mtext("D", adj=0, line=0.5)
  abline(v=focdist, col='red', lwd=2, lty=2)
  
  boxplot(lapply(conf_stats_ses_list, function(x) x$Specificity), 
          ylab="Specificity (SES)", 
          xlab="Radius (m)", axes=F, col=cols)
  polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
  boxplot(lapply(conf_stats_ses_list, function(x) x$Specificity), 
          axes=F, col=cols, add=T)
  axis(1, labels=names(stem.23$abund), at=1:20); axis(2); graphics::box()
  mtext("E", adj=0, line=0.5)
  abline(v=focdist, col='red', lwd=2, lty=2)
  
  boxplot(lapply(conf_stats_ses_list, function(x) x$MCC), 
          ylab="MCC (SES)", 
          xlab="Radius (m)", axes=F, col=cols)
  polygon(x=c(-1,200,200,-1), y=c(-1.96, -1.96, 1.96, 1.96), lty=0, col='grey')
  boxplot(lapply(conf_stats_ses_list, function(x) x$MCC), 
          axes=F, col=cols, add=T)
  axis(1, labels=names(stem.23$abund), at=1:20); axis(2); graphics::box()
  mtext("F", adj=0, line=0.5)
  abline(v=focdist, col='red', lwd=2, lty=2)
  
  dev.off()
  
}
