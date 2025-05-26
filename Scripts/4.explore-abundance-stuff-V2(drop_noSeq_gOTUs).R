########################################################
### Explore false presences and abundance
### This version drops species without sequence data
########################################################

library(viridis)

### Color palette for plotting
cp <- rev(viridis::viridis(20))

### Make labels for saving plots
labels <- c("lenient_10k",
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

### Select data file to load (for testing)
# data_selector <- 1

### Select data file to load (for full loop)
for(data_selector in seq_along(labels)){
  
  file <- paste0("Processed_data/stem-soil-40pt-data-", 
                 labels[data_selector], "-drop_noSeq_gOTUs-20250406.RDA")

  ### Load data
  dnamat <- readRDS(file)[[1]]
  stem.23 <- readRDS(file)[[2]]
  
  ### Count number of valid eDNA sample points (out of 40)
  (npts <- nrow(stem.23$abund[[1]]) - length(grep('random', rownames(stem.23$abund[[1]]))))
  
  ### gOTU full plot summaries
  lfdp23 <- readRDS("Raw_data/LFDP2023-extract-v2-20240427-gOTUs.RDA")
  
  ### Make presence absence matrices
  dnamat.pa <- 1*(dnamat>0)
  stem.pa.list <- lapply(stem.23$abund, function(x) 1*(x>0))
  
  ### Make a matrix of false positives at 5 m
  fp <- 1 * (dnamat.pa==1 & stem.pa.list$`5`[1:npts,]==0)
  
  ### Make a matrix of false absences at 5 m
  fa <- 1 * (dnamat.pa==0 & stem.pa.list$`5`[1:npts,]==1)
  
  ### Make a matrix of true presences at 5 m
  tp <- 1 * (dnamat.pa==1 & stem.pa.list$`5`[1:npts,]==1)
  
  ### Make a matrix of false absences at 5 m
  ta <- 1 * (dnamat.pa==0 & stem.pa.list$`5`[1:npts,]==0)
  
  ### Get proportion of sites where a gOTU was a false positive at 5 m
  fp_prop <- colSums(fp)/npts
  fa_prop <- colSums(fa)/npts
  tp_prop <- colSums(tp)/npts
  ta_prop <- colSums(ta)/npts
  
### Plot!

  pdf(paste0("Figures/compare_filtering/Fig4.Abundance-conf_mat_",
             labels[data_selector],
             "-drop_noSeq_gOTUs.pdf"), width = 6, height = 7)
  
  par(mfcol=c(2,2), mar=c(4,4,2,1))
  
  cp <- viridis::viridis_pal(option = "A")(20)[cut(log10((16*lfdp23$total_ba)+0.0001), 20)]
  
  plot(lfdp23$total_abund, tp_prop, log='x',
       xlab="Total Abundance", ylab="Prop sites where TP @5m",
       pch=21, bg=cp, ylim=c(0,1), cex=1.5)
  mtext("A", adj=0, line=0.5)
  
  legend_image <- as.raster(matrix(rev(viridis_pal(option = "A")(50)), ncol=1))
  rasterImage(legend_image, 1, 0.5, 2, 0.95)                       ## the gradient
  polygon(x=c(1,2,2,1), y=c(0.5,0.5,0.95,0.95))
  text(2, 0.925, "High BA", adj=-0.1, cex=0.75)
  text(2, 0.525, "Low BA", adj=-0.1, cex=0.75)
  
  plot(lfdp23$total_abund, ta_prop, log='x',
       xlab="Total Abundance", ylab="Prop sites where TA @5m",
       pch=21, bg=cp, ylim=c(0,1), cex=1.5)
  mtext("B", adj=0, line=0.5)
  
  plot(lfdp23$total_abund, fp_prop, log='x',
       xlab="Total Abundance", ylab="Prop sites where FP @ 5m",
       pch=21, bg=cp, ylim=c(0,1), cex=1.5)
  mtext("C", adj=0, line=0.5)
  
  plot(lfdp23$total_abund, fa_prop, log='x',
       xlab="Total Abundance", ylab="Prop sites where FA @5m",
       pch=21, bg=cp, ylim=c(0,1), cex=1.5)
  mtext("D", adj=0, line=0.5)
  
  dev.off()
  
}


