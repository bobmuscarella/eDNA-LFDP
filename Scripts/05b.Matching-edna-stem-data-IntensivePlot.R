########################################################
### PREPARE THE GOTU COLLAPSED DATASET OF THE INTENSIVE PLOT SAMPLES
########################################################

library(phyloseq)
library(EDIutils)

### Get intensive plot sample points phyloseq object (list of 12, for each bioinfo method)
glist <- readRDS("Processed_data/PR_eDNA-for-analysis_intenplot_2025-04-01.RData")

for(g in seq_along(glist)){
  
  message(paste("Now working on iteration", g, names(glist)[g]))
  
  gdat <- glist[[g]]
  
  ### Load the LFDP 2023 census extract data
  tree <- readRDS("Raw_data/LFDP2023-extract-20240510.RDA")
  codes <- readxl::read_xlsx("Raw_data/LFDP-SPcodes.xlsx")
  
  ### Load the sample point coordinates
  sample_xy <- read.csv("Raw_data/LFDP-eDNA-xy-V2.csv")
  
  ### Load link between OTUs and RefIDs
  otupotuR1 <- readRDS("Processed_data/OTU-to-RefIDs-List_R1.rds")[[6]]
  otupotuR2 <- readRDS("Processed_data/OTU-to-RefIDs-List_R2.rds")[[6]]
  
  ### Determine which otu.potu mapping to use
  if(!grepl("repfiltered", names(glist)[g])) {
    otu.potu.link <- otupotuR1
  } else {
    otu.potu.link <- otupotuR2
  }
  
  # Get coordinates from sample points included
  xy <- data.frame(X=140, Y=260)
  
  ### Make the rows of the stem data match the rows of the eDNA otu-table data
  ### Create stem data for 2023
  stem_abund <- lapply(tree$abund, function(x) {
    rbind(x[match(paste(xy$X, xy$Y), paste(sample_xy$PX, sample_xy$PY)),],
          x[88:1087,])})
  for(i in seq_along(stem_abund)){
    rownames(stem_abund[[i]]) <- c("g", paste0('random', 1:1000))
  }
  
  stem_ba <- lapply(tree$ba, function(x) {
    rbind(x[match(paste(xy$X, xy$Y), paste(sample_xy$PX, sample_xy$PY)),],
          x[88:1087,])})
  for(i in seq_along(stem_ba)){
    rownames(stem_ba[[i]]) <- c("g", paste0('random', 1:1000))
  }
  
  stem_nn <- rbind(tree$nearest_sp[match(paste(xy$X, xy$Y),
                                         paste(sample_xy$PX, sample_xy$PY)),],
                   tree$nearest_sp[88:1087,])
  rownames(stem_nn) <- c("g", paste0('random', 1:1000))
  
  stem <- list(abund=stem_abund, ba=stem_ba, nearest_sp=stem_nn)
  
  ##### THIS BLOCK COLLAPSES THE DNA DATA TO THE gOTU CLUSTERS
  ### Identify all POTUs/OTUs we need to keep because they are assigned to a gOTU
  potus.keep <- codes$POTU[!is.na(codes$gOTU)]
  otus.keep <- otu.potu.link$OTU[match(potus.keep, otu.potu.link$abundance)]
  gotus.keep <- codes$gOTU[!is.na(codes$gOTU)]
  codes.keep <- codes$`SPECIES CODE`[!is.na(codes$gOTU)]
  
  # POTUs that we keep and to be collapsed
  codes.collapse.list <- split(codes.keep, gotus.keep)
  otus.collapse.list <- split(otus.keep, gotus.keep)
  
  # loop to collapse taxa in phyloseq objects
  for(i in seq_along(otus.collapse.list)){
    if(length(unique(otus.collapse.list[[i]][!is.na(otus.collapse.list[[i]])])) > 1){
      taxcollapse <- unique(otus.collapse.list[[i]][!is.na(otus.collapse.list[[i]])])
      taxcollapse <- taxcollapse[taxcollapse %in% rownames(gdat@tax_table)]
      gdat <- merge_taxa(gdat, taxcollapse, 1)
    }
  }
  
  ##### THIS BLOCK COLLAPSES THE STEM DATA TO THE gOTU CLUSTERS
  stem.otu <- list()
  for(r in seq_along(stem$abund)){
    
    # Initialize a matrix to hold collapsed results for a given radius
    tmat_abund <- matrix(ncol=length(codes.collapse.list), nrow=nrow(xy)+1000)
    tmat_ba <- matrix(ncol=length(codes.collapse.list), nrow=nrow(xy)+1000)
    tmat_nn <- matrix(ncol=length(codes.collapse.list), nrow=nrow(xy)+1000)
    
    # loop to collapse taxa in phyloseq objects
    for(i in seq_along(codes.collapse.list)){
      tmp_abund <- stem$abund[[r]][,which(colnames(stem$abund[[1]]) %in% codes.collapse.list[[i]])]
      tmp_ba <- stem$ba[[r]][,which(colnames(stem$ba[[1]]) %in% codes.collapse.list[[i]])]
      tmp_nn <- stem$nearest_sp[,which(colnames(stem$nearest_sp) %in% codes.collapse.list[[i]])]
      
      if(is.matrix(tmp_abund)){
        tmat_abund[,i] <- rowSums(tmp_abund)
        tmat_ba[,i] <- rowSums(tmp_ba)
        tmat_nn[,i] <- rowSums(tmp_nn)
      } else {
        tmat_abund[,i] <- tmp_abund
        tmat_ba[,i] <- tmp_ba
        tmat_nn[,i] <- tmp_nn
      }
    }
    
    tmat_abund <- as.data.frame(tmat_abund)
    tmat_ba <- as.data.frame(tmat_ba)
    tmat_nn <- as.data.frame(tmat_nn)
    
    names(tmat_nn) <- names(tmat_ba) <- names(tmat_abund) <- names(codes.collapse.list)
    rownames(tmat_nn) <- rownames(tmat_ba) <- rownames(tmat_abund) <- rownames(stem$abund[[1]])
    
    stem.otu$abund[[r]] <- tmat_abund
    stem.otu$ba[[r]] <- tmat_ba
    stem.otu$nn <- tmat_nn
    
    names(stem.otu$ba)[r] <- names(stem.otu$abund)[r] <- names(stem$abund)[r]
  }
  
  
  ### To make the columns match, and add zero columns to OTU data
  
  # The new phyloseq object should have these OTUs as colnames
  collapsed.otus <- unlist(lapply(otus.collapse.list, function(x) x[!is.na(x)][1]))
  
  # Confirm the stem data has these names exactly
  if(!all(names(collapsed.otus) == names(stem.otu$abund[[1]]))){
    message("Warning: all(names(collapsed.otus) != names(stem.otu$abund[[1]]))")
  }
  
  # Some gOTU names are in the soil data but not the collapsed list because they aren't trees
  qs <- colnames(gdat@otu_table)[!colnames(gdat@otu_table) %in% collapsed.otus]
  otu.potu.link[otu.potu.link$OTU %in% qs,]
  
  # Drop these from the OTU Table
  gdat <- prune_taxa(colnames(gdat@otu_table)[colnames(gdat@otu_table) %in% collapsed.otus], gdat)
  
  # Some gOTU names are not in the soil data (not detected).  For these we add a zero column.
  (absent.gOTUs <- collapsed.otus[!collapsed.otus %in% colnames(gdat@otu_table)])
  
  # Make the tables comparable
  dnamat <- gdat@otu_table
  colnames(dnamat) <- names(collapsed.otus)[match(colnames(dnamat), collapsed.otus)]
  dnamat <- dnamat[,match(colnames(stem.otu$abund[[1]]), colnames(dnamat))]
  colnames(dnamat) <- colnames(stem.otu$abund[[1]])
  dnamat <- as.data.frame(dnamat)
  dnamat[is.na(dnamat)] <- 0
  
  #################################
  ### SAVE DATA FOR DOWNSTREAM ANALYSES
  #################################
  
  dnamat.pa <- 1 * (dnamat > 0)
  dnamat.pa.G <- 1 * (colSums(dnamat.pa) > 0)
  
  saveRDS(list(dnamat.pa.G=dnamat.pa.G,
               dnamat.pa=dnamat.pa,
               dnamat=dnamat,
               stem.otu=stem.otu),
          paste0("Processed_data/stem-soil-intensive-data-", names(glist)[g], "-20250406.RDA")
  )
}  
