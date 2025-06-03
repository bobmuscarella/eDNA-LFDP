### Generating upset plots of shared taxa between pooled / unpooled big grid samples
### Heatmap plot of intensive grid sample, and rarefaction / hill numbers plots Figures S5 - S9

library(microViz)
library(dada2)
library(MiscMetabar)

subsample_etc<- readRDS("Processed_data/subsample-smallplot.RData")

subsample_etc_gotus<- readRDS("Processed_data/subsample-smallplot-gotus.RData")

## making sure Zingiberales order is swapped out with HELCAR - no banana recorded in plot but is easily misidentified with Heliconia caribaea
str(tax_table(subsample_etc$ssdat))

# Define the correct taxonomy for Heliconia caribaea
heliconia_tax <- c(
  phylum       = "Streptophyta",
  class        = "Liliopsida",
  order        = "Zingiberales",
  family       = "Heliconiaceae",
  genus        = "Heliconia",
  species      = "Heliconia caribaea",
  speciesOTU   = "Heliconia caribaea_1"
)

# Function to update taxonomy for gOTU76
update_gOTU76_tax <- function(physeq_obj) {
  tax <- tax_table(physeq_obj)
  if ("gOTU76" %in% rownames(tax)) {
    tax["gOTU76", ] <- heliconia_tax
    tax_table(physeq_obj) <- tax
  }
  return(physeq_obj)
}

# Apply the function to each element in your named list
subsample_etc<- lapply(subsample_etc, update_gOTU76_tax)

subsample_etc_gotus<- lapply(subsample_etc_gotus, update_gOTU76_tax)


subsample_etc$ssdat
subsample_etc$repfil_ssdat
#### Making some upset plots for sub-sample & mixed samples - c12
tab <- subsample_etc_gotus$ssdat
c12 <- subset_samples(tab, realsample1 == "C12")
c12 <- subset_samples(c12, DNAtreat == "clean")
c12<- phyloseq_validate(c12, remove_undetected = TRUE)
k1 <- upset_pq(c12, fact = "sample_from_sheet", name='sample/subsample', set_sizes=FALSE) + 
  ggtitle("All data, libs >10k, no other filtering, \nmixed sample 12 & subsamples") +
  theme(plot.title = element_text(size = 10), axis.title.x=element_blank())

tab <- subsample_etc_gotus$rare_ssdat
c12 <- subset_samples(tab, realsample1 == "C12")
c12 <- subset_samples(c12, DNAtreat == "clean")
c12<- phyloseq_validate(c12, remove_undetected = TRUE)
k2 <- upset_pq(c12, fact = "sample_from_sheet", name='sample/subsample', set_sizes=FALSE) + 
  ggtitle("All data, libs >10k, rarefied, \nmixed sample 12 & subsamples") +
  theme(plot.title = element_text(size = 10), axis.title.x=element_blank())

tab <- subsample_etc_gotus$repfil_ssdat
c12 <- subset_samples(tab, realsample1 == "C12")
c12 <- subset_samples(c12, DNAtreat == "clean")
c12<- phyloseq_validate(c12, remove_undetected = TRUE)
k3 <- upset_pq(c12, fact = "sample_from_sheet", name='sample/subsample', set_sizes=FALSE) + 
  ggtitle("Libs >10k, OTUs in => 2 PCRs, \nmixed sample 12 & subsamples") +
  theme(plot.title = element_text(size = 10), axis.title.x=element_blank())


tab <- subsample_etc_gotus$repfil_ssdattf1
c12 <- subset_samples(tab, realsample1 == "C12")
c12 <- subset_samples(c12, DNAtreat == "clean")
c12<- phyloseq_validate(c12, remove_undetected = TRUE)
k4 <- upset_pq(c12, fact = "sample_from_sheet", name='sample/subsample', set_sizes=FALSE) + 
  ggtitle("Libs >10k, OTUs in => 2 PCRs,\nOTUs >1 read per 10k, mixed sample 12 & subsamples") +
  theme(plot.title = element_text(size = 10), axis.title.x=element_blank())
library(ggpubr)
figure <- ggarrange(k1, k2, k3, k4,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)
ggsave(plot = figure, file="Figures/FigureS5.pdf", width = 212, height = 210, units = "mm")

#### Making some upset plots for sub-sample & mixed samples - c18
tab <- subsample_etc_gotus$ssdat
c18 <- subset_samples(tab, realsample1 == "C18")
c18 <- subset_samples(c18, DNAtreat == "clean")
c18<- phyloseq_validate(c18, remove_undetected = TRUE)
k1 <- upset_pq(c18, fact = "sample_from_sheet", name='sample/subsample', set_sizes=FALSE) + 
  ggtitle("All data, libs >10k, no other filtering, \nmixed sample 18 & subsamples") +
  theme(plot.title = element_text(size = 10), axis.title.x=element_blank())

tab <- subsample_etc_gotus$rare_ssdat
c18 <- subset_samples(tab, realsample1 == "C18")
c18 <- subset_samples(c18, DNAtreat == "clean")
c18<- phyloseq_validate(c18, remove_undetected = TRUE)
k2 <- upset_pq(c18, fact = "sample_from_sheet", name='sample/subsample', set_sizes=FALSE) + 
  ggtitle("All data, libs >10k, rarefied, \nmixed sample 18 & subsamples") +
  theme(plot.title = element_text(size = 10), axis.title.x=element_blank())

tab <- subsample_etc_gotus$repfil_ssdat
c18 <- subset_samples(tab, realsample1 == "C18")
c18 <- subset_samples(c18, DNAtreat == "clean")
c18<- phyloseq_validate(c18, remove_undetected = TRUE)
k3 <- upset_pq(c18, fact = "sample_from_sheet", name='sample/subsample', set_sizes=FALSE) + 
  ggtitle("Libs >10k, OTUs in => 2 PCRs, \nmixed sample 18 & subsamples") +
  theme(plot.title = element_text(size = 10), axis.title.x=element_blank())


tab <- subsample_etc_gotus$repfil_ssdattf1
c18 <- subset_samples(tab, realsample1 == "C18")
c18 <- subset_samples(c18, DNAtreat == "clean")
c18<- phyloseq_validate(c18, remove_undetected = TRUE)
k4 <- upset_pq(c18, fact = "sample_from_sheet", name='sample/subsample', set_sizes=FALSE) + 
  ggtitle("Libs >10k, OTUs in => 2 PCRs,\nOTUs >1 read per 10k, mixed sample 18 & subsamples") +
  theme(plot.title = element_text(size = 10), axis.title.x=element_blank())
library(ggpubr)
figure <- ggarrange(k1, k2, k3, k4,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)
ggsave(plot = figure, file="Figures/FigureS6.pdf", width = 210, height = 210, units = "mm")


#### Making some upset plots for sub-sample & mixed samples - c8
tab <- subsample_etc_gotus$ssdat
c8 <- subset_samples(tab, realsample1 == "C8")
c8 <- subset_samples(c8, DNAtreat == "clean")
c8<- phyloseq_validate(c8, remove_undetected = TRUE)
k1 <- upset_pq(c8, fact = "sample_from_sheet", name='sample/subsample', set_sizes=FALSE) + 
  ggtitle("All data, libs >10k, no other filtering, \nmixed sample 8 & subsamples") +
  theme(plot.title = element_text(size = 10), axis.title.x=element_blank())

tab <- subsample_etc_gotus$rare_ssdat
c8 <- subset_samples(tab, realsample1 == "C8")
c8 <- subset_samples(c8, DNAtreat == "clean")
c8<- phyloseq_validate(c8, remove_undetected = TRUE)
k2 <- upset_pq(c8, fact = "sample_from_sheet", name='sample/subsample', set_sizes=FALSE) + 
  ggtitle("All data, libs >10k, rarefied, \nmixed sample 8 & subsamples") +
  theme(plot.title = element_text(size = 10), axis.title.x=element_blank())

tab <- subsample_etc_gotus$repfil_ssdat
c8 <- subset_samples(tab, realsample1 == "C8")
c8 <- subset_samples(c8, DNAtreat == "clean")
c8<- phyloseq_validate(c8, remove_undetected = TRUE)
k3 <- upset_pq(c8, fact = "sample_from_sheet", name='sample/subsample', set_sizes=FALSE) + 
  ggtitle("Libs >10k, OTUs in => 2 PCRs, \nmixed sample 8 & subsamples") +
  theme(plot.title = element_text(size = 10), axis.title.x=element_blank())


tab <- subsample_etc_gotus$repfil_ssdattf1
c8 <- subset_samples(tab, realsample1 == "C8")
c8 <- subset_samples(c8, DNAtreat == "clean")
c8<- phyloseq_validate(c8, remove_undetected = TRUE)
k4 <- upset_pq(c8, fact = "sample_from_sheet", name='sample/subsample', set_sizes=FALSE) + 
  ggtitle("Libs >10k, OTUs in => 2 PCRs,\nOTUs >1 read per 10k, mixed sample 8 & subsamples") +
  theme(plot.title = element_text(size = 10), axis.title.x=element_blank())
library(ggpubr)
figure <- ggarrange(k1, k2, k3, k4,
                    labels = c("A", "B", "C", "D"),
                    ncol = 2, nrow = 2)
ggsave(plot = figure, file="Figures/FigureS7.pdf", width = 210, height = 210, units = "mm")

#### taking a quick look at the same with the small grid samples...
tab <- subsample_etc_gotus$ssdat
tabinc1 <- subset_samples(tab, experiment == "bigplot/smallss")
tabinc1 <- phyloseq_validate(tabinc1, remove_undetected = TRUE)
tab <- subsample_etc_gotus$rare_ssdat
tabinc2 <- subset_samples(tab, experiment == "bigplot/smallss")
tabinc2 <- phyloseq_validate(tabinc2, remove_undetected = TRUE)
tab <- subsample_etc_gotus$repfil_ssdat
tabinc3 <- subset_samples(tab, experiment == "bigplot/smallss")
tabinc3 <- phyloseq_validate(tabinc3, remove_undetected = TRUE)
tab <- subsample_etc_gotus$repfil_ssdattf1
tabinc4 <- subset_samples(tab, experiment == "bigplot/smallss")
tabinc4 <- phyloseq_validate(tabinc4, remove_undetected = TRUE)

library(grid)
library(circlize)
offset = 1e-10
col_fun = colorRamp2(c(0, 0.00015, 1, 1+offset),
                     c("white", "greenyellow", "darkgreen","darkolivegreen"))


k1 <- tabinc1 %>%
  tax_transform("compositional", rank = "speciesOTU") %>%
  comp_heatmap(sample_names_show = TRUE,col = col_fun, tax_seriation= "Identity", show_heatmap_legend = FALSE)
k1g <- grid.grabExpr(ComplexHeatmap::draw(k1))
k2 <- tabinc2 %>%
  tax_transform("compositional", rank = "speciesOTU") %>%
  comp_heatmap(sample_names_show = TRUE,col = col_fun, tax_seriation= "Identity")
k2g <- grid.grabExpr(ComplexHeatmap::draw(k2))
k3 <- tabinc3 %>%
  tax_transform("compositional", rank = "speciesOTU") %>%
  comp_heatmap(sample_names_show = TRUE,col = col_fun, tax_seriation= "Identity", show_heatmap_legend = FALSE)
k3g <- grid.grabExpr(ComplexHeatmap::draw(k3))
k4 <- tabinc4 %>%
  tax_transform("compositional", rank = "speciesOTU") %>%
  comp_heatmap(sample_names_show = TRUE,col = col_fun, tax_seriation= "Identity")
k4g <- grid.grabExpr(ComplexHeatmap::draw(k4))
figure1 <- ggarrange(k1g,k2g,k4g,k4g,labels = c("A", "B", "C", "D"),
                     ncol = 2, nrow = 2)
figure1
ggsave(plot = figure1, file="Figures/FigureS8.pdf", width = 210, height = 210, units = "mm")


library(vegan)
## adequate sampling on a sequence library size per sample basis.. both big grid and intensive plot samples (need ASVs in columns, samples in rows)
ttab <- t(otu_table(tabinc1))
class(ttab) <- "matrix" 
seqrare_smallgrid <- rarecurve(ttab, step = 10, tidy = TRUE, col = "blue", label = FALSE, ylab = "OTU richness", xlab = "Sequence Count", xlim=c(0,200000), frame.plot = FALSE)

## getting biggrid data in
libfilts<- readRDS("Processed_data/Biggrid_libsize_filter.RDA")
ttab <- otu_table(libfilts$biggrid_10k)
class(ttab) <- "matrix" 
seqrare_biggrid <- rarecurve(ttab, step = 10, tidy = TRUE, col = "blue", label = FALSE, ylab = "OTU richness", xlab = "Sequence Count", xlim=c(0,200000), frame.plot = FALSE)

## species accumulation (incidence) on the entire grid basis - hill numbers:
## first need to convert data to raw incidence data and require ASVs in rows and samples in columns

tabinc1.incidence <- (as(otu_table(tabinc1), "matrix"))
tabinc1.incidence <- ifelse(tabinc1.incidence > 0 , 1, 0)

tabinc2.incidence <- (as(otu_table(tabinc2), "matrix"))
tabinc2.incidence <- ifelse(tabinc2.incidence > 0 , 1, 0)

tabinc3.incidence <- (as(otu_table(tabinc3), "matrix"))
tabinc3.incidence <- ifelse(tabinc3.incidence > 0 , 1, 0)

tabinc4.incidence <- (as(otu_table(tabinc4), "matrix"))
tabinc4.incidence <- ifelse(tabinc4.incidence > 0 , 1, 0)

library(iNEXT.3D)
output_TD_inci1 <- iNEXT3D(tabinc1.incidence, diversity = 'TD', q = c(0, 1, 2), endpoint = 100,
                           datatype = "incidence_raw")
tdx1 <- ggiNEXT3D(output_TD_inci1, type = 1, facet.var = "Assemblage")

output_TD_inci2 <- iNEXT3D(tabinc2.incidence, diversity = 'TD', q = c(0, 1, 2),endpoint = 100,
                           datatype = "incidence_raw")
tdx2 <- ggiNEXT3D(output_TD_inci2, type = 1, facet.var = "Assemblage")

output_TD_inci3 <- iNEXT3D(tabinc3.incidence, diversity = 'TD', q = c(0, 1, 2),endpoint = 100,
                           datatype = "incidence_raw")
tdx3 <- ggiNEXT3D(output_TD_inci3, type = 1, facet.var = "Assemblage")

output_TD_inci4 <- iNEXT3D(tabinc4.incidence, diversity = 'TD', q = c(0, 1, 2),endpoint = 100,
                           datatype = "incidence_raw")
tdx4 <- ggiNEXT3D(output_TD_inci4, type = 1, facet.var = "Assemblage")

p1 <- output_TD_inci1$TDiNextEst$coverage_based[c(3,12,19),]
p1$data <- rep("RD", 3)
p1$mt1  <- p1$mT-0.5
p2 <- output_TD_inci2$TDiNextEst$coverage_based[c(3,12,19),]
p2$data <- rep("RD, Rarefied", 3)
p2$mt1  <- p2$mT-0.25
p3 <- output_TD_inci3$TDiNextEst$coverage_based[c(3,12,19),]
p3$data <- rep("PCRrep&Abund.Filtered", 3)
p3$mt1  <- p3$mT+0.25
p4 <- output_TD_inci4$TDiNextEst$coverage_based[c(3,12,19),]
p4$data <- rep("PCRrep&Abund.Filtered, Rarefied", 3)
p4$mt1  <- p4$mT+0.5

allnext <- rbind(p1,p2,p3,p4)

pd <- position_dodge(preserve = "total") # move them .05 to the left and right

d <- allnext %>% dplyr::filter(data %in% c('10k RD', '10k RD, Rarefied', 'PCRrep.Filtered', 'PCRrep&Abund.Filtered'))
d$data <- factor(d$data, levels=c('10k RD', '10k RD, Rarefied', 'PCRrep. Filtered', 'PCRrep&Abund.Filtered'))

t1 <- subset(tdx1$data, Order.q == 0)
t1$data <- rep("10k, RD", nrow(t1))
t2 <- subset(tdx2$data, Order.q == 0)
t2$data <- rep("10k RD, Rarefied", nrow(t2))
t3 <- subset(tdx3$data, Order.q == 0)
t3$data <- rep("PCRrep. Filtered", nrow(t3))
t4 <- subset(tdx4$data, Order.q == 0)
t4$data <- rep("PCRrep&Abund.Filtered", nrow(t4))

tall <- rbind(t1, t2, t3, t4)

d <- tall %>% dplyr::filter(data %in% c('10k, RD', '10k RD, Rarefied', 'PCRrep. Filtered', 'PCRrep&Abund.Filtered'))
d$data <- factor(tall$data, levels=c('10k, RD', '10k RD, Rarefied', 'PCRrep. Filtered', 'PCRrep&Abund.Filtered'))

# First, plot the larger dataset in light gray with a legend
gg1 <- ggplot() + 
  # First, plot the larger dataset in light gray with a legend
  geom_line(data = seqrare_biggrid, aes(x = Sample, y = Species, group = Site, color = "Big Grid")) +
  # Then, overlay the smaller dataset in black with a legend
  geom_line(data = seqrare_smallgrid, aes(x = Sample, y = Species, group = Site, color = "Small Grid")) +
  # Set axis limits
  coord_cartesian(ylim = c(0, 25), xlim = c(0, 200000)) +
  # Axis labels
  xlab("Read Count") + 
  ylab("Reference Library ASVs Diversity") +
  # Manually define legend labels and colors
  scale_color_manual(values = c("Big Grid" = "lightgray", "Small Grid" = "black")) +
  # Remove legend title and adjust theme
  guides(color = guide_legend(title = NULL)) +  # Removes the legend title
  theme(
    panel.background = element_blank(), # Removes the grey background
    panel.grid.major = element_blank(), # Removes major gridlines
    panel.grid.minor = element_blank(), # Removes minor gridlines
    axis.line = element_line(color = "black"), # Adds black axes lines
    legend.position = c(0.85, 0.15), # Position legend in bottom right
    legend.key = element_blank() # Removes legend key background
  )

# Print the plot
print(gg1)

# Print the plot to take a look 
print(gg1)


gg2 <- ggplot(data=d, aes(x=x, y=y, ymin=LCL, ymax=UCL, fill=data, lty=x>12)) + 
  coord_cartesian(ylim = c(0, 60), xlim = c(0, 90)) +
  geom_line(aes(color=data)) + 
  guides(lty="none")+
  geom_ribbon(alpha=0.5) +
  theme(
    panel.background = element_blank(), # Removes the grey background
    panel.grid.major = element_blank(), # Removes major gridlines
    panel.grid.minor = element_blank(), # Removes minor gridlines
    axis.line = element_line(color = "black"), # Adds black axes lines
    legend.justification = c(0, 1), 
    legend.position = c(0, 1), 
    legend.key = element_rect(colour = NA, fill = NA),
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.box.background = element_blank()
  ) +
  xlab("Number of samples from intensive plot") + 
  ylab("Reference Library ASVs") +
  
  # Light grey dashed vertical line from y = 0 to y = 40 at x = 3
  geom_segment(aes(x = 3, xend = 3, y = 0, yend = 40), linetype = "dashed", color = "grey70", size = 0.7) +
  
  # Black dashed vertical line from y = 0 to y = 40 at x = 12
  geom_segment(aes(x = 12.5, xend = 12.5, y = 0, yend = 40), linetype = "dashed", color = "black", size = 0.7) +
  
  # Other data points
  geom_point(aes(x = 45, y = 36.668733), colour= "#F8766D", shape = 16, show.legend = FALSE, size = 3) +
  geom_point(aes(x = 86, y = 37.688344), colour= "#7CAE00", shape = 16, show.legend = FALSE, size = 3) +
  geom_point(aes(x = 86, y = 35.688344), colour= "#00BFC4", shape = 16, show.legend = FALSE, size = 3) +
  geom_point(aes(x = 63, y = 25.573465), colour= "#C77CFF", shape = 16, show.legend = FALSE, size = 3)

# Print the plot to take a look 
print(gg2)

figure2 <- ggarrange(gg1, gg2, labels = c("A", "B"),
                     ncol = 2, nrow = 1)
figure2
ggsave(plot = figure2, file="Figures/FigureS9.pdf", width = 210, height = 105, units = "mm")

#