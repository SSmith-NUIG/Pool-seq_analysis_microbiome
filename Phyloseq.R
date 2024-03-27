library(tidyverse)
library(broom)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(janitor)
library("dada2")
library(DAtest)
library(vegan)
#BiocManager::install("phyloseq")
library(phyloseq)
library(plyr)
library(microbiome)

#devtools::install_github('schuyler-smith/phylosmith')
library(phylosmith)

#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)



otutab <- read.csv("/home/stephen/Documents/Thesis/Microbiome/Results/kaiju_filtered.otu.tab", sep = "\t", row.names = 1, header = TRUE, check.names=FALSE)
otutab = data.frame(otutab)
names(otutab) <- sub("^\\d+_", "", names(otutab))
OTU = otu_table(otutab, taxa_are_rows = TRUE)

taxtab <- read.csv("/home/stephen/Documents/Thesis/Microbiome/Results/kaiju_filtered.tax.tab", sep = "\t", row.names = 1, header = TRUE, 
                   stringsAsFactors = FALSE, na.strings=c(NA,"NA"," NA"))
taxtab = data.frame(taxtab)
taxtab = taxtab[complete.cases(taxtab[ , 6]),]
taxmat = as.matrix(taxtab)
TAX = tax_table(taxmat)

metadata = read.csv("/home/stephen/Documents/Thesis/Microbiome/Results/metadata_added.csv", sep = ",", skip=1, header=FALSE)
metadata <- plyr::rename(metadata, c("V1"="ColonyName", "V2"="Subspecies", "V3"="ColonyType", "V4"="Treatment", "V5"="Type")) #renaming the columns

META = sample_data(metadata)
rownames(META) <-metadata$ColonyName
physeq = phyloseq(OTU, TAX, META)

saveRDS(physeq, "/home/stephen/Documents/Thesis/Microbiome/Results/microbiome.phyloseq.rds")

# can take quite a while
physeq_df <- psmelt(physeq)

# basic abundance plot
rawkaijubarplot <- ggplot(physeq_df, aes(x = Sample, y = Abundance, fill = Domain)) + 
  theme_bw() +   geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle=90, size=8)) +
  facet_wrap(~ColonyType, scale="free_x") +
  ggtitle("Raw read numbers per sample")
rawkaijubarplot

physeq_bac <- subset_taxa(physeq, Domain == "Bacteria")
saveRDS(physeq_bac, "/home/stephen/Documents/Thesis/Microbiome/Results/microbiome.phyloseq_bac.rds")
physeq_bac = readRDS("/home/stephen/Documents/Thesis/Microbiome/Results/microbiome.phyloseq_bac.rds")

physeq_bac_cutoff = preDA(physeq_bac, min.abundance = 0.0050)
physeq_bac_cutoff
saveRDS(physeq_bac_cutoff, "/home/stephen/Documents/Thesis/Microbiome/Results/microbiome.phyloseq_bac_cutoff.rds")
physeq_bac_cutoff = readRDS("/home/stephen/Documents/Thesis/Microbiome/Results/microbiome.phyloseq_bac_cutoff.rds")

# more detailed abundance barplot with phylosmith
phylogeny_profile(physeq_bac_cutoff, classification = 'Species', 
                  treatment = "Type", merge = TRUE, 
                  relative_abundance = TRUE) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(legend.title = element_text(size = 15), 
  legend.text = element_text(size = 15))

# barplot of relative abundances 
physeq_relat_abund <- transform_sample_counts(physeq_bac_cutoff, function(x){x / sum(x)})
taxa_abundance_bars(physeq_relat_abund, classification = 'Species',
                    treatment = "Type", transformation = "mean") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=18)) +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size=15))

# pairwise adonis testing 

# zero correction and log transformation
physeq_zc = transform_sample_counts(physeq_bac_cutoff, function(y) sapply(y, function(x) ifelse(x==0, 1, (1-(sum(y==0)*1)/sum(y))*x)))
                                                                          
physeq_clr = transform_sample_counts(physeq_zc, function(x) log(x/exp(mean(log(x)))))
# Distance matrix creation
clr_dist_matrix = phyloseq::distance(physeq_clr, method = "euclidean")
# pairwise adonis testing
pairwise.adonis(clr_dist_matrix, sample_data(physeq_clr)$Type, sim.method = "eucledian", p.adjust.m = "bonferroni")

# diversity plot
plot_richness(physeq_bac_cutoff, "Type",  measures=c("Shannon")) + 
  geom_violin(fill='#A4A4A4',aes(colour=ColonyType)) + 
  geom_boxplot(aes(colour=ColonyType),width=0.1) + 
  geom_jitter(shape=16, position=position_jitter(0.01)) + 
  stat_summary(fun=median, geom="point", size=1, color="red") +
  theme(plot.title = element_text(size=30)) +
  theme(axis.text=element_text(size=17), axis.title=element_text(size=18)) + 
  theme(legend.position="none") + 
  scale_x_discrete(labels=c('Managed Treated', 'Managed Untreated', 'Free Living'))

# MICROVIZ BOXPLOTS
## BAR PLOTS 

my_comparisons <- rev(list(c("Wild","Managed_Untreated"),c("Wild","Managed_Treated"),
                           c("Managed_Untreated","Managed_Treated")))
p <- list()
for (i in as.data.frame(physeq_bac_cutoff@tax_table)$Species){
  if (!is.na(i)){
  plot_data <- physeq_bac_cutoff %>%
    tax_fix() %>%
    tax_transform("compositional", rank = "Species") %>%
    tax_transform("log2", zero_replace = 0.00001, chain = TRUE) %>%
    ps_get() %>%
    ps_otu2samdat(i) %>% # adds i as sample data!
    samdat_tbl()
  j = as.name(i)

  v = ggplot(plot_data, aes_string(x = "Type", j,group="Type")) +
    geom_boxplot(aes(col=Type),width = 0.5, colour = "grey35",outlier.shape=NA) +
    geom_jitter(aes(col=Type),width = 0.2, alpha = 0.5) +
    scale_y_continuous(
    ) +
    stat_compare_means(aes(label=..p.adj..), 
                       comparisons=my_comparisons, 
                       label = "p.format", 
                       method = "wilcox.test", 
                       p.adjust.method = "fdr",
                       symnum.args = list(cutpoints = c(0, 0.0001, 
                                                        0.001, 
                                                        0.01, 
                                                        0.05, 
                                                        Inf), 
                                          symbols = c("****", "***", "**", "*", "ns"))) +
    guides(color = FALSE, size = FALSE) +
    theme_bw()
  plot(v)
  p[[i]] = v
  }
}
plot_data <- physeq_bac %>%
  tax_fix() %>%
  tax_transform("compositional", rank = "Species") %>%
  tax_transform("log2", zero_replace = 0.00001, chain = TRUE) %>%
  ps_get() %>%
  ps_otu2samdat(" Photorhabdus thracensis") %>% # adds Parabacteroides as sample data!
  samdat_tbl()

v=ggplot(plot_data, aes(x = Type, y= ` Photorhabdus thracensis`, group=Type)) +
  geom_boxplot(aes(col=Type),width = 0.5, colour = "grey35", outlier.shape=NA) +
  geom_jitter(aes(col=Type),width = 0.2, alpha = 0.5) +
  ylim(15,5) +
  scale_y_continuous() +
  stat_compare_means(aes(label=..p.adj..), 
                     comparisons=my_comparisons, 
                     label = "p.format",
                     method = "wilcox.test", 
                     p.adjust.method = "fdr",
                     symnum.args = list(cutpoints = c(0, 0.0001, 
                                                      0.001, 
                                                      0.01, 
                                                      0.05, 
                                                      Inf), 
                                        symbols = c("****", "***", "**", "*", "ns"))) +
  guides(color = FALSE, size = FALSE)+
  theme_bw()
v
p[[" Photorhabdus thracensis"]] = v
p_core = list(p$` Snodgrassella alvi`, 
              p$` Lactobacillus mellis`,
              p$` Lactobacillus apis`,
              p$` Gilliamella apicola`,
              p$` Bartonella apis`,
              p$` Bifidobacterium asteroides`,
              p$` Frischella perrara`,
              p$` Arsenophonus nasoniae`,
              p$` Photorhabdus thracensis`)

grid.arrange(grobs = p_core, ncol = 3, labels = c("A","B","C","D","E","F","G","H","   I")) ## display plot

ggarrange(p_core[[1]] +rremove("xlab") +theme(axis.text.x=element_blank()) +ylim(-7,0.5), 
          p_core[[2]] +rremove("xlab") + theme(axis.text.x=element_blank()) +ylim(-13,1), 
          p_core[[3]] +rremove("xlab") +theme(axis.text.x=element_blank()) +ylim(-10,0.5),
          p_core[[4]] +rremove("xlab") +theme(axis.text.x=element_blank()) + ylim(-6,0.5), 
          p_core[[5]] +rremove("xlab") + theme(axis.text.x=element_blank()) +ylim(-8,2), 
          p_core[[6]] +rremove("xlab") +theme(axis.text.x=element_blank()) +ylim(-12,0.5),
          p_core[[7]] +scale_x_discrete(labels=c('Managed Treated', 'Managed Untreated', 'Free Living'))+ylim(-8,0.5), 
          p_core[[8]] +scale_x_discrete(labels=c('Managed Treated', 'Managed Untreated', 'Free Living'))+ylim(-15,5), 
          p_core[[9]] +scale_x_discrete(labels=c('Managed Treated', 'Managed Untreated', 'Free Living'))+ylim(-16,-3),
          ncol = 3, nrow=3,
          labels = c("A","B"," C",
                     "D","E","F",
                     "G","H"," I"),
          label.x=0.15, label.y =.95)
