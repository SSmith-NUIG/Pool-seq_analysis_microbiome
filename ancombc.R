physeq_bac= readRDS("/data2/ssmith/microbiome/microbiome.phyloseq_bac.rds")


library(ANCOMBC)
output_all = ancombc2(data = physeq_bac, tax_level = "Species",

                  fix_formula = "Type", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "Type", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(-1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2, 2, 1),
                                       solver = "ECOS",
                                       B = 100))



res_pair = output_all$res_pair
res_pair$taxon[res_pair$taxon == "OTU00092"] = "Apibacter mensalis *"
res_pair$taxon[res_pair$taxon == "OTU00044"] = "Unknown Lactobacillus *"
res_pair$taxon[res_pair$taxon == "OTU00034"] = "Lactobacillus apis *"
res_pair$taxon[res_pair$taxon == "OTU00031"] = "Lactobacillus kimbladii *"
res_pair$taxon[res_pair$taxon == "OTU00026"] = "Lactobacillus kullabergensis *"
res_pair$taxon[res_pair$taxon == "OTU00002"] = "Snodgrassella alvi"
res_pair$taxon[res_pair$taxon == "OTU00003"] = "Lactobacillus mellis"
res_pair$taxon[res_pair$taxon == "OTU00004"] = "Bartonella apis"
res_pair$taxon[res_pair$taxon == "OTU00006"] = "Unknown Gilliamella"
res_pair$taxon[res_pair$taxon == "OTU00008"] = "Gilliamella apicola"
res_pair$taxon[res_pair$taxon == "OTU00013"] = "Bifidobacterium coryneforme"
res_pair$taxon[res_pair$taxon == "OTU00015"] = "Frischella perrara"
res_pair$taxon[res_pair$taxon == "OTU00021"] = "Unknown bifidobacterium"
res_pair$taxon[res_pair$taxon == "OTU00022"] = "Bifidobacterium asteroides"
res_pair$taxon[res_pair$taxon == "OTU00032"] = "Lactobacillus sp. wkB8"
res_pair$taxon[res_pair$taxon == "OTU00033"] = "Commensalibactera sp. MX01"
res_pair$taxon[res_pair$taxon == "OTU00035"] = "Gilliamella sp. N-G2"
res_pair$taxon[res_pair$taxon == "OTU00038"] = "Lactobacillus mellifer"
res_pair$taxon[res_pair$taxon == "OTU00063"] = "Lactobacillus melliventris"
res_pair$taxon[res_pair$taxon == "OTU00080"] = "Commensalibactera intestini"
res_pair$taxon[res_pair$taxon == "OTU01039"] = "Arsenophonus nasoniae"
res_pair$taxon[res_pair$taxon == "OTU06567"] = "Candidatus Arsenophonus triatominarum"
# * = signficant between Untreated and Treated
# + = signficiant between Wild and Treated
# ^ significant between wild and untreated
# "diff_TypeWild", "diff_TypeWild_TypeManaged_Untreated"
res_pair[,c("taxon","diff_TypeWild_TypeManaged_Untreated")]
res_pair$taxon[res_pair$taxon == "Lactobacillus kullabergensis *"] = "	Lactobacillus kullabergensis *+"
res_pair$taxon[res_pair$taxon == "Lactobacillus kimbladii *"] = "Lactobacillus kimbladii *+"
res_pair$taxon[res_pair$taxon == "Lactobacillus sp. wkB8 *"] = "Lactobacillus sp. wkB8 +^"
res_pair$taxon[res_pair$taxon == "Lactobacillus apis *"] = "Lactobacillus apis *+"
res_pair$taxon[res_pair$taxon == "Genus Lactobacillus *"] = "Unknown Lactobacillus +^"
res_pair$taxon[res_pair$taxon == "Lactobacillus melliventris *"] = "Lactobacillus melliventris +^"
res_pair$taxon[res_pair$taxon == "Apibacter mensalis *"] = "Apibacter mensalis *+"
res_pair$taxon[res_pair$taxon == "Genus Bifidobacterium *"] = "Unknown Bifidobacterium +"
res_pair$taxon[res_pair$taxon == "Genus Gilliamella"] = "Unknown Gilliamella"
# lean = wild
# overwight = untreated
df_fig_pair1 = res_pair %>%
  dplyr::filter(diff_TypeManaged_Untreated == 1 |
                  diff_TypeWild == 1 | 
                  diff_TypeWild_TypeManaged_Untreated == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_TypeManaged_Untreated == 1, 
                              round(lfc_TypeManaged_Untreated, 2), 0),
                lfc2 = ifelse(diff_TypeWild == 1, 
                              round(lfc_TypeWild, 2), 0),
                lfc3 = ifelse(diff_TypeWild_TypeManaged_Untreated == 1, 
                              round(lfc_TypeWild_TypeManaged_Untreated, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)


df_fig_pair1 = res_pair %>%
  dplyr::filter(diff_TypeManaged_Untreated %in% c(1,0) |
                  diff_TypeWild %in% c(1,0) | 
                  diff_TypeWild_TypeManaged_Untreated %in% c(1,0)) %>%
  dplyr::mutate(lfc1 = ifelse(diff_TypeManaged_Untreated %in% c(1,0), 
                              round(lfc_TypeManaged_Untreated, 2), 0),
                lfc2 = ifelse(diff_TypeWild %in% c(1,0), 
                              round(lfc_TypeWild, 2), 0),
                lfc3 = ifelse(diff_TypeWild_TypeManaged_Untreated %in% c(1,0), 
                              round(lfc_TypeWild_TypeManaged_Untreated, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(q_TypeWild)


df_fig_pair1$group = recode(df_fig_pair1$group, 
                           `lfc1` = "Managed_Untreated - Managed_Treated",
                           `lfc2` = "Wild - Managed_Treated",
                           `lfc3` = "Wild - Managed_Untreated")
df_fig_pair1$group = factor(df_fig_pair1$group, 
                           levels = c("Managed_Untreated - Managed_Treated",
                                      "Wild - Managed_Treated", 
                                      "Wild - Managed_Untreated"))

lo = floor(min(df_fig_pair1$value))
up = ceiling(max(df_fig_pair1$value))
mid = (lo + up)/2
fig_pair = df_fig_pair1 %>%
  ggplot(aes(x = group, y = reorder(taxon,q_TypeWild), fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = "Log fold change") +
  geom_text(aes(group, taxon, label = value, color = "black"), size = 5) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to Managed_Treated colonies") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 11))
png(filename="/home/stephen/ancombc1.png",  width=700, height=350)
fig_pair
dev.off()
