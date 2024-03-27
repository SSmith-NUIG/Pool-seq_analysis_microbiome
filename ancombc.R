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


saveRDS(output_all, "/data2/ssmith/microbiome/ancom_output_all_bac.rds")
