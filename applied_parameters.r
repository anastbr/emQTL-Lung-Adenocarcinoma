
source("emQTL_functions.r")
source("one_click_emQTL_function.r")


tic()

one_click_emQTL(
  path_import = getwd(),
  path_save = getwd(),
  data_df1.rds ="TCGA-LUAD.htseq_fpkm-uq.rds",
  data_df2.rds = "TCGA-LUAD.methylation450_processed_imp.rds",
  biopsy_type_str = "primary",
  seq_biotype_str = "protein_coding",
  rm_gender_effect_yesORno_str = "yes",
  IQR_level_filter = 0.1,
  signif_level_bonferroni = 0.05,
  num_assoc_low_thresh_discovery = 1,
  num_assoc_in_validation = 5,
  save_ID = "manuscript",
  seed = 42 #or NULL if the seed should not be defined
)

toc()
