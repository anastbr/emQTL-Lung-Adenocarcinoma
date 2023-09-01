## -----------------------------------------------------------------------------
source("emQTL_all_functions.r")
library(qs) #save data as qs



one_click_emQTL <- function(path_import,path_save,data_df1.rds,data_df2.rds,biopsy_type_str,seq_biotype_str,rm_gender_effect_yesORno_str,IQR_level_filter,signif_level_bonferroni,num_assoc_low_thresh_discovery,num_assoc_in_validation, save_ID, seed){
  tic()
  #SEQ data preprosessing
  seq_data <- import_rds_data(path_import, data_df1.rds)
  seq_data$Ensembl_ID <- sapply(strsplit(seq_data$Ensembl_ID,".",fixed = T), getElement, 1)
  msigdb_genes <- import_transcript_info_msigdb()
  #filtered_msig <- filter_duplicates_msigdb(msigdb_genes) #you can also save and load the data to save some time
  msigdb_genes <- msigdb_genes[!duplicated(msigdb_genes), ]
  
  seq_data <- ENS_to_symbol_converter(msigdb_genes,seq_data)
  seq_data <- group_by_median(seq_data) #merging gene expression    
  seq_data <- convert_point_to_hyphen(seq_data)

  #seq_data <- seq_data[1:1000,]
  seq_data <- select_sample_type(seq_data, biopsy_type_str) #"primary"
  seq_data <- merge_aliquots(seq_data) ###wait with this one (comes too early in compared to the original script)
  seq_data <- filter_biotypes(seq_data,seq_biotype_str) #protein_coding
  seq_data <- remove_zero_rows(seq_data) #Meybe wait with this step!!!!
qsave(seq_data, paste0(path_save,"/preprocessed_seq_",dim(seq_data)[1],"x",dim(seq_data)[2],"_",save_ID,".qs"))
 
 ##### OR ##### use trancsripts annotated with older msigdbr version (january 2022)
  
  # seq_data <- readRDS("/open/tmp/Anastasia/LUAD/Gender_Filtered_LUAD/SexFiltered_emQTL/MsigDB_anotated_seq_38048x585.rds") #originally saved preprocessed data
  # seq_data <- select_sample_type(seq_data, biopsy_type_str) #"primary"
  # seq_data <- merge_aliquots(seq_data) ###wait with this one (comes too early in compared to the original script)
  # seq_data <- filter_biotypes(seq_data,seq_biotype_str) #protein_coding
  # seq_data <- remove_zero_rows(seq_data) #This is how it is in the original script

 #######Methylation data preprocessing  

  meth_data <- import_rds_data(path_import,data_df2.rds)
  #meth_data <- meth_data[1:2000,] #testing data
  
  meth_data <- rm_gender_effect_updated(meth_data,rm_gender_effect_yesORno_str)
  meth_data <- merge_aliquots(meth_data) ###wait with this one (comes too early in compared to the original script)
  meth_data <- select_sample_type(meth_data, biopsy_type_str) #"primary"
qsave(meth_data, paste0(path_save,"/preprocessed_meth_",dim(meth_data)[1],"x",dim(meth_data)[2],"_",save_ID))

  
  #splitting into discovery and validation
  set.seed(seed)
  splittet_samp <- split_discovery_validation(seq_data, meth_data)
  discovery_samp <- splittet_samp$discovery
  validation_samp <- splittet_samp$validation
  
  #correlation in discovery cohort
  seq_discovery <- seq_data[,colnames(seq_data) %in% discovery_samp]  
  meth_discovery <- meth_data[,colnames(meth_data) %in% discovery_samp] 
  
  meth_discovery <- filter_IQR(meth_discovery,IQR_level_filter)	
  meth_discovery <- match_column_order(seq_discovery, meth_discovery)
 
  seq_discovery <- remove_zero_rows(seq_discovery)
    # Expr: remove genes with var==0
  seq_discovery <- seq_discovery[apply(seq_discovery,1,var)>0 , ]

  #My method
  identical(colnames(meth_discovery),colnames(seq_discovery))	
  cor_discovery <- pearsons_correaltion_original(seq_discovery,meth_discovery)
  p_discovery <- as.data.frame(cor_discovery$p)
  r_discovery <- as.data.frame(cor_discovery$r)

  bonferroni_corect <- p_val_correction_bonferroni_regular(p_discovery,signif_level_bonferroni,num_assoc_low_thresh_discovery)
  p_filt_discovery <- bonferroni_corect$p_filt
  p_cut_discovery <- bonferroni_corect$pcut	

  # Correlation in validation cohort
  # Selecting genes and CpGs based on results from correlation in discovery cohort.
  
  seq_validation <- seq_data[,colnames(seq_data) %in% validation_samp]
  seq_validation <- seq_validation[rownames(seq_validation) %in% rownames(p_filt_discovery),]
  seq_validation <- remove_zero_rows(seq_validation)
    # Expr: remove genes with var==0
  seq_validation <- seq_validation[apply(seq_validation,1,var)>0 , ]

  meth_validation <- meth_data[,colnames(meth_data) %in% validation_samp]
  meth_validation <- meth_validation[rownames(meth_validation) %in% colnames(p_filt_discovery),]
  meth_validation <- match_column_order(seq_validation, meth_validation)
	
   # My Method
	identical(colnames(seq_validation),colnames(meth_validation))
  validation_corr <- pearsons_correaltion_original(seq_validation,meth_validation)
  p_validation <- validation_corr$p 
  
  cat("Is the dimention in p_filt_dicovery and p_validation still the same after removing zero expression?\n")  
  print(identical(colnames(p_filt_discovery),colnames(p_validation)))
  print(identical(rownames(p_filt_discovery),rownames(p_validation)))
  
  
	#Because we remove zero rows in validation cohort, common gene should be further used.
	cat("Is the order of genes in discovery and validation the same?\n")
	print(identical(rownames(p_filt_discovery),rownames(p_validation)))
  p_filt_discovery <- p_filt_discovery[rownames(p_filt_discovery) %in% rownames(p_validation),]
  p_filt_discovery <- match_rows_order(p_validation,p_filt_discovery)
	print(identical(rownames(p_filt_discovery),rownames(p_validation)))

	
  p_bool_discovery <- p_filt_discovery < p_cut_discovery  
  num_signif_discovery <- sum(p_bool_discovery)
  cat(paste0("Number of significant association in discovery cohort: ",num_signif_discovery,"\n"))
	
	
  p_cut_validation <- p_val_correction_bonferroni_based_on_signif(signif_level_bonferroni,num_signif_discovery)
  p_bool_validation <- p_validation < p_cut_validation

  # filtering based on signif values both in discovery and validation data.
  cat("Is the order of columns and rows the same in matrices with p-values in discovery and validation cohorts?\n")
  print(identical(rownames(p_filt_discovery),rownames(p_validation)))
  print(identical(colnames(p_filt_discovery),colnames(p_validation)))
  bool_signif_p <- p_bool_discovery & p_bool_validation
    
  #Preparing for Spectral Coclustering
  num_signif_validation <- sum(bool_signif_p)
  cat(paste("Number of significant association after validation is: ",num_signif_validation,"\n")) 
  cat(paste0("Percentage of significant association that validated is:  ",sprintf("%.2f",sum(bool_signif_p)/sum(p_bool_discovery)*100)," %","\n"))

  rowsums <- rowSums(bool_signif_p == TRUE) >= num_assoc_in_validation 
  colsums <- colSums(bool_signif_p == TRUE) >= num_assoc_in_validation
  
  p_final <- p_filt_discovery[rowsums,colsums]    

  r_emQTL_matrix <- r_discovery[rownames(r_discovery) %in% rownames(p_final), colnames(r_discovery) %in% colnames(p_final)]
  r_emQTL_matrix <- match_column_order(p_final, r_emQTL_matrix)
  r_emQTL_matrix <- match_rows_order(p_final, r_emQTL_matrix)  
  
  # save results as matrix with r and p-values
  qsave(r_emQTL_matrix,paste0(path_save,"/emQTLs_mat_r_val_",save_ID,".qs"))
  qsave(p_final,paste0(path_save,"/emQTLs_mat_p_val_",save_ID,".qs"))
    
  #bool_signif_p <- p_final < p_cut_discovery & p_final < p_cut_validation
  bool_signif_p <- bool_signif_p[rownames(bool_signif_p) %in% rownames(p_final), colnames(bool_signif_p) %in% colnames(p_final)]
  cat("Final number of emQTLs:\n")
  print(sum(bool_signif_p))
  
  bool_signif_p <- match_column_order(p_final, bool_signif_p)
  bool_signif_p <- match_rows_order(p_final, bool_signif_p)
  qsave(bool_signif_p, paste0(path_save,"/emQTLs_mat_bool_signif_p_val_",save_ID,".qs"))
  
    # Create a list with significant genes and CpGs
  boolean_list <- reshape2::melt(as.matrix(bool_signif_p), value.name="Validated", varnames=c("Gene", "Probe"))
  p_val_list <- reshape2::melt(as.matrix(p_final), value.name="p_values", varnames=c("Gene", "Probe"))
  r_val_list <- reshape2::melt(as.matrix(r_emQTL_matrix), value.name="r_values", varnames=c("Gene", "Probe"))  
  list_all <- cbind(p_val_list,r_val_list,boolean_list) #be carefull, use only if data has the same order
  #emQTLs <- list_all[list_all$boolean_values == TRUE,] #save full version not only sifnificant and then filter out during condor or else clustering.
  emQTLs <- list_all[,c("Gene", "Probe", "p_values","r_values", "Validated")]
  qsave(emQTLs, paste0(path_save,"/emQTLs_complete_list_",save_ID,".qs"))
  
  # record the running time
  timer <- toc()
  
  percentage <- num_signif_validation/sum(p_bool_discovery)*100
  input_parameters <- data.frame(matrix(ncol = 2, nrow = 27))
  colnames(input_parameters) <- c("Input_parameters", "Value")
  packages <-as.data.frame( installed.packages()[c("msigdbr","dplyr","TCGAutils","matrixStats","EnsDb.Hsapiens.v79"),c("Package","Version")])
  input_parameters$Input_parameters <- c("path_import","path_save","data_df1.rds","data_df2.rds","biopsy_type_str","seq_biotype_str","rm_gender_effect_yesORno_str","IQR_level_filter","signif_level_bonferroni","num_assoc_low_thresh_discovery","num_assoc_in_validation","","Results","Discovery_pcut","Validation_pcut","number_assoc_discovery","number_assoc_validated", "number validated (%)","final_n_emQTLs", "time elapsed","","Packages used","msigdbr","dplyr","TCGAutils","matrixStats","EnsDb.Hsapiens.v79")
  input_parameters$Value <- c(path_import,path_save,data_df1.rds,data_df2.rds,biopsy_type_str,seq_biotype_str,rm_gender_effect_yesORno_str,IQR_level_filter,signif_level_bonferroni,num_assoc_low_thresh_discovery,num_assoc_in_validation,"","Values",p_cut_discovery,p_cut_validation,num_signif_discovery,num_signif_validation,percentage,sum(bool_signif_p),timer$callback_msg,"","Version",packages$Version[1], packages$Version[2],packages$Version[3],packages$Version[4],packages$Version[5])
  
  
  write.table(input_parameters, paste0(path_save,"/Summary_",save_ID,".csv"), col.names = TRUE,row.names = FALSE, quote = FALSE,sep =";" )

  cat("---------FINISHED-------\n")
  cat(paste0("Dimention of the final correlation matrix includes ",dim(r_emQTL_matrix)[1]," genes and ", dim(r_emQTL_matrix)[2]," CpGs","\n"))
  cat("discovery_p_cut\n")
  cat(paste0(p_cut_discovery,"\n"))
  cat("validation p_cut\n")
  cat(paste0(p_cut_validation, "\n"))
  cat("Number of emQTLs:\n")
  cat(sum(bool_signif_p))
}


## -----------------------------------------------------------------------------
