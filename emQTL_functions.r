
#Import required libraries
library(msigdbr)
library(dplyr)
library(TCGAutils) #used to identify biopsy origin in TCGA samples
library(tictoc)# taking time of the run
library(matrixStats) # IQR 
library(EnsDb.Hsapiens.v79) #biotype identification



# Function for calculating the mean of CpGs across all samples.
#input path as string or getwd()
#input file name as string with .rds at the end
import_rds_data <- function (path, file_name_str.rds) {
  
data <- readRDS(paste0(path,"/",file_name_str.rds))

cat("running function import_rds_data","\n")
cat(paste0("The imported data is: ", file_name_str.rds,"\n"))
cat("Dimention of the imported data is:")
print(dim(data))
return(data)
}


#Obtain info on transcripts by IDs in Homo Sapiens using msigdb library
import_transcript_info_msigdb <- function(){
  genesets <- msigdbr(species = "Homo sapiens")
  print("Running the function:  import_transcript_info_msigdb")
  print("Selected annotations: ensembl and gene symbols")
  genesets <- genesets[,c("human_ensembl_gene","human_gene_symbol")]
  print(head(genesets))
return(genesets)
}



filter_duplicates_msigdb <- function(msigdbr_transcripts_df){
  print("running function filter_duplicates")
  print("Original dimention:")
  print(dim(msigdbr_transcripts_df))

  msigdbr_transcripts_df <- msigdbr_transcripts_df %>%   group_by(human_ensembl_gene,human_gene_symbol) %>% dplyr::filter(row_number() == 1)
  print("After removing duplicates in msigdb:")
  print(dim(msigdbr_transcripts_df))
  print(head(msigdbr_transcripts_df))
return(msigdbr_transcripts_df)
}



#from Ensembl ID to gene symbol
ENS_to_symbol_converter <- function(msigdb,seq){
  print("Dimention of seq:")
  print(dim(seq))
  seq_merged <- merge(seq,msigdb, by.x = "Ensembl_ID", by.y = "human_ensembl_gene")
  seq_merged$Ensembl_ID = NULL #Removing Ensembl_ID column
  print("Dimention of preprocessed seq:")
  print(dim(seq_merged))
  return(seq_merged)
}



# Grouping by median
group_by_median <- function(seq_df){
  if(sum(duplicated(seq_df$human_gene_symbol))==0){return(subset(seq_df,select=-human_gene_symbol))}
  cat("Running function group_by_median","\n")
  cat("Dimention of the input data:","\n")
  print(dim(seq_df))
  
  dupe_genes = seq_df$human_gene_symbol[duplicated(seq_df$human_gene_symbol)]
  singles = subset(seq_df,  !(human_gene_symbol %in% dupe_genes))
  dupes   = subset(seq_df ,   human_gene_symbol %in% dupe_genes) 
  
  rownames(singles) = singles$human_gene_symbol
  singles$human_gene_symbol <- NULL
  
  seq_median <- dupes %>% group_by(human_gene_symbol) %>% dplyr::summarise_all(median)
  seq_median <- as.data.frame(seq_median)
  
  rownames(seq_median) <- seq_median$human_gene_symbol # Set symbols as rownames
  seq_median$human_gene_symbol <- NULL
  
  #check if the order of colnames are the same before combining data by rbind.
  print(all(colnames(singles)==colnames(seq_median))) #method 1

  seq_median = rbind(singles, seq_median) 
  
  cat("Dimention after merging by median:","\n")
  print(dim(seq_median))
  return(seq_median)
}



#filtering biotypes
#biotype_str = (LRG_gene, miRNA, miscRNA, protein_coding, or all)
#"3prime_overlapping_ncrna","antisense","IG_C_gene","IG_C_pseudogene","IG_D_gene","IG_J_gene","IG_J_pseudogene","IG_V_gene","IG_V_pseudogene","lincRNA","LRG_gene","macro_lncRNA","miRNA","misc_RNA","Mt_rRNA","Mt_tRNA","non_coding","polymorphic_pseudogene","processed_pseudogene","processed_transcript","protein_coding","pseudogene","ribozyme","rRNA","scaRNA","sense_intronic","sense_overlapping","snoRNA","snRNA","TEC","TR_C_gene","TR_D_gene","TR_J_gene","TR_J_pseudogene","TR_V_gene","TR_V_pseudogene","transcribed_processed_pseudogene","transcribed_unitary_pseudogene","transcribed_unprocessed_pseudogene","translated_processed_pseudogene","translated_processed_pseudogene","translated_unprocessed_pseudogene","unitary_pseudogene","unprocessed_pseudogene","vaultRNA".
filter_biotypes <- function(seq_data_df,seq_biotype_str){
  my_genes <- rownames(seq_data_df)
  biotypes <- ensembldb::select(EnsDb.Hsapiens.v79, keys = my_genes, columns=c("SYMBOL","GENEBIOTYPE"), keytype="SYMBOL")
  my_genes_grouped <- biotypes %>% group_by(GENEBIOTYPE) %>% tally() #count biotype distribution
  gene_biotypes <- my_genes_grouped$GENEBIOTYPE
  if(seq_biotype_str == "all"){return(seq_data_df)}
  if(seq_biotype_str %in% gene_biotypes == FALSE){stop("This biotype does not exist. Choose one on the list")}
#Group by biotypes and count number of genes within each group i the input data frame
  selected_biotypes <- biotypes[biotypes$GENEBIOTYPE == seq_biotype_str,]
  selected_seq_df <- seq_data_df[intersect(rownames(seq_data_df),selected_biotypes$SYMBOL),]
  cat(paste0("Number of transcrips in the input data: ",dim(seq_data_df)[1],"\n"))
  cat(paste0("Number of transcrips annotated as ",seq_biotype_str," in the output data: ",dim(selected_seq_df)[1],"\n"))
  return(selected_seq_df)
}



#Matching column names (sample IDs in SEQ with columns names in methylation data
convert_point_to_hyphen <- function(data_frame){
  print("Running convert_point_to_hyphen")
  print("Before conversion")
  print(head(colnames(data_frame)))
  colnames(data_frame) <- gsub(".","-",colnames(data_frame),fixed = TRUE)
  print("After conversion")
  print(head(colnames(data_frame)))
  return(data_frame)
}



#Extract any type of biopsy#########
#primary <- TCGAsampleSelect(my_samples,01)
#metastatic <- TCGAsampleSelect(my_samples,06)
#normal <- TCGAsampleSelect(my_samples,11)
select_sample_type <- function(data_df,sample_type_str){ #"normal","metastatic","primary"
  if(is.character(sample_type_str)==FALSE){stop("sample_type_str should be string")}
  if(sample_type_str == "all"){ return(data_df)}
  my_samples <- colnames(data_df)
  sample_type <- ifelse(sample_type_str=="normal",11,ifelse(sample_type_str=="metastatic",06,ifelse(sample_type_str=="primary",01,                                                    stop("ERROR: WRONG sample_type_str. Use primary/normal/metastatic/all"))))
  
  selected_samp_type <- my_samples[TCGAsampleSelect(my_samples,sample_type)]
  print("Number of samples with selecting biopsy type:")
  print(length(selected_samp_type))
  data_selected_samp <- data_df[,colnames(data_df) %in% selected_samp_type]
  print("New dimention:")
  print(dim(data_selected_samp))
  return(data_selected_samp)
}



#Merge aliquots
merge_aliquots <- function(data_df) {
  data_df <- as.data.frame(data_df)
  cat("Dimention of the input data:","\n")
  print(dim(data_df))
  cat("Running function merge_aliquots","\n")
  cat("Original sampleIDs","\n")
  print(head(colnames(data_df)))
  tcga_pattern = "^(.*-[0-9][0-9])[A-Za-z]*$" 
  if(!all(grepl(tcga_pattern, colnames(data_df)))) {
      wrong = colnames(data_df)[
        !grepl(tcga_pattern, colnames(data_df))]
      stop(sprintf("Wrong sampleID format in %d columns : %s ...", 
                      length(wrong), paste(head(wrong), collapse=",")))
  }
    
  colnames(data_df) <- gsub(tcga_pattern, "\\1", colnames(data_df))
  
  ifelse(sum(duplicated(colnames(data_df)))!= 0,(data_df <- sapply(unique(colnames(data_df)), function(sampleid)
  rowMeans(data_df[grepl(sampleid,colnames(data_df))]))), data_df)
  
  cat("After merging aliquotes sampleIDs","\n")
  print(head(colnames(data_df)))
  cat("Dimention after merging:","\n")
  print(dim(data_df))
  data_df <- as.data.frame(data_df)
  return(data_df)
}



#Matching sample IDs in seq and meth
get_common_samples <- function(data_df1, data_df2){  
  common_samples <- intersect(colnames(data_df1), colnames(data_df2))
  print("number of common samples:")
  print(length(common_samples))
  return(common_samples)
}


# Splitting data into discovery and validation
split_discovery_validation <- function(data_df1, data_df2){
  common_samples <- get_common_samples(data_df1, data_df2)
  discovery <- sample(common_samples,length(common_samples)/2)
  validation <- common_samples[!common_samples %in% discovery]
 sampleID_split_list = list(
  discovery = discovery,
  validation = validation 
 )
  
  return(sampleID_split_list)
}

#!!!!!!!!!!!probably don't have use for this function anymore
match_samplIDs_common <- function(data_df,common_samples){
  print("Original dimention")
  print(dim(data_df))
  matched_df <- data_df[,colnames(data_df) %in% common_samples]
  print("After matching:")
  print(dim(matched_df))
  return(matched_df)
}

#call this function inside the correlation function
match_column_order <- function(data_df1, data_df2_to_match){
  print("Does order of sampleIDs match in both data frames?")
  print(all(colnames(data_df1)==colnames(data_df2_to_match))) #method 1
  print(identical(colnames(data_df1),colnames(data_df2_to_match))) #method 2

  print("data_df1")
  print(head(colnames(data_df1)))
  indx <- match(colnames(data_df1),colnames(data_df2_to_match))  
  #Match the colnames in data_df1 with                                                         colnames in data_df1
	#The second input will change while the first one stays the same.
  data_df2_ordered <- data_df2_to_match[,indx]
  print("data_df2")
  print(head(colnames(data_df2_ordered)))
  print("Does order of sampleIDs match in both data frames now?")
  print(all(colnames(data_df1)==colnames(data_df2_ordered))) #method 1
  print(identical(colnames(data_df1),colnames(data_df2_ordered))) #method 2
  return(data_df2_ordered)
}

#call this function inside the correlation function
match_rows_order <- function(data_df1, data_df2_to_match){
  print("Does order of rows match in both data frames?")
  print(all(rownames(data_df1)==rownames(data_df2_to_match))) #method 1
  print(identical(rownames(data_df1),rownames(data_df2_to_match))) #method 2

  print("data_df1")
  print(head(rownames(data_df1)))
  indx <- match(rownames(data_df1),rownames(data_df2_to_match))  
  #Match the colnames in data_df1 with                                                         colnames in data_df1
	#The second input will change while the first one stays the same.
  data_df2_ordered <- data_df2_to_match[indx,]
  print("data_df2")
  print(head(rownames(data_df2_ordered)))
  print("Does order of sampleIDs match in both data frames now?")
  print(all(rownames(data_df1)==rownames(data_df2_ordered))) #method 1
  print(identical(rownames(data_df1),rownames(data_df2_ordered))) #method 2
  return(data_df2_ordered)
}

#call this function inside correlation function
remove_zero_rows <- function(data_df){
  print("Running function remove_zero_rows")
  print("Dimention of input data:")
  print(dim(data_df))
  #data_df <- data.matrix(data_df, rownames.force = TRUE)
  data_df <- data_df[rowSums(data_df) != 0,] #remove rows with all zeros
  print("Dimention after remove zero rows:")
  print(dim(data_df))
  return(data_df)
}

#remove CpGs with low variation
filter_IQR <- function(data_df, IQR_level_filter){
  data_df <- as.matrix(data_df) # input has to be a matrix (latest update)
  #IQR_level_filter <- numeric(IQR_level_filter)
  data_df <- data_df[rowIQRs(data_df)>IQR_level_filter,]
  data_df <- as.data.frame(data_df)
  return(data_df)
}

#filter out CpGs on X nad Y chromosome
rm_gender_effect_updated <- function(data_df,rm_gender_effect_yesORno_str){
  if(rm_gender_effect_yesORno_str=="no"){return(data_df)}
  print("number of CpGs in the input data:")
  print(dim(data_df)[1])
  #CpG_info <- read.table("/data2/thomas/breast450k/NewProbeinfo.txt", header = TRUE, sep ="\t")
  CpG_info <- read.table("/open/work/Anastasia/one_click_emQTL/NewProbeinfo.txt", header = TRUE, sep ="\t")
  CpG_X_Y <- subset(CpG_info, Chr == "Y" | Chr=="X")
  data_df_noXY <- data_df[!rownames(data_df) %in% CpG_X_Y$Probe,]
  print("number of CpGs ion X and Y chromosomes removed:")
  print(dim(data_df)[1]-dim(data_df_noXY)[1])
  print("Final number of remaining CpGs")
  print(dim(data_df_noXY)[1])
  data_df_noXY <- as.data.frame(data_df_noXY)
    return(data_df_noXY)
}

rm_gender_effect <- function(data_df,rm_gender_effect_yesORno_str){
CpG_info <- read.table("/data2/thomas/breast450k/NewProbeinfo.txt", header = TRUE, sep ="\t")

only_cpg <- as.data.frame(rownames(data_df))
colnames(only_cpg) = "CpG"

complete <- merge(CpG_info,only_cpg, by.x = "Probe", by.y = "CpG" )

Y <- subset(complete, Chr == "Y")
X <- subset(complete, Chr == "X")

data_df <- data_df[setdiff(rownames(data_df), X$Probe),]
data_df <- data_df[setdiff(rownames(data_df), Y$Probe),]

    return(data_df)
 }


#Correlation emQTL
pearsons_correaltion <- function(data_df1_seq, data_df2_meth){
  if("matrix" %in% class(data_df1_seq)) {
    r <- cor(t(data_df1_seq),t(data_df2_meth))
    n <- ncol(data_df1_seq)
  } else {
    r <- cor(data_df1_seq, data_df2_meth)
    n = length(data_df1_seq)
  }
  
  df <- n-2
  t <- (r*sqrt(df))/(sqrt(1-r^2))
  p <- 2*pt(abs(t),df,lower.tail=FALSE) # pt gives the distribution function
  res <- list(
    r = r,
    p = p
  )
  return(res)
}

	pearsons_correaltion_original <- function(data_df1_seq, data_df2_meth){
	
	r <- cor(t(data_df1_seq),t(data_df2_meth))
	n <- ncol(data_df1_seq)
	df <- n-2
	t <- (r*sqrt(df))/(sqrt(1-r^2))
	p <- 2*pt(abs(t),df,lower.tail=FALSE)# pt gives the distribution function
	res <- list(
	r = r,
	p = p
	)
	return(res)
	}



p_val_correction_bonferroni_regular <- function( p_df,signif_level,num_assoc_low_thresh){
  pcut <- signif_level/(as.numeric(nrow(p_df))*as.numeric(ncol(p_df))) 
  rowsums <- rowSums(p_df < pcut) >= num_assoc_low_thresh 
  colsums <- colSums(p_df < pcut) >= num_assoc_low_thresh
  p_filt <- p_df[rowsums,colsums]
  res <- list(
    pcut = pcut,
    p_filt = p_filt
  )
  return(res)
}


p_val_correction_bonferroni_based_on_signif <- function(signif_level,num_signif){
  pcut <- signif_level/num_signif
  return(pcut)
}

