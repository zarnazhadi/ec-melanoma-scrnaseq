# Load required libraries - quick change
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)
library(scater)
library(stringr)
library(dplyr)

# Define a list of patient IDs
patients = c("P01", "P02", "P03", 
             "P04", "P06", "P07", "P11", 
             "P12", "P13", "P14", 
             "P15", "P18", "P20", 
             "P21", "P23", "P24",
             "P25", "P26", "P27", 
             "P28", "P29", "P30", 
             "P31", "P33", "P35")

# Define the base path for OneDrive
work_dir = "/Users/bs21zh/OneDrive - University of Leeds/Sade Feldman et al/RAW/"
setwd(work_dir)

# Iterate over each patient
for (patient_id in patients) {
  
  #-----------------------------IMPORT—————————————————————----------------
  
  # Import data
  counts <- read.table(paste0(patient_id, ".tsv"), header = TRUE, sep = ",")
  
  #-----------------------------PSEUDOGENES—————————————————————------------
  # Identify and remove pseudogenes
  par_genes <- counts %>%
                dplyr::filter(str_detect(gene_id, "_PAR_Y")) %>%
                pull(gene_id)
  counts <- counts %>%
              dplyr::filter(!gene_id %in% par_genes)
  
  #-----------------------------GENE SYMBOLS—————————————————————-----------
  
  # Identify symbols for transcript IDs
  gene_ids <- data.frame(gene_id = counts$gene_id) %>%
                mutate(id = strsplit(gene_id, "[.]") %>% sapply(function(x) x[1]))
  
  # Check available key types
  master_gene_table <- mapIds(org.Hs.eg.db,
                              keys = gene_ids$id,
                              keytype = "ENSEMBL",
                              column = "SYMBOL")
  
  #-----------------------------MITOCHONDRIAL GENES—————————————————————-----
  
  # Identify mitochondrial genes
  mt_genes <- as.data.frame(genes(EnsDb.Hsapiens.v86, filter = ~ seq_name == "MT"))
  mt_genes <- mt_genes %>% dplyr::select(gene_name)
  master_gene_table[rownames(mt_genes)] <- mt_genes$gene_name
  
  #-----------------------------DUP AND N/A————-----—————————————————-----
  # Identify NA rows and replace NA with transcript ID
  length(which(is.na(master_gene_table))) # 29855 duplicated or NA genes
  
  gene_ids <- gene_ids %>%
    mutate(symbol = as.vector(master_gene_table))
  
  na_idx <- which(is.na(gene_ids$symbol))
  
  gene_ids <- gene_ids %>%
    mutate(symbol = ifelse(is.na(symbol), id, symbol))
  
  # Identify duplicates 
  dup_idx <- which(duplicated(gene_ids$symbol))
  duplicated <- gene_ids[dup_idx,] %>% na.omit
  duplicated_symbol <- duplicated$symbol
  
  counts$gene_symbol <- gene_ids$symbol
  
  all_dup <- counts[counts$gene_symbol %in% duplicated_symbol,]
  
  colMeans <- all_dup %>%
    group_by(gene_symbol) %>%
    dplyr::select(-gene_id, -gene_symbol) %>%
    summarise(across(everything(), mean, na.rm = TRUE))
  
  # Remove duplicates and add new values
  new_counts <-
    counts %>% dplyr::select(gene_symbol, everything(), -gene_id) %>%
    dplyr::filter(!gene_symbol %in% colMeans$gene_symbol) %>% rbind(colMeans)
  
  # Print previous and new dimensions
  cat(paste(patient_id), "\n","Previous dimensions:", paste(dim(counts), collapse = "x"), "\n",
      "New dimensions:", paste(dim(new_counts), collapse = "x"), "\n")
  
  #----------------------------------SAVE————————————------—————————-----
  # Save tables
  outdir = "/Users/bs21zh/OneDrive - University of Leeds/Sade Feldman et al/FILTERED/"
  write.table(new_counts$gene_symbol,
              paste0(outdir, patient_id, "_gene_metadata.tsv"))
  write.table(new_counts, paste0(outdir, patient_id, "_clean.tsv"))
}

#-----------------------------Comments—————————————————————----------------

# Define the source and destination paths
#source_path <- paste0(onedrive_path, patient_id, ".tsv")
#destination_path <- paste0(work_dir, "/", patient_id, ".tsv")

# Copy the file from the source to the destination
#file.copy(from = source_path, to = destination_path)

# Check if the copy was successful
#if (file.exists(destination_path)) {
#cat("File copied successfully.\n")

# If the copy was successful, you can choose to remove the original file
#file.remove(source_path)
#cat("Original file removed.\n")
#} else {
# cat("File copy failed.\n")
#}
