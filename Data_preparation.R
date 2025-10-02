# Load required libraries
library("Seurat")
library(ggplot2)
library(dplyr)

# Define a list of patient IDs
# patients = c("P01", "P02", "P03", 
#              "P04", "P06", "P11", 
#              "P12", "P13", "P14", 
#              "P15", "P18", "P20", 
#              "P21", "P23", "P24",
#              "P25", "P26", "P27", 
#              "P28", "P29", "P30", 
#              "P31", "P33", "P35")

patients = c("P01", "P02", "P03", 
             "P04", "P06", "P07",
             "P12", "P15", "P20", 
             "P24", "P25", "P26", 
             "P27", "P28", "P29",
             "P31", "P33", "P35")

workdir = "/Users/bs21zh/OneDrive - University of Leeds/Sade Feldman et al/FILTERED/"
outdir = "/Users/bs21zh/OneDrive - University of Leeds/Sade Feldman et al/FIGURES"
setwd(workdir)

metadata = read.csv("/Users/bs21zh/Library/CloudStorage/OneDrive-UniversityofLeeds/Sade Feldman et al/mergedmetadata.csv")
metadata = metadata %>% filter(Treatment == "Untreated")

# Create an empty data frame to store summary information
pre_QC_summary_df <- data.frame(
  Patient_ID = character(),
  Min_RNA_Count = numeric(),
  Max_RNA_Count = numeric(),
  Median_RNA_Count = numeric(),
  Lower_Quartile_RNA_Count = numeric(),
  Upper_Quartile_RNA_Count = numeric(),
  Min_Features = numeric(),
  Max_Features = numeric(),
  Pearsons_Coef_Features_Counts = numeric(),
  Min_Mito_Genes = numeric(),
  Max_Mito_Genes = numeric(),
  Median_Mito_Genes = numeric(),
  Lower_Quartile_Mito_Genes = numeric(),
  Upper_Quartile_Mito_Genes = numeric(),
  Number_of_Cells = numeric(),
  stringsAsFactors = FALSE
)

# Create an empty data frame to store summary information
post_QC_summary_df <- data.frame(
  Patient_ID = character(),
  Min_RNA_Count = numeric(),
  Max_RNA_Count = numeric(),
  Median_RNA_Count = numeric(),
  Lower_Quartile_RNA_Count = numeric(),
  Upper_Quartile_RNA_Count = numeric(),
  Min_Features = numeric(),
  Max_Features = numeric(),
  Pearsons_Coef_Features_Counts = numeric(),
  Min_Mito_Genes = numeric(),
  Max_Mito_Genes = numeric(),
  Median_Mito_Genes = numeric(),
  Lower_Quartile_Mito_Genes = numeric(),
  Upper_Quartile_Mito_Genes = numeric(),
  Number_of_Cells = numeric(),
  stringsAsFactors = FALSE
)

mito_summary <- data.frame(
  Patient_ID = character(),
  Raw_Cell_Count = numeric(),
  Cell_Count_10 = numeric(),
  Cell_Count_20 = numeric(),
  Cell_Count_30 = numeric(),
  Cell_Count_40 = numeric(),
  stringsAsFactors = FALSE
)

summary <- data.frame(
  Patient_ID = character(),
  Raw_Cell_Count = numeric(),
  Min_Features = numeric(),
  Mito_20 = numeric(),
  Doublet_Filter = numeric(),
  #Untreated_Cells = numeric(),
  stringsAsFactors = FALSE
)


#patient_id = "P07"
for (patient_id in patients) {
  
  #-------------------------------RAW DATA PREP--------------------------------------------------
  
  # Read and prepare files
  counts <- read.table(paste0(patient_id, "_clean.tsv"),
               header = T)

  # convert to matrix
  counts.matrix <- as.matrix(counts[, -1])
  rownames(counts.matrix) <- counts$gene_symbol
  
  counts.matrix = counts.matrix[, colnames(counts.matrix) %in% metadata$Cell_ID]
  n_raw = length(colnames(counts))
  
  #----------------------------------QC--------------------------------------------------------
  
  # Initialize the Seurat object with raw (non-normalized) data
  # Set exclusion criteria: remove cells with less than 1000 genes
  sample <- CreateSeuratObject(
    counts = counts.matrix,
    min.features = 1000,
    project = patient_id,
    assay = "RNA"
  )
  meta <- sample@meta.data
  
  # Summarise median RNA counts per cell
  median_RNA_count <- summary(meta$nCount_RNA)[["Median"]]
  max_RNA_count <- summary(meta$nCount_RNA)[["Max."]]
  min_RNA_count <- summary(meta$nCount_RNA)[["Min."]]
  upper_RNA_count <- summary(meta$nCount_RNA)[["3rd Qu."]]
  lower_RNA_count <- summary(meta$nCount_RNA)[["1st Qu."]]
  
  # Summarise the number of detected features (genes) per cell
  min_features <- summary(meta$nFeature_RNA)[["Min."]]
  max_features <- summary(meta$nFeature_RNA)[["Max."]]
  
  # Count the number of cells in the sample
  n_cells <- length(sample$nFeature_RNA)
  
  # # Calculate the percentage of mitochondrial genes for each cell
  sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-", assay = "RNA")

  # Summarise the number of detected features (genes) per cell
  # min_mito <- summary(sample$percent.mt)[["Min."]]
  # max_mito <- summary(sample$percent.mt)[["Max."]]
  # median_mito <- summary(sample$percent.mt)[["Median"]]
  # upper_mito <- summary(sample$percent.mt)[["3rd Qu."]]
  # lower_mito <- summary(sample$percent.mt)[["1st Qu."]]
  
  # Pearson's correlation coefficient
  #pcc <- round(cor(meta$nFeature_RNA, meta$nCount_RNA,  method = "pearson", use="complete"), digits=2)
  
  # Add summary information to the data frame
  # pre_QC_summary_df <- rbind(
  #   pre_QC_summary_df,
  #   data.frame(
  #     Patient_ID = patient_id,
  #     Min_RNA_Count = min_RNA_count,
  #     Max_RNA_Count = max_RNA_count,
  #     Median_RNA_Count = median_RNA_count,
  #     Lower_Quartile_RNA_Count = lower_RNA_count,
  #     Upper_Quartile_RNA_Count = upper_RNA_count,
  #     Min_Features = min_features,
  #     Max_Features = max_features,
  #     Pearsons_Coef_Features_Counts = pcc,
  #     Min_Mito_Genes = min_mito,
  #     Max_Mito_Genes = max_mito,
  #     Median_Mito_Genes = median_mito,
  #     Lower_Quartile_Mito_Genes = lower_mito,
  #     Upper_Quartile_Mito_Genes = upper_mito,
  #     Number_of_Cells = n_cells
  #   )
  # )
  
  #-----------------------------QC Plots—————————————————————----------------

  pdf(paste0("QC_plots/", patient_id, "_QC_plots.pdf"))

  # violin plot of features and counts
  vln_plot <- VlnPlot(
    sample,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0.1
  ) &
    theme(plot.title = element_text(size = 10))

  print(vln_plot)

  print(vln_plot + scale_y_log10())

  # pearson's correlation coefficient between features and counts
  feature_plot <- FeatureScatter(sample,
                                 feature1 = "nFeature_RNA",
                                 feature2 = "nCount_RNA") +
    scale_y_log10() + scale_x_log10()

  print(feature_plot)

  FeatureScatter(sample,
                 feature1 = "nFeature_RNA",
                 feature2 = "percent.mt")

  dev.off()
  
  #-----------------------------FILTER—————————————————————----------------

  #filtered_sample_5 <- subset(sample, subset = percent.mt < 5)
  # filtered_sample_10 <- subset(sample, subset = percent.mt < 10)
  filtered_sample_20 <- subset(sample, subset = percent.mt < 20)
  # filtered_sample_30 <- subset(sample, subset = percent.mt < 30)
  # filtered_sample_40 <- subset(sample, subset = percent.mt < 40)

  post_mito_filter = length(filtered_sample_20$nFeature_RNA) # get number of cells after this filter

  # filter out potential doublets
  filtered_sample_20 <- subset(filtered_sample_20, subset = nFeature_RNA < 10000)
  post_doublet_filter = length(filtered_sample_20$nFeature_RNA) # get number of cells after this filter

  # filter out untreated cells
  #filtered_sample_20[["CellName"]] <- colnames(filtered_sample_20)
  #filtered_sample_20 = subset(filtered_sample_20, subset = CellName %in% metadata$Cell_ID)
  #post_untreated_filter = length(filtered_sample_20$nFeature_RNA) # get number of cells after this filter

  meta <- filtered_sample_20@meta.data

  # Summarise median RNA counts per cell
  # median_RNA_count <- summary(meta$nCount_RNA)[["Median"]]
  # max_RNA_count <- summary(meta$nCount_RNA)[["Max."]]
  # min_RNA_count <- summary(meta$nCount_RNA)[["Min."]]
  # upper_RNA_count <- summary(meta$nCount_RNA)[["3rd Qu."]]
  # lower_RNA_count <- summary(meta$nCount_RNA)[["1st Qu."]]


  # Summarise the number of detected features (genes) per cell
  # min_features <- summary(meta$nFeature_RNA)[["Min."]]
  # max_features <- summary(meta$nFeature_RNA)[["Max."]]

  # Count the number of cells in the sample
  #n_cells_5 <- length(filtered_sample_5$nFeature_RNA)
  #n_cells_10 <- length(filtered_sample_10$nFeature_RNA)
  #n_cells_20 <- length(filtered_sample_20$nFeature_RNA)
  #n_cells_30 <- length(filtered_sample_30$nFeature_RNA)
  #n_cells_40 <- length(filtered_sample_40$nFeature_RNA)

  # mito_summary <- rbind(
  #   mito_summary,
  #   data.frame(
  #     Patient_ID = patient_id,
  #     Raw_Cell_Count = n_cells,
  #     #Cell_Count_5 = n_cells_5,
  #     Cell_Count_10 = n_cells_10,
  #     Cell_Count_20 = n_cells_20,
  #     Cell_Count_30 = n_cells_30,
  #     Cell_Count_40 = n_cells_40,
  #     stringsAsFactors = FALSE
  #   )
  # )
#
  # Summarise the number of detected features (genes) per cell
  # min_mito <- summary(meta$percent.mt)[["Min."]]
  # max_mito <- summary(meta$percent.mt)[["Max."]]
  # median_mito <- summary(meta$percent.mt)[["Median"]]
  # upper_mito <- summary(meta$percent.mt)[["3rd Qu."]]
  # lower_mito <- summary(meta$percent.mt)[["1st Qu."]]

  # Pearson's correlation coefficient
  # pcc <- round(cor(meta$nFeature_RNA, meta$nCount_RNA,  method = "pearson", use="complete"), digits=2)

  # Add summary information to the data frame
  # post_QC_summary_df <- rbind(
  #   post_QC_summary_df,
  #   data.frame(
  #     Patient_ID = patient_id,
  #     Min_RNA_Count = min_RNA_count,
  #     Max_RNA_Count = max_RNA_count,
  #     Median_RNA_Count = median_RNA_count,
  #     Lower_Quartile_RNA_Count = lower_RNA_count,
  #     Upper_Quartile_RNA_Count = upper_RNA_count,
  #     Min_Features = min_features,
  #     Max_Features = max_features,
  #     Pearsons_Coef_Features_Counts = pcc,
  #     Min_Mito_Genes = min_mito,
  #     Max_Mito_Genes = max_mito,
  #     Median_Mito_Genes = median_mito,
  #     Lower_Quartile_Mito_Genes = lower_mito,
  #     Upper_Quartile_Mito_Genes = upper_mito,
  #     Number_of_Cells = n_cells
  #   )
  # )
  
  summary <- rbind(
    summary,
    data.frame(
    Patient_ID = patient_id,
    Raw_Cell_Count = n_raw,
    Min_Features = n_cells,
    Mito_20 = post_mito_filter,
    Doublet_Filter = post_doublet_filter
    #Untreated_Cells = post_untreated_filter
  )
  )
  #-----------------------------SAVE—————————————————————----------------

  #saveRDS(filtered_sample_20, file = paste0(patient_id, "_untreated.rds"))
}

write.table(qc_summary, paste0(workdir, "QC_summary.csv"))

library(tidyr)

long <- summary %>% 
  pivot_longer(
    cols = 2:5, 
    names_to = "Filters",
    values_to = "Cell_Count"
  )

long$Filters = factor(long$Filters, levels = c("Raw_Cell_Count", "Min_Features", "Mito_20" , "Doublet_Filter"
                                               #"Untreated_Cells"
                                               ))

tiff(file.path(outdir, "1.cell_filters.tiff"), units = "in", width = 10, height = 7, res = 900)

ggplot(long, aes(x = Patient_ID, y = Cell_Count, fill = Filters)) +
  geom_bar(stat = "identity", width = 0.5, position = position_dodge()) +
  scale_fill_brewer(palette = "Set2") + labs(x = 'Sample', y = 'Number of Cells')

dev.off()
