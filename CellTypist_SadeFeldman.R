library("dplyr")
library("Seurat")
library("SeuratObject")
library(ggplot2)

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

r_path = "/Users/bs21zh/OneDrive - University of Leeds/Sade Feldman et al/FILTERED/"
outdir = "/Users/bs21zh/OneDrive - University of Leeds/Sade Feldman et al/FIGURES"

setwd(r_path)

sample_list <- list()

# Loop through samples MEL001 to MEL033
for (patient_id in patients) {
  # Read the RDS file for the current sample
  file_path <- paste0(patient_id, "_untreated.rds")
  mel_data <- readRDS(file_path)
  
  # Add the data frame to the list
  sample_list[[patient_id]] <- mel_data
}

sample <- merge(sample_list[["P01"]], y = c(sample_list[[patients[2]]], sample_list[[patients[3]]],
                                               sample_list[[patients[4]]], sample_list[[patients[5]]],
                                               sample_list[[patients[6]]], sample_list[[patients[7]]],
                                               sample_list[[patients[8]]], sample_list[[patients[9]]],
                                               sample_list[[patients[10]]], sample_list[[patients[11]]],
                                               sample_list[[patients[12]]], sample_list[[patients[13]]],
                                               sample_list[[patients[14]]], sample_list[[patients[15]]],
                                               sample_list[[patients[16]]], sample_list[[patients[17]]],
                                               sample_list[[patients[18]]]), 
                add.cell.ids = patients, project = "Sade Feldman")

#-----------------------------Normalisation—————————————————————----------------
sample <- NormalizeData(sample)
sample <- JoinLayers(sample)

workdir <- "/Users/bs21zh/OneDrive - University of Leeds/Sade Feldman et al/INTEGRATED/"
setwd(workdir)

sample = FindVariableFeatures(sample, selection.method = "vst", nfeatures = 10000)

top10 = head(VariableFeatures(sample), 10)
top10

# scale
all.genes = rownames(sample)
scaled_sample = ScaleData(sample, features = all.genes)

# pca
scaled_sample <- RunPCA(scaled_sample)
saveRDS(scaled_sample, "sample_scaled.RDS")

tiff(file.path(outdir, "2.elbow_plot.tiff"), units = "in", width = 10, height = 7, res = 900)

ElbowPlot(scaled_sample, ndims = 30)

dev.off()

tiff("pca_sample_origin.tiff", units = "in", width = 10, height = 7, res = 900)

DimPlot(scaled_sample, reduction = "pca") 

dev.off()

scaled_umap = scaled_sample %>% 
  FindNeighbors(dims = 1:7) %>% 
  FindClusters(resolution = 0.8) %>% 
  RunUMAP(dims = 1:7)

tiff("umap_unannotated.tiff", units = "in", width = 10, height = 7, res = 900)

d1 = DimPlot(scaled_umap, label.size = 4, repel = T, label = T)
d2 = DimPlot(scaled_umap, group.by = "orig.ident", label.size = 4) + ggtitle("Grouped by Sample ID")
d1 | d2

dev.off()

saveRDS(scaled_umap, "sample_umap.rds")

library("SeuratDisk")
library("SCopeLoomR")
SeuratDisk::as.loom(scaled_umap, "sample.loom")

r_path = "/Users/bs21zh/OneDrive - University of Leeds/Sade Feldman et al/INTEGRATED/"
setwd(r_path)

scaled_umap = readRDS("sample_umap.rds")

# Analyse CellTypist results
annot = read.csv("celltypist_annotations_low.csv", header = T, sep = "\t")

# Add cell annotations to Seurat UMAP
all(annot$CellID == rownames(scaled_umap@meta.data))
scaled_umap@meta.data = cbind(scaled_umap@meta.data, annot$majority_voting)
colnames(scaled_umap@meta.data)[7] = "majority_voting"

tiff("umap_annotated.tiff", units = "in", width = 10, height = 7, res = 900)

DimPlot(scaled_umap, reduction = "umap", group.by = "majority_voting", label = T, repel = T) + 
  ggtitle("UMAP grouped by CellTypist Annotation")

dev.off()

library(patchwork)

# Define the cell groups to highlight
all_cells = dim(scaled_umap@meta.data)[1]
nk_cells <- rownames(scaled_umap@meta.data[scaled_umap@meta.data$majority_voting == "NK cells",])
cd16_neg_nk_cells <- rownames(scaled_umap@meta.data[scaled_umap@meta.data$majority_voting == "CD16- NK cells",])
cd16_pos_nk_cells <- rownames(scaled_umap@meta.data[scaled_umap@meta.data$majority_voting == "CD16+ NK cells",])

# Create individual DimPlots
p1 <- DimPlot(scaled_umap, reduction = "pca", 
              cells.highlight = c(nk_cells, cd16_neg_nk_cells, cd16_pos_nk_cells)) + 
  NoLegend() + ggtitle("CD16+ NK cells, CD16- NK cells and NK cells")

p2 <- DimPlot(scaled_umap, reduction = "pca", 
              cells.highlight = nk_cells) + 
  labs(subtitle = paste0("n = ", length(nk_cells), " (", signif(length(nk_cells)/all_cells *100, digits = 3), "%)")) +
  NoLegend() + ggtitle("NK Cells")

p3 <- DimPlot(scaled_umap, reduction = "pca", 
              cells.highlight = cd16_neg_nk_cells) + 
  labs(subtitle = paste0("n = ", length(cd16_neg_nk_cells), " (", signif(length(cd16_neg_nk_cells)/all_cells *100, digits = 3), "%)")) +
  NoLegend() + ggtitle("CD16- NK cells")

p4 <- DimPlot(scaled_umap, reduction = "pca", 
              cells.highlight = cd16_pos_nk_cells) +
  labs(subtitle = paste0("n = ", length(cd16_pos_nk_cells), " (", signif(length(cd16_pos_nk_cells)/all_cells *100, digits = 3), "%)")) +
  NoLegend() + ggtitle("CD16+ NK cells")

# Combine all plots into one figure using patchwork
combined_plot <- p1 + p2 + p3 + p4 + plot_layout(ncol = 2)  # Arrange in 2x2 grid


tiff("pca_NK.tiff", units = "in", width = 10, height = 7, res = 900)

combined_plot

dev.off()


# Feature Plots
pos_sig = list(c("NCR1", "NCAM1", "GNLY", "FCGR3A", "KLRB1"))
neg_sig = list(c("CD3D", "CD3E", "CD3G", "GZMK"))

# CellTypist signatures
cd16_neg_sig = list(c("TYROBP", "LDB2", "CCL3", "CLIC3", "NKG7", "GSTP1",
                      "FCER1G", "IRF8", "CXCL3", "ADGRG3"))
cd16_pos_sig = list(c("FCER1G", "GNLY", "TYROBP", "IGFBP7", "GZMB", "PTGDS",
                      "MYOM2", "FCGR3A", "SPON2", "SH2D1B"))
nk_cell_sig = list(c("NKG7", "GNLY", "CST3", "MALAT1", "SPINK2",
                     "KLRD1", "HMOX1", "IGFBP4", "TRBV10-2", "CD3E"))

scaled_umap = AddModuleScore(scaled_umap, features = pos_sig, name = "Positive_Signature")
scaled_umap = AddModuleScore(scaled_umap, features = neg_sig, name = "Negative_Signature")
scaled_umap = AddModuleScore(scaled_umap, features = cd16_neg_sig, name = "CD16_Negative_Signature")
scaled_umap = AddModuleScore(scaled_umap, features = cd16_pos_sig, name = "CD16_Positive_Signature")
scaled_umap = AddModuleScore(scaled_umap, features = nk_cell_sig, name = "NK_Cell_Signature")

library(viridis)

f1 = FeaturePlot(scaled_umap, pt.size = 0.1,
                 features = "Positive_Signature1", order = TRUE) &
  scale_colour_viridis(option="rocket", direction = -1) & ggtitle("Positive NK Cell Signature") &
  labs(subtitle = "NCR1, NCAM1, GNLY, FCGR3A, KLRB1")
f2 = FeaturePlot(scaled_umap, pt.size = 0.1,
                 features = "Negative_Signature1", order = TRUE) &
  scale_colour_viridis(option="rocket", direction = -1) & ggtitle("Negative NK Cell Signature") &
  labs(subtitle = "CD3D, CD3E, CD3G")

combined_plot <- f1 + f2 + plot_layout(ncol = 2)

tiff("umap_signature.tiff", units = "in", width = 10, height = 7, res = 900)

print(combined_plot)

dev.off()

tiff("dotplot_signature.tiff", units = "in", width = 10, height = 7, res = 900)

DotPlot(scaled_umap,
        features = c(unlist(pos_sig), unlist(neg_sig)),
        group.by = "majority_voting") +
  ylab("CellTypist Annotation") +
  scale_colour_viridis(option = "rocket", direction = -1) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))

dev.off()

tiff("featureplot_signature.tiff", units = "in", width = 10, height = 7, res = 900)

FeaturePlot(scaled_umap, features = unlist(c(pos_sig, neg_sig)), cols = c("lightgrey", "red3"))

dev.off()

# Scatter plot
# Generate a color palette using Set2
library(RColorBrewer)
mycolors2 = c("#66C2A5", "#FC8D62", rep("grey", 10),  rep("#8DA0CB", 2))

tiff("scatter.tiff", units = "in", width = 13, height = 7, res = 900)
ggplot(scaled_umap@meta.data, aes(x = Negative_Signature1, y = Positive_Signature1, color = majority_voting)) + 
  geom_point() + scale_color_manual(values = mycolors2) +
  labs(y = "Positive Signature Score", x = "Negative Signature Score", color = "CellTypist Annotation") +
  geom_vline(data=scaled_umap@meta.data, 
             mapping=aes(xintercept=0), 
             linetype="dashed", 
             color="black") + 
  geom_hline(data=scaled_umap@meta.data, 
             mapping=aes(yintercept=0), 
             linetype="dashed", 
             color="black") 

dev.off()

#------------------------------ VENN DIAGRAM ----------------------------------

high_pos_sig_threshold1 = rownames(scaled_umap@meta.data[scaled_umap@meta.data$Positive_Signature1 > 0.5 & scaled_umap@meta.data$Negative_Signature1 < 0,])
high_pos_sig_threshold2 = rownames(scaled_umap@meta.data[scaled_umap@meta.data$Positive_Signature1 > 0.5 & scaled_umap@meta.data$Negative_Signature1 < 1.5,])
high_pos_sig = rownames(scaled_umap@meta.data[scaled_umap@meta.data$Positive_Signature1 > 0 & scaled_umap@meta.data$Negative_Signature1 < 0,])

library(VennDiagram)

# 3 circles
myCol <- brewer.pal(3, "Set2")
# Chart
venn.diagram(
  x = list(cd16_pos_nk_cells, high_pos_sig, cd16_neg_nk_cells),
  category.names = c("CD16+ NK Cells" ,
                     "Positive Score > 0\nNegative Score < 0",
                     "CD16- NK Cells"),
  filename = 'venn1.png',
  output = T,
  
  # Output features
  imagetype = "png" ,
  height = 800 ,
  width = 900 ,
  resolution = 500,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 0.3,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = .3,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 180, 0),
  cat.dist = c(0.05, 0.05, 0.025),
  cat.fontfamily = "sans"
)

venn.diagram(
  x = list(cd16_pos_nk_cells, high_pos_sig_threshold1, cd16_neg_nk_cells),
  category.names = c(
    "CD16+ NK Cells" ,
    "Positive Score > 0.5\nNegative Score < 0",
    "CD16- NK Cells"
  ), filename = 'venn3.png',
  output = T,
  
  # Output features
  imagetype = "png" ,
  height = 800 ,
  width = 900 ,
  resolution = 500,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 0.3,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = .3,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 180, 0),
  cat.dist = c(0.05, 0.05, 0.025),
  cat.fontfamily = "sans",
  rotation = 1
)

venn.diagram(
  x = list(cd16_pos_nk_cells, high_pos_sig_threshold2, cd16_neg_nk_cells),
  category.names = c(
    "CD16+ NK Cells" ,
    "Positive Score > 0.5\nNegative Score < 1.5",
    "CD16- NK Cells"
  ), filename = 'venn2.png',
  output = T,
  
  # Output features
  imagetype = "png" ,
  height = 800 ,
  width = 900 ,
  resolution = 500,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = 0.3,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = .3,
  #cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(0, 180, 0),
  cat.dist = c(0.05, 0.05, 0.025),
  cat.fontfamily = "sans"
)

FeaturePlot(scaled_umap, features = c("CLIC3", "SH2D1B", "MYOM2", "ADGRG3"), cols = c("lightgrey", "red2"))
