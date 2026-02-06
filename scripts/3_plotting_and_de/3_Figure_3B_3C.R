###########################################################
##   Script to generate Figures 3B and 3C                ##
##   Input:                                              ##
##      master rds object: cd4_tcells_hpv_abseq.rds      ##
###########################################################

#### Set up ####
library(Seurat)
library(tidyverse)
library(scCustomize)
library(ggpubr)
library(grDevices)
set.seed(42)

#setwd("Set main working directory")
source("utils/3_utils.R")
source("scripts/plotting_setup.R")

out_fold <- "output/script3"

#### Prep seurat objects ####
scObj.cd4 <- readRDS(file.path(out_fold, "../script2/cd4_tcells_hpv_abseq.rds"))
DefaultAssay(scObj.cd4) <- "integrated"

########################################
##   UMAP of cell type                ##
##   Figure 3B                        ##
########################################
output_dimplot(scObj.cd4, out_fold, "Fig3B_dimplot_cd4_celltype.tif",
               group = "celltype", width = 3250, height = 2500, cols = IMC_colors)

########################################
##   Cluster marker bubble plot       ##
##   Figure 3C                        ##
########################################
# Custom markers request
markers_abfilt <- c("CD62L", "CCR5", "CD197-CCR7", "CD11b", "CD25-IL2RA", "CD31",
                    "CD45RA", "CD45RO", "CD56", "CD183-CXCR3", "CXCR5", "IL7R",
                    "KLRD1", "CD279-PDCD1", "CustomAbSeq32", "CD21", "CD27")

markers_filt <- c("ARG2", "BACH2", "CCL20", "CCR4", "CCR6", "CCR7", "CSF1", "CXCR5",
                  "DUSP1", "EGR1", "FOS", "FOXP3", "GATA3", "GNLY", "GZMB", "GZMH",
                  "HSPA1A", "IFIT1", "IKZF2", "IL2RA", "IL4R", "IL6ST",
                  "JUN", "KIT", "KLRK1", "LEF1", "MAF", "MX1", "NKG7", "NKRF",
                  "NR4A1", "OAS1", "OASL", "PASK", "PRF1", "PTGDR2", "PTTG1", "PTTG2",
                  "SELL", "TNF", "TIGIT", "TNF", "TNFSF9", "TRDC", "TRGC2", "XCL2")

# Check that these are appropriate search terms and amend above if needed
for (i in 1:length(markers_abfilt)) {
  print(i)
  x <- markers_abfilt[i]
  print(rownames(scObj.cd4@assays$IAB$data)[grep(x, rownames(scObj.cd4@assays$IAB$data))])
}

# Fix the matching issues and check again
markers_abfilt <- c("CD62L", "CCR5", "CD197-CCR7", "CD11b", "CD25-IL2RA", "CD31:",
                    "CD45RA", "CD45RO", "CD56", "CD183-CXCR3", "CXCR5", "IL7R", "KLRD1",
                    "CD279-PDCD1", "CustomAbSeq32", "CD21:", "CD27-")

for (i in 1:length(markers_abfilt)) {
  print(i)
  x <- markers_abfilt[i]
  print(rownames(scObj.cd4@assays$IAB$data)[grep(x, rownames(scObj.cd4@assays$IAB$data))])
}

# Tidy up search terms to use for features
markers_abfilt <- paste0(markers_abfilt, collapse = "|")
markers_abfilt <- rownames(scObj.cd4@assays$IAB$data)[grep(markers_abfilt, rownames(scObj.cd4@assays$IAB$data))]

# Get hierarchical clustering to manually order x-axis as it's missing from this function
# Abseq
psd_exp <- Seurat::AverageExpression(scObj.cd4, assays = "IAB",
                             features = markers_abfilt,
                             group.by = "celltype")
expr_matrix <- data.matrix(psd_exp[["IAB"]])
# Compute distance between columns (x-axis variables)
distance_matrix <- base::get("dist", asNamespace("stats"))(expr_matrix, method = "euclidean")
hc <- hclust(distance_matrix, method = "ward.D2")
ordered_x <- rownames(expr_matrix)[hc$order]
# Reorder short-names based on NEW plot order
renamed_markers_abfilt <- c("KLRD1","IL2RA","TIGIT", "CXCR5","SELL","CCR7","CD27","CCR5","IL7R", "CR2",
                            "CXCR3","PDCD1", "CD45RO", "PECAM1", "ITGAM","CD45RA", "NCAM1")

# RNA
psd_exp <- Seurat::AverageExpression(scObj.cd4, assays = "RNA",
                                     features = markers_filt,
                                     group.by = "celltype")
expr_matrix <- data.matrix(psd_exp[["RNA"]])
# Compute distance between columns (x-axis variables)
distance_matrix <- base::get("dist", asNamespace("stats"))(expr_matrix, method = "euclidean") 
hc <- hclust(distance_matrix, method = "ward.D2")
ordered_x_RNA <- rownames(expr_matrix)[hc$order]

# Plot
rna_plot <- scCustomize::DotPlot_scCustom(scObj.cd4,
                                          features = unique(factor(ordered_x_RNA,
                                                            levels = ordered_x_RNA)),
                                          remove_axis_titles = FALSE,
                                          cluster.idents = TRUE) +
  scale_x_discrete() +
  RotatedAxis() +
  labs(x = "mRNA", y = "") +
  scale_colour_gradient2(low = "#003B7F", mid = "#ffffff", high = "#CA3433")

# Enforce the order of the idents to match RNA plotting
groups_order <- rev(c("CD4+ Th17", "CD4+ Th2", "CD4+ CXCR5+ pTfh",
                      "CD4+ T-regs", "CD4+ VD1 γδ", "CD4+ Transitional Naive",
                      "CD4+ Naive", "CD4+ Mitotic", "CD4+ Early Activated",
                      "CD4+ AP-1+ Naive", "CD4+ Activated Naive",
                      "ILC2", "CD4+ Type 1\nIFN responsive", "CD4+ Cytotoxic",
                      "CD4+ VD2 γδ", "CD4+ iNKT"))
# Reset the idents
Idents(scObj.cd4) <- factor(scObj.cd4$celltype, levels = groups_order)

ab_plot <- scCustomize::DotPlot_scCustom(scObj.cd4,
                                         features = unique(factor(ordered_x,
                                                           levels = ordered_x)),
                                         assay = "IAB", remove_axis_titles = FALSE,
                                         cluster.idents = FALSE) +
  labs(x = "Protein (AbSeq)",
       y = "Cluster") +
  scale_x_discrete(labels = unique(renamed_markers_abfilt)) +
  RotatedAxis() +
  scale_colour_gradient2(low =  "#003B7F", mid = "#ffffff", high = "#CA3433") +
  theme(legend.position = "none")

dot_plots <- ggpubr::ggarrange(ab_plot, rna_plot, widths = c(0.5, 1))

tiff(file.path(out_fold, paste0("Fig3C_marker_bubble_Custom_Markers",
                                "Red_Blue",".tif")),
     width = 4000, height = 1000,
     res = 300,
     compression = "lzw")
print(dot_plots)
dev.off()