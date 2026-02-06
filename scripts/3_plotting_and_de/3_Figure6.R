###########################################################
##   Script to plot                                      ##
##   clustered seurat CD4 objects                        ##
##   Input:                                              ##
##      master rds object: cd4_tcells_hpv_abseq.rds      ##
###########################################################

#### Set up ####
library(Seurat)
library(tidyverse)
library(ggpubr)
library(grDevices)
set.seed(42)

source("utils/3_utils.R")
source("scripts/plotting_setup.R")

#setwd("Set main working directory")
out_fold <- "output/script3"

#### Prep seurat objects ####
scObj.cd4 <- readRDS(file.path(out_fold, "../script2/cd4_tcells_hpv_abseq.rds"))
DefaultAssay(scObj.cd4) <- "integrated"


########################################
##   Tcell Fitness signature          ##
##   TESPA1, GPR183, CD28             ##
########################################
genes <- c("TESPA1", "GPR183", "CD28")
geneSave <- paste(genes, collapse = "_")

# For RNA  assay
# Collect data and shape in to a format where you can add to metadata
TCFscore <- as.data.frame(t(as.data.frame(scObj.cd4[["RNA"]]$data[genes,])))
# Create TCF score (average of the 3 genes)
TCFscore <- cbind(TCFscore, "TCFscore" = rowSums(TCFscore[,1:length(genes)])/length(genes))
# Add to metadata
scObj.cd4 <- AddMetaData(scObj.cd4, TCFscore$TCFscore, col.name = "TCFscore")

########################################
##   UMAP Feature Plot - TCF          ##
##   FIGURE 6A                        ##
########################################
tiff(file.path(out_fold, sprintf("%s.tif", "Fig6A_UMAP_TCFscore")),
     width = 3300, height = 3000,
     res = 300,
     compression = "lzw")
plot <- Seurat::FeaturePlot(scObj.cd4,
                            features = "TCFscore",
                            reduction =  "wnn.umap",
                            pt.size = 1)  +
  consistentTheme +
  labs(x = "UMAP1",
       y = "UMAP2",
       title = "")
print(plot)
dev.off()

########################################
##   UMAP Feature Plot - TCF          ##
##   FIGURE 6B                        ##
########################################
tmp <- Seurat::AggregateExpression(scObj.cd4, assays = "RNA",
                                   slot = "counts",
                                   group.by = c("Participant"),
                                   normalization.method = "LogNormalize") %>%
  as.data.frame()
# subset for genes
tmp <- tmp[genes, ]

# Get names for patients for C1 and C2
scObj.cd4$Groups <- scObj.cd4$Cohorts
C1names <- unique(scObj.cd4$Participant[scObj.cd4$Groups == "Group 1"])
C2names <- unique(scObj.cd4$Participant[scObj.cd4$Groups == "Group 2"])

# Take mean of the values
tmp <- rbind(tmp, TCFscore = colSums(tmp)/length(genes))
colnames(tmp) <- stringr::str_replace(colnames(tmp), "RNA\\.", "")
tmp <- as.data.frame(t(as.data.frame(tmp)))
tmp <- cbind(tmp, "Groups" = ifelse(rownames(tmp) %in% C1names, "Group 1", "Group 2"))
tmp <- tmp[, c("TCFscore", "Groups")]

# Wilcox test
pVal <- wilcox.test(tmp$TCFscore[tmp$Groups == "Group 1"],
                    tmp$TCFscore[tmp$Groups == "Group 2"], alternative = "greater")
pVal <- pVal$p.value

# Set colours
col <- c("Group 1" = alpha("#7AD0FF", 0.5), "Group 2" = alpha("#1857F4", 0.7))
tiff(file.path(out_fold, sprintf("%s.tif", "Fig6B_Boxplot_TCFscore_Group")),
     width = 3000, height = 2800,
     res = 300,
     compression = "lzw")
ggplot(tmp, aes(x = Groups, y = TCFscore, fill = Groups)) +
  geom_boxplot(outlier.shape = NA) +
  ylim(c(0, 1000)) +
  scale_fill_manual(values = col) +
  theme_minimal() +
  consistentTheme +
  geom_point(position = position_jitterdodge(dodge.width = 0.78,
                                             jitter.width = 0.00,
                                             seed = 1),
             alpha = 0.6, size = 5) +
  labs(x = "",
       y = "TCF Score") +
  ggpubr::stat_compare_means(comparisons =  list(c("Group 1", "Group 2")),
                             method = "wilcox.test", label = "p.signif",
                             size = 14)
dev.off()


#########################################################
##   Statistical summary groups per cell type          ##
##   FIGURE 6C                                         ##
#########################################################
# Aggregate expression as a form of pseudobulk, done per cluster, and patient
tmp <- Seurat::AggregateExpression(scObj.cd4, assays = "RNA",
                                   slot = "counts",
                                   group.by = c("celltype", "Participant"),
                                   normalization.method = "LogNormalize") %>%
  as.data.frame()
# subset for genes
tmp <- tmp[genes, ]

# Get names for patients for C1 and C2
C1names <- unique(scObj.cd4$Participant[scObj.cd4$Cohorts == "Group 1"])
C2names <- unique(scObj.cd4$Participant[scObj.cd4$Cohorts == "Group 2"])
# Check there are no matches
sum(C1names %in% C2names) + sum(C2names %in% C1names) # 0 - no matches

# Take mean of the values
tmp <- rbind(tmp, mean = colSums(tmp)/length(genes))

# Take Apart C1 and C2 and make summary table of values for comparing, replicates as row
CohortNameList <- list("C1names" = C1names, "C2names" = C2names)
CohortMeanTCF <- list()

for (c in 1:length(CohortNameList)){
  cohort_patients <- CohortNameList[[c]]
  # Search for patient IDs in colnames to subset by cohort
  pattern <- paste(cohort_patients, collapse = "|")
  tmpdf <- as.data.frame(t(tmp[, grep(pattern, colnames(tmp))]))
  print("How many clusters missing/have no cells in this cohort?")
  print((length(unique(scObj.cd4$celltype)) * 8) - nrow(tmpdf))
  # Simplify colnames for easier subsetting
  rownames(tmpdf) <- stringr::str_replace_all(rownames(tmpdf), "\\.\\.", ".")
  rownames(tmpdf) <- stringr::str_replace(rownames(tmpdf), "RNA\\.", "")
  rownames(tmpdf) <- stringr::str_replace(rownames(tmpdf), "CD4\\.", "")
  # Make patient and celltype columns for easier table manipulation
  tmpdf <- cbind(tmpdf, "celltype" = stringr::str_extract(rownames(tmpdf), "^[^_]+"))
  tmpdf <- cbind(tmpdf, "Subject" = stringr::str_extract(rownames(tmpdf), "(?<=_)[^_]*$"))
  # Spread the data according to variables
  tmpdf <- tmpdf[, c("mean", "celltype", "Subject")] %>%
    dplyr::mutate(row = Subject) %>% # Turn pt to row
    tidyr::spread(key = celltype, value = mean) # spread the mean data across the celltypes
  
  # Re-map colnames so they look tidier
  map <-  c("CD4+ Activated Naive", "CD4+ AP-1+ Naive", "CD4+ CXCR5+ pTfh",
            "CD4+ Cytotoxic", "CD4+ Early Activated", "ILC2", "CD4+ iNKT",
            "CD4+ Mitotic", "CD4+ Naive", "CD4+ T-regs", "CD4+ Th17",
            "CD4+ Th2", "CD4+ Transitional Naive",  "CD4+ Type 1\nIFN responsive",
            "CD4+ VD1 γδ", "CD4+ VD2 γδ")
  names(map) <-  colnames(tmpdf)[3:ncol(tmpdf)]
  colnames(tmpdf)[3:ncol(tmpdf)] <- map[colnames(tmpdf)[3:ncol(tmpdf)]]
  
  # Assign to list
  CohortMeanTCF[[c]] <- tmpdf
  
  # Tidy up
  rm(cohort_patients, pattern, tmpdf)
}
names(CohortMeanTCF) <- c("Group 1", "Group 2")

# Do wilcox test on a loop for all the cell types
iterateOver <- colnames(CohortMeanTCF$`Group 1`)[3:ncol(CohortMeanTCF$`Group 1`)]
mwuList <- list()
for (i in 1:length(iterateOver)){
  tmp <- iterateOver[i]
  mwuList[[i]] <- wilcox.test(CohortMeanTCF$`Group 1`[, tmp],
                              CohortMeanTCF$`Group 2`[, tmp],
                              alternative = "greater")
  rm(tmp)
}
names(mwuList) <- iterateOver

# Extract pVals for correction
pVal <- as.numeric()
for (i in 1:length(iterateOver)){
  tmp <- mwuList[[i]]$p.value
  pVal <- c(pVal, tmp)
  rm(tmp)
}
names(pVal) <- unique(scObj.cd4$celltype)
pValDF <- data.frame(pVal = unlist(pVal))

# Correct by two different methods
pValDF$pBonf <- signif(p.adjust(pVal, method = "bonferroni"), 2)
sum(pValDF$pBonf < 0.05) # with FWER, 1 significant
which(pValDF$pBonf < 0.05)
pValDF$pBH <- signif(p.adjust(pVal, method = "BH"), 2)
sum(pValDF$pBH < 0.05) # with FDR correction, 8 significant
which(pValDF$pBH  < 0.05)
## Use FDR as it's more common, although there is risk of false positives
# One step line to replace no sig with gaps and to paste FDR to sig
pValDF$pBH <- ifelse(pValDF$pBH <= 0.05, pValDF$pBH, NA)
CohortMeanTCF_agg <- dplyr::bind_rows(CohortMeanTCF, .id = "Groups")
CohortMeanTCF_agg <- CohortMeanTCF_agg[c(1, 4:ncol(CohortMeanTCF_agg))]

# Reshape
CohortMeanTCF_agg <- tidyr::pivot_longer(CohortMeanTCF_agg,
                                         cols = colnames(CohortMeanTCF_agg)[2:ncol(CohortMeanTCF_agg)],
                                         names_to = "Cluster")

# Define plotting order so it looks neat
order <- CohortMeanTCF_agg %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarise(avg = mean(value)) %>%
  dplyr::arrange(desc(avg)) %>%
  dplyr::select(Cluster) %>%
  as.data.frame() %>%
  unlist()
CohortMeanTCF_agg$Cluster <- factor(CohortMeanTCF_agg$Cluster, levels = order)

# Save the Figure 6C
col <- c("Group 1" = alpha("#7AD0FF", 0.5), "Group 2" = alpha("#1857F4", 0.7))
tiff(file.path(out_fold,  sprintf("%s.tif", "Fig_6C_Boxplot_FDR_TCFscore_Pt_Celltype")),
     width = 3330, height = 2000,
     res = 300,
     compression = "lzw")
ggplot(CohortMeanTCF_agg, aes(x = Cluster, y = value, fill = Groups)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = col) +
  ylim(c(0, 265)) +
  theme_minimal() +
  consistentTheme +
  geom_point(position = position_jitterdodge(
             dodge.width = 0.78,
             jitter.width = 0.00, seed = 1),
             alpha = 0.6) +
  geom_label(data = pValDF, aes(x = rownames(pValDF),
                                y = 250,
                                label = pBH,
                                fill = celltype),
             na.rm = TRUE,  vjust = -0.5, fill = alpha("green", 0.4)) +
  labs(x = "",
       y = "TCF Score")
dev.off()