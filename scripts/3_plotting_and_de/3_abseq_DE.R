###########################################################
##   Script to explore and plot                          ##
##   clustered seurat CD4 objects                        ##
##   Input:                                              ##
##      master rds object: cd4_tcells_hpv_abseq.rds      ##
###########################################################

#### Set up ####
library(Seurat)
library(tidyverse)
library(scCustomize)
library(ggpubr)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(ReactomePA)
library(org.Hs.eg.db)
library(grDevices)
library(ggarchery)
set.seed(42)

#setwd("Set main working directory")
source("utils/3_utils.R")
source("scripts/plotting_setup.R")

out_fold <- "output/script3"
dir.create(out_fold, showWarnings = FALSE)

#### Prep seurat objects ####
scObj.cd4 <- readRDS(file.path(out_fold, "../script2/cd4_tcells_hpv_abseq.rds"))
DefaultAssay(scObj.cd4) <- "integrated"

#################################################
##   Cells from Group1 or Group2 per cluster   ##
#################################################
# Get table of numbers per cluster
m <- scObj.cd4@meta.data
m1 <- m %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(n_celltype = n())
m2 <- m %>%
  dplyr::group_by(celltype, hpv_clearance) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::left_join(m1) %>%
  dplyr::mutate(prop = n / n_celltype)
m2pos <- m2 %>%
  dplyr::filter(hpv_clearance == "HPV detected") %>%
  dplyr::ungroup() %>%
  dplyr::arrange(prop)

# Add whole dataset numbers
m3 <- m %>% 
  dplyr::group_by(hpv_clearance) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate("Cluster" = "All") %>%
  dplyr::rename(c("celltype" = "Cluster")) %>%
  dplyr::mutate(n_celltype = sum(n)) %>%
  dplyr::mutate(prop = n / n_celltype)
m_all <- rbind(m3, m2) %>%
  dplyr::select(celltype, hpv_clearance, n, n_celltype, prop) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(prop)

# Save table
write.table(m_all, file.path(out_fold, "table_cd4_dnapos_v_dnaneg_per_cluster.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#######################################
##   Cells per patient per cluster   ##
##   Input for Sup Fig 6             ##
#######################################
# Get table of numbers per cluster
pat <- m %>% 
  dplyr::group_by(celltype, Participant) %>%
  dplyr::summarise(n_participant = n()) %>%
  dplyr::left_join(m1) %>%
  dplyr::mutate(prop = n_participant / n_celltype)

# Add whole dataset numbers
pat_total <- m %>% 
  dplyr::group_by(Participant) %>%
  dplyr::summarise(n_participant = n()) %>%
  dplyr::mutate("Cluster" = "All") %>%
  dplyr::rename(c("celltype" = "Cluster")) %>%
  dplyr::mutate(n_celltype = sum(n_participant)) %>%
  dplyr::mutate(prop = n_participant / n_celltype)
pat_all <- rbind(pat_total, pat) %>%
  dplyr::select(celltype, Participant, n_participant, n_celltype, prop) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(prop)

# Save table
write.table(pat_all, file.path(out_fold, "table_cd4_patients_per_cluster.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

################################################
## Celltype vs HPV clearance barplot          ##
## Stacked barplot                            ##
## FIGURE 4A                                  ##
################################################
m_all$celltype <- factor(m_all$celltype,
                         levels = c(as.character(m2pos$celltype), "All"))

bar_hpv <- ggplot(m_all, aes(fill = hpv_clearance,
                             x = celltype,
                             y = prop)) +
  geom_bar(position = "stack",
           stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        panel.grid.major.y = element_line(color = "black")) +
  scale_fill_manual(values = unname(IMC_colors)[2:1],
                    name = "HPV clearance",
                    labels = c("HPV detected", "HPV non-detected")) +
  xlab("Cluster") +
  ylab("Proportion in cluster")
bar_hpv
# Save
grDevices::cairo_pdf(filename = file.path(out_fold,
                     "stackbar_cd4_dnapos_v_dnaneg_per_cluster.pdf"),
                     width = 6, height = 3)
bar_hpv
dev.off()

################################################
## Diff Exp cell type and HPV clearance       ##
## in CD4+ AP-1+ naive cell population        ##
################################################
# Add metadata col for cell type + clearance status
scObj.cd4$celltype_hpv <- paste(scObj.cd4$celltype,
                                scObj.cd4$hpv_clearance,
                                sep = "_")

# Run differential expression
ap1_dif_exp <- clust_de(scObj.cd4, colname = "celltype_hpv",
                        id1 = "CD4+ AP-1+ Naive_HPV detected",
                        id2 = "CD4+ AP-1+ Naive_HPV non-detected",
                        out_fold = out_fold,
                        filename = "table_cd4_diff_exp_ap1.naive.cd4_dnapos_v_dnaneg.txt")


######################
## Volcano plot     ##
## FIGURE 4D        ##
######################
plot_volc(ap1_dif_exp, title = "AP-1+ naive CD4",
          subtitle = "DNA positive vs DNA cleared",
          height = 10, width = 8, out_fold = out_fold,
          filename = "volc_cd4_ap1.naive.cd4_dnapos_v_dnaneg.pdf",
          arrows = TRUE)


################################################
## Diff Exp cell type and HPV clearance       ##
## across all cells                           ##
## SUP DATA 2                                 ##
################################################
de_hpv <- clust_de(scObj.cd4, colname = "hpv_clearance",
                   id1 = "HPV detected", id2 = "HPV non-detected",
                   out_fold = out_fold,
                   filename = "table_cd4_diff_exp_dnapos_v_dnaneg.txt")

# Volcano plot
plot_volc(de_hpv, title = "DNA positive vs DNA cleared", subtitle = "All cells",
          height = 10, width = 8, out_fold = out_fold,
          filename = "volc_cd4_dnapos_v_dnaneg.pdf", arrows = TRUE)

# Get mean expression for each gene
Idents(scObj.cd4) <- "avg_exp"
avg_exp <- as.data.frame(Seurat::AverageExpression(scObj.cd4)$integrated)
avg_exp_ab <- as.data.frame(Seurat::AverageExpression(scObj.cd4)$IAB)
avg_exp$gene <- rownames(avg_exp)
avg_exp_ab$gene <- rownames(avg_exp_ab)
avg_exp_all <- rbind(avg_exp, avg_exp_ab)
de_hpv_ma <- dplyr::left_join(de_hpv, avg_exp_all) %>%
  dplyr::rename(c("baseMean" = "all",
                  "log2FoldChange" = "avg_log2FC",
                  "padj" = "p_val_adj"))
write.table(de_hpv_ma,
            file = file.path(out_fold, "table_cd4_dnapos_v_dnaneg_avg_exp.txt"), 
            quote = FALSE, row.names = FALSE, sep = "\t")


############################
## DEseq plots - MA plot  ##
## FIGURE 3E              ##
############################
maplot <- ggpubr::ggmaplot(de_hpv_ma, fdr = 0.05, fc = 1.41421,
                           genenames = de_hpv_ma$gene,
                           top = 20,
                           size = 2) +
  consistentTheme +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1))

tiff(file.path(out_fold, "Fig3E_maplot_cd4_dnapos_v_dnaneg_newFC_top20ANN.tif"),
     width = 2250, height = 1500,
     res = 300,
     compression = "lzw")
print(maplot)
dev.off()


### Over-representation analysis for de genes ###
Idents(scObj.cd4) <- scObj.cd4$celltype

# Remove abseq genes, get sig p vals
de_hpv_filt2 <- de_hpv %>%
  dplyr::filter(!grepl("-pAbO", gene)) %>%
  dplyr::filter(p_val_adj <= 0.05)

# Prep list of lfcs
# returns: list(l2fc, l2fc_up, l2fc_down)
l2fc_list <- prep_lfcs(de_hpv_filt2, 0.5)

# Convert to entrez ID
# returns: enz_uni, enz_up, enz_down, l2fc
enz_list <- entrez_convert(rownames(scObj.cd4@assays$integrated), l2fc_list)

# clusterProfiler
# returns: list(go_red, go_gse_red, reactome, reactome_gse))
# No upregulated genes so using down only
gsea_list <- run_gsea_go(enz_list[3][[1]]$SYMBOL,
                         rownames(scObj.cd4@assays$integrated),
                         enz_list[4][[1]],
                         enz_list[3][[1]]$ENTREZID,
                         enz_list[1][[1]]$ENTREZID)

# Save tables of results
write.table(gsea_list[1][[1]],
            file = file.path(out_fold, "table_cd4_overrep_go_dnapos_v_neg_down_lfc0.5.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gsea_list[2][[1]],
            file = file.path(out_fold, "table_cd4_gsea_go_dnapos_v_neg_lfc0.5.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#############################
## ORA - Gene Ontology     ##
## FIGURE 3E               ##
#############################
plot_gsea_go(gsea_list[2][[1]],
             gsea_list[1][[1]],
             gsea_list[4][[1]],
             gsea_list[3][[1]],
             filename_end = "cd4_dnapos_v_neg_down_lfc0.5.pdf")