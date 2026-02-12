###############################################################################
##   Script to generate Figure 5G                                            ##
##   Smoothed gene expression trajectories for diff exp genes                ##
##   identified as significantly differentially expressed                    ##
##   in the lineage 3 condition-specific association test                    ## 
##   Input:                                                                  ##
##      sce_cell: slingshot_traj_analysis_celltype_CD4_naive.RDS             ##
##      sceGAM: fit_gam_hpv_clearance_knot5.RDS                              ##
###############################################################################

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(tradeSeq)
library(SingleCellExperiment)
library(showtext)


showtext::showtext_auto(TRUE)      # enable showtext
showtext::showtext_opts(dpi = 300) # rescale txt properly

set.seed(42)

source("utils/4_utils.R")
source("scripts/plotting_setup.R")

#setwd("Set main working directory")
in_fold <- "output/script4"
out_fold <- "output/script4"

# Loading from RDS obj
# Trajectory - slingshot
sce_cell <- readRDS(file.path(in_fold, "slingshot_traj_analysis_celltype_CD4_naive.RDS"))
# DE - TradeSeq
sceGAM <- readRDS(file.path(in_fold, "fit_gam_hpv_clearance_knot5.RDS"))

#################################     5_diff_exp_gene_along_pseudotime_per_lineage
fname <- "6_diff_exp_gene_along_condition"
dir.create(file.path(out_fold, fname))

rowData(sceGAM)$condTest <- tradeSeq::conditionTest(models = sceGAM,
                                                    lineages = TRUE,
                                                    global = TRUE,
                                                    pairwise = TRUE)
conditionTest <- rowData(sceGAM)$condTest

sigGenes_cond <- conditionTest %>%
  dplyr::mutate(padj = p.adjust(pvalue, "fdr")) %>%
  dplyr::filter(padj <= 0.05) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::select(Gene, "waldStat", "df", "pvalue", "padj")

sigGenes_cond1 <- conditionTest %>%
  dplyr::mutate(padj_lineage1 = p.adjust(pvalue_lineage1, "fdr")) %>%
  dplyr::filter(padj_lineage1 <= 0.05) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::select(Gene, paste0(c("waldStat", "df", "pvalue", "padj"), "_lineage1"))

sigGenes_cond2 <- conditionTest %>%
  dplyr::mutate(padj_lineage2 = p.adjust(pvalue_lineage2, "fdr")) %>%
  dplyr::filter(padj_lineage2 <= 0.05) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::select(Gene, paste0(c("waldStat", "df", "pvalue", "padj"), "_lineage2"))


sigGenes_cond3 <- conditionTest %>%
  dplyr::mutate(padj_lineage3 = p.adjust(pvalue_lineage3, "fdr")) %>%
  dplyr::filter(padj_lineage3 <= 0.05) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::select(Gene, paste0(c("waldStat", "df", "pvalue", "padj"), "_lineage3"))


sigGenes_cond4 <- conditionTest %>% 
  dplyr::mutate(padj_lineage4 = p.adjust(pvalue_lineage4, "fdr")) %>%
  dplyr::filter(padj_lineage4 <= 0.05) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::select(Gene, paste0(c("waldStat", "df", "pvalue", "padj"), "_lineage4"))

sigGenes_cond5 <- conditionTest %>% 
  dplyr::mutate(padj_lineage5 = p.adjust(pvalue_lineage5, "fdr")) %>%
  dplyr::filter(padj_lineage5 <= 0.05) %>%
  tibble::rownames_to_column("Gene")

# Save all significant genes (and all results before filtering)
write.table(conditionTest,
  file.path(out_fold, "6_diff_exp_gene_along_condition", "CD4_naive_tradeseq_diff_exp_along_pseudotime_CD4_naive_conditionTest.txt"),
            sep = "\t", quote = FALSE, row.names = TRUE)
write.table(sigGenes_cond,
  file.path(out_fold,  "6_diff_exp_gene_along_condition", "CD4_naive_tradeseq_diff_exp_along_pseudotime_global_padj0.05_sel_columns_conditionTest.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sigGenes_cond1,
  file.path(out_fold, "6_diff_exp_gene_along_condition", "CD4_naive_tradeseq_diff_exp_along_pseudotime_lin1_padj0.05_sel_columns_conditionTest.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sigGenes_cond2,
  file.path(out_fold, "6_diff_exp_gene_along_condition", "CD4_naive_tradeseq_diff_exp_along_pseudotime_lin2_padj0.05_sel_columns_conditionTest.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sigGenes_cond3,
  file.path(out_fold, "6_diff_exp_gene_along_condition", "CD4_naive_tradeseq_diff_exp_along_pseudotime_lin3_padj0.05_sel_columns_conditionTest.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sigGenes_cond4,
  file.path(out_fold, "6_diff_exp_gene_along_condition", "CD4_naive_tradeseq_diff_exp_along_pseudotime_lin4_padj0.05_sel_columns_conditionTest.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sigGenes_cond5,
  file.path(out_fold, "6_diff_exp_gene_along_condition", "CD4_naive_tradeseq_diff_exp_along_pseudotime_lin5_padj0.05_sel_columns_conditionTest.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

####################################################

curves_colors <- list(
  "lin1" = c("#FDBD13", "#1857F4", rep("#FFFFFF00", 8)),
  "lin2" = c(rep("#FFFFFF00", 2), "#FDBD13", "#1857F4", rep("#FFFFFF00", 6)),
  "lin3" = c(rep("#FFFFFF00", 4), "#FDBD13", "#1857F4", rep("#FFFFFF00", 4)),
  "lin4" = c(rep("#FFFFFF00", 6), "#FDBD13", "#1857F4", rep("#FFFFFF00", 2)),
  "lin5" = c(rep("#FFFFFF00", 8), "#FDBD13", "#1857F4"))

################################################
## Pseudotime values against gene expression  ##
## per condition across lineage 3             ##
## FIGURE 5G                                  ##
## #############################################
cols <- curves_colors[[3]]
unique_pos <- unique(sigGenes_cond3$Gene)

plt_pos <- list()
for (gene in unique_pos) {
  stats <- sigGenes_cond3 %>% dplyr::filter(Gene == gene) %>% convert_padj_wald_annotation()
  plt_pos[[gene]] <- get_smoothers_per_condition_for_one_lineage_with_cond_stats(model = sceGAM,
                                                                                 gene = gene,
                                                                                 cols = cols,
                                                                                 lineage = 3,
                                                                                 stat = stats)
}

a <- ggpubr::ggarrange(
      plt_pos[[1]],
      plt_pos[[2]],
      plt_pos[[3]],
      plt_pos[[4]],
      plt_pos[[5]],
      plt_pos[[6]],
      plt_pos[[7]],
      common.legend = TRUE, ncol = 4, nrow = 2)
b <- ggpubr::ggarrange(
      plt_pos[[8]],
      plt_pos[[9]],
      plt_pos[[10]],
      plt_pos[[11]],
      plt_pos[[12]],
      plt_pos[[13]],
      plt_pos[[14]],
      common.legend = TRUE, ncol = 4, nrow = 2)

saveRDS(a,
  file.path(out_fold, fname, sprintf("Figure_5G_CD4_naive_lineage%s_per_condition_diff_expressed_genes_Oct2025_4x3_part1_conditionTest.rds", "3")))
tiff(file.path(out_fold, fname, sprintf("Figure_5G_CD4_naive_lineage%s_per_condition_diff_expressed_genes_Oct2025_4x3_part1_conditionTest.tiff", "3")),
     width = 12, height = 12,
     units = "in",
     res = 300,
     compression = "lzw")
print(a)
dev.off()

saveRDS(b,
  file.path(out_fold, fname, sprintf("Figure_5G_CD4_naive_lineage%s_per_condition_diff_expressed_genes_Oct2025_4x3_part2_conditionTest.rds", "3")))
tiff(file.path(out_fold, fname, sprintf("Figure_5G_CD4_naive_lineage%s_per_condition_diff_expressed_genes_Oct2025_4x3_part2_conditionTest.tiff", "3")),
     width = 12, height = 12,
     units = "in",
     res = 300,
     compression = "lzw")
print(b)
dev.off()
