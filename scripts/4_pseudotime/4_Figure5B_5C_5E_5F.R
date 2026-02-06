###############################################################################
##   Script to generate Figures 5B, 5E and 5F                                ##
##   Smoothed gene expression profiles of top12 genes from global,           ##
##   lineage 3 and lineage 5 association tests                               ##
##   Input:                                                                  ##
##      master obj: cd4_tcells_hpv_abseq.rds                                 ##
##      sce_cell: slingshot_traj_analysis_celltype_CD4_naive.RDS             ##
##      sceGAM: fit_gam_no_cond_knot5.RDS                                    ##
###############################################################################

library(Seurat)
library(tidyverse)
library(ggplot2)
library(slingshot)
library(tradeSeq)
library(showtext)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

showtext::showtext_auto(TRUE)      # enable showtext
showtext::showtext_opts(dpi = 300) # rescale txt properly

set.seed(42)

source("utils/4_utils.R")
source("scripts/plotting_setup.R")

#setwd("Set main working directory")
in_fold <- "output/script4"
out_fold <- "output/script4"

# Loading from RDS obj
# integrated object
scObj.cd4 <- readRDS(file.path(out_fold, "../script2/cd4_tcells_hpv_abseq.rds"))
DefaultAssay(scObj.cd4) <- "integrated"
# Set to cell type
Idents(scObj.cd4) <- "celltype"
# Remove ILC2 and gd clusters
scObj.cd4.filt <- subset(scObj.cd4, idents = c("ILC2", "CD4+ VD2 γδ", "CD4+ VD1 γδ"), invert = TRUE)

# Trajectory - slingshot
sce_cell <- readRDS(file.path(in_fold, "slingshot_traj_analysis_celltype_CD4_naive.RDS"))
# DE - TradeSeq
sceGAM <- readRDS(file.path(in_fold, "fit_gam_no_cond_knot5.RDS"))

###############################################
##  GO GSEA for global results from          ##
##  association tests                        ##
##  FIGURE 5C                                ##
## ############################################

#################################     5_diff_exp_gene_along_pseudotime_per_lineage
fname <- "5_diff_exp_gene_along_pseudotime_per_lineage"
dir.create(file.path(out_fold, fname))

# General trends of gene expression across pseudotime
# This can be interpreted as testing whether the average gene expression is significantly changing along pseudotime.
# lineages=TRUE argument calculates results for each lineage separately, as well as global test
rowData(sceGAM)$assocRes <- tradeSeq::associationTest(sceGAM,
                                                      lineages = TRUE,
                                                      contrastType = "consecutive")
assocRes <- rowData(sceGAM)$assocRes

# Add padj values and filter for significant genes per lineage and global using padj/fdr <= 0.05
sigGenes <- assocRes %>%
  dplyr::mutate(padj = p.adjust(pvalue, "fdr")) %>%
  dplyr::filter(padj <= 0.05) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::select(Gene, waldStat, df, pvalue, padj, meanLogFC)

sigGenes1 <- assocRes %>%
  dplyr::mutate(padj_1 = p.adjust(pvalue_1, "fdr")) %>%
  dplyr::filter(padj_1 <= 0.05) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::select(Gene, paste0(c("waldStat", "df", "pvalue", "padj"), "_1"), meanLogFC)

sigGenes2 <- assocRes %>%
  dplyr::mutate(padj_2 = p.adjust(pvalue_2, "fdr")) %>%
  dplyr::filter(padj_2 <= 0.05) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::select(Gene, paste0(c("waldStat", "df", "pvalue", "padj"), "_2"), meanLogFC)

sigGenes3 <- assocRes %>%
  dplyr::mutate(padj_3 = p.adjust(pvalue_3, "fdr")) %>%
  dplyr::filter(padj_3 <= 0.05) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::select(Gene, paste0(c("waldStat", "df", "pvalue", "padj"), "_3"), meanLogFC)

sigGenes4 <- assocRes %>%
  dplyr::mutate(padj_4 = p.adjust(pvalue_4, "fdr")) %>%
  dplyr::filter(padj_4 <= 0.05) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::select(Gene, paste0(c("waldStat", "df", "pvalue", "padj"), "_4"), meanLogFC)

sigGenes5 <- assocRes %>% 
  dplyr::mutate(padj_5 = p.adjust(pvalue_5, "fdr")) %>%
  dplyr::filter(padj_5 <= 0.05) %>%
  tibble::rownames_to_column("Gene") %>%
  dplyr::select(Gene, paste0(c("waldStat", "df", "pvalue", "padj"), "_5"), meanLogFC)

# Save all significant genes (and all results before filtering)
write.table(assocRes, file.path(out_fold, fname, "CD4_naive_tradeseq_diff_exp_along_pseudotime_CD4_naive.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sigGenes,
            file.path(out_fold, fname, "CD4_naive_tradeseq_diff_exp_along_pseudotime_global_padj0.05_sel_columns.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sigGenes1,
            file.path(out_fold, fname, "CD4_naive_tradeseq_diff_exp_along_pseudotime_lin1_padj0.05_sel_columns.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sigGenes2,
            file.path(out_fold, fname, "CD4_naive_tradeseq_diff_exp_along_pseudotime_lin2_padj0.05_sel_columns.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sigGenes3,
            file.path(out_fold, fname, "CD4_naive_tradeseq_diff_exp_along_pseudotime_lin3_padj0.05_sel_columns.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sigGenes4,
            file.path(out_fold, fname, "CD4_naive_tradeseq_diff_exp_along_pseudotime_lin4_padj0.05_sel_columns.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sigGenes5,
            file.path(out_fold, fname, "CD4_naive_tradeseq_diff_exp_along_pseudotime_lin5_padj0.05_sel_columns.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Using results from associationTest (above)
for (lin in c("global")) {
  print(lin)
  var <- ifelse(lin == "global", "sigGenes", paste0("sigGenes", lin))
  # GO
  for (cat in c("CC", "BP", "MF")) {
    print(cat)
    enrichment_go(stats = get(var),
                  lineage = lin,
                  all_genes = rownames(sceGAM),
                  catname = cat,
                  dir = file.path(out_fold, fname),
                  type = "ORA")
  }
}

#################################################
##  Pseudotime values against gene expression  ##
##  per lineage (significant genes (top 12)    ##
##  FIGURE 5B (All lineage)                    ##
##  FIGURE 5E (Lineage 3)                      ##
##  FIGURE 5F (Lineage 5)                      ##
## ##############################################

sigGenes <- sigGenes %>%
  dplyr::rename("waldStat_0" = "waldStat",
                "padj_0" = "padj")
list_genes <- list(sigGenes, sigGenes1, sigGenes2, sigGenes3, sigGenes4, sigGenes5)
names(list_genes) <- c("All lineages", "Lineage1", "Lineage2",  "Lineage3", "Lineage4", "Lineage5")

for (i in 0:length(list_genes)) {
  # Get name of lineage and list of significant genes - filtered for top 12
  lin_name <- names(list_genes)[[i + 1]]
  fil <- list_genes[[i + 1]] %>%
    dplyr::arrange(desc(!!sym(paste0("waldStat_", i)))) %>%
    head(12)
  genes <- fil$Gene
  
  if (length(genes) == 0) {
    next
  }

  # Get wald stat and adj p
  wald_padj <- list_genes[[i + 1]] %>%
    dplyr::select(Gene, "waldStat" = paste0("waldStat_",i), "padj" = paste0("padj_", i)) %>%
    dplyr::mutate(padj = formatC(padj, format = "e", digits = 1),
                  waldStat = round(waldStat, 1)) %>%
    dplyr::filter(Gene %in% genes)
  
  # Run plotting
  plot_multi_gene(genes, lin_name, wald_padj, sce_cell, scObj.cd4.filt, type = "top_12_by_wald_sig_0.05", out_fold)
}
