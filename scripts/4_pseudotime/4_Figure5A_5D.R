###############################################################################
##   Script to generate Figures 5A and 5D                                    ##
##   Input:                                                                  ##
##      sce_cell: slingshot_traj_analysis_celltype_CD4_naive.RDS             ##
##      sceGAM: fit_gam_no_cond_knot5.RDS                                    ##
###############################################################################

library(tidyverse)
library(RColorBrewer)
library(grDevices)
library(slingshot)
library(showtext)

showtext::showtext_auto(TRUE)      # enable showtext
showtext::showtext_opts(dpi = 300) # rescale txt properly

set.seed(42)

source("scripts/plotting_setup.R")
out_fold <- "output/script4"
fname <- "figures"
dir.create(file.path(out_fold, fname), showWarnings = FALSE)

# Loading from RDS obj
# Trajectory - slingshot
sce_cell <- readRDS(file.path(out_fold, "slingshot_traj_analysis_celltype_CD4_naive.RDS"))
# DE - TradeSeq
sceGAM <- readRDS(file.path(out_fold, "fit_gam_no_cond_knot5.RDS"))

###############################################
## Plots from celltype traj analysis         ##
## FIGURE 5A                                 ##
## ############################################

# For each cell: take the minimum non-NA pseudotime and assign a colour
pt <- slingPseudotime(sce_cell)
global_pt <- apply(pt, 1, min, na.rm = TRUE)
global_pt[is.infinite(global_pt)] <- NA

colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral")[-6])(100)
plotcol <- colors[cut(global_pt, breaks = 100)] 

tiff(file.path(out_fold, fname, "Figure_5A_pseudotime_cd4_celltype_root_naive_CD4_annotated.tiff"),
  width = 8, height = 8, units = "in", res = 300, compression = "lzw")
plot(reducedDims(sce_cell)$WNN.UMAP, col = plotcol,
     pch = 16, asp = 1, cex.lab = 1.5, cex.axis = 2,
     xlab = "UMAP1", ylab = "UMAP2")
lines(slingshot::SlingshotDataSet(sce_cell),
      linInd = 1, type = "lineage", lty = "dotted")
text(x = 10, y = 0.5, labels = "L1", cex = 2, col = "black")
lines(slingshot::SlingshotDataSet(sce_cell),
      linInd = 2, type = "lineage", lty = "dashed")
text(x = 4.7, y = -3.5, labels = "L2", cex = 2, col = "black")
lines(slingshot::SlingshotDataSet(sce_cell),
      linInd = 3, type = "lineage", col = "black")
text(x = 2, y = 3, labels = "L3", cex = 2, col = "black")
lines(slingshot::SlingshotDataSet(sce_cell),
      linInd = 4, type = "lineage", lty = "twodash", col = "black")
text(x = -1.2, y = -5.6, labels = "L4", cex = 2, col = "black")
lines(slingshot::SlingshotDataSet(sce_cell),
      linInd = 5, type = "lineage", lty = "dotdash", col = "black")
text(x = -4.4, y = 1.9, labels = "L5", cex = 2, col = "black")
p <- recordPlot()
saveRDS(p, file.path(out_fold, fname, "Figure_5A_pseudotime_cd4_celltype_root_naive_CD4_annotated.rds"))
dev.off()


###############################################
##  Cell proportions across pseudotime       ##
##  bins                                     ##
##  FIGURE 5D                                ##
## ############################################
pseudo_tab_cell <- as.data.frame(slingPseudotime(sce_cell)) %>%
  tibble::rownames_to_column("cell")
pseudo_tab_cell$celltype <- colData(sce_cell)$celltype

pseudo_tab_cell$hpv_clearance <- colData(sce_cell)$hpv_clearance

cols_dna <- IMC_colors[1:2]
names(cols_dna) <- c("Group 1", "Group 2")

pseudo_cell_hpv_long <- tidyr::pivot_longer(pseudo_tab_cell,
                                            cols = starts_with("Lineage"),
                                            names_to = "trajectory",
                                            values_to = "pseudotime") %>%
                        dplyr::filter(!is.na(pseudotime))

scaling <- pseudo_cell_hpv_long %>%
  dplyr::group_by(hpv_clearance, trajectory) %>%
  dplyr::summarise(n = n())
scaling2 <- pseudo_cell_hpv_long %>%
  dplyr::group_by(trajectory) %>%
  dplyr::summarise(n_traj = n())
scaling <- dplyr::left_join(scaling, scaling2) %>%
  dplyr::mutate(prop = ((n_traj - n) / n_traj) * 2)

num_bins <- 15

pseudo_cell_hpv_long <- pseudo_cell_hpv_long %>%
  dplyr::mutate(bins = cut(pseudotime,
                breaks = num_bins,
                include.lowest = TRUE,
                labels = as.character(1:num_bins)))
# scale based on total number of cells in HPV category by lineage
pseudo_binned_cellty_hpv <- pseudo_cell_hpv_long %>%
  dplyr::group_by(trajectory, hpv_clearance, celltype, bins) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::arrange(as.numeric(bins)) %>%
  dplyr::mutate(x_lab = paste0(bins, " - ", hpv_clearance)) %>%
  dplyr::left_join(scaling %>%
  dplyr::select(hpv_clearance, trajectory, prop)) %>%
  dplyr::mutate(n_scale = n * prop)

bar_celltype_bin2 <- ggplot2::ggplot(pseudo_binned_cellty_hpv,
                                     aes(x = factor(x_lab, unique(pseudo_binned_cellty_hpv$x_lab)),
                                         y = n_scale,
                                         fill = celltype)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = IMC_colors) +
  ggplot2::facet_grid(trajectory ~ .) +
  gg_theme_specs(fs = 16) +
  ggplot2::scale_x_discrete(guide = guide_axis(angle = 45))

# Create blanks
blanks <- pseudo_binned_cellty_hpv %>%
  dplyr::ungroup() %>%
  dplyr::select(-c(n_scale, x_lab, hpv_clearance)) %>%
  dplyr::mutate(n = 0) %>%
  dplyr::mutate(x_lab = bins) %>%
  dplyr::mutate(hpv_clearance = "") %>%
  unique()

pseudo_binned_cellty_hpv2_L3_5 <- as.data.frame(rbind(pseudo_binned_cellty_hpv, blanks)) %>%
  dplyr::mutate(alpha = dplyr::case_when(hpv_clearance == "Group 2" ~ 1,
                                         hpv_clearance == "Group 1" ~ 1,
                                         .default = 0)) %>%
  dplyr::arrange(as.numeric(bins), trajectory, celltype) %>%
  dplyr::filter(trajectory %in% c("Lineage3", "Lineage5"))

# order of x labs
ord1 <- as.numeric(unique(pseudo_binned_cellty_hpv2_L3_5$bins))
ord2 <- unique(pseudo_binned_cellty_hpv2_L3_5$hpv_clearance)
t <- expand.grid(ord1 = ord1, ord2 = ord2)
t <- t[order(t$ord1), ]
ord <- paste(t$ord1, t$ord2, sep = " - ")
ord <- stringr::str_replace_all(ord, " - $", "")

label_naming <- pseudo_binned_cellty_hpv2_L3_5 %>%
  dplyr::select(x_lab, bins) %>%
  dplyr::group_by(bins) %>%
  dplyr::mutate(x_lab = gsub("^[0-9]+$", " ", x_lab)) %>%
  dplyr::distinct()

# Plot
bar_celltype_bin3 <- ggplot2::ggplot(pseudo_binned_cellty_hpv2_L3_5,
                                     aes(x = factor(x_lab, ord),
                                         y = n_scale,
                                         fill = celltype),
                                     alpha = alpha) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::scale_fill_manual(values = IMC_colors) +
  ggplot2::facet_grid(trajectory ~ .) +
  gg_theme_specs(fs = 16) +
  ggplot2::scale_x_discrete(guide = guide_axis(angle = 45),
                            labels = label_naming$x_lab) +
  ggplot2::xlab("Pseudotime bin - category") +
  ggplot2::ylab("Scaled number of cells") +
  ggplot2::labs(fill = "")

saveRDS(bar_celltype_bin3,
  file.path(out_fold, fname, "Figure_5D_pseudotime_binned_boxplot_cd4_together_root_naive_CD4_only_L3_L5.rds"))
tiff(file.path(out_fold, fname, "Figure_5D_pseudotime_binned_boxplot_cd4_together_root_naive_CD4_only_L3_L5.tiff"),
     width = 12,
     height = 8,
     units = "in",
     res = 300,
     compression = "lzw")
print(bar_celltype_bin3)
dev.off()
