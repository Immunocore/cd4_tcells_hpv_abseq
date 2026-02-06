###############################################################
##   Script to carry out trajectory analysis on CD4+ object  ##
##   with slinghot and tradeSeq                              ##
##   Input:                                                  ##
##      master rds object: cd4_tcells_hpv_abseq.rds          ##
###############################################################

# Load packages and required objects and rename cell clusters
library(scran)
library(Seurat)
library(tidyverse)
library(slingshot)
library(tradeSeq)
library(Matrix) # for drop0 command
library(BiocParallel)
library(showtext)

showtext::showtext_auto(TRUE)      # enable showtext
showtext::showtext_opts(dpi = 300) # rescale txt properly

set.seed(42)

#setwd("Set main working directory")
source("utils/4_utils.R")

in_fold <- "output/script2"
out_fold <- "output/script4"

# Load master object with filtered and renamed cells and metadata
scObj.cd4 <- readRDS(file.path(in_fold, 'cd4_tcells_hpv_abseq.rds'))

# Set to cell type
Idents(scObj.cd4) <- "celltype"

# Remove ILC2 and gd clusters
scObj.cd4.filt <- subset(scObj.cd4, idents = c("ILC2", "CD4+ VD2 γδ", "CD4+ VD1 γδ"), invert = TRUE)

# Convert data to singlecellexperiment
sce <- as.SingleCellExperiment(scObj.cd4.filt, assay = "integrated")

#### Trajectory analysis on celltypes ####
if (!file.exists(file.path(out_fold, "slingshot_traj_analysis_celltype_CD4_naive.RDS"))) {
  sce_cell <- sce
  sce_cell <- slingshot::slingshot(sce_cell, reducedDim = 'WNN.UMAP', clusterLabels = colData(sce_cell)$celltype,
                        approx_points = 150, start.clus = 'CD4+ Naive')
  saveRDS(sce_cell, file.path(out_fold, "slingshot_traj_analysis_celltype_CD4_naive.RDS"))
} else {
  message("Loading from RDS obj...")
  sce_cell <- readRDS(file.path(out_fold, "slingshot_traj_analysis_celltype_CD4_naive.RDS"))

}

#################################################################################

# Select cells and genes to use - based on the integrated assay
sel_cells <- split(Cells(scObj.cd4.filt), Idents(scObj.cd4.filt))

# sampling 50% of each cell type population - might still be too big?
sel_cells <- unlist(lapply(sel_cells, function(x) {
  set.seed(1)
  return(sample(x, size = ceiling(length(x)*0.5)))
}))

# Get gene variance using RNA assay - select top 500
gv <- as.data.frame(na.omit(scran::modelGeneVar(scObj.cd4.filt@assays$RNA@data[ ,sel_cells])))
gv <- gv[order(gv$bio, decreasing = TRUE), ]
sel_genes <- sort(rownames(gv)[1:500])

###################################
#### Conditions - hpv clearance ###
###################################
if (!file.exists(file.path(out_fold, "evaluatek_500genes_3to8_results_hpv_clearance.RDS"))) {
  kn_cond <- tradeSeq::evaluateK(
    counts = drop0(scObj.cd4.filt@assays$RNA@counts[sel_genes, sel_cells]),
    # Pseudotime and cellWeights from the previous slingshot analysis
    pseudotime = slingshot::slingPseudotime(sce_cell, na = FALSE)[sel_cells, ],
    cellWeights = slingshot::slingCurveWeights(sce_cell)[sel_cells, ],
    conditions = factor(scObj.cd4.filt[, sel_cells]@meta.data %>% dplyr::select(hpv_clearance) %>% pull()),
    nGenes = 500,
    BPPARAM = BiocParallel::MulticoreParam(workers = 16),
    k = 3:8
  )
  saveRDS(kn_cond, file.path(out_fold, "evaluatek_500genes_3to8_results_hpv_clearance.RDS"))
} else {
  message(sprintf("evaluateK executed"))
}

# Execute fitGAM with different number of knots
for (k in c(5, 6, 7, 8)) {
 if (!file.exists(file.path(out_fold, sprintf("fit_gam_hpv_clearance_knot%d.RDS", k)))) {
    print(sprintf("Executing nknots=%d...", k))
    time_st <- Sys.time()
    sceGAM_cond <- tradeSeq::fitGAM(
      # Counts from RNA assay - using raw counts, because the package needs raw
      counts = drop0(scObj.cd4.filt@assays$RNA@counts[sel_genes, sel_cells]),
      # Pseudotime and cellWeights from the previous slingshot analysis
      pseudotime = slingshot::slingPseudotime(sce_cell, na = FALSE)[sel_cells, ],
      cellWeights = slingshot::slingCurveWeights(sce_cell)[sel_cells, ],
      conditions = factor(scObj.cd4.filt[, sel_cells]@meta.data %>% dplyr::select(hpv_clearance) %>% pull()),
      # Should we check the nknots parameter?
      nknots = k, verbose = TRUE, parallel = TRUE, sce = TRUE,
      BPPARAM = BiocParallel::MulticoreParam(workers=16)
    )
    saveRDS(sceGAM_cond, file.path(out_fold, sprintf("fit_gam_hpv_clearance_knot%d.RDS", k)))
    time_end = Sys.time()
    print(time_end-time_st)
  } else {
    message(sprintf("fitGAM with nknots %d executed", k))
  } 
}

###################################
#### No Conditions - cell types ###
###################################
if (!file.exists(file.path(out_fold, "evaluatek_500genes_3to8_results_no_cond.RDS"))) {
  kn <- tradeSeq::evaluateK(
    counts = drop0(scObj.cd4.filt@assays$RNA@counts[sel_genes, sel_cells]),
    # Pseudotime and cellWeights from the previous slingshot analysis
    pseudotime = slingshot::slingPseudotime(sce_cell, na = FALSE)[sel_cells, ],
    cellWeights = slingshot::slingCurveWeights(sce_cell)[sel_cells, ],
    nGenes = 500,
    BPPARAM = BiocParallel::MulticoreParam(workers=16),
    k = 3:8
  )
  saveRDS(kn, file.path(out_fold, "evaluatek_500genes_3to8_results_no_cond.RDS"))
} else {
  message(sprintf("evaluateK executed"))
}

# Fit model
if (!file.exists(file.path(out_fold, "fit_gam_no_cond_knot5.RDS"))) {
  sceGAM <- tradeSeq::fitGAM(
    # Counts from RNA assay - using raw counts, because the package needs raw
    counts = drop0(scObj.cd4.filt@assays$RNA@counts[sel_genes, sel_cells]),
    # Pseudotime and cellWeights from the previous slingshot analysis
    pseudotime = slingshot::slingPseudotime(sce_cell, na = FALSE)[sel_cells, ],
    cellWeights = slingshot::slingCurveWeights(sce_cell)[sel_cells, ],
    # Should we check the nknots parameter?
    nknots = 5, verbose = TRUE, parallel = TRUE, sce = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(workers=16)
  )
  # Save as rds - because it takes so long to run
  saveRDS(sceGAM, file.path(out_fold, "fit_gam_no_cond_knot5.RDS"))
} else {
  message("No condition processed")
}
