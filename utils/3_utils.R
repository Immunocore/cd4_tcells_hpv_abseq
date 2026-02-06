###################################################
##   Utils script with                           ##
##   helper functions for section 3              ##
###################################################

#' Output UMAP dimension plot to TIFF file
#'
#' @param scobj Seurat object
#' @param out_fold Output folder path
#' @param filename Output filename
#' @param group Column name in metadata to group by
#' @param reduction Dimensionality reduction to use
#' @param title Plot title
#' @param label_box Logical, draw label boxes
#' @param a Alpha transparency (not used in current implementation)
#' @param label Logical, show labels
#' @param label_size Label font size
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param cols Colour vector for groups
#'
#' @return Prints plot, saves TIFF file, returns nothing
output_dimplot <- function(scobj, out_fold, filename, group,
                           reduction = "wnn.umap", title = "",
                           label_box = TRUE, a = 1, label = TRUE, label_size = 5,
                           width = 11, height = 10,
                           cols = unname(IMC_colors)) {
  
  tiff(file.path(out_fold, filename),
       width = width, height = height,
       res = 300,
       compression = "lzw")
  plot <- Seurat::DimPlot(scobj, reduction = reduction, label = label,
                          group.by = group, label.size = label_size,
                          label.box = label_box, repel = TRUE,
                          cols = cols) +
    labs(x = "UMAP1",
         y = "UMAP2",
         title = title)
  print(plot)
  dev.off()
  
}

#' Create marker gene bubble plot for RNA and AbSeq assays
#'
#' @param scobj Seurat object
#' @param out_fold Output folder path
#' @param table_filename Filename for marker table output
#' @param plot_filename Filename for plot output
#' @param wid Plot width in inches
#' @param hei Plot height in inches
#'
#' @return Combined dot plot object
marker_bubble <- function(scobj, out_fold, table_filename,
                          plot_filename, wid = 20, hei = 8) {
  # Get cluster marker genes
  
  markers <- Seurat::FindAllMarkers(scobj, only.pos = TRUE)
  markers_filt <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::filter(p_val_adj <= 0.05) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    dplyr::top_n(n = 2, avg_log2FC) %>%
    dplyr::mutate("type" = "RNA") %>%
    dplyr::mutate(gene_short = gene) %>%
    dplyr::mutate(gene_short = stringr::str_replace(gene_short, "-refseq", ""))
  
  markers_ab <- Seurat::FindAllMarkers(scobj, only.pos = TRUE, assay = "IAB")
  markers_abfilt <- markers_ab %>%
    dplyr::group_by(cluster) %>%
    dplyr::filter(p_val_adj <= 0.05) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    dplyr::top_n(n = 2, avg_log2FC) %>%
    dplyr::mutate("type" = "Ab") %>%
    dplyr::mutate(gene2 = stringr::str_replace(gene, "HLA-DR", "HLADR")) %>%
    dplyr::mutate(gene2 = stringr::str_replace(gene2, "M1-70", "M170")) %>%
    tidyr::separate(gene2, c(NA, "gene_short", NA), sep = "-", remove = TRUE)
  
  all_markers <- rbind(markers_abfilt, markers_filt)
  write.table(all_markers, file.path(out_fold, table_filename),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  ab_plot <- scCustomize::DotPlot_scCustom(scobj,
                                           features = unique(markers_abfilt$gene),
                                           assay = "IAB",
                                           remove_axis_titles = FALSE) +
    ggplot2::xlab("Protein(AbSeq)") +
    ggplot2::ylab("Cluster") +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_x_discrete(labels = unique(markers_abfilt$gene_short)) +
    RotatedAxis() +
    ggplot2::scale_colour_gradient(low = "#ffffff", high = "#206ba3")

  rna_plot <- scCustomize::DotPlot_scCustom(scobj,
                                            features = unique(markers_filt$gene),
                                            remove_axis_titles = FALSE) +
    ggplot2::scale_x_discrete(labels = unique(markers_filt$gene_short)) +
    RotatedAxis() +
    ggplot2::xlab("mRNA") +
    ggplot2::ylab("") +
    ggplot2::scale_colour_gradient(low = "#ffffff", high = "#206ba3")
  
  dot_plots <- ggpubr::ggarrange(ab_plot, rna_plot, widths = c(0.5, 1))
  
  cairo_pdf(file.path(out_fold, plot_filename), width = wid, height = hei)
  print(dot_plots)
  dev.off()
  
  return(dot_plots)
}

#' Perform differential expression analysis between two groups within clusters
#'
#' @param scobj Seurat object
#' @param colname Column name in metadata to set as identity
#' @param id1 First group identifier
#' @param id2 Second group identifier
#' @param out_fold Output folder path
#' @param filename Output filename for results table
#'
#' @return Data frame with differential expression results for RNA and AbSeq assays
clust_de <- function(scobj, colname, id1, id2, out_fold, filename) {
  
  Idents(scobj) <- colname
  
  # RNA
  de.out <- Seurat::FindMarkers(scobj, ident.1 = id1, ident.2 = id2,
                                verbose = FALSE,
                                logfc.threshold = 0.25, min.pct = 0.1)
  de.out["gene"] <- rownames(de.out)
  # Ab
  de.out.ab <- Seurat::FindMarkers(scobj, ident.1 = id1, ident.2 = id2,
                                   verbose = FALSE, assay = "IAB",
                                   logfc.threshold = 0.25, min.pct = 0.1)
  de.out.ab["gene"] <- rownames(de.out.ab)
  
  # join
  de.out.all <- rbind(de.out, de.out.ab)
  
  write.table(de.out.all, file.path(out_fold, filename), quote = FALSE, row.names = FALSE, sep = "\t")
  
  return(de.out.all)
  
}

#' Create volcano plot of differential expression results
#'
#' @param table Data frame with differential expression results (expects 'gene', 'avg_log2FC', 'p_val_adj' columns)
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param height Plot height in inches
#' @param width Plot width in inches
#' @param out_fold Output folder path
#' @param filename Output filename
#' @param arrows Logical, draw connectors to significant points
#'
#' @return Saves PDF file, returns nothing
plot_volc <- function(table, title, subtitle, height, width,
                      out_fold, filename, arrows = TRUE) {
  
  pdf(file.path(out_fold, filename), width = width, height = height)
  print(EnhancedVolcano::EnhancedVolcano(table, lab = table$gene,
                                         x = "avg_log2FC", y = "p_val_adj",
                                         title = title, subtitle = subtitle,
                                         drawConnectors = arrows,
                                         widthConnectors = 0.75,
                                         col = c("black", "black",
                                                 "black", "#1857F4")))
  dev.off()
}


##########################################
#### Gene overrepresentation analysis ####
##########################################

#' Prepare log fold change vectors for GSEA analysis
#'
#' @param datatab Data frame with differential expression results 
#'                (must have 'avg_log2FC' column)
#' @param lfc Log fold change threshold
#'
#' @return List: all LFCs (sorted), up LFCs, down LFCs
prep_lfcs <- function(datatab, lfc) {
  
  l2fc <- datatab$avg_log2FC
  names(l2fc) <- row.names(datatab)
  l2fc <- sort(l2fc, decreasing = TRUE)
  l2fc_up <- l2fc[l2fc > lfc]
  l2fc_down <- l2fc[l2fc < -lfc]
  
  return(list(l2fc, l2fc_up, l2fc_down))
}

#' Convert gene symbols to Entrez IDs for GSEA analysis
#'
#' @param backgr_in Background gene list (symbols)
#' @param prep_lfcs_out Output from prep_lfcs function
#' @param from Source ID type (default: "SYMBOL")
#' @param to Target ID type (default: "ENTREZID")
#' @param org Organism database (default: "org.Hs.eg.db")
#'
#' @return List: background Entrez IDs, upregulated Entrez IDs,
#'         downregulated Entrez IDs, named LFC vector with Entrez IDs
entrez_convert <- function(backgr_in, prep_lfcs_out,
                           from = "SYMBOL", to = "ENTREZID",
                           org = "org.Hs.eg.db") {
  
  l2fc <- prep_lfcs_out[1][[1]]
  
  # background  
  enz_uni <- clusterProfiler::bitr(backgr_in,
                                   fromType = from, toType = to, OrgDb = org)
  
  # Convert to entrez ID
  enz_up <- clusterProfiler::bitr(names(prep_lfcs_out[2][[1]]),
                                  fromType = from, toType = to, OrgDb = org)
  enz_down <- clusterProfiler::bitr(names(prep_lfcs_out[3][[1]]),
                                    fromType = from, toType = to, OrgDb = org)
  
  enz_both <- clusterProfiler::bitr(names(l2fc),
                                    fromType = from, toType = to, OrgDb = org)
  names(l2fc) <- enz_both$ENTREZID
  l2fc <- l2fc[!is.na(names(l2fc))]
  
  return(list(enz_uni, enz_up, enz_down, l2fc))
}

#' Run GSEA and GO enrichment analysis
#'
#' @param genelist_sym Gene list as symbols
#' @param backgr_sym Background gene list as symbols
#' @param lfc_enz Named log fold change vector with Entrez IDs
#' @param genelist_enz Gene list as Entrez IDs
#' @param backgr_enz Background gene list as Entrez IDs
#' @param org Organism database object
#' @param ont GO ontology (default: "BP")
#' @param qcut Q-value cutoff (default: 0.05)
#' @param pcut P-value cutoff (default: 0.05)
#'
#' @return List: GO ORA; GO GSEA, Reactome ORA, Reactome GSEA
run_gsea_go <- function(genelist_sym, backgr_sym, lfc_enz,
                        genelist_enz, backgr_enz,
                        org = org.Hs.eg.db, ont = "BP",
                        qcut = 0.05, pcut = 0.05) {
  
  go <- clusterProfiler::enrichGO(genelist_sym, OrgDb = org, keyType = "SYMBOL",
                                  ont = ont, pvalueCutoff = pcut,
                                  qvalueCutoff = qcut,
                                  universe = backgr_sym)
  go_gse <- clusterProfiler::gseGO(lfc_enz, ont = ont, keyType = "ENTREZID",
                                   OrgDb = org, pvalueCutoff = pcut)
  
  # Reduce term redundancy
  go_red <- clusterProfiler::simplify(go)
  go_gse_red <- clusterProfiler::simplify(go_gse)
  
  # ReactomePA
  reactome <- ReactomePA::enrichPathway(genelist_enz, organism = "human",
                                        pvalueCutoff = pcut,
                                        qvalueCutoff = qcut,
                                        universe = backgr_enz)
  reactome_gse <- ReactomePA::gsePathway(lfc_enz, pvalueCutoff = pcut)
  
  return(list(go_red, go_gse_red, reactome, reactome_gse))
}

#' Plot GSEA and GO enrichment results
#'
#' @param gsea_go GSEA GO enrichment result object
#' @param overrep_go ORA GO enrichment result object
#' @param gsea_react GSEA Reactome enrichment result object
#' @param overrep_react ORA Reactome enrichment result object
#' @param filename_end Suffix for output filenames
#'
#' @return Saves PDF files, returns nothing
plot_gsea_go <- function(gsea_go, overrep_go,
                         gsea_react, overrep_react, filename_end) {
  
  # Plot results
  if (nrow(as.data.frame(gsea_go)) > 0) {
    dot_go_gse <- enrichplot::dotplot(gsea_go, showCategory = 15, split = ".sign")
    go_graph2 <- enrichplot::goplot(gsea_go)
    go_graph3 <- enrichplot::goplot(gsea_go)
    
    pdf(file.path(out_fold, paste0("dotplot_gsea_go_", filename_end)))
    print(dot_go_gse)
    dev.off()
    pdf(file.path(out_fold, paste0("graph_gsea_go_", filename_end)))
    print(go_graph2)
    dev.off()
    pdf(file.path(out_fold, paste0("graph_gsea_go2_", filename_end)),
        width = 15)
    print(go_graph3)
    dev.off()
  }
  
  if (nrow(as.data.frame(overrep_go)) > 0) {
    dot_go <- enrichplot::dotplot(overrep_go, showCategory = 15)
    
    pdf(file.path(out_fold, paste0("dotplot_overrep_go_", filename_end)))
    print(dot_go)
    dev.off()
  }
  
  if (nrow(as.data.frame(gsea_react)) > 0) {
    dot_react_gse <- enrichplot::dotplot(gsea_react, showCategory = 15,
                                         split = ".sign")
    
    pdf(file.path(out_fold, paste0("dotplot_gsea_react_", filename_end)))
    print(dot_react_gse)
    dev.off()
  }
  
  if (nrow(as.data.frame(overrep_react)) > 0) {
    dot_react <- enrichplot::dotplot(overrep_react, showCategory = 15)
    
    pdf(file.path(out_fold, paste0("dotplot_overrep_react_", filename_end)))
    print(dot_react)
    dev.off()
  }
}


#' Create data frame with cell type proportions per patient
#'
#' @param scObj Seurat object with 'Subject' and cell type identities
#'
#' @return Data frame with columns: Subject, celltype, count
create_patient_proportions_data <- function(scObj) {
  DefaultAssay(scObj) <- "integrated"
  celltypes <- levels(Idents(scObj))
  patients <- unique(scObj$Subject)
  
  data <- list()
  for (cell in celltypes) {
    sub <- subset(scObj, idents = cell)
    
    for (pat in patients) {
      
      count <- sub@meta.data %>%
        dplyr::filter(Subject == pat) %>%
        nrow()
      
      d <- data.frame(Subject = pat, celltype = cell, count = count)
      data <- c(data, list(d))
    }
    
  }
  return(purrr::list_rbind(data))
}

#' Run GSEA and GO analysis per cluster with up/downregulated genes
#'
#' @param dif_exp Data frame with differential expression results
#' @param adj_p Adjusted p-value threshold for filtering
#' @param lfc Log fold change threshold
#' @param back_enz Background gene list as Entrez IDs
#' @param back_sym Background gene list as symbols
#' @param name Name prefix for output files
#'
#' @return Saves PDF files via plot_gsea_go, returns nothing
per_clust_gsea_go <- function(dif_exp, adj_p, lfc, back_enz, back_sym, name) {
  
  dif_exp_f <- dif_exp %>%
    dplyr::filter(p_val_adj <= adj_p) %>%
    dplyr::arrange(desc(avg_log2FC))
  
  # Get entrez ids
  enz_both <- clusterProfiler::bitr(dif_exp_f$gene,
                                    fromType = "SYMBOL", toType = "ENTREZID",
                                    OrgDb = "org.Hs.eg.db")
  l2fc <- dif_exp_f$avg_log2FC
  names(l2fc) <- enz_both$ENTREZID
  l2fc <- l2fc[!is.na(names(l2fc))]
  
  # upregulated genes only
  up <- dif_exp_f %>%
    dplyr::filter(avg_log2FC > lfc)
  enz_up <- clusterProfiler::bitr(up$gene,
                                  fromType = "SYMBOL", toType = "ENTREZID",
                                  OrgDb = "org.Hs.eg.db")
  
  # downregulated genes only
  down <- dif_exp_f %>%
    dplyr::filter(avg_log2FC < lfc)
  enz_down <- clusterProfiler::bitr(down$gene,
                                    fromType = "SYMBOL", toType = "ENTREZID",
                                    OrgDb = "org.Hs.eg.db")
  
  # Run GSEA and GO
  gsea_list_up <- run_gsea_go(up$gene, back_sym,
                              l2fc, enz_up$ENTREZID, back_enz)

  gsea_list_down <- run_gsea_go(down$gene, back_sym,
                                l2fc, enz_down$ENTREZID, back_enz)

  # Plots
  plot_gsea_go(gsea_list_up[2][[1]], gsea_list_up[1][[1]],
               gsea_list_up[4][[1]], gsea_list_up[3][[1]],
               filename_end = paste0(name, "_lfc0.5.pdf"))
  # Plots
  plot_gsea_go(gsea_list_down[2][[1]], gsea_list_down[1][[1]],
               gsea_list_down[4][[1]], gsea_list_down[3][[1]],
               filename_end = paste0(name, "_lfc-0.5.pdf"))
  
}