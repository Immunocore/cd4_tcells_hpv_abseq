###################################################
##   Utils script with                           ##
##   helper functions for section 4              ##
###################################################

source("scripts/plotting_setup.R")

#' Plot multiple genes along pseudotime for specified lineage
#'
#' @param genes Character vector of gene names to plot
#' @param lin_name Character string, name of the lineage to plot
#' @param wald_padj Data frame with wald statistics and adjusted p-values
#' @param sce_cell SingleCellExperiment object with pseudotime information
#' @param scObj.cd4.filt Filtered Seurat object with CD4 cells
#' @param type Character string, type identifier for output filename
#' @param out_fold Character string, output folder path
#'
#' @return Saves RDS and TIFF files to disk, returns nothing
plot_multi_gene <- function(genes, lin_name, wald_padj, sce_cell, scObj.cd4.filt, type, out_fold) {

  name <- lin_name

  # Get data frame for 1 lineage, all significant genes - with expression value, pseudotime value and cluster
  pseudotime <- slingshot::slingPseudotime(sce_cell)
  expression <- t(scObj.cd4.filt@assays$RNA@data[genes, ]) #data = log normed counts?
  cluster <- scObj.cd4.filt$celltype
  plot_data <- data.frame(Pseudotime = pseudotime,
                          Expression = expression,
                          Cluster = cluster)

  # All together if its the all lineages one
  if (name == "All lineages") {
    plot_data <- plot_data %>%
      tidyr::pivot_longer(cols = dplyr::starts_with("Expression"),
                          names_to = "Gene",
                          values_to = "Expression") %>%
      dplyr::mutate(Gene = stringr::str_replace_all(Gene, "HLA\\.", "HLA-")) %>%
      dplyr::mutate(Gene = stringr::str_replace(Gene, "Expression\\.", "")) %>%
      tidyr::pivot_longer(cols = dplyr::starts_with("Pseudotime"),
                          names_to = "Lineage",
                          values_to = "Pseudotime") %>%
      dplyr::mutate(Lineage = stringr::str_replace(Lineage, "Pseudotime\\.", "")) %>%
      dplyr::filter(!is.na(Pseudotime))

    frac <- 0.1
  } else {
    plot_data <- plot_data %>%
      dplyr::select(c(paste0("Pseudotime.", name), Cluster, dplyr::starts_with("Expression"))) %>%
      dplyr::rename("Pseudotime" = paste0("Pseudotime.", name)) %>%
      tidyr::pivot_longer(cols = dplyr::starts_with("Expression"),
                          names_to = "Gene", 
                          values_to = "Expression") %>%
      dplyr::mutate(Gene = stringr::str_replace_all(Gene, "HLA\\.", "HLA-")) %>%
      dplyr::mutate(Gene = stringr::str_replace(Gene, "Expression\\.", "")) %>%
      dplyr::filter(!is.na(Pseudotime)) %>%
      dplyr::filter(Gene %in% genes)
    frac <- 0.5
  }

  # Subset randomly so less points to plot
  plot_data_s <- plot_data %>%
    dplyr::sample_frac(frac)

  # Plot it
  gg_genes <- ggplot2::ggplot(plot_data_s,
                              aes(x = Pseudotime,
                                  y = Expression,
                                  color = Cluster),
                              alpha = 0.3) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "loess", color = "black") +
    ggplot2::facet_wrap(~ Gene, ncol = 6) +
    ggplot2::labs(title = name,
                  y = "Normalised RNA expression") +
    ggplot2::geom_text(data = wald_padj,
                       aes(x = max(plot_data_s$Pseudotime) / 2 + 0.5,
                       y = max(plot_data_s$Expression) - 0.05,
                       label = paste0("padj = ", padj, "\nWald = ", waldStat)),
                       inherit.aes = FALSE,
                       size = 3.5) +
    ggplot2::scale_color_manual(values = IMC_colors) +
    gg_theme_specs(fs = 16) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 4)))

  filp <- paste0("CD4_naive_gene_expression_by_pseudotime_", type, "_", name, ".tiff")
  saveRDS(gg_genes, file.path(out_fold, gsub("tiff", "rds", filp)))
  tiff(file.path(out_fold, filp),
       width = 16, height = 8,
       units = "in",
       res = 300,
       compression = "lzw")
  print(gg_genes)
  dev.off()
}

###########################################
#### Gene expression along pseudotime #####
###########################################

#' Rank genes by Wald statistic for GSEA
#'
#' @param stats Data frame with gene statistics
#' @param lineage Lineage name (global, 1-5)
#'
#' @return Named numeric vector of ranked genes (sorted by Wald statistic, decreasing)
rank_genes <- function(stats, lineage = NULL) {
  col_name <- ifelse(lineage == "global", "waldStat", paste0("waldStat_", lineage))
  gene_list <- stats %>%
    dplyr::select(dplyr::contains(col_name)) %>%
    dplyr::pull()
  
  names(gene_list) <- stats$Gene
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list
}


#' Run GO enrichment analysis (ORA or GSEA)
#'
#' @param stats Data frame with gene statistics
#' @param lineage Lineage name (global, 1-5)
#' @param all_genes Character vector of all genes (background)
#' @param catname GO category (default: "CC")
#' @param dir Output directory path
#' @param type Analysis type: "ORA" or "GSEA" (default: "ORA")
#' @param suffix Suffix for output filenames (default: "")
#'
#' @return List with enrichment results and plot
enrichment_go <- function(stats, lineage, all_genes,
                          catname = "CC",
                          dir,
                          type = "ORA",
                          suffix = "") {
  title_name <- ifelse(catname == "CC", "Cellular component",
                       ifelse(catname == "MF", "Molecular function", "Biological processes"))

  if (type == "ORA") {
    genes <- stats %>%
      dplyr::select("Gene") %>%
      dplyr::pull()
    p <- clusterProfiler::enrichGO(gene = genes,
                                   keyType = "SYMBOL",
                                   OrgDb = org.Hs.eg.db,
                                   ont = catname,
                                   minGSSize = 5,
                                   pvalueCutoff = 0.05)
  } else if (type == "GSEA") {
    genes <- rank_genes(stats, lineage)
    p <- clusterProfiler::gseGO(gene = genes,
                                keyType = "SYMBOL",
                                OrgDb = org.Hs.eg.db,
                                ont = catname,
                                minGSSize = 5,
                                pvalueCutoff = 0.05)
  }

  if (nrow(as.data.frame(p)) > 0) {
    plt <- enrichplot::dotplot(p,
                               showCategory = 20,
                               font.size = 10) +
              gg_theme_specs(fs = 16) +
              enrichplot::ggtitle(sprintf("GO: %s", title_name))
    dir.create(file.path(dir, type))
    saveRDS(plt,
            file.path(dir, type,
            sprintf("Figure_5C_GO_%s_de_gene_along_pseudotime_category_%s_%s_%s.rds", type, catname, ifelse(is.null(lineage), "global", lineage), suffix)))
    tiff(file.path(dir, type, sprintf("Figure_5C_GO_%s_de_gene_along_pseudotime_category_%s_%s_%s.tiff", type, catname, ifelse(is.null(lineage), "global", lineage), suffix)),
         width = 7, height = 11,
         units = "in",
         res = 300,
         compression = "lzw")
    print(plt)
    dev.off()
    write.table(as.data.frame(p),
                file.path(dir, type, sprintf("GO_%s_de_gene_along_pseudotime_category_%s_%s_%s.tsv", type, catname, ifelse(is.null(lineage), "global", lineage), suffix)),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    return(list(res = p, plt = plt))
  } else {
    message("Nothing to return")
  }
}

#' Plot smoothers per condition for selected lineage with 
#' optional statistics display from conditionTest
#'
#' @param model GAM model object
#' @param gene Character string, gene name to plot
#' @param cols Colour vector for conditions
#' @param lineage Character string, lineage name
#' @param stat Character string, statistics annotation (optional)
#'
#' @return ggplot object
get_smoothers_per_condition_for_one_lineage_with_cond_stats <- function(model,
                                                                        gene,
                                                                        cols,
                                                                        lineage,
                                                                        stat = NULL) {  
    p <- tradeSeq::plotSmoothers(model,
                                 SummarizedExperiment::assays(model)[["counts"]],
                                 gene = gene,
                                 size = 2,
                                 alpha = 1,
                                 border = FALSE,
                                 curvesCols = as.vector(cols)) +
                  ggplot2::scale_colour_manual(values = cols,
                                               labels = ifelse(cols == "#FFFFFF00", "", ifelse(cols == "#FDBD13", "Group 1", "Group 2"))) +
                  ggplot2::labs(title = gene, color = "", y = "log(expr + 1)") +
                  ggplot2::guides(colour = guide_legend(nrow = 1, override.aes = list(size = 4))) +
                  gg_theme_specs(fs = 16)
    ymax <- max(log(p$data$gene_count + 1))
    if (!is.null(stat)) {
      p <- p +
        ggplot2::annotate("text", x = 3.5, y = ymax - 0.10, label = stat, color = "black", size = 4.5, fontface = "bold")
    }
    return(p)
}

#' Convert padj and Wald statistics to annotation string
#' for plotSmoothers with statistics
#'
#' @param df Data frame with padj and waldStat columns
#'
#' @return Character string with formatted annotation, or NULL if df is empty
convert_padj_wald_annotation <- function(df) {
  if (nrow(df) == 0) {
    res <- NULL
  } else {
    padj_val <- df %>% dplyr::select(dplyr::matches("padj")) %>% dplyr::pull()
    wald_val <- df %>% dplyr::select(dplyr::matches("waldStat")) %>% dplyr::pull()
    res <- paste0("padj = ", formatC(padj_val, format = "e", digits = 1), "\n", "Wald = ", round(wald_val, 1))
  }
  return(res)
}
