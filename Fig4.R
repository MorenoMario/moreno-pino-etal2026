# Fig4.R
# Description: Generates publication-ready functional heatmaps for marine microorganism genomic data.
# Purpose: Visualizes the presence/absence and quantitative abundance of biogeochemical cycling genes
#          (nitrogen, sulfur, photosynthesis, carbon fixation) across MAGs. Produces both transposed
#          and standard heatmap orientations, focused panels per metabolic pathway, composite multi-panel
#          figures (A-B-C), and supplementary tables. Designed for mBio publication standards (600 DPI).

input_file         <- "heatmap.txt"
output_prefix      <- "heatmap_publication"

enable_intelligent_filtering <- FALSE
focus_mode                  <- "all"
min_presence_threshold      <- 0
max_presence_percentage     <- 100

color_palette_option <- 4

dist_rows          <- "binary"
dist_cols          <- "binary"
cluster_rows_flag  <- TRUE
cluster_cols_flag  <- FALSE
border_color       <- "#CCCCCC"
cellwidth          <- 8
cellheight         <- 10
angle_col          <- 45
angle_row          <- 0

fontsize_row_transposed_quant   <- 7
fontsize_col_transposed_quant   <- 6
fontsize_row_transposed_binary  <- 8
fontsize_col_transposed_binary  <- 6

dist_rows_normal          <- "binary"
dist_cols_normal          <- "binary"
cluster_rows_flag_normal  <- FALSE
cluster_cols_flag_normal  <- TRUE
border_color_normal       <- "#CCCCCC"
cellwidth_normal          <- 8
cellheight_normal         <- 6
angle_col_normal          <- 45
angle_row_normal          <- 0

fontsize_row_normal_quant   <- 6
fontsize_col_normal_quant   <- 8
fontsize_row_normal_binary  <- 5
fontsize_col_normal_binary  <- 8

generate_transposed <- TRUE
generate_normal     <- TRUE
generate_both_types <- TRUE
generate_focused    <- TRUE
output_dpi          <- 600
output_format       <- "both"

presence_colors <- c("#fff5eb", "#EC9078")
quantitative_colors <- c("#FFFFFF", "#EFF3FF", "#C6DBEF", "#9ECAE1",
                         "#4292C6", "#2171B5")

group_colors <- c("Widespread" = "#1F77B4", "Restricted" = "#FF7F0E")

functional_colors_option1 <- c(
  "Assimilatory nitrate reduction"    = "#2166AC",
  "Assimilatory nitrite reduction"    = "#4393C3",
  "Periplasmic nitrate reduction"     = "#92C5DE",
  "Respiratory nitrate reduction"     = "#053061",
  "Nitric oxide reduction"            = "#67A9CF",
  "Respiratory nitrite reduction"     = "#D1E5F0",
  "Nitrous oxide reduction"           = "#5BA3D0",
  "DNRA"                             = "#1E5A73",
  "Sulfate reduction (Apr)"          = "#D73027",
  "Sulfate reduction (Dsr)"          = "#F46D43",
  "Sulfate activation"               = "#FDAE61",
  "Sulfur oxidation"                 = "#FEE090",
  "Anoxygenic photosynthesis"        = "#238B45",
  "Bacteriochlorophyll synthesis"    = "#74C476",
  "Carbon fixation (other)"          = "#7B3F99",
  "Calvin-Benson-Bassham cycle"      = "#AE7BAE",
  "Photosystem I"                    = "#41AB5D",
  "Photosystem II"                   = "#A1D99B",
  "Other functions"                  = "#737373"
)

functional_colors_option2 <- c(
  "Assimilatory nitrate reduction"    = "#440154",
  "Assimilatory nitrite reduction"    = "#482677",
  "Periplasmic nitrate reduction"     = "#3F4A8A",
  "Respiratory nitrate reduction"     = "#31678E",
  "Nitric oxide reduction"            = "#26838F",
  "Respiratory nitrite reduction"     = "#1F9D8A",
  "Nitrous oxide reduction"           = "#6CCE59",
  "DNRA"                             = "#B6DE2B",
  "Sulfate reduction (Apr)"          = "#FEE825",
  "Sulfate reduction (Dsr)"          = "#FDAE61",
  "Sulfate activation"               = "#F46D43",
  "Sulfur oxidation"                 = "#D73027",
  "Anoxygenic photosynthesis"        = "#A50F15",
  "Bacteriochlorophyll synthesis"    = "#CB181D",
  "Carbon fixation (other)"          = "#EF3B2C",
  "Calvin-Benson-Bassham cycle"      = "#FB6A4A",
  "Photosystem I"                    = "#FC9272",
  "Photosystem II"                   = "#FCBBA1",
  "Other functions"                  = "#969696"
)

functional_colors_option3 <- c(
  "Assimilatory nitrate reduction"    = "#8DD3C7",
  "Assimilatory nitrite reduction"    = "#FFFFB3",
  "Periplasmic nitrate reduction"     = "#BEBADA",
  "Respiratory nitrate reduction"     = "#FB8072",
  "Nitric oxide reduction"            = "#80B1D3",
  "Respiratory nitrite reduction"     = "#FDB462",
  "Nitrous oxide reduction"           = "#B3DE69",
  "DNRA"                             = "#FCCDE5",
  "Sulfate reduction (Apr)"          = "#D9D9D9",
  "Sulfate reduction (Dsr)"          = "#BC80BD",
  "Sulfate activation"               = "#CCEBC5",
  "Sulfur oxidation"                 = "#FFED6F",
  "Anoxygenic photosynthesis"        = "#E78AC3",
  "Bacteriochlorophyll synthesis"    = "#A6D854",
  "Carbon fixation (other)"          = "#FFD92F",
  "Calvin-Benson-Bassham cycle"      = "#E5C494",
  "Photosystem I"                    = "#B3B3B3",
  "Photosystem II"                   = "#8FA3D3",
  "Other functions"                  = "#999999"
)

functional_colors_option4 <- c(
  "Assimilatory nitrate reduction"    = "#08519C",
  "Assimilatory nitrite reduction"    = "#2171B5",
  "Periplasmic nitrate reduction"     = "#4292C6",
  "Respiratory nitrate reduction"     = "#6BAED6",
  "Nitric oxide reduction"            = "#9ECAE1",
  "Respiratory nitrite reduction"     = "#C6DBEF",
  "Nitrous oxide reduction"           = "#3182BD",
  "DNRA"                             = "#084594",
  "Sulfate reduction (Apr)"          = "#A63603",
  "Sulfate reduction (Dsr)"          = "#CC4C02",
  "Sulfate activation"               = "#EC7014",
  "Sulfur oxidation"                 = "#FE9929",
  "Anoxygenic photosynthesis"        = "#00441B",
  "Bacteriochlorophyll synthesis"    = "#238B45",
  "Photosystem I"                    = "#41AB5D",
  "Photosystem II"                   = "#74C476",
  "Carbon fixation (other)"          = "#54278F",
  "Calvin-Benson-Bassham cycle"      = "#756BB1",
  "Other functions"                  = "#525252"
)

functional_colors_option5 <- c(
  "Assimilatory nitrate reduction"    = "#35978f",
  "Assimilatory nitrite reduction"    = "#35978f",
  "Periplasmic nitrate reduction"     = "#35978f",
  "Respiratory nitrate reduction"     = "#35978f",
  "Nitric oxide reduction"            = "#35978f",
  "Respiratory nitrite reduction"     = "#35978f",
  "Nitrous oxide reduction"           = "#35978f",
  "DNRA"                             = "#35978f",
  "Sulfate reduction (Apr)"          = "#5ab4ac",
  "Sulfate reduction (Dsr)"          = "#5ab4ac",
  "Sulfate activation"               = "#5ab4ac",
  "Sulfur oxidation"                 = "#5ab4ac",
  "Anoxygenic photosynthesis"        = "#f6e8c3",
  "Bacteriochlorophyll synthesis"    = "#f6e8c3",
  "Photosystem I"                    = "#f6e8c3",
  "Photosystem II"                   = "#f6e8c3",
  "Carbon fixation (other)"          = "#fdae61",
  "Calvin-Benson-Bassham cycle"      = "#fdae61",
  "Other functions"                  = "#595959"
)

functional_colors_option6 <- c(

  "Assimilatory nitrate reduction"    = "#08306B",
  "Assimilatory nitrite reduction"    = "#2171B5",
  "Periplasmic nitrate reduction"     = "#4292C6",
  "Respiratory nitrate reduction"     = "#6BAED6",
  "Nitric oxide reduction"            = "#9ECAE1",
  "Respiratory nitrite reduction"     = "#C6DBEF",
  "Nitrous oxide reduction"           = "#3182BD",
  "DNRA"                             = "#08519C",

  "Sulfate reduction (Apr)"          = "#7F2704",
  "Sulfate reduction (Dsr)"          = "#A63603",
  "Sulfate activation"               = "#CC4C02",
  "Sulfur oxidation"                 = "#EC7014",

  "Anoxygenic photosynthesis"        = "#00441B",
  "Bacteriochlorophyll synthesis"    = "#238B45",
  "Photosystem I"                    = "#41AB5D",
  "Photosystem II"                   = "#74C476",

  "Carbon fixation (other)"          = "#4A1486",
  "Calvin-Benson-Bassham cycle"      = "#6A51A3",

  "Other functions"                  = "#737373"
)

functional_colors <- switch(as.character(color_palette_option),
                            "1" = functional_colors_option1,
                            "2" = functional_colors_option2,
                            "3" = functional_colors_option3,
                            "4" = functional_colors_option4,
                            "5" = functional_colors_option5,
                            "6" = functional_colors_option6,
                            functional_colors_option4
)

cat("Using color palette option:", color_palette_option, "\n")

required_packages <- c("pheatmap", "RColorBrewer", "ggplot2", "dplyr",
                       "viridis", "grid", "gridExtra", "scales", "cowplot", "stringr")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

intelligent_filtering <- function(data, min_presence = 2, max_presence_pct = 90) {
  cat("Aplicando filtrado inteligente de genes...\n")

  func_cols <- 5:ncol(data)
  gene_presence <- colSums(data[, func_cols] > 0)
  total_organisms <- nrow(data)

  min_filter <- gene_presence >= min_presence
  max_filter <- gene_presence <= (total_organisms * max_presence_pct / 100)

  keep_genes <- min_filter & max_filter

  cat("FILTRADO ESTADÍSTICO:\n")
  cat("- Genes muy raros (< ", min_presence, " organismos): ", sum(!min_filter), "\n")
  cat("- Genes muy comunes (> ", max_presence_pct, "%): ", sum(!max_filter), "\n")
  cat("- Genes informativos conservados: ", sum(keep_genes), " de ", length(keep_genes), "\n")

  if (sum(keep_genes) > 0) {
    filtered_data <- data[, c(1:4, which(keep_genes) + 4)]
    return(filtered_data)
  } else {
    cat("ADVERTENCIA: Filtrado demasiado estricto. Usando datos originales.\n")
    return(data)
  }
}

focus_by_process <- function(data, process = "all") {
  if (process == "all") {
    return(list(data = data, suffix = ""))
  }

  func_cols <- 5:ncol(data)
  gene_names <- colnames(data)[func_cols]

  if (process == "nitrogen") {
    nitrogen_patterns <- c("nasc", "narb", "nira", "nirb", "nird",
                           "napa", "napb", "narg", "narh", "nary", "narz", "nxra", "nxrb",
                           "norb", "norc", "nirs", "nosz", "nrfh", "nrfa")
    keep_genes <- sapply(gene_names, function(x) {
      any(sapply(nitrogen_patterns, function(p) grepl(p, tolower(x))))
    })
    suffix <- "- Nitrogen cycling genes"

  } else if (process == "sulfur") {
    sulfur_patterns <- c("apra", "aprb", "dsra", "dsrb", "sat", "met3",
                         "soxx", "soxb", "soxc", "soxy", "soxz")
    keep_genes <- sapply(gene_names, function(x) {
      any(sapply(sulfur_patterns, function(p) grepl(p, tolower(x))))
    })
    suffix <- "- Sulfur cycling genes"

  } else if (process == "photosynthesis") {
    photo_patterns <- c("pufa", "pufb", "puca", "pufl", "pufm", "bche",
                        "chlb", "chll", "chln", "psaa", "psab", "psad",
                        "psaf", "psae", "psbc", "psbb", "psbe", "psbf",
                        "psba", "psbd", "rbcl", "cbbl", "rbcs", "cbbs", "acsf", "chle")
    keep_genes <- sapply(gene_names, function(x) {
      any(sapply(photo_patterns, function(p) grepl(p, tolower(x))))
    })
    suffix <- "- Photosynthesis and carbon fixation genes"

  } else {
    return(list(data = data, suffix = ""))
  }

  if (sum(keep_genes) > 0) {
    focused_data <- data[, c(1:4, which(keep_genes) + 4)]
    cat("FOCALIZACIÓN EN", toupper(process), ":\n")
    cat("- Genes seleccionados:", sum(keep_genes), "\n")
    return(list(data = focused_data, suffix = suffix))
  } else {
    cat("ADVERTENCIA: No se encontraron genes para", process, ". Usando todos.\n")
    return(list(data = data, suffix = ""))
  }
}

categorize_function_enhanced <- function(gene_name) {
  gene_lower <- tolower(gene_name)

  if (grepl("nasc|nasa", gene_lower)) return("Assimilatory nitrate reduction")
  if (grepl("narb", gene_lower)) return("Assimilatory nitrite reduction")
  if (grepl("nira|nirb|nird", gene_lower)) return("Assimilatory nitrite reduction")
  if (grepl("napa|napb", gene_lower)) return("Periplasmic nitrate reduction")
  if (grepl("narg|narh|nary|narz|nxra|nxrb", gene_lower)) return("Respiratory nitrate reduction")
  if (grepl("norb|norc", gene_lower)) return("Nitric oxide reduction")
  if (grepl("nirs", gene_lower)) return("Respiratory nitrite reduction")
  if (grepl("nosz", gene_lower)) return("Nitrous oxide reduction")
  if (grepl("nrfh|nrfa", gene_lower)) return("DNRA")

  if (grepl("apra|aprb", gene_lower)) return("Sulfate reduction (Apr)")
  if (grepl("dsra|dsrb", gene_lower)) return("Sulfate reduction (Dsr)")
  if (grepl("sat|met3", gene_lower)) return("Sulfate activation")
  if (grepl("soxx|soxb|soxc|soxy|soxz", gene_lower)) return("Sulfur oxidation")

  if (grepl("pufa|pufb|puca|pufl|pufm", gene_lower)) return("Anoxygenic photosynthesis")
  if (grepl("bche|chlb|chll|chln", gene_lower)) return("Bacteriochlorophyll synthesis")

  if (grepl("acsf|chle", gene_lower)) return("Carbon fixation (other)")
  if (grepl("rbcl|cbbl|rbcs|cbbs", gene_lower)) return("Calvin-Benson-Bassham cycle")

  if (grepl("psaa|psab|psad|psaf|psae", gene_lower)) return("Photosystem I")
  if (grepl("psbc|psbb|psbe|psbf|psba|psbd", gene_lower)) return("Photosystem II")

  return("Other functions")
}

debug_categorization <- function(data) {
  func_cols <- 5:ncol(data)
  clean_names <- colnames(data)[func_cols]
  clean_names <- gsub(" K[0-9]+.*$", "", clean_names)
  clean_names <- gsub("^[^,]+, ", "", clean_names)

  func_categories <- sapply(clean_names, categorize_function_enhanced)

  cat("ANÁLISIS DE GENES FUNCIONALES\n")
  cat(strrep("=", 70), "\n")

  for(i in 1:length(clean_names)) {
    cat(sprintf("%-30s -> %s\n", clean_names[i], func_categories[i]))
  }

  cat("\n", strrep("=", 70), "\n")
  cat("RESUMEN POR CATEGORÍA FUNCIONAL:\n")
  cat(strrep("=", 70), "\n")

  category_counts <- table(func_categories)
  for(cat_name in names(category_counts)) {
    cat(sprintf("%-40s: %2d genes\n", cat_name, category_counts[cat_name]))
  }

  return(func_categories)
}

create_functional_diversity_report <- function(data) {
  func_cols <- 5:ncol(data)

  org_diversity <- data.frame(
    Organism = data$name,
    Distribution = data$final.clasification.group,
    Total_genes = rowSums(data[, func_cols] > 0),
    Nitrogen_genes = rowSums(data[, func_cols[grepl("nasc|narb|nir|nap|nar|nor|nos|nrf",
                                                    colnames(data)[func_cols], ignore.case = TRUE)]] > 0),
    Sulfur_genes = rowSums(data[, func_cols[grepl("apr|dsr|sat|met3|sox",
                                                  colnames(data)[func_cols], ignore.case = TRUE)]] > 0),
    Photosynthesis_genes = rowSums(data[, func_cols[grepl("puf|puc|bch|chl|psa|psb|rbc|cbb|acs",
                                                          colnames(data)[func_cols], ignore.case = TRUE)]] > 0),
    stringsAsFactors = FALSE
  )

  return(org_diversity)
}

perform_basic_stats <- function(data) {
  func_cols <- 5:ncol(data)

  stats_results <- list()

  for (i in func_cols) {
    gene_name <- colnames(data)[i]

    contingency_table <- table(data$final.clasification.group, data[, i] > 0)

    if (all(dim(contingency_table) == c(2, 2))) {
      fisher_test <- fisher.test(contingency_table)
      stats_results[[gene_name]] <- list(
        p_value = fisher_test$p.value,
        odds_ratio = fisher_test$estimate,
        significant = fisher_test$p.value < 0.05
      )
    }
  }

  return(stats_results)
}

create_both_combined_figures <- function(data_filtered, output_prefix) {

  cat("\n", strrep("=", 70), "\n")
  cat("GENERANDO FIGURAS COMPUESTAS CON PANELES A-B-C\n")
  cat("(Nitrogen, Sulfur, Photosynthesis - AMBAS ORIENTACIONES)\n")
  cat(strrep("=", 70), "\n")

  suppressWarnings(suppressMessages({
    if (!require("gridExtra", quietly = TRUE)) {
      install.packages("gridExtra", repos = "https://cloud.r-project.org")
      library(gridExtra)
    }
    if (!require("cowplot", quietly = TRUE)) {
      install.packages("cowplot", repos = "https://cloud.r-project.org")
      library(cowplot)
    }
  }))

  processes <- c("nitrogen", "sulfur", "photosynthesis")
  panels_transposed <- list()
  panels_normal <- list()

  for (process in processes) {
    process_result <- focus_by_process(data_filtered, process)
    process_data <- process_result$data
    process_suffix <- process_result$suffix

    if (ncol(process_data) > 4) {
      cat(" - Creando paneles para", process, "...\n")

      hm_transposed <- create_publication_heatmap(
        process_data,
        type = "binary",
        transposed = TRUE,
        title_suffix = paste("in marine microorganisms", process_suffix)
      )
      panels_transposed[[process]] <- hm_transposed$gtable

      hm_normal <- create_publication_heatmap(
        process_data,
        type = "binary",
        transposed = FALSE,
        title_suffix = paste("in marine microorganisms", process_suffix)
      )
      panels_normal[[process]] <- hm_normal$gtable

    } else {
      cat(" - Saltando", process, "(genes insuficientes)\n")
    }
  }

  panel_labels <- c("A", "B", "C")

  if (length(panels_transposed) > 0) {
    pdf(paste0(output_prefix, "_panels_ABC_transposed.pdf"), width = 18, height = 12)
    plot_transposed <- plot_grid(
      plotlist = panels_transposed,
      labels = panel_labels[1:length(panels_transposed)],
      ncol = 3,
      label_size = 16,
      label_fontface = "bold"
    )
    print(plot_transposed)
    dev.off()

    png(paste0(output_prefix, "_panels_ABC_transposed.png"),
        width = 18 * output_dpi, height = 12 * output_dpi, res = output_dpi)
    print(plot_transposed)
    dev.off()

    cat("✓ Figura TRANSPUESTA guardada:", paste0(output_prefix, "_panels_ABC_transposed.pdf/png"), "\n")
  }

  if (length(panels_normal) > 0) {
    pdf(paste0(output_prefix, "_panels_ABC_normal.pdf"), width = 20, height = 8)
    plot_normal <- plot_grid(
      plotlist = panels_normal,
      labels = panel_labels[1:length(panels_normal)],
      ncol = length(panels_normal),
      label_size = 16,
      label_fontface = "bold"
    )
    print(plot_normal)
    dev.off()

    png(paste0(output_prefix, "_panels_ABC_normal.png"),
        width = 20 * output_dpi, height = 8 * output_dpi, res = output_dpi)
    print(plot_normal)
    dev.off()

    cat("✓ Figura NORMAL guardada:", paste0(output_prefix, "_panels_ABC_normal.pdf/png"), "\n")
  }

  return(list(transposed = panels_transposed, normal = panels_normal))
}

create_publication_heatmap <- function(data, type = "binary", transposed = FALSE, title_suffix = "") {

  func_cols <- 5:ncol(data)
  mat <- as.matrix(data[, func_cols])
  rownames(mat) <- data$name

  clean_names <- colnames(mat)
  clean_names <- gsub(" K[0-9]+.*$", "", clean_names)
  clean_names <- gsub("^[^,]+, ", "", clean_names)

  clean_names <- paste0(toupper(substring(clean_names, 1, 1)),
                        substring(clean_names, 2))
  colnames(mat) <- clean_names

  if (transposed) {
    mat <- t(mat)
    title_suffix <- paste("(transposed)", title_suffix)
  }

  if (transposed) {
    func_categories <- sapply(rownames(mat), categorize_function_enhanced)
    annotation_row <- data.frame(
      Function = factor(func_categories),
      row.names = rownames(mat)
    )
    annotation_col <- data.frame(
      Distribution = factor(data$final.clasification.group,
                            levels = c("wide", "rare"),
                            labels = c("Widespread", "Restricted")),
      row.names = data$name
    )
    annotation_col <- annotation_col[colnames(mat), , drop = FALSE]
  } else {
    func_categories <- sapply(colnames(mat), categorize_function_enhanced)
    annotation_col <- data.frame(
      Function = factor(func_categories),
      row.names = colnames(mat)
    )
    annotation_row <- data.frame(
      Distribution = factor(data$final.clasification.group,
                            levels = c("wide", "rare"),
                            labels = c("Widespread", "Restricted")),
      row.names = rownames(mat)
    )
  }

  if (transposed) {
    if (type == "binary") {
      fontsize_row_param <- fontsize_row_transposed_binary
      fontsize_col_param <- fontsize_col_transposed_binary
    } else {
      fontsize_row_param <- fontsize_row_transposed_quant
      fontsize_col_param <- fontsize_col_transposed_quant
    }
    cluster_rows_param <- cluster_rows_flag
    cluster_cols_param <- cluster_cols_flag
    cellwidth_param <- cellwidth
    cellheight_param <- cellheight
    border_param <- border_color
    distance_rows_param <- dist_rows
    distance_cols_param <- dist_cols
  } else {
    if (type == "binary") {
      fontsize_row_param <- fontsize_row_normal_binary
      fontsize_col_param <- fontsize_col_normal_binary
    } else {
      fontsize_row_param <- fontsize_row_normal_quant
      fontsize_col_param <- fontsize_col_normal_quant
    }
    cluster_rows_param <- cluster_rows_flag_normal
    cluster_cols_param <- cluster_cols_flag_normal
    cellwidth_param <- cellwidth_normal
    cellheight_param <- cellheight_normal
    border_param <- border_color_normal
    distance_rows_param <- dist_rows_normal
    distance_cols_param <- dist_cols_normal
  }

  if (type == "binary") {
    mat[mat > 0] <- 1
    colors <- presence_colors
    main_title <- paste("Functional gene presence/absence", title_suffix)
    breaks_param <- c(-0.5, 0.5, 1.5)
    legend_breaks_param <- c(0, 1)
    legend_labels_param <- c("Absent", "Present")
  } else {
    colors <- colorRampPalette(quantitative_colors)(50)
    main_title <- paste("Functional gene abundance", title_suffix)
    breaks_param <- NULL
    legend_breaks_param <- NULL
    legend_labels_param <- NULL
  }

  annotation_colors <- list(
    Function = functional_colors,
    Distribution = group_colors
  )

  p <- pheatmap(
    mat,
    color = colors,
    breaks = breaks_param,
    cluster_rows = cluster_rows_param,
    cluster_cols = cluster_cols_param,
    clustering_distance_rows = distance_rows_param,
    clustering_distance_cols = distance_cols_param,
    clustering_method = "ward.D2",

    annotation_row = if (exists("annotation_row")) annotation_row else NULL,
    annotation_col = if (exists("annotation_col")) annotation_col else NULL,
    annotation_colors = annotation_colors,

    cellwidth = cellwidth_param,
    cellheight = cellheight_param,
    fontsize_row = fontsize_row_param,
    fontsize_col = fontsize_col_param,

    border_color = border_param,
    main = main_title,
    show_rownames = TRUE,
    show_colnames = TRUE,

    treeheight_row = if (cluster_rows_param) 30 else 0,
    treeheight_col = if (cluster_cols_param) 30 else 0,

    legend = TRUE,
    legend_breaks = legend_breaks_param,
    legend_labels = legend_labels_param
  )

  return(p)
}

save_heatmap <- function(heatmap_obj, filename_base, orientation, type, process = "") {
  if (process != "") {
    filename_base <- paste0(filename_base, "_", process)
  }

  if (output_format %in% c("pdf", "both")) {
    pdf(paste0(filename_base, "_", orientation, "_", type, ".pdf"),
        width = if(orientation == "transposed") 14 else 12,
        height = if(orientation == "transposed") 10 else 12)
    print(heatmap_obj)
    dev.off()
    cat("Saved:", paste0(filename_base, "_", orientation, "_", type, ".pdf"), "\n")
  }

  if (output_format %in% c("png", "both")) {
    png(paste0(filename_base, "_", orientation, "_", type, ".png"),
        width = (if(orientation == "transposed") 14 else 12) * output_dpi,
        height = (if(orientation == "transposed") 10 else 12) * output_dpi,
        res = output_dpi)
    print(heatmap_obj)
    dev.off()
    cat("Saved:", paste0(filename_base, "_", orientation, "_", type, ".png"), "\n")
  }
}

preview_color_palette <- function(colors_list, palette_name) {
  library(ggplot2)

  color_data <- data.frame(
    Category = names(colors_list),
    Color = as.character(colors_list),
    y = seq_along(colors_list)
  )

  p <- ggplot(color_data, aes(x = 1, y = y, fill = Color)) +
    geom_tile(height = 0.9, width = 0.9) +
    scale_fill_identity() +
    geom_text(aes(label = Category), x = 1.5, hjust = 0, size = 3) +
    theme_void() +
    theme(legend.position = "none",
          plot.margin = margin(1,6,1,1, "cm")) +
    labs(title = paste("Functional Colors -", palette_name)) +
    xlim(0, 4)

  print(p)
  return(p)
}

cat("GENERADOR DE HEATMAPS PUBLICATION-READY\n")
cat(strrep("=", 70), "\n")

if (file.exists(input_file)) {
  data <- read.delim(input_file, sep = "\t", header = TRUE,
                     stringsAsFactors = FALSE, check.names = FALSE)

  cat("Datos cargados exitosamente:\n")
  cat("- Organismos:", nrow(data), "\n")
  cat("- Genes funcionales:", ncol(data) - 4, "\n")

  data_filtered <- data[rowSums(data[, 5:ncol(data)]) > 0, ]
  cat("- Organismos con genes funcionales:", nrow(data_filtered), "\n")

  if (enable_intelligent_filtering) {
    data_filtered <- intelligent_filtering(data_filtered,
                                           min_presence_threshold,
                                           max_presence_percentage)
  }

  focused_result <- focus_by_process(data_filtered, focus_mode)
  data_focused <- focused_result$data
  title_suffix <- focused_result$suffix

  palette_names <- c("Biochemical Cycles", "Viridis Modified", "Soft Professional",
                     "Scientific Journal Style", "Minimalist", "Nature Optimized")
  cat("- Paleta de colores:", color_palette_option, "-", palette_names[min(color_palette_option, 6)], "\n")
  cat("- Modo de focalización:", toupper(focus_mode), "\n")
  cat("- Filtrado inteligente:", ifelse(enable_intelligent_filtering, "ACTIVADO", "DESACTIVADO"), "\n")

  cat("\n")
  debug_result <- debug_categorization(data_focused)

  cat("\n", strrep("=", 70), "\n")
  cat("GENERANDO HEATMAPS ESTÁNDAR\n")
  cat(strrep("=", 70), "\n")

  if (generate_transposed) {
    cat("\nCreando heatmaps TRANSPUESTOS (genes como filas, organismos como columnas)...\n")

    if (generate_both_types) {
      cat("- Heatmap binario transpuesto...\n")
      heatmap_trans_binary <- create_publication_heatmap(data_focused, type = "binary",
                                                         transposed = TRUE,
                                                         title_suffix = paste("in marine microorganisms", title_suffix))
      save_heatmap(heatmap_trans_binary, output_prefix, "transposed", "binary")
    }

    if (generate_both_types && max(data_focused[, 5:ncol(data_focused)]) > 1) {
      cat("- Heatmap cuantitativo transpuesto...\n")
      heatmap_trans_quant <- create_publication_heatmap(data_focused, type = "quantitative",
                                                        transposed = TRUE,
                                                        title_suffix = paste("in marine microorganisms", title_suffix))
      save_heatmap(heatmap_trans_quant, output_prefix, "transposed", "quantitative")
    }
  }

  if (generate_normal) {
    cat("\nCreando heatmaps NORMALES (organismos como filas, genes como columnas)...\n")

    if (generate_both_types) {
      cat("- Heatmap binario normal...\n")
      heatmap_norm_binary <- create_publication_heatmap(data_focused, type = "binary",
                                                        transposed = FALSE,
                                                        title_suffix = paste("in marine microorganisms", title_suffix))
      save_heatmap(heatmap_norm_binary, output_prefix, "normal", "binary")
    }

    if (generate_both_types && max(data_focused[, 5:ncol(data_focused)]) > 1) {
      cat("- Heatmap cuantitativo normal...\n")
      heatmap_norm_quant <- create_publication_heatmap(data_focused, type = "quantitative",
                                                       transposed = FALSE,
                                                       title_suffix = paste("in marine microorganisms", title_suffix))
      save_heatmap(heatmap_norm_quant, output_prefix, "normal", "quantitative")
    }
  }

  if (generate_focused && focus_mode == "all") {
    cat("\n", strrep("=", 70), "\n")
    cat("GENERANDO HEATMAPS FOCALIZADOS POR PROCESO METABÓLICO\n")
    cat(strrep("=", 70), "\n")

    processes <- c("nitrogen", "sulfur", "photosynthesis")

    for (process in processes) {
      cat("\n--- PROCESO:", toupper(process), "---\n")

      process_result <- focus_by_process(data_filtered, process)
      process_data <- process_result$data
      process_suffix <- process_result$suffix

      if (ncol(process_data) > 4) {

        cat("Creando heatmap transpuesto focalizado para", process, "...\n")
        heatmap_process_transposed <- create_publication_heatmap(
          process_data,
          type = "binary",
          transposed = TRUE,
          title_suffix = paste("in marine microorganisms", process_suffix)
        )
        save_heatmap(heatmap_process_transposed, output_prefix, "transposed", "binary", process)

        cat("Creando heatmap normal focalizado para", process, "...\n")
        heatmap_process_normal <- create_publication_heatmap(
          process_data,
          type = "binary",
          transposed = FALSE,
          title_suffix = paste("in marine microorganisms", process_suffix)
        )
        save_heatmap(heatmap_process_normal, output_prefix, "normal", "binary", process)

      } else {
        cat("No se encontraron genes suficientes para", process, "\n")
      }
    }
  }

  cat("\n", strrep("=", 70), "\n")
  cat("EXPORTANDO DATOS SUPLEMENTARIOS\n")
  cat(strrep("=", 70), "\n")

  write.csv(data_focused, paste0(output_prefix, "_data_filtered.csv"), row.names = FALSE)
  cat("- Datos filtrados:", paste0(output_prefix, "_data_filtered.csv"), "\n")

  func_summary <- data.frame(
    Gene = colnames(data_focused)[5:ncol(data_focused)],
    Functional_Category = sapply(colnames(data_focused)[5:ncol(data_focused)],
                                 categorize_function_enhanced),
    Present_in_Organisms = colSums(data_focused[, 5:ncol(data_focused)] > 0),
    Percentage_Present = round(colSums(data_focused[, 5:ncol(data_focused)] > 0) /
                                 nrow(data_focused) * 100, 1),
    stringsAsFactors = FALSE
  )

  write.csv(func_summary, paste0(output_prefix, "_functional_summary.csv"), row.names = FALSE)
  cat("- Resumen funcional:", paste0(output_prefix, "_functional_summary.csv"), "\n")

  binary_matrix <- data_focused[, 5:ncol(data_focused)]
  binary_matrix[binary_matrix > 0] <- 1
  cooccurrence_matrix <- t(as.matrix(binary_matrix)) %*% as.matrix(binary_matrix)
  write.csv(cooccurrence_matrix, paste0(output_prefix, "_cooccurrence_matrix.csv"))
  cat("- Matriz de co-ocurrencia:", paste0(output_prefix, "_cooccurrence_matrix.csv"), "\n")

  cat("\n", strrep("=", 70), "\n")
  cat("ANÁLISIS ESTADÍSTICO ADICIONAL\n")
  cat(strrep("=", 70), "\n")

  diversity_report <- create_functional_diversity_report(data_focused)
  write.csv(diversity_report, paste0(output_prefix, "_diversity_report.csv"), row.names = FALSE)
  cat("- Reporte de diversidad:", paste0(output_prefix, "_diversity_report.csv"), "\n")

  cat("Realizando tests estadísticos...\n")
  tryCatch({
    stats_results <- perform_basic_stats(data_focused)

    significant_genes <- sum(sapply(stats_results, function(x) x$significant), na.rm = TRUE)
    total_tested <- length(stats_results)

    cat("- Genes con asociación significativa (p<0.05):", significant_genes, "de", total_tested, "\n")

    if (length(stats_results) > 0) {
      stats_df <- data.frame(
        Gene = names(stats_results),
        P_value = sapply(stats_results, function(x) x$p_value),
        Odds_Ratio = sapply(stats_results, function(x) as.numeric(x$odds_ratio)),
        Significant = sapply(stats_results, function(x) x$significant),
        stringsAsFactors = FALSE
      )
      write.csv(stats_df, paste0(output_prefix, "_statistical_tests.csv"), row.names = FALSE)
      cat("- Resultados estadísticos:", paste0(output_prefix, "_statistical_tests.csv"), "\n")
    }

  }, error = function(e) {
    cat("Advertencia: No se pudieron realizar todos los tests estadísticos\n")
  })

  if (exists("data_filtered")) {

    cat("\n", strrep("=", 70), "\n")
    cat("CREANDO FIGURAS COMPUESTAS FINALES\n")
    cat(strrep("=", 70), "\n")

    combined_figures <- create_both_combined_figures(data_filtered, output_prefix)

    cat("\n", strrep("=", 70), "\n")
    cat("FIGURAS COMPUESTAS GENERADAS:\n")
    cat(strrep("=", 70), "\n")

    compound_files <- list.files(pattern = paste0("^", output_prefix, "_panels_ABC"), full.names = FALSE)
    if (length(compound_files) > 0) {
      for (f in compound_files) {
        cat("✓ ", f, "\n")
      }
    } else {
      cat("No se generaron figuras compuestas\n")
    }

    cat("\nRECOMENDACIONES FINALES:\n")
    cat("- Para publicación: Use la versión NORMAL (organismos como filas)\n")
    cat("- Para análisis detallado: Use la versión TRANSPUESTA (genes como filas)\n")
    cat("- Ambas versiones están optimizadas para publicaciones científicas\n")

  }

  cat("\n", strrep("=", 70), "\n")
  cat("HEATMAPS PUBLICATION-READY GENERADOS EXITOSAMENTE!\n")
  cat(strrep("=", 70), "\n")

  cat("CONFIGURACIÓN UTILIZADA:\n")
  cat("- Paleta de colores: ", color_palette_option, " (", palette_names[min(color_palette_option, 6)], ")\n", sep="")
  cat("- Focalización: ", toupper(focus_mode), "\n")
  cat("- Filtrado inteligente: ", ifelse(enable_intelligent_filtering, "SÍ", "NO"), "\n")
  cat("- Umbral mínimo presencia: ", min_presence_threshold, " organismos\n")
  cat("- Umbral máximo presencia: ", max_presence_percentage, "%\n")
  cat("- Versión transpuesta: ", ifelse(generate_transposed, "SÍ", "NO"), "\n")
  cat("- Versión normal: ", ifelse(generate_normal, "SÍ", "NO"), "\n")
  cat("- Versiones focalizadas: ", ifelse(generate_focused, "SÍ", "NO"), "\n")
  cat("- Formato de salida: ", output_format, " (", output_dpi, " DPI)\n")

  cat("\nARCHIVOS GENERADOS:\n")
  files_created <- list.files(pattern = paste0("^", output_prefix), full.names = FALSE)
  for (f in files_created) {
    cat("- ", f, "\n")
  }

  cat("\n", strrep("=", 70), "\n")
  cat("RESUMEN FINAL DE ANÁLISIS\n")
  cat(strrep("=", 70), "\n")

  cat("DATOS PROCESADOS:\n")
  cat("- Organismos totales:", nrow(data), "\n")
  cat("- Organismos con genes funcionales:", nrow(data_filtered), "\n")
  cat("- Organismos en análisis final:", nrow(data_focused), "\n")
  cat("- Genes funcionales originales:", ncol(data) - 4, "\n")
  cat("- Genes en análisis final:", ncol(data_focused) - 4, "\n")

  dist_summary <- table(data_focused$final.clasification.group)
  cat("\nDISTRIBUCIÓN DE ORGANISMOS:\n")
  for (dist in names(dist_summary)) {
    cat("- ", stringr::str_to_title(dist), ": ", dist_summary[dist], " organismos\n", sep="")
  }

  cat("\nCATEGORÍAS FUNCIONALES EN ANÁLISIS:\n")
  if (ncol(data_focused) > 4) {
    final_categories <- table(sapply(colnames(data_focused)[5:ncol(data_focused)],
                                     categorize_function_enhanced))
    for (cat_name in names(final_categories)) {
      cat("- ", cat_name, ": ", final_categories[cat_name], " genes\n", sep="")
    }
  }

  cat("\n", strrep("=", 70), "\n")
  cat("RECOMENDACIONES PARA PUBLICACIÓN:\n")
  cat(strrep("=", 70), "\n")
  cat("1. FIGURA PRINCIPAL: Use heatmap focalizado (nitrogen, sulfur, o photosynthesis)\n")
  cat("2. MATERIAL SUPLEMENTARIO: Heatmap completo como referencia\n")
  cat("3. PALETA RECOMENDADA: Opción 4 o 6 para revistas Nature/Science\n")
  cat("4. FORMATO: PDF vectorial para mejor calidad en publicación\n")
  cat("5. RESOLUCIÓN: 600 DPI configurado automáticamente\n")
  cat("6. FIGURAS COMPUESTAS: Use versión normal para comparación de organismos\n")

  if (enable_intelligent_filtering) {
    cat("\nNOTA: Filtrado inteligente activado - genes muy raros o muy comunes excluidos\n")
    cat("      Esto mejora la interpretabilidad pero puede omitir genes de interés específico\n")
  }

  cat("\n", strrep("=", 70), "\n")
  cat("OPCIONES DE PALETAS DISPONIBLES:\n")
  cat(strrep("=", 70), "\n")
  cat("Cambie 'color_palette_option' para diferentes paletas:\n")
  cat("1. Biochemical Cycles - Colores por proceso metabólico\n")
  cat("2. Viridis Modified - Amigable para daltónicos\n")
  cat("3. Soft Professional - Colores suaves y pasteles\n")
  cat("4. Scientific Journal Style - Estilo Nature/Science (RECOMENDADO)\n")
  cat("5. Minimalist - Solo 4 colores principales, muy limpio\n")
  cat("6. Nature Optimized - NUEVO! Optimizado para publicaciones de alto impacto\n")

  cat("\nPara previsualizar paletas, descomente:\n")
  cat("# preview_color_palette(functional_colors_option4, 'Opción 4')\n")

} else {
  cat("ERROR: Archivo de entrada '", input_file, "' no encontrado!\n")
  cat("Asegúrese de que el archivo existe en el directorio actual.\n")
}

cat("\n", strrep("=", 70), "\n")
cat("SCRIPT COMPLETADO EXITOSAMENTE\n")
cat("Timestamp:", Sys.time(), "\n")
cat("Colores wide/rare mantenidos: wide =", group_colors["wide"], ", rare =", group_colors["rare"], "\n")
cat("Configuración optimizada para publicaciones de alto impacto\n")
cat("NUEVA FUNCIONALIDAD: Heatmaps normales y transpuestos + Figuras compuestas A-B-C\n")
cat(strrep("=", 70), "\n")

cat("\n¡Análisis completado! Revise los archivos generados.\n")
cat("NUEVAS FIGURAS DISPONIBLES:\n")
cat("- Heatmaps focalizados en AMBAS orientaciones (normal + transpuesto)\n")
cat("- Figuras compuestas con 3 paneles A-B-C para cada orientación\n")
cat("- Versión normal recomendada para publicación principal\n")
cat("- Versión transpuesta para análisis detallado de genes\n\n")
