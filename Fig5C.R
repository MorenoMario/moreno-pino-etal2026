# Fig5C.R
# Description: Generates compact antimicrobial resistance (AMR) heatmaps and bar plots from genus-level data.
# Purpose: Visualizes the distribution of antimicrobial resistance gene families across bacterial genera,
#          annotated by taxonomic family and biogeographic distribution (widespread vs. restricted).
#          Outputs include a clustered heatmap and a stacked bar plot summarizing resistance by group.

options(stringsAsFactors = FALSE)
suppressWarnings(suppressMessages({
  library(utils)
}))

cat("=== AMR HEATMAP GENERATOR (FORMATO COMPACTO) ===\n")
cat("Adaptado para datos de resistencia antimicrobiana\n")
cat("Directorio actual:", getwd(), "\n")

config <- list(
  input_file = "heatmap_AMR.txt",
  output_prefix = "amr_analysis_compact",
  abundance_type = "count",
  abundance_transform = "none",
  min_prevalence = 0.0,
  top_n_genera = 50,
  distance_method = "jaccard",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "average",

  cellwidth = 8,
  cellheight = 8,
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 9,
  angle_col = 90,
  border_color = NA,

  heatmap_colors = c("#FFFFFF", "#FFF3CD", "#FFE69C", "#FFD93D", "#FF8F00", "#E65100", "#BF360C"),
  group_colors = c("wide" = "#1976D2", "rare" = "#FF8A65"),
  dpi = 300,
  width_heatmap = 10,
  height_base = 0.15,
  max_height = 20,

  order_by_family = TRUE
)

args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
  if (grepl("^--top=", arg)) {
    val <- sub("^--top=", "", arg)
    config$top_n_genera <- if (tolower(val) %in% c("all", "inf", "na")) Inf else as.numeric(val)
  }
  if (grepl("^--transform=", arg)) {
    val <- tolower(sub("^--transform=", "", arg))
    allowed <- c("none", "log2", "log10", "sqrt")
    if (val %in% allowed) config$abundance_transform <- val
  }
  if (grepl("^--order-by-family=", arg)) {
    val <- tolower(sub("^--order-by-family=", "", arg))
    config$order_by_family <- val %in% c("1","true","t","yes","y")
  }
}

required_packages <- c(
  "pheatmap", "RColorBrewer", "ggplot2", "dplyr",
  "tidyr", "vegan", "scales", "dendextend",
  "stringr", "reshape2", "viridis", "tibble"
)
for (pkg in required_packages) {
  if (!suppressWarnings(require(pkg, character.only = TRUE))) {
    message(sprintf("[INFO] Instalando paquete '%s'...", pkg))
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
      library(pkg, character.only = TRUE)
    }, error = function(e) {
      warning(sprintf("[WARN] No se pudo instalar/cargar '%s': %s", pkg, e$message))
    })
  }
}

if (!suppressWarnings(require("Cairo", character.only = TRUE))) {
  message("[INFO] (opcional) Instala 'Cairo' para PDFs más compatibles en macOS: install.packages('Cairo')")
}

extract_genus_amr <- function(tax_string) {
  if (grepl("g__", tax_string)) {
    genus <- stringr::str_extract(tax_string, "g__[^;]+")
    genus <- gsub("g__", "", genus)
    if (!is.na(genus) && nzchar(genus) && !grepl("^\\s*$", genus)) return(genus)
  }

  if (grepl("f__", tax_string)) {
    family <- stringr::str_extract(tax_string, "f__[^;]+")
    family <- gsub("f__", "", family)
    if (!is.na(family) && nzchar(family) && !grepl("^\\s*$", family)) {
      return(paste0("Unclassified_", family))
    }
  }
  return("Unclassified")
}

extract_family_amr <- function(tax_string) {
  if (grepl("f__", tax_string)) {
    fam <- stringr::str_extract(tax_string, "f__[^;]+")
    fam <- gsub("f__", "", fam)
    if (!is.na(fam) && nzchar(fam) && !grepl("^\\s*$", fam)) return(fam)
  }
  return("Unknown")
}

apply_transform_amr <- function(mat, method) {
  method <- tolower(method)
  if (method == "none") return(mat)
  if (method == "log2")  return(log2(mat + 1))
  if (method == "log10") return(log10(mat + 1))
  if (method == "sqrt")  return(sqrt(mat))
  warning(sprintf("[WARN] Transformación no soportada '%s'", method))
  mat
}

load_and_process_amr_data <- function(filename, config) {
  cat("\n=== CARGANDO DATOS AMR ===\n")
  data <- read.delim(filename, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  cat("Dimensiones originales:", nrow(data), "filas ×", ncol(data), "columnas\n")

  name_col  <- which(names(data) == "name")
  tax_col   <- which(names(data) == "tax")
  group_col <- which(names(data) == "group")
  if (!length(name_col) || !length(tax_col) || !length(group_col)) {
    stop("El archivo debe contener las columnas: 'name', 'tax', 'group'.")
  }

  resistance_cols  <- setdiff(seq_len(ncol(data)), c(name_col, tax_col, group_col))
  resistance_names <- names(data)[resistance_cols]
  cat("Columnas de resistencia encontradas:", length(resistance_names), "\n")
  print(resistance_names)

  sample_names <- data[[name_col]]
  taxonomies   <- data[[tax_col]]
  group_info   <- data[[group_col]]
  names(group_info) <- sample_names

  resistance_matrix <- suppressWarnings(apply(as.matrix(data[, resistance_cols, drop = FALSE]), c(1, 2), as.numeric))
  rownames(resistance_matrix) <- sample_names
  colnames(resistance_matrix) <- resistance_names

  genera    <- sapply(taxonomies, extract_genus_amr)
  families  <- sapply(taxonomies, extract_family_amr)
  names(genera)   <- sample_names
  names(families) <- sample_names

  cat("\n=== AGREGANDO POR GÉNERO ===\n")
  df_for_agg <- as.data.frame(resistance_matrix)
  df_for_agg$Genus  <- genera
  df_for_agg$Group  <- group_info
  df_for_agg$Family <- families

  genus_resistance <- aggregate(df_for_agg[, resistance_names],
                                by = list(Genus = df_for_agg$Genus),
                                FUN = sum, na.rm = TRUE)
  rownames(genus_resistance) <- genus_resistance$Genus
  genus_matrix <- as.matrix(genus_resistance[, resistance_names, drop = FALSE])
  cat("Géneros únicos:", nrow(genus_matrix), "\n")

  genus_groups <- df_for_agg |>
    dplyr::group_by(Genus) |>
    dplyr::summarise(
      n_wide = sum(Group == "wide", na.rm = TRUE),
      n_rare = sum(Group == "rare", na.rm = TRUE),
      .groups = "drop"
    )
  genus_groups$Group <- ifelse(genus_groups$n_wide >= genus_groups$n_rare, "wide", "rare")
  genus_group_info <- genus_groups$Group
  names(genus_group_info) <- genus_groups$Genus

  fam_by_genus <- sapply(split(df_for_agg$Family, df_for_agg$Genus), function(fv) {
    if (length(fv) == 0) return("Unknown")
    tt <- sort(table(fv), decreasing = TRUE)
    names(tt)[1]
  })
  fam_by_genus <- as.character(fam_by_genus)

  list(
    genus_matrix     = genus_matrix,
    genus_groups     = genus_group_info,
    resistance_types = resistance_names,
    original_data    = data,
    genera_info      = genera,
    sample_groups    = group_info,
    genera_family    = fam_by_genus
  )
}

create_amr_heatmap_compact <- function(data_list, config) {
  cat("\n=== CREANDO HEATMAP AMR COMPACTO ===\n")

  genus_matrix <- data_list$genus_matrix

  heatmap_matrix <- apply_transform_amr(genus_matrix, config$abundance_transform)

  genus_totals <- rowSums(heatmap_matrix)
  top_n <- config$top_n_genera
  if (is.infinite(top_n)) {
    top_genera <- names(genus_totals)
  } else {
    top_n_selected <- min(as.integer(top_n), length(genus_totals))
    top_genera <- names(sort(genus_totals, decreasing = TRUE))[seq_len(top_n_selected)]
  }
  heatmap_matrix_filtered <- heatmap_matrix[top_genera, , drop = FALSE]
  cat(sprintf("Mostrando %d géneros con mayor resistencia\n", length(top_genera)))

  row_genera <- rownames(heatmap_matrix_filtered)
  dist_vec <- data_list$genus_groups[row_genera]
  fam_all  <- data_list$genera_family
  fam_vec  <- fam_all[row_genera]

  dist_vec[is.na(dist_vec) | dist_vec == ""] <- "unknown"
  fam_vec[is.na(fam_vec) | fam_vec == ""] <- "Unknown"

  annotation_row <- data.frame(
    Distribution = dist_vec,
    Family = fam_vec,
    row.names = row_genera,
    stringsAsFactors = FALSE
  )

  families_present <- sort(unique(annotation_row$Family))
  if (length(families_present) == 0) families_present <- "Unknown"
  kfam <- length(families_present)
  if (kfam <= 12) {
    fam_cols <- RColorBrewer::brewer.pal(max(3, kfam), "Set3")[seq_len(kfam)]
  } else {
    fam_cols <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(kfam)
  }
  names(fam_cols) <- families_present

  dist_present <- sort(unique(annotation_row$Distribution))
  dist_cols <- config$group_colors
  missing_dist <- setdiff(dist_present, names(dist_cols))
  if (length(missing_dist) > 0) {
    dist_cols <- c(dist_cols, stats::setNames(rep("#B0BEC5", length(missing_dist)), missing_dist))
  }

  ann_colors <- list(
    Distribution = dist_cols,
    Family = fam_cols
  )
  colors <- colorRampPalette(config$heatmap_colors)(50)

  if (isTRUE(config$order_by_family)) {
    ord <- order(annotation_row$Family)
    heatmap_matrix_filtered <- heatmap_matrix_filtered[ord, , drop = FALSE]
    annotation_row <- annotation_row[ord, , drop = FALSE]
  }

  n_genera <- nrow(heatmap_matrix_filtered)
  height_inches <- min(config$max_height, max(8, n_genera * config$height_base + 3))
  width_inches  <- config$width_heatmap

  filename_png <- paste0(config$output_prefix, "_heatmap.png")
  filename_pdf <- paste0(config$output_prefix, "_heatmap.pdf")

  generate_amr_heatmap <- function(filename, format = "png") {
    open_dev <- function() {
      if (format == "png") {
        png(filename,
            width = width_inches * config$dpi,
            height = height_inches * config$dpi,
            res = config$dpi)
      } else {
        if (capabilities("cairo")) {
          grDevices::cairo_pdf(filename, width = width_inches, height = height_inches, onefile = FALSE)
        } else {
          grDevices::pdf(filename, width = width_inches, height = height_inches, useDingbats = FALSE)
        }
      }
    }

    open_dev()
    on.exit({ try(grDevices::dev.off(), silent = TRUE) }, add = TRUE)

    cluster_cols_obj <- FALSE
    cluster_rows_obj <- FALSE
    if (ncol(heatmap_matrix_filtered) > 2 && isTRUE(config$cluster_cols)) {
      try({
        dist_cols <- stats::dist(t(heatmap_matrix_filtered))
        cluster_cols_obj <- stats::hclust(dist_cols, method = config$clustering_method)
      }, silent = TRUE)
    }
    if (nrow(heatmap_matrix_filtered) > 2 && isTRUE(config$cluster_rows)) {
      try({
        dist_rows <- stats::dist(heatmap_matrix_filtered)
        cluster_rows_obj <- stats::hclust(dist_rows, method = config$clustering_method)
      }, silent = TRUE)
    }

    ok <- TRUE
    tryCatch({
      pheatmap::pheatmap(
        heatmap_matrix_filtered,
        color = colors,
        cluster_rows = cluster_rows_obj,
        cluster_cols = cluster_cols_obj,
        annotation_row = annotation_row,
        annotation_colors = ann_colors,
        cellwidth = config$cellwidth,
        cellheight = config$cellheight,
        fontsize = config$fontsize,
        fontsize_row = config$fontsize_row,
        fontsize_col = config$fontsize_col,
        angle_col = config$angle_col,
        border_color = config$border_color,
        show_rownames = TRUE,
        show_colnames = TRUE,
        main = "Antimicrobial Resistance by Genus",
        treeheight_row = 30,
        treeheight_col = 30,
        annotation_legend = TRUE,
        annotation_names_row = FALSE,
        drop_levels = TRUE
      )
    }, error = function(e) {
      message("[WARN] pheatmap falló, usando fallback stats::heatmap(): ", e$message)
      ok <<- FALSE
    })

    if (!ok) {

      stats::heatmap(
        heatmap_matrix_filtered,
        scale = "none",
        margins = c(12, 15),
        cexRow = 0.6,
        cexCol = 0.8,
        main = "Antimicrobial Resistance by Genus"
      )
    }

    try({
      if (file.exists(filename)) {
        sz <- file.info(filename)$size
        if (is.finite(sz) && sz < 2048) {
          message(sprintf(
            "[WARN] El archivo quedó muy pequeño (%d bytes). Si Preview no abre, prueba con Adobe Reader o reinstala Cairo.",
            sz
          ))
        }
      }
    }, silent = TRUE)

    cat("✓ Guardado:", filename, "\n")
  }

  generate_amr_heatmap(filename_png, "png")
  generate_amr_heatmap(filename_pdf, "pdf")

  list(
    heatmap_matrix = heatmap_matrix_filtered,
    annotation_row = annotation_row
  )
}

create_resistance_barplot <- function(data_list, config) {
  cat("\n=== CREANDO GRÁFICO DE BARRAS DE RESISTENCIA ===\n")

  genus_matrix <- data_list$genus_matrix
  genus_groups <- data_list$genus_groups

  resistance_long <- genus_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Genus") %>%
    tidyr::pivot_longer(cols = -Genus, names_to = "Resistance_Type", values_to = "Count") %>%
    dplyr::mutate(Distribution = genus_groups[Genus])

  resistance_summary <- resistance_long %>%
    dplyr::group_by(Resistance_Type, Distribution) %>%
    dplyr::summarise(Total_Count = sum(Count, na.rm = TRUE), .groups = "drop")

  p_resistance <- ggplot2::ggplot(resistance_summary, ggplot2::aes(x = Resistance_Type, y = Total_Count, fill = Distribution)) +
    ggplot2::geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = config$group_colors) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      legend.position = "top"
    ) +
    ggplot2::labs(
      title = "Total Antimicrobial Resistance by Type and Distribution",
      x = "Resistance Type",
      y = "Total Count",
      fill = "Distribution"
    )

  filename_bar <- paste0(config$output_prefix, "_resistance_barplot.png")
  ggplot2::ggsave(filename_bar, p_resistance, width = 10, height = 6, dpi = config$dpi)
  cat("✓ Gráfico de barras guardado:", filename_bar, "\n")

  list(plot = p_resistance, data = resistance_summary)
}

run_amr_analysis <- function(config) {
  cat("\n", paste(rep("=", 70), collapse = ""), "\n")
  cat("ANÁLISIS AMR COMPACTO\n")
  cat(paste(rep("=", 70), collapse = ""), "\n")

  data_list <- load_and_process_amr_data(config$input_file, config)

  heatmap_result <- create_amr_heatmap_compact(data_list, config)
  barplot_result <- create_resistance_barplot(data_list, config)

  write.csv(data_list$genus_matrix, paste0(config$output_prefix, "_genus_resistance_matrix.csv"))
  write.csv(barplot_result$data, paste0(config$output_prefix, "_resistance_summary.csv"))

  summary_stats <- data.frame(
    Metric = c("Total_Genera", "Total_Resistance_Types", "Total_Samples",
               "Wide_Distribution_Genera", "Rare_Distribution_Genera"),
    Value = c(
      nrow(data_list$genus_matrix),
      ncol(data_list$genus_matrix),
      length(data_list$genera_info),
      sum(data_list$genus_groups == "wide", na.rm = TRUE),
      sum(data_list$genus_groups == "rare", na.rm = TRUE)
    )
  )
  write.csv(summary_stats, paste0(config$output_prefix, "_summary_stats.csv"), row.names = FALSE)

  cat("\n✔️ ANÁLISIS AMR COMPLETADO\n")
  cat("Archivos generados:\n")
  files <- list.files(pattern = paste0("^", config$output_prefix))
  for (f in files) cat(sprintf("- %s\n", f))

  invisible(list(data = data_list, heatmap = heatmap_result, barplot = barplot_result))
}

cat("\nConfiguración AMR:\n")
cat("- Archivo de entrada:", config$input_file, "\n")
cat("- Tipo de abundancia:", config$abundance_type, "\n")
cat("- Transformación:", config$abundance_transform, "\n")
cat("- Top géneros:", config$top_n_genera, "\n")
cat("- Ordenar por familia:", config$order_by_family, "\n\n")

if (!file.exists(config$input_file)) {
  cat("[ERROR] Archivo no encontrado:", config$input_file, "\n")
  cat("Asegúrate de que el archivo esté en el directorio actual.\n")
  quit(status = 1)
}

run_amr_analysis(config)

cat("\n🎯 Listo: heatmap con colores por FAMILIA + distribución wide/rare\n")
