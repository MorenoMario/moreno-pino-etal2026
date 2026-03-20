# Fig_s1.R
# Description: Generates taxonomic composition bar plots and abundance heatmaps for Chilean marine
#              microbiome samples with a fixed geographic sample order.
# Purpose: Aggregates metagenomic data at a configurable taxonomic level (default: genus), applies
#          log2 or z-score transformations, and produces: (1) a stacked bar plot of relative abundance
#          (Fig 2B), (2) a pheatmap heatmap with row clustering but fixed column order (Fig 2D), and
#          (3) an informational dendrogram of sample clustering (Fig 2C). Designed for mBio.

config <- list(
  input_file = "heatmap_revisado.txt",
  output_prefix = "test",

  taxonomic_level = "genus",

  ranking_method = "mean_abundance",
  min_prevalence = 0.0,

  top_n_figure_2b = 10,
  top_n_figure_2d = 90,

  abundance_type = "raw",
  abundance_transform = "log2",
  barplot_style = "nature_discrete",
  nature_color_style = "nature_discrete",
  barplot_order = "descending",

  distance_method = "bray",
  cluster_rows = TRUE,
  cluster_cols = TRUE,

  zscore_distance_rows = "euclidean",
  zscore_distance_cols = "euclidean",
  zscore_clustering_rows = "average",
  zscore_clustering_cols = "average",

  abundance_distance_rows = "bray",
  abundance_distance_cols = "bray",
  abundance_clustering_rows = "average",
  abundance_clustering_cols = "average",

  binary_distance = "jaccard",
  binary_clustering = "average",

  clip_quantile = 0.01,
  zscore_symmetric_limits = TRUE,
  zscore_manual_limits = c(-5, 5),
  zscore_legend_breaks = c(-5, -2, 0, 2, 5),

  add_region_annotation = TRUE,
  add_group_gaps = TRUE,

  cellwidth = 10,
  cellheight = 9,
  fontsize = 9,
  fontsize_row = 8,
  fontsize_col = 10,
  angle_col = 45,
  border_color = "white",
  dpi = 300,
  output_format = "pdf",
  width_heatmap_compact = 8,
  height_base = 0.18,
  max_height = 18,

  heatmap_colors_continuous = c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272",
                                "#FB6A4A", "#EF3B2C", "#CB181D", "#A50F15", "#67000D"),
  heatmap_colors_binary = c("#FFFFFF", "#2171B5"),
  group_colors = c("wide" = "#1F77B4", "rare" = "#FF7F0E"),
  wide_color_start = "#E8F4FD",
  wide_color_end = "#1F77B4",
  rare_color_start = "#FFF2CC",
  rare_color_end = "#FF7F0E"
)

SAMPLE_ORDER_FIXED <- c("Las Docas", "Algarrobo", "Navidad", "Topocalma",
                        "Pargua", "Los Chonos", "Ilque", "San Antonio")

cat("=== ANÁLISIS DE MICROBIOMA MARINO - COSTA CHILENA ===\n")
cat("Nivel taxonómico:", toupper(config$taxonomic_level), "\n")
cat("Transformación:", config$abundance_transform, "\n")
cat("Orden de muestras:", paste(SAMPLE_ORDER_FIXED, collapse = " → "), "\n")

required_packages <- c(
  "pheatmap", "RColorBrewer", "ggplot2", "dplyr",
  "tidyr", "vegan", "scales", "stringr", "viridis"
)

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Instalando", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

sitios_info <- data.frame(
  Sitio = SAMPLE_ORDER_FIXED,
  Latitud_aprox = c(-33.48, -33.35, -33.97, -34.15,
                    -41.88, -45.87, -42.46, -33.58),
  Region = c("North", "North", "Center", "Center",
             "South Exposed", "South Exposed", "South Protected", "South Protected"),
  Order_Position = 1:8,
  stringsAsFactors = FALSE
)

cat("✓ Información de sitios cargada con orden fijo establecido\n")
cat("- North:", paste(sitios_info$Sitio[sitios_info$Region == "North"], collapse = ", "), "\n")
cat("- Center:", paste(sitios_info$Sitio[sitios_info$Region == "Center"], collapse = ", "), "\n")
cat("- South Exposed:", paste(sitios_info$Sitio[sitios_info$Region == "South Exposed"], collapse = ", "), "\n")
cat("- South Protected:", paste(sitios_info$Sitio[sitios_info$Region == "South Protected"], collapse = ", "), "\n")

extract_taxonomy <- function(tax_string, level = "genus") {
  level_prefixes <- c(
    "kingdom" = "k__", "phylum" = "p__", "class" = "c__",
    "order" = "o__", "family" = "f__", "genus" = "g__", "species" = "s__"
  )

  if (!level %in% names(level_prefixes)) {
    stop("Nivel taxonómico no válido. Opciones: ", paste(names(level_prefixes), collapse = ", "))
  }

  prefix <- level_prefixes[[level]]

  if (grepl(prefix, tax_string)) {
    taxon <- stringr::str_extract(tax_string, paste0(prefix, "[^;]+"))
    taxon <- gsub(prefix, "", taxon)

    if (taxon == "" || is.na(taxon)) {
      level_order <- names(level_prefixes)
      current_index <- which(level_order == level)

      for (i in (current_index - 1):1) {
        if (i > 0) {
          upper_level <- level_order[i]
          upper_prefix <- level_prefixes[[upper_level]]

          if (grepl(upper_prefix, tax_string)) {
            upper_taxon <- stringr::str_extract(tax_string, paste0(upper_prefix, "[^;]+"))
            upper_taxon <- gsub(upper_prefix, "", upper_taxon)

            if (upper_taxon != "" && !is.na(upper_taxon)) {
              return(paste0("Unclassified_", upper_taxon))
            }
          }
        }
      }
      return(paste0("Unclassified_", level))
    }
    return(taxon)
  }
  return(paste0("Unclassified_", level))
}

get_level_name <- function(level) {
  level_names <- c(
    "kingdom" = "Reino", "phylum" = "Filo", "class" = "Clase", "order" = "Orden",
    "family" = "Familia", "genus" = "Género", "species" = "Especie"
  )
  return(level_names[[level]])
}

get_level_label_en <- function(level) {
  level_labels <- c(
    "kingdom" = "Kingdom", "phylum" = "Phylum", "class" = "Class", "order" = "Order",
    "family" = "Family", "genus" = "Genus", "species" = "Species"
  )
  return(level_labels[[level]])
}

get_output_params <- function(format) {
  switch(tolower(format),
         "pdf" = list(ext = ".pdf", width_unit = "in", height_unit = "in"),
         "png" = list(ext = ".png", width_unit = "px", height_unit = "px"),
         "both" = list(ext = c(".pdf", ".png"), width_unit = c("in", "px"), height_unit = c("in", "px")),
         list(ext = ".pdf", width_unit = "in", height_unit = "in"))
}

save_plot <- function(plot_obj, filename_base, width_val, height_val, dpi_val, format) {
  if (format == "both") {
    ggplot2::ggsave(paste0(filename_base, ".pdf"), plot_obj,
                    width = width_val, height = height_val,
                    units = "in", dpi = dpi_val, device = "pdf")
    cat("✓ Guardado PDF:", paste0(filename_base, ".pdf"), "\n")

    ggplot2::ggsave(paste0(filename_base, ".png"), plot_obj,
                    width = width_val, height = height_val,
                    units = "in", dpi = dpi_val, device = "png")
    cat("✓ Guardado PNG:", paste0(filename_base, ".png"), "\n")

  } else if (format == "pdf") {
    ggplot2::ggsave(paste0(filename_base, ".pdf"), plot_obj,
                    width = width_val, height = height_val,
                    units = "in", dpi = dpi_val, device = "pdf")
    cat("✓ Guardado PDF:", paste0(filename_base, ".pdf"), "\n")

  } else {
    ggplot2::ggsave(paste0(filename_base, ".png"), plot_obj,
                    width = width_val, height = height_val,
                    units = "in", dpi = dpi_val, device = "png")
    cat("✓ Guardado PNG:", paste0(filename_base, ".png"), "\n")
  }
}

create_graphics_device <- function(filename_base, width_val, height_val, dpi_val, format) {
  if (format == "both") {
    filename <- paste0(filename_base, ".pdf")
    pdf(filename, width = width_val, height = height_val)
    return(list(filename = filename, second_format = "png"))
  } else if (format == "pdf") {
    filename <- paste0(filename_base, ".pdf")
    pdf(filename, width = width_val, height = height_val)
    return(list(filename = filename, second_format = NULL))
  } else {
    filename <- paste0(filename_base, ".png")
    png(filename, width = width_val * dpi_val, height = height_val * dpi_val, res = dpi_val)
    return(list(filename = filename, second_format = NULL))
  }
}

rank_taxa_by_abundance <- function(taxa_matrix, method = "mean_abundance") {
  if (method == "mean_abundance") {
    taxa_scores <- rowMeans(taxa_matrix)
    rank_desc <- "Mean abundance across samples"
  } else if (method == "max_abundance") {
    taxa_scores <- apply(taxa_matrix, 1, max)
    rank_desc <- "Maximum abundance in any sample"
  } else if (method == "median_abundance") {
    taxa_scores <- apply(taxa_matrix, 1, median)
    rank_desc <- "Median abundance across samples"
  } else if (method == "prevalence_weighted") {
    prevalence <- rowSums(taxa_matrix > 0) / ncol(taxa_matrix)
    mean_abundance <- rowMeans(taxa_matrix)
    taxa_scores <- mean_abundance * prevalence
    rank_desc <- "Mean abundance weighted by prevalence"
  } else if (method == "variance") {
    taxa_scores <- apply(taxa_matrix, 1, var)
    rank_desc <- "Variance across samples"
  } else {
    taxa_scores <- rowSums(taxa_matrix > 0) / ncol(taxa_matrix)
    rank_desc <- "Prevalence"
  }

  taxa_ranking <- sort(taxa_scores, decreasing = TRUE)
  return(list(
    scores = taxa_scores,
    ranking = taxa_ranking,
    method = method,
    description = rank_desc
  ))
}

get_top_taxa <- function(taxa_matrix, n, ranking_method = "mean_abundance") {
  if (is.infinite(n) || n >= nrow(taxa_matrix)) {
    return(rownames(taxa_matrix))
  }

  ranking_result <- rank_taxa_by_abundance(taxa_matrix, ranking_method)
  top_taxa <- names(ranking_result$ranking)[1:min(n, length(ranking_result$ranking))]
  return(top_taxa)
}

generate_nature_colors <- function(n, style = "nature_discrete") {
  if (n <= 0) return(character(0))

  if (style == "nature_discrete") {
    nature_base <- c(
      "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
      "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85",
      "#F7DC6F", "#BB8FCE", "#85C1E9", "#F8C471", "#82E0AA",
      "#F1948A", "#D7BDE2", "#A9DFBF", "#F9E79F", "#AED6F1"
    )

    if (n <= length(nature_base)) {
      return(nature_base[1:n])
    } else {
      return(colorRampPalette(nature_base)(n))
    }
  } else {
    return(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n))
  }
}

generate_nature_gradient_colors <- function(abundances, distributions, config) {
  n <- length(abundances)
  colors <- character(n)

  if (max(abundances) > min(abundances)) {
    norm_abundances <- (abundances - min(abundances)) / (max(abundances) - min(abundances))
  } else {
    norm_abundances <- rep(0.5, n)
  }

  for (i in 1:n) {
    if (distributions[i] == "wide") {
      colors[i] <- colorRampPalette(c(config$wide_color_start, config$wide_color_end))(100)[
        max(1, min(100, round(norm_abundances[i] * 99 + 1)))
      ]
    } else {
      colors[i] <- colorRampPalette(c(config$rare_color_start, config$rare_color_end))(100)[
        max(1, min(100, round(norm_abundances[i] * 99 + 1)))
      ]
    }
  }
  return(colors)
}

generate_safe_colors <- function(n) {
  if (n <= 0) return(character(0))
  if (n == 1) return("#1F78B4")
  if (n == 2) return(c("#E31A1C", "#1F78B4"))
  if (n <= 8) return(RColorBrewer::brewer.pal(max(3, n), "Dark2")[1:n])
  if (n <= 11) return(RColorBrewer::brewer.pal(n, "Spectral"))

  base_colors <- c(
    RColorBrewer::brewer.pal(11, "Spectral"),
    RColorBrewer::brewer.pal(8, "Dark2")
  )

  if (n <= 19) return(base_colors[1:n])
  return(colorRampPalette(base_colors)(n))
}

apply_abundance_transform <- function(mat, method) {
  method_lower <- tolower(method)

  if (method_lower == "none") {
    return(mat)
  } else if (method_lower == "log2") {
    return(log2(mat + 1))
  } else if (method_lower == "log10") {
    return(log10(mat + 1))
  } else if (method_lower == "sqrt") {
    return(sqrt(mat))
  } else if (method_lower == "asinh") {
    return(asinh(mat))
  } else if (method_lower == "zscore") {
    zscore_mat <- t(apply(mat, 1, function(row) {
      if (sd(row) == 0) return(rep(0, length(row)))
      (row - mean(row)) / sd(row)
    }))
    colnames(zscore_mat) <- colnames(mat)
    rownames(zscore_mat) <- rownames(mat)
    return(zscore_mat)
  }
  return(mat)
}

get_diverging_palette <- function(n = 101) {

  return(colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(n))
}

compute_distance_matrix <- function(matrix_data, margin = "rows", transform_method = "log2") {

  mat_for_dist <- if (margin == "rows") matrix_data else t(matrix_data)

  if (tolower(config$abundance_type) == "binary") {

    return(vegan::vegdist(mat_for_dist, method = config$binary_distance))

  } else if (transform_method == "zscore") {

    if (margin == "rows") {
      distance_method <- config$zscore_distance_rows
    } else {
      distance_method <- config$zscore_distance_cols
    }

    if (distance_method == "correlation") {

      cor_mat <- suppressWarnings(cor(t(mat_for_dist), method = "spearman", use = "pairwise.complete.obs"))
      cor_mat[is.na(cor_mat)] <- 0
      return(as.dist(1 - cor_mat))
    } else {

      return(stats::dist(mat_for_dist, method = "euclidean"))
    }

  } else {

    if (margin == "rows") {
      distance_method <- config$abundance_distance_rows
    } else {
      distance_method <- config$abundance_distance_cols
    }

    if (distance_method %in% c("bray", "jaccard", "sorensen")) {
      return(vegan::vegdist(mat_for_dist, method = distance_method))
    } else {
      return(stats::dist(mat_for_dist, method = distance_method))
    }
  }
}

get_clustering_method <- function(transform_method = "log2", margin = "rows") {
  method <- ""

  if (tolower(config$abundance_type) == "binary") {
    method <- config$binary_clustering

  } else if (transform_method == "zscore") {
    if (margin == "rows") {
      method <- config$zscore_clustering_rows
    } else {
      method <- config$zscore_clustering_cols
    }

  } else {

    if (margin == "rows") {
      method <- config$abundance_clustering_rows
    } else {
      method <- config$abundance_clustering_cols
    }
  }

  if (method == "ward.D2") {
    method <- "ward"
  }

  return(method)
}

validate_clustering_method <- function(method) {
  valid_methods <- c("ward.D2", "ward", "single", "complete", "average", "mcquitty", "median", "centroid")

  if (!method %in% valid_methods) {
    warning("Método de clustering no válido: ", method, ". Usando 'complete' por defecto.")
    return("complete")
  }

  return(method)
}

configure_zscore_range <- function(min_val = -4, max_val = 4, legend_breaks = NULL) {
  cat("\n=== CONFIGURANDO RANGO DE Z-SCORE ===\n")

  if (min_val >= max_val) {
    stop("El valor mínimo debe ser menor que el máximo")
  }

  config$zscore_manual_limits <<- c(min_val, max_val)

  if (is.null(legend_breaks)) {

    if (abs(min_val) == abs(max_val)) {

      step <- (max_val - min_val) / 4
      config$zscore_legend_breaks <<- c(min_val, min_val + step, 0, max_val - step, max_val)
    } else {

      config$zscore_legend_breaks <<- seq(min_val, max_val, length.out = 5)
    }
  } else {
    config$zscore_legend_breaks <<- legend_breaks
  }

  cat("✓ Rango Z-score configurado:", min_val, "a", max_val, "\n")
  cat("✓ Breaks de leyenda:", paste(config$zscore_legend_breaks, collapse = ", "), "\n")
  cat("✓ Para aplicar los cambios, ejecuta el análisis nuevamente\n")
}

reset_zscore_auto <- function() {
  cat("\n=== RESTABLECIENDO Z-SCORE AUTOMÁTICO ===\n")

  config$zscore_manual_limits <<- NULL
  config$zscore_legend_breaks <<- NULL

  cat("✓ Z-score configurado en modo automático\n")
  cat("✓ Los límites se calcularán basándose en los datos (quantile", config$clip_quantile, ")\n")
}

change_clustering_config <- function(transform = NULL,
                                     zscore_dist_rows = NULL, zscore_dist_cols = NULL,
                                     zscore_clust_rows = NULL, zscore_clust_cols = NULL,
                                     abundance_dist_rows = NULL, abundance_dist_cols = NULL,
                                     abundance_clust_rows = NULL, abundance_clust_cols = NULL) {

  cat("\n=== CONFIGURANDO MÉTODOS DE CLUSTERING ===\n")

  if (!is.null(transform)) {
    valid_transforms <- c("none", "log2", "log10", "sqrt", "asinh", "zscore")
    if (transform %in% valid_transforms) {
      config$abundance_transform <<- transform
      cat("✓ Transformación cambiada a:", transform, "\n")
    } else {
      cat("⚠️ Transformación no válida. Opciones:", paste(valid_transforms, collapse = ", "), "\n")
    }
  }

  if (!is.null(zscore_dist_rows)) {
    if (zscore_dist_rows %in% c("euclidean", "correlation")) {
      config$zscore_distance_rows <<- zscore_dist_rows
      cat("✓ Distancia Z-score filas:", zscore_dist_rows, "\n")
    } else {
      cat("⚠️ Distancia Z-score no válida. Opciones: euclidean, correlation\n")
    }
  }

  if (!is.null(zscore_dist_cols)) {
    if (zscore_dist_cols %in% c("euclidean", "correlation")) {
      config$zscore_distance_cols <<- zscore_dist_cols
      cat("✓ Distancia Z-score columnas:", zscore_dist_cols, "\n")
    } else {
      cat("⚠️ Distancia Z-score no válida. Opciones: euclidean, correlation\n")
    }
  }

  if (!is.null(zscore_clust_rows)) {
    zscore_clust_rows <- validate_clustering_method(zscore_clust_rows)
    config$zscore_clustering_rows <<- zscore_clust_rows
    cat("✓ Clustering Z-score filas:", zscore_clust_rows, "\n")
  }

  if (!is.null(zscore_clust_cols)) {
    zscore_clust_cols <- validate_clustering_method(zscore_clust_cols)
    config$zscore_clustering_cols <<- zscore_clust_cols
    cat("✓ Clustering Z-score columnas:", zscore_clust_cols, "\n")
  }

  if (!is.null(abundance_dist_rows)) {
    if (abundance_dist_rows %in% c("bray", "euclidean", "manhattan", "jaccard", "sorensen")) {
      config$abundance_distance_rows <<- abundance_dist_rows
      cat("✓ Distancia abundancia filas:", abundance_dist_rows, "\n")
    } else {
      cat("⚠️ Distancia abundancia no válida. Opciones: bray, euclidean, manhattan, jaccard, sorensen\n")
    }
  }

  if (!is.null(abundance_dist_cols)) {
    if (abundance_dist_cols %in% c("bray", "euclidean", "manhattan", "jaccard", "sorensen")) {
      config$abundance_distance_cols <<- abundance_dist_cols
      cat("✓ Distancia abundancia columnas:", abundance_dist_cols, "\n")
    } else {
      cat("⚠️ Distancia abundancia no válida. Opciones: bray, euclidean, manhattan, jaccard, sorensen\n")
    }
  }

  if (!is.null(abundance_clust_rows)) {
    abundance_clust_rows <- validate_clustering_method(abundance_clust_rows)
    config$abundance_clustering_rows <<- abundance_clust_rows
    cat("✓ Clustering abundancia filas:", abundance_clust_rows, "\n")
  }

  if (!is.null(abundance_clust_cols)) {
    abundance_clust_cols <- validate_clustering_method(abundance_clust_cols)
    config$abundance_clustering_cols <<- abundance_clust_cols
    cat("✓ Clustering abundancia columnas:", abundance_clust_cols, "\n")
  }

  cat("\n📊 Configuración actual:\n")
  if (config$abundance_transform == "zscore") {
    cat("- Transformación: Z-score\n")
    cat("- Distancia filas:", config$zscore_distance_rows, "| Clustering filas:", config$zscore_clustering_rows, "\n")
    cat("- Distancia columnas:", config$zscore_distance_cols, "| Clustering columnas:", config$zscore_clustering_cols, "\n")
  } else {
    cat("- Transformación:", config$abundance_transform, "\n")
    cat("- Distancia filas:", config$abundance_distance_rows, "| Clustering filas:", config$abundance_clustering_rows, "\n")
    cat("- Distancia columnas:", config$abundance_distance_cols, "| Clustering columnas:", config$abundance_clustering_cols, "\n")
  }
}

show_clustering_config <- function() {
  cat("\n📊 CONFIGURACIÓN ACTUAL DE CLUSTERING:\n")
  cat("═══════════════════════════════════════════\n")

  cat("🧬 Transformación actual:", config$abundance_transform, "\n\n")

  cat("📈 Para Z-SCORE:\n")
  cat("  Filas (taxones):\n")
  cat("    - Distancia:", config$zscore_distance_rows, "\n")
  cat("    - Clustering:", config$zscore_clustering_rows, "\n")
  cat("  Columnas (muestras):\n")
  cat("    - Distancia:", config$zscore_distance_cols, "\n")
  cat("    - Clustering:", config$zscore_clustering_cols, "\n\n")

  cat("📊 Para ABUNDANCIAS RELATIVAS (log2, log10, etc.):\n")
  cat("  Filas (taxones):\n")
  cat("    - Distancia:", config$abundance_distance_rows, "\n")
  cat("    - Clustering:", config$abundance_clustering_rows, "\n")
  cat("  Columnas (muestras):\n")
  cat("    - Distancia:", config$abundance_distance_cols, "\n")
  cat("    - Clustering:", config$abundance_clustering_cols, "\n\n")

  cat("🔳 Para DATOS BINARIOS:\n")
  cat("  - Distancia:", config$binary_distance, "\n")
  cat("  - Clustering:", config$binary_clustering, "\n\n")

  cat("💡 Opciones disponibles:\n")
  cat("  Distancias Z-score: euclidean, correlation\n")
  cat("  Distancias abundancia: bray, euclidean, manhattan\n")
  cat("  Métodos clustering: average, complete, ward.D2\n")
}

cat("\n=== CARGANDO Y PROCESANDO DATOS ===\n")

data <- read.delim(config$input_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
cat("Dimensiones originales:", nrow(data), "filas ×", ncol(data), "columnas\n")

tax_col <- 1
group_col <- ncol(data)
sample_cols <- 2:(ncol(data) - 1)
sample_names <- trimws(names(data)[sample_cols])

missing_samples <- setdiff(SAMPLE_ORDER_FIXED, sample_names)
extra_samples <- setdiff(sample_names, SAMPLE_ORDER_FIXED)

if (length(missing_samples) > 0) {
  cat("⚠️ ADVERTENCIA: Muestras faltantes en el orden fijo:", paste(missing_samples, collapse = ", "), "\n")
}
if (length(extra_samples) > 0) {
  cat("⚠️ ADVERTENCIA: Muestras extra no incluidas en el orden fijo:", paste(extra_samples, collapse = ", "), "\n")
}

group_info <- data[[group_col]]
names(group_info) <- data[[tax_col]]

abundance_matrix <- as.matrix(data[, sample_cols, drop = FALSE])
rownames(abundance_matrix) <- data[[tax_col]]
colnames(abundance_matrix) <- sample_names
abundance_matrix <- apply(abundance_matrix, c(1,2), as.numeric)

samples_present <- intersect(SAMPLE_ORDER_FIXED, colnames(abundance_matrix))
abundance_matrix <- abundance_matrix[, samples_present, drop = FALSE]

cat("✓ Matriz de abundancias creada:", nrow(abundance_matrix), "×", ncol(abundance_matrix), "\n")
cat("✓ Orden de muestras aplicado:", paste(colnames(abundance_matrix), collapse = " → "), "\n")

cat("\n=== AGREGANDO POR NIVEL TAXONÓMICO:", toupper(config$taxonomic_level), "===\n")

taxa_at_level <- sapply(rownames(abundance_matrix), function(x) {
  extract_taxonomy(x, config$taxonomic_level)
})
names(taxa_at_level) <- rownames(abundance_matrix)

df_for_agg <- as.data.frame(abundance_matrix)
df_for_agg$Taxon <- taxa_at_level

taxon_matrix <- aggregate(. ~ Taxon, data = df_for_agg, FUN = sum)
rownames(taxon_matrix) <- taxon_matrix$Taxon
taxon_matrix <- as.matrix(taxon_matrix[, -1, drop = FALSE])

taxon_matrix <- taxon_matrix[, samples_present, drop = FALSE]

if (tolower(config$abundance_type) == "binary") {
  taxon_matrix[taxon_matrix > 0] <- 1
}

cat("Taxones únicos agregados:", nrow(taxon_matrix), "\n")

if (config$min_prevalence > 0) {
  prevalence <- rowSums(taxon_matrix > 0) / ncol(taxon_matrix)
  keep_taxa_basic <- prevalence >= config$min_prevalence
  taxon_matrix <- taxon_matrix[keep_taxa_basic, , drop = FALSE]
  cat("Taxones después del filtro de prevalencia:", nrow(taxon_matrix), "\n")
}

ranking_result <- rank_taxa_by_abundance(taxon_matrix, config$ranking_method)
cat("✓ Ranking generado. Top 5:", paste(names(ranking_result$ranking)[1:5], collapse = ", "), "\n")

cat("\n=== CREANDO FIGURA 2B: BARPLOT CON ORDEN FIJO DE MUESTRAS ===\n")

top_taxa_2b <- get_top_taxa(taxon_matrix, config$top_n_figure_2b, config$ranking_method)

taxon_distribution_type <- sapply(unique(rownames(taxon_matrix)), function(taxon) {
  taxa_of_level <- names(taxa_at_level)[taxa_at_level == taxon]
  distributions <- group_info[taxa_of_level]
  n_wide <- sum(distributions == "wide", na.rm = TRUE)
  n_rare <- sum(distributions == "rare", na.rm = TRUE)
  ifelse(n_wide >= n_rare, "wide", "rare")
}, USE.NAMES = TRUE)

taxon_df <- as.data.frame(taxon_matrix)
taxon_df$Taxon <- rownames(taxon_df)

taxon_long <- taxon_df %>%
  tidyr::pivot_longer(cols = -Taxon, names_to = "Sample", values_to = "Presence") %>%
  dplyr::mutate(Distribution = taxon_distribution_type[Taxon])

taxon_composition <- taxon_long %>%
  dplyr::group_by(Sample) %>%
  dplyr::mutate(Total = sum(Presence, na.rm = TRUE)) %>%
  dplyr::mutate(RelativeAbundance = ifelse(Total > 0, (Presence / Total) * 100, 0)) %>%
  dplyr::ungroup()

if (config$barplot_order == "ascending") {
  taxon_order <- names(sort(ranking_result$scores, decreasing = FALSE))
} else if (config$barplot_order == "descending") {
  taxon_order <- names(sort(ranking_result$scores, decreasing = TRUE))
} else {
  taxon_total_abundance <- taxon_composition %>%
    dplyr::group_by(Taxon, Distribution) %>%
    dplyr::summarise(TotalAbundance = sum(RelativeAbundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(Distribution, dplyr::desc(TotalAbundance))
  taxon_order <- taxon_total_abundance$Taxon
}

taxon_composition_grouped <- taxon_composition %>%
  dplyr::mutate(Taxon_grouped = ifelse(Taxon %in% top_taxa_2b, Taxon, "Others")) %>%
  dplyr::group_by(Sample, Taxon_grouped, Distribution) %>%
  dplyr::summarise(RelativeAbundance = sum(RelativeAbundance, na.rm = TRUE), .groups = "drop")

taxon_composition_grouped$Sample <- factor(taxon_composition_grouped$Sample, levels = samples_present)

unique_taxa <- unique(taxon_composition_grouped$Taxon_grouped)
taxon_levels <- intersect(taxon_order, unique_taxa)
if ("Others" %in% unique_taxa) {
  taxon_levels <- c(taxon_levels, "Others")
}
taxon_composition_grouped$Taxon_grouped <- factor(taxon_composition_grouped$Taxon_grouped, levels = taxon_levels)

if (config$barplot_style == "nature_discrete") {
  n_colors <- length(taxon_levels)
  colors <- generate_nature_colors(n_colors, config$nature_color_style)
  if ("Others" %in% taxon_levels) {
    colors[taxon_levels == "Others"] <- "#CCCCCC"
  }
  names(colors) <- taxon_levels
} else if (config$barplot_style == "nature_gradient") {
  selected_scores <- ranking_result$scores[intersect(taxon_levels, names(ranking_result$scores))]
  selected_distributions <- taxon_distribution_type[names(selected_scores)]
  colors <- generate_nature_gradient_colors(selected_scores, selected_distributions, config)
  if ("Others" %in% taxon_levels) colors <- c(colors, "#CCCCCC")
  names(colors) <- taxon_levels
} else {
  colors <- generate_safe_colors(length(taxon_levels))
  if ("Others" %in% taxon_levels) colors[taxon_levels == "Others"] <- "gray75"
  names(colors) <- taxon_levels
}

level_name_en <- get_level_label_en(config$taxonomic_level)

p_2b <- ggplot2::ggplot(taxon_composition_grouped, aes(x = Sample, y = RelativeAbundance, fill = Taxon_grouped)) +
  ggplot2::geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.1) +
  ggplot2::scale_fill_manual(values = colors, name = level_name_en) +
  ggplot2::scale_x_discrete(drop = FALSE) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::labs(
    title = paste("Taxonomic composition of marine microbiome -", level_name_en, "level"),
    subtitle = sprintf("Top %d %s (ranked by %s, %s style) - Fixed sample order",
                       length(top_taxa_2b), tolower(level_name_en), config$ranking_method, config$barplot_style),
    x = "Sampling sites (ordered: Las Docas → Algarrobo → Navidad → Topocalma → Pargua → Los Chonos → Ilque → San Antonio)",
    y = "Relative abundance (%)"
  ) +
  ggplot2::theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 10, color = "black", face = "bold"),
    axis.title.x = element_text(size = 8),
    legend.text = element_text(size = 8, face = "italic", color = "black"),
    legend.title = element_text(size = 9, face = "bold", color = "black"),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(size = 12, face = "bold", color = "black"),
    plot.subtitle = element_text(size = 9, color = "gray30"),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank()
  ) +
  ggplot2::scale_y_continuous(expand = c(0, 0), labels = function(x) paste0(x, "%"))

filename_2b_base <- paste0(config$output_prefix, "_figure_2B_", config$taxonomic_level, "_",
                           config$ranking_method, "_", config$barplot_style, "_FIXED_ORDER")
save_plot(p_2b, filename_2b_base, width_val = 14, height_val = 6, dpi_val = config$dpi, format = config$output_format)

cat("\n=== CREANDO FIGURA 2D: HEATMAP CON ORDEN FIJO DE MUESTRAS ===\n")

top_taxa_2d <- get_top_taxa(taxon_matrix, config$top_n_figure_2d, config$ranking_method)

if (tolower(config$abundance_type) == "binary") {
  heatmap_matrix <- taxon_matrix[top_taxa_2d, , drop = FALSE]
} else {
  heatmap_matrix <- apply_abundance_transform(taxon_matrix[top_taxa_2d, , drop = FALSE], config$abundance_transform)
}

heatmap_matrix <- heatmap_matrix[, samples_present, drop = FALSE]

taxon_distribution_heatmap <- data.frame(
  Distribution = taxon_distribution_type[top_taxa_2d],
  row.names = top_taxa_2d,
  stringsAsFactors = FALSE
)

annotation_col <- NULL
if (config$add_region_annotation) {
  annotation_col <- data.frame(
    Region = factor(
      sitios_info$Region[match(colnames(heatmap_matrix), sitios_info$Sitio)],
      levels = c("North", "Center", "South Exposed", "South Protected")
    )
  )
  rownames(annotation_col) <- colnames(heatmap_matrix)
}

ann_colors <- list(Distribution = config$group_colors)
if (!is.null(annotation_col)) {
  ann_colors$Region <- c(
    "North" = "#E64B35",
    "Center" = "#4DBBD5",
    "South Exposed" = "#00A087",
    "South Protected" = "#3C5488"
  )
}

breaks_vec <- NULL
main_title <- ""

if (tolower(config$abundance_type) == "binary") {
  colors_heatmap <- config$heatmap_colors_binary
  main_title <- sprintf("%s presence/absence (top %d by %s) - Fixed sample order",
                        level_name_en, length(top_taxa_2d), config$ranking_method)
} else if (config$abundance_transform == "zscore") {

  z_vals <- as.numeric(heatmap_matrix)
  z_vals <- z_vals[is.finite(z_vals)]

  if (!is.null(config$zscore_manual_limits) && length(config$zscore_manual_limits) == 2) {

    lim_min <- config$zscore_manual_limits[1]
    lim_max <- config$zscore_manual_limits[2]
    cat("✓ Usando límites manuales Z-score:", lim_min, "a", lim_max, "\n")

  } else {

    lim <- stats::quantile(abs(z_vals), probs = config$clip_quantile, na.rm = TRUE)
    lim <- max(lim, 1)

    if (config$zscore_symmetric_limits) {
      lim_min <- -lim
      lim_max <- lim
    } else {
      lim_min <- min(z_vals, na.rm = TRUE)
      lim_max <- max(z_vals, na.rm = TRUE)
    }
    cat("✓ Usando límites automáticos Z-score:", round(lim_min, 2), "a", round(lim_max, 2), "\n")
  }

  heatmap_for_colors <- pmin(pmax(heatmap_matrix, lim_min), lim_max)

  n_breaks <- 101
  breaks_vec <- seq(lim_min, lim_max, length.out = n_breaks)
  colors_heatmap <- get_diverging_palette(length(breaks_vec) - 1)

  main_title <- sprintf("%s Z-score (top %d by %s) — range [%.1f, %.1f] - Fixed sample order",
                        level_name_en, length(top_taxa_2d), config$ranking_method, lim_min, lim_max)
} else {
  colors_heatmap <- colorRampPalette(config$heatmap_colors_continuous)(50)

  transform_label <- switch(config$abundance_transform,
                            "log2" = "log2(x+1)",
                            "log10" = "log10(x+1)",
                            "sqrt" = "√x",
                            "asinh" = "asinh(x)",
                            "none" = "raw")

  main_title <- sprintf("%s abundance %s (top %d by %s) - Fixed sample order",
                        level_name_en, transform_label, length(top_taxa_2d), config$ranking_method)
}

cluster_cols_obj <- FALSE
cluster_rows_obj <- FALSE

if (nrow(heatmap_matrix) > 2 && config$cluster_rows) {
  tryCatch({
    dist_rows <- compute_distance_matrix(heatmap_matrix, margin = "rows", transform_method = config$abundance_transform)
    clustering_method_rows <- get_clustering_method(config$abundance_transform, margin = "rows")
    clustering_method_rows <- validate_clustering_method(clustering_method_rows)

    cluster_rows_obj <- hclust(dist_rows, method = clustering_method_rows)

    cat("✓ Clustering filas - Distancia:",
        if (config$abundance_transform == "zscore") config$zscore_distance_rows else
          if (tolower(config$abundance_type) == "binary") config$binary_distance else config$abundance_distance_rows,
        ", Método:", clustering_method_rows, "\n")

  }, error = function(e) {
    cat("[WARN] No se pudo hacer clustering de filas:", e$message, "\n")
    cluster_rows_obj <- FALSE
  })
}

cat("✓ Clustering de columnas DESACTIVADO para mantener orden fijo de muestras\n")

gaps_row <- NULL
if (config$add_group_gaps) {

  dist_order <- order(taxon_distribution_heatmap$Distribution, rownames(taxon_distribution_heatmap))
  heatmap_matrix <- heatmap_matrix[dist_order, , drop = FALSE]
  if (exists("heatmap_for_colors")) {
    heatmap_for_colors <- heatmap_for_colors[dist_order, , drop = FALSE]
  }
  taxon_distribution_heatmap <- taxon_distribution_heatmap[dist_order, , drop = FALSE]

  group_vector <- as.character(taxon_distribution_heatmap$Distribution)
  change_positions <- which(group_vector[-1] != group_vector[-length(group_vector)])
  if (length(change_positions) > 0) {
    gaps_row <- change_positions
  }
}

filename_2d_base <- paste0(config$output_prefix, "_figure_2D_", config$taxonomic_level, "_", config$ranking_method, "_FIXED_ORDER")
width_inches_heatmap <- config$width_heatmap_compact
height_inches_heatmap <- min(config$max_height, max(8, nrow(heatmap_matrix) * config$height_base + 3))

device_info <- create_graphics_device(filename_2d_base, width_inches_heatmap, height_inches_heatmap,
                                      config$dpi, config$output_format)

tryCatch({

  matrix_for_plot <- if (exists("heatmap_for_colors")) heatmap_for_colors else heatmap_matrix

  pheatmap::pheatmap(
    matrix_for_plot,
    color = colors_heatmap,
    breaks = breaks_vec,
    cluster_rows = cluster_rows_obj,
    cluster_cols = cluster_cols_obj,
    annotation_row = taxon_distribution_heatmap,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    cellwidth = config$cellwidth,
    cellheight = config$cellheight,
    fontsize = config$fontsize,
    fontsize_row = if (nrow(heatmap_matrix) > 40) 7 else config$fontsize_row,
    fontsize_col = config$fontsize_col,
    angle_col = config$angle_col,
    border_color = config$border_color,
    show_rownames = TRUE,
    show_colnames = TRUE,
    main = main_title,
    treeheight_row = 30,
    treeheight_col = 0,
    annotation_legend = TRUE,
    annotation_names_row = FALSE,
    annotation_names_col = FALSE,
    gaps_row = gaps_row,
    legend_breaks = if (config$abundance_transform == "zscore" && !is.null(config$zscore_legend_breaks)) config$zscore_legend_breaks else NULL,
    legend_labels = if (config$abundance_transform == "zscore" && !is.null(config$zscore_legend_breaks)) as.character(config$zscore_legend_breaks) else NULL
  )
}, error = function(e) {
  cat("[ERROR] Problema con pheatmap:", e$message, "\n")
  cat("Creando heatmap básico...\n")

  matrix_for_plot <- if (exists("heatmap_for_colors")) heatmap_for_colors else heatmap_matrix
  heatmap(matrix_for_plot, scale = "none", margins = c(8, 12), cexRow = 0.6, cexCol = 0.8, Colv = NA)
  title(main_title, cex.main = 1.2)
})

dev.off()
cat("✓ Figura 2D guardada:", device_info$filename, "\n")

if (!is.null(device_info$second_format)) {
  filename_png <- paste0(filename_2d_base, ".png")
  png(filename_png,
      width = width_inches_heatmap * config$dpi,
      height = height_inches_heatmap * config$dpi,
      res = config$dpi)

  tryCatch({
    matrix_for_plot <- if (exists("heatmap_for_colors")) heatmap_for_colors else heatmap_matrix

    pheatmap::pheatmap(
      matrix_for_plot,
      color = colors_heatmap,
      breaks = breaks_vec,
      cluster_rows = cluster_rows_obj,
      cluster_cols = cluster_cols_obj,
      annotation_row = taxon_distribution_heatmap,
      annotation_col = annotation_col,
      annotation_colors = ann_colors,
      cellwidth = config$cellwidth,
      cellheight = config$cellheight,
      fontsize = config$fontsize,
      fontsize_row = if (nrow(heatmap_matrix) > 40) 7 else config$fontsize_row,
      fontsize_col = config$fontsize_col,
      angle_col = config$angle_col,
      border_color = config$border_color,
      show_rownames = TRUE,
      show_colnames = TRUE,
      main = main_title,
      treeheight_row = 30,
      treeheight_col = 0,
      annotation_legend = TRUE,
      annotation_names_row = FALSE,
      annotation_names_col = FALSE,
      gaps_row = gaps_row,
      legend_breaks = if (config$abundance_transform == "zscore" && !is.null(config$zscore_legend_breaks)) config$zscore_legend_breaks else NULL,
      legend_labels = if (config$abundance_transform == "zscore" && !is.null(config$zscore_legend_breaks)) as.character(config$zscore_legend_breaks) else NULL
    )
  }, error = function(e) {
    cat("[ERROR] Problema con pheatmap PNG:", e$message, "\n")
    matrix_for_plot <- if (exists("heatmap_for_colors")) heatmap_for_colors else heatmap_matrix
    heatmap(matrix_for_plot, scale = "none", margins = c(8, 12), cexRow = 0.6, cexCol = 0.8, Colv = NA)
    title(main_title, cex.main = 1.2)
  })

  dev.off()
  cat("✓ Figura 2D guardada:", filename_png, "\n")
}

cat("\n=== CREANDO FIGURA 2C: DENDROGRAMA DE MUESTRAS (INFORMACIÓN SOLAMENTE) ===\n")

if (tolower(config$abundance_type) == "binary") {
  dist_matrix <- vegan::vegdist(t(taxon_matrix), method = config$binary_distance)
  clustering_method_samples <- validate_clustering_method(config$binary_clustering)
} else if (config$abundance_transform == "zscore") {

  if (config$zscore_distance_cols == "correlation") {
    cor_mat <- suppressWarnings(cor(taxon_matrix, method = "spearman", use = "pairwise.complete.obs"))
    cor_mat[is.na(cor_mat)] <- 0
    dist_matrix <- as.dist(1 - cor_mat)
  } else {
    dist_matrix <- stats::dist(t(taxon_matrix), method = "euclidean")
  }
  clustering_method_samples <- validate_clustering_method(config$zscore_clustering_cols)
} else {

  if (config$abundance_distance_cols %in% c("bray", "jaccard", "sorensen")) {
    dist_matrix <- vegan::vegdist(t(taxon_matrix), method = config$abundance_distance_cols)
  } else {
    dist_matrix <- stats::dist(t(taxon_matrix), method = config$abundance_distance_cols)
  }
  clustering_method_samples <- validate_clustering_method(config$abundance_clustering_cols)
}

hclust_result <- hclust(dist_matrix, method = clustering_method_samples)

cat("✓ Dendrograma - Distancia:",
    if (config$abundance_transform == "zscore") config$zscore_distance_cols else
      if (tolower(config$abundance_type) == "binary") config$binary_distance else config$abundance_distance_cols,
    ", Método:", clustering_method_samples, "\n")
cat("ℹ️ NOTA: Este dendrograma es solo informativo. En las figuras principales se mantiene el orden fijo.\n")

filename_2c_base <- paste0(config$output_prefix, "_figure_2C_", config$taxonomic_level, "_INFO_ONLY")
width_inches_dendro <- 8
height_inches_dendro <- 5

device_info <- create_graphics_device(filename_2c_base, width_inches_dendro, height_inches_dendro,
                                      config$dpi, config$output_format)

plot(hclust_result,
     main = paste("Sample clustering -", level_name_en, "level (INFORMATIVO - no usado en figuras)"),
     xlab = "Orden real en figuras: Las Docas → Algarrobo → Navidad → Topocalma → Pargua → Los Chonos → Ilque → San Antonio",
     sub = "",
     cex.main = 0.9)
dev.off()
cat("✓ Figura 2C guardada:", device_info$filename, "\n")

if (!is.null(device_info$second_format)) {
  filename_png <- paste0(filename_2c_base, ".png")
  png(filename_png,
      width = width_inches_dendro * config$dpi,
      height = height_inches_dendro * config$dpi,
      res = config$dpi)

  plot(hclust_result,
       main = paste("Sample clustering -", level_name_en, "level (INFORMATIVO - no usado en figuras)"),
       xlab = "Orden real en figuras: Las Docas → Algarrobo → Navidad → Topocalma → Pargua → Los Chonos → Ilque → San Antonio",
       sub = "",
       cex.main = 0.9)
  dev.off()
  cat("✓ Figura 2C guardada:", filename_png, "\n")
}

cat("\n=== EXPORTANDO DATOS Y REPORTES ===\n")

write.csv(taxon_matrix, paste0(config$output_prefix, "_", config$taxonomic_level, "_matrix_full_FIXED_ORDER.csv"))
cat("✓ Matriz completa guardada con orden fijo\n")

write.csv(taxon_composition_grouped,
          paste0(config$output_prefix, "_composition_", config$taxonomic_level, "_", config$ranking_method, "_FIXED_ORDER.csv"))
cat("✓ Datos de composición guardados con orden fijo\n")

taxon_stats <- data.frame(
  Taxon = names(ranking_result$scores),
  Level = config$taxonomic_level,
  Primary_Score = ranking_result$scores[names(ranking_result$scores)],
  Mean_Abundance = rowMeans(taxon_matrix),
  Max_Abundance = apply(taxon_matrix, 1, max),
  Median_Abundance = apply(taxon_matrix, 1, median),
  Std_Dev = apply(taxon_matrix, 1, sd),
  Prevalence = rowSums(taxon_matrix > 0) / ncol(taxon_matrix),
  Distribution = taxon_distribution_type[names(ranking_result$scores)],
  Rank = rank(-ranking_result$scores, ties.method = "min"),
  stringsAsFactors = FALSE
)
taxon_stats <- taxon_stats[order(taxon_stats$Rank), ]
write.csv(taxon_stats,
          paste0(config$output_prefix, "_", config$taxonomic_level, "_ranking_", config$ranking_method, ".csv"),
          row.names = FALSE)
cat("✓ Estadísticas de taxones guardadas\n")

overlap_taxa <- intersect(top_taxa_2b, top_taxa_2d)
consistency_report <- data.frame(
  Figure = c("2B_Barplot", "2D_Heatmap"),
  Taxonomic_Level = config$taxonomic_level,
  Top_N_Used = c(length(top_taxa_2b), length(top_taxa_2d)),
  Ranking_Method = config$ranking_method,
  Color_Style = config$barplot_style,
  Taxa_Overlap = length(overlap_taxa),
  Percentage_Overlap = round(length(overlap_taxa) / max(length(top_taxa_2b), length(top_taxa_2d)) * 100, 1),
  Sample_Order_Fixed = TRUE,
  Sample_Order = paste(samples_present, collapse = " → "),
  stringsAsFactors = FALSE
)
write.csv(consistency_report,
          paste0(config$output_prefix, "_consistency_", config$taxonomic_level, "_FIXED_ORDER.csv"),
          row.names = FALSE)
cat("✓ Reporte de consistencia guardado con información de orden fijo\n")

sitios_info_ordered <- sitios_info[sitios_info$Sitio %in% samples_present, ]
sitios_info_ordered <- sitios_info_ordered[match(samples_present, sitios_info_ordered$Sitio), ]
write.csv(sitios_info_ordered, paste0(config$output_prefix, "_sitios_muestreo_FIXED_ORDER.csv"), row.names = FALSE)
cat("✓ Información de sitios guardada con orden fijo aplicado\n")

generate_taxonomic_comparison <- function(levels_to_compare = c("phylum", "class", "family", "genus")) {
  cat("\n=== GENERANDO COMPARACIÓN DE NIVELES TAXONÓMICOS ===\n")

  comparison_results <- list()
  original_level <- config$taxonomic_level

  for (level in levels_to_compare) {
    cat("- Procesando nivel:", level, "\n")

    tryCatch({

      taxa_level <- sapply(rownames(abundance_matrix), function(x) extract_taxonomy(x, level))

      df_temp <- as.data.frame(abundance_matrix)
      df_temp$Taxon <- taxa_level

      taxon_temp <- aggregate(. ~ Taxon, data = df_temp, FUN = sum)
      rownames(taxon_temp) <- taxon_temp$Taxon
      taxon_temp <- as.matrix(taxon_temp[, -1, drop = FALSE])

      taxon_temp <- taxon_temp[, samples_present, drop = FALSE]

      if (config$min_prevalence > 0) {
        prevalence_temp <- rowSums(taxon_temp > 0) / ncol(taxon_temp)
        keep_temp <- prevalence_temp >= config$min_prevalence
        taxon_temp <- taxon_temp[keep_temp, , drop = FALSE]
      }

      taxon_dist_temp <- sapply(unique(rownames(taxon_temp)), function(taxon) {
        taxa_of_level <- names(taxa_level)[taxa_level == taxon]
        distributions <- group_info[taxa_of_level]
        n_wide <- sum(distributions == "wide", na.rm = TRUE)
        n_rare <- sum(distributions == "rare", na.rm = TRUE)
        ifelse(n_wide >= n_rare, "wide", "rare")
      }, USE.NAMES = TRUE)

      comparison_results[[level]] <- list(
        n_taxa = nrow(taxon_temp),
        n_wide = sum(taxon_dist_temp == "wide", na.rm = TRUE),
        n_rare = sum(taxon_dist_temp == "rare", na.rm = TRUE)
      )

    }, error = function(e) {
      cat("   [ERROR] en nivel", level, ":", e$message, "\n")
    })
  }

  comparison_summary <- data.frame(
    Level = names(comparison_results),
    N_Taxa = sapply(comparison_results, function(x) x$n_taxa),
    N_Wide = sapply(comparison_results, function(x) x$n_wide),
    N_Rare = sapply(comparison_results, function(x) x$n_rare),
    Sample_Order_Fixed = TRUE,
    stringsAsFactors = FALSE
  )

  write.csv(comparison_summary,
            paste0(config$output_prefix, "_taxonomic_levels_summary_FIXED_ORDER.csv"),
            row.names = FALSE)

  cat("✓ Comparación de niveles taxonómicos completada con orden fijo\n")
  print(comparison_summary)

  return(comparison_results)
}

change_taxonomic_level <- function(new_level) {
  valid_levels <- c("kingdom", "phylum", "class", "order", "family", "genus", "species")
  if (!new_level %in% valid_levels) {
    stop("Nivel no válido. Opciones: ", paste(valid_levels, collapse = ", "))
  }

  cat("Cambiando nivel taxonómico de", config$taxonomic_level, "a", new_level, "\n")
  config$taxonomic_level <<- new_level
  cat("✓ Configuración actualizada. Ejecuta el script nuevamente.\n")
  cat("ℹ️ El orden fijo de muestras se mantendrá:", paste(SAMPLE_ORDER_FIXED, collapse = " → "), "\n")
}

cat("\n✔️ ANÁLISIS DE MICROBIOMA MARINO COMPLETADO CON ORDEN FIJO!\n")
cat("\n📊 RESUMEN ESPECÍFICO PARA TU DATASET:\n")
cat("- Nivel taxonómico analizado:", toupper(config$taxonomic_level), "(", get_level_name(config$taxonomic_level), ")\n")
cat("- Taxones totales disponibles:", nrow(taxon_matrix), "\n")
cat("- Taxones clasificados como 'wide':", sum(taxon_distribution_type == "wide", na.rm = TRUE), "\n")
cat("- Taxones clasificados como 'rare':", sum(taxon_distribution_type == "rare", na.rm = TRUE), "\n")
cat("- Top taxones en barplot 2B:", length(top_taxa_2b), "\n")
cat("- Top taxones en heatmap 2D:", length(top_taxa_2d), "\n")
cat("- Taxones compartidos entre figuras:", length(overlap_taxa), "\n")
cat("- Porcentaje de solapamiento:", round(length(overlap_taxa)/max(length(top_taxa_2b), length(top_taxa_2d))*100, 1), "%\n")

cat("\n🔄 ORDEN FIJO DE MUESTRAS APLICADO:\n")
cat("- Orden establecido:", paste(samples_present, collapse = " → "), "\n")
cat("- ✓ Barplot: Orden fijo mantenido\n")
cat("- ✓ Heatmap: Sin clustering de columnas, orden fijo preservado\n")
cat("- ✓ Dendrograma: Solo informativo, no afecta el orden en las figuras\n")

cat("\n🎨 CONFIGURACIÓN APLICADA:\n")
cat("- Nivel taxonómico:", config$taxonomic_level, "\n")
cat("- Transformación:", config$abundance_transform, "\n")
cat("- Método de ranking:", config$ranking_method, "\n")
cat("- Estilo de colores:", config$barplot_style, "\n")
cat("- Formato de salida:", config$output_format, "\n")
cat("- Clustering de filas (taxones):", ifelse(config$cluster_rows, "ACTIVADO", "DESACTIVADO"), "\n")
cat("- Clustering de columnas (muestras): DESACTIVADO (orden fijo)", "\n")

cat("\n🌊 SITIOS ANALIZADOS POR REGIÓN (EN ORDEN FIJO):\n")
for (site in samples_present) {
  region <- sitios_info$Region[sitios_info$Sitio == site]
  position <- which(samples_present == site)
  cat(sprintf("  %d. %s (%s)\n", position, site, region))
}

cat("\n📁 ARCHIVOS GENERADOS CON ORDEN FIJO:\n")
files_generated <- list.files(pattern = paste0("^", config$output_prefix), full.names = FALSE)
level_files <- files_generated[grepl(config$taxonomic_level, files_generated)]
fixed_order_files <- level_files[grepl("FIXED_ORDER", level_files)]
for (f in fixed_order_files) {
  cat("- ", f, "\n")
}

cat("\n💡 NUEVAS FUNCIONES DISPONIBLES (CON ORDEN FIJO):\n")
cat("- configure_zscore_range(): Configurar rango personalizado de Z-score\n")
cat("- reset_zscore_auto(): Restablecer a modo automático\n")
cat("- show_clustering_config(): Mostrar configuración actual de clustering\n")
cat("- change_clustering_config(): Cambiar métodos de clustering\n")
cat("- change_taxonomic_level(): Cambiar nivel taxonómico (mantiene orden fijo)\n")
cat("- generate_taxonomic_comparison(): Comparar múltiples niveles (con orden fijo)\n")

cat("\n🎨 EJEMPLOS DE CONFIGURACIÓN DE Z-SCORE:\n")
cat("# Configurar rango -4 a +4 (como solicitas):\n")
cat("configure_zscore_range(-4, 4)\n\n")

cat("# Configurar rango personalizado con breaks específicos:\n")
cat("configure_zscore_range(-3, 3, legend_breaks = c(-3, -1.5, 0, 1.5, 3))\n\n")

cat("# Rango asimétrico:\n")
cat("configure_zscore_range(-2, 5)\n\n")

cat("# Volver a modo automático:\n")
cat("reset_zscore_auto()\n\n")

cat("\n🔧 EJEMPLOS DE CONFIGURACIÓN DE CLUSTERING:\n")
cat("# Para Z-score con euclidiana y complete:\n")
cat("change_clustering_config(transform = 'zscore',\n")
cat("                        zscore_dist_rows = 'euclidean',\n")
cat("                        zscore_clust_rows = 'complete')\n\n")

cat("# Para abundancias relativas con ward:\n")
cat("change_clustering_config(transform = 'log2',\n")
cat("                        abundance_clust_rows = 'ward.D2')\n\n")

cat("# Ver configuración actual:\n")
cat("show_clustering_config()\n")

show_clustering_config()

cat("\nTop 10 taxones más abundantes (", config$taxonomic_level, "):\n")
top_abundant <- head(names(ranking_result$ranking), 10)
for (i in 1:length(top_abundant)) {
  abundance_val <- round(ranking_result$ranking[i], 2)
  distribution <- taxon_distribution_type[top_abundant[i]]
  cat(sprintf("  %2d. %s (%.2f, %s)\n", i, top_abundant[i], abundance_val, distribution))
}

cat("\n📈 Mostrando gráfico en RStudio...\n")
print(p_2b)

cat("\n🚀 ¡Análisis completo con orden fijo de muestras implementado!\n")
cat("🔬 Revisa los archivos generados y el gráfico en el panel de Plots.\n")
cat("📊 Todas las figuras mantienen el orden: Las Docas → Algarrobo → Navidad → Topocalma → Pargua → Los Chonos → Ilque → San Antonio\n")
cat("🧬 Los heatmaps incluyen clustering mejorado de taxones pero NO de muestras (orden fijo).\n")
cat("📋 Los archivos CSV incluyen '_FIXED_ORDER' en el nombre para identificarlos.\n")
