# Fig_3B.R
# Description: Generates combined ComplexHeatmap figures of quorum sensing (QS) gene distributions
#              across MAGs and genera, using Bray-Curtis clustering and functional category annotations.
# Purpose: Visualizes QS gene abundance by functional category and individual gene (name2) for both
#          MAG-level and genus-level aggregations. MAG row labels include taxonomic annotations
#          (order and genus). Outputs PDF figures and supplementary CSV tables for mBio submission.

rm(list = ls())

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(viridis)
  library(grid)
  library(vegan)
  library(pheatmap)
  library(stringr)
})

USE_CLUSTERING    <- TRUE
CLUSTERING_METHOD <- "average"
DISTANCE_METHOD   <- "bray"

COLOR_SCHEME   <- "Reds"
SHOW_NUMBERS   <- FALSE
MIN_PREVALENCE <- 2

FIG_WIDTH  <- 180
FIG_HEIGHT <- 200

group_colors <- c("wide" = "#1976D2", "rare" = "#FF8A65")

extract_genus_from_tax <- function(tax_string) {
  if (grepl("g__", tax_string)) {
    genus <- stringr::str_extract(tax_string, "g__[^;]+")
    genus <- gsub("g__", "", genus)
    if (!is.na(genus) && nzchar(genus) && !grepl("^\\s*$", genus)) return(genus)
  }

  if (grepl("o__", tax_string)) {
    order <- stringr::str_extract(tax_string, "o__[^;]+")
    order <- gsub("o__", "", order)
    if (!is.na(order) && nzchar(order) && !grepl("^\\s*$", order)) {
      return(paste0("Unclassified_", order))
    }
  }
  return("Unclassified")
}

create_mag_names_with_taxonomy <- function(mag_ids, nodes_data) {
  tax_info <- nodes_data %>%
    filter(MAG_ID %in% mag_ids) %>%
    mutate(

      order_clean = case_when(
        is.na(Order) | Order == "" | Order == "Unknown" ~ "Unknown",
        TRUE ~ Order
      ),

      genus_clean = case_when(
        is.na(Genus) | Genus == "" | Genus == "Unknown" ~ "Unknown",
        TRUE ~ Genus
      ),

      tax_label = paste0("(o__", order_clean, "; g__", genus_clean, ")"),

      mag_name_with_tax = paste0(MAG_ID, " ", tax_label)
    ) %>%
    select(MAG_ID, mag_name_with_tax)

  name_map <- setNames(tax_info$mag_name_with_tax, tax_info$MAG_ID)

  name_map[mag_ids]
}

drop_empty_rc <- function(mat, tag = "matrix") {
  mat[is.na(mat)] <- 0
  rs <- rowSums(mat)
  cs <- colSums(mat)
  keep_r <- rs > 0
  keep_c <- cs > 0
  removed_r <- sum(!keep_r); removed_c <- sum(!keep_c)
  if (removed_r > 0) cat(sprintf("• [%s] Filas todo-cero eliminadas: %d\n", tag, removed_r))
  if (removed_c > 0) cat(sprintf("• [%s] Columnas todo-cero eliminadas: %d\n", tag, removed_c))
  mat[keep_r, keep_c, drop = FALSE]
}

safe_dist <- function(mat, method = "bray", transpose = FALSE) {
  if (transpose) mat <- t(mat)
  if (nrow(mat) < 2) return(NULL)
  if (method == "bray") {
    d <- suppressWarnings(vegdist(mat, method = "bray"))
  } else {
    d <- dist(mat, method = method)
  }
  if (any(!is.finite(d))) return(NULL)
  d
}

safe_hclust <- function(dist_obj, method = "average") {
  if (is.null(dist_obj)) return(FALSE)
  hclust(dist_obj, method = method)
}

make_col_fun <- function(mat, scheme = "viridis") {
  rng <- range(mat, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2])) rng <- c(0, 1)
  if (rng[1] == rng[2]) rng <- c(rng[1], rng[1] + 1)

  cols <- switch(
    scheme,
    "viridis" = viridis(9),
    "RdBu"    = rev(brewer.pal(11, "RdBu")),
    "YlOrRd"  = brewer.pal(9, "YlOrRd"),
    "Reds"    = brewer.pal(9, "Reds"),
    "custom"  = c("white", "#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26", "#A50F15"),
    viridis(9)
  )
  colorRamp2(seq(rng[1], rng[2], length.out = length(cols)), cols)
}

cat("=== CARGANDO DATOS ===\n")

qs_matrix  <- read_tsv("./Matriz_qs_short.tsv", show_col_types = FALSE)
nodes      <- read_tsv("./global_node_metrics.tsv", show_col_types = FALSE)
categories <- read_tsv("categories.txt", show_col_types = FALSE)

cat("\n=== VERIFICACIÓN DE INTEGRIDAD ===\n")
genes_in_matrix     <- names(qs_matrix)[-1]
genes_in_categories <- categories$accession
missing_categories  <- setdiff(genes_in_matrix, genes_in_categories)

if (length(missing_categories) > 0) {
  cat("⚠️ Genes en matriz sin categoría (no aparecerán en heatmaps):\n")
  cat(paste("  -", head(missing_categories, 5)), sep = "\n")
  rem <- length(missing_categories) - min(5, length(missing_categories))
  if (rem > 0) cat(sprintf("  ... y %d más\n", rem))
}

missing_in_matrix <- setdiff(genes_in_categories, colnames(qs_matrix))
if (length(missing_in_matrix) > 0) {
  cat("\n⚠️ Genes en categorías pero no en matriz (ignorados):\n")
  cat(paste("  -", head(missing_in_matrix, 3)), sep = "\n")
}

cat(sprintf("\n📊 Resumen:\n"))
cat(sprintf("  - Genes totales en matriz: %d\n", length(genes_in_matrix)))
cat(sprintf("  - Genes con categoría: %d (%.1f%%)\n",
            sum(genes_in_matrix %in% genes_in_categories),
            100 * sum(genes_in_matrix %in% genes_in_categories) / length(genes_in_matrix)))
cat(sprintf("  - Categorías únicas: %d\n", dplyr::n_distinct(categories$Category)))

mags_in_network <- nodes$MAG_ID
qs_filtered <- qs_matrix %>%
  select(accession, any_of(mags_in_network))

qs_mat <- qs_filtered %>%
  column_to_rownames("accession") %>%
  as.matrix() %>%
  t()

qs_by_category <- qs_mat %>%
  as.data.frame() %>%
  rownames_to_column("MAG_ID") %>%
  pivot_longer(-MAG_ID, names_to = "accession", values_to = "count") %>%
  left_join(categories %>% select(accession, Category), by = "accession") %>%
  filter(!is.na(Category)) %>%
  group_by(MAG_ID, Category) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = Category, values_from = total, values_fill = 0) %>%
  column_to_rownames("MAG_ID") %>%
  as.matrix()

qs_by_name <- qs_mat %>%
  as.data.frame() %>%
  rownames_to_column("MAG_ID") %>%
  pivot_longer(-MAG_ID, names_to = "accession", values_to = "count") %>%
  left_join(categories %>% select(accession, name2), by = "accession") %>%
  filter(!is.na(name2)) %>%
  group_by(MAG_ID, name2) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  pivot_wider(names_from = name2, values_from = total, values_fill = 0) %>%
  column_to_rownames("MAG_ID") %>%
  as.matrix()

cat("\n=== ANÁLISIS POR GÉNERO ===\n")

tax_info <- nodes %>%
  mutate(
    Genus_clean = ifelse(is.na(Genus) | Genus == "" | Genus == "Unknown",
                         paste0("Unclassified_", Order),
                         Genus)
  ) %>%
  filter(MAG_ID %in% rownames(qs_by_category))

qs_genus_category <- qs_by_category %>%
  as.data.frame() %>%
  rownames_to_column("MAG_ID") %>%
  left_join(tax_info %>% select(MAG_ID, Genus_clean, Group), by = "MAG_ID") %>%
  group_by(Genus_clean) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  column_to_rownames("Genus_clean") %>%
  as.matrix()

qs_genus_genes <- qs_by_name %>%
  as.data.frame() %>%
  rownames_to_column("MAG_ID") %>%
  left_join(tax_info %>% select(MAG_ID, Genus_clean, Group), by = "MAG_ID") %>%
  group_by(Genus_clean) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop") %>%
  column_to_rownames("Genus_clean") %>%
  as.matrix()

genus_distribution <- tax_info %>%
  group_by(Genus_clean) %>%
  summarise(
    n_wide = sum(Group == "wide", na.rm = TRUE),
    n_rare = sum(Group == "rare", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Group = ifelse(n_wide >= n_rare, "wide", "rare"))

genus_group_info <- genus_distribution$Group
names(genus_group_info) <- genus_distribution$Genus_clean

cat("\n=== LIMPIEZA DE MATRICES ===\n")
qs_by_category <- drop_empty_rc(qs_by_category, tag = "Categorías")
qs_by_name     <- drop_empty_rc(qs_by_name,     tag = "Genes (name2)")
qs_genus_category <- drop_empty_rc(qs_genus_category, tag = "Géneros-Categorías")
qs_genus_genes <- drop_empty_rc(qs_genus_genes, tag = "Géneros-Genes")

gene_prev  <- colSums(qs_by_name > 0)
genes_keep <- names(gene_prev[gene_prev >= MIN_PREVALENCE])
if (length(genes_keep) == 0) {
  stop("No quedan genes tras aplicar MIN_PREVALENCE. Reduce MIN_PREVALENCE o revisa los datos.")
}
qs_name_filt <- qs_by_name[, genes_keep, drop = FALSE]
qs_name_filt <- drop_empty_rc(qs_name_filt, tag = "Genes filtrados")

common_genes <- intersect(colnames(qs_genus_genes), genes_keep)
qs_genus_genes_filt <- qs_genus_genes[, common_genes, drop = FALSE]
qs_genus_genes_filt <- drop_empty_rc(qs_genus_genes_filt, tag = "Géneros-Genes filtrados")

common_mags <- intersect(rownames(qs_by_category), rownames(qs_name_filt))
qs_by_category <- qs_by_category[common_mags, , drop = FALSE]
qs_name_filt   <- qs_name_filt[common_mags, ,   drop = FALSE]

cat("\n=== APLICANDO NOMBRES CON TAXONOMÍA ===\n")

mag_names_with_tax <- create_mag_names_with_taxonomy(rownames(qs_by_category), nodes)

rownames(qs_by_category) <- mag_names_with_tax[rownames(qs_by_category)]
rownames(qs_name_filt) <- mag_names_with_tax[rownames(qs_name_filt)]

cat(sprintf("✓ Aplicados nombres taxonómicos a %d MAGs\n", length(mag_names_with_tax)))

cat("Ejemplos de nombres generados:\n")
sample_names <- head(mag_names_with_tax, 3)
for(i in seq_along(sample_names)) {
  cat(sprintf("  %s\n", sample_names[i]))
}

row_clust <- FALSE
col_clust_cat  <- FALSE
col_clust_name <- FALSE
row_clust_genus <- FALSE
col_clust_genus_cat <- FALSE
col_clust_genus_genes <- FALSE

if (USE_CLUSTERING) {
  cat("\n=== APLICANDO CLUSTERING ===\n")
  cat(sprintf("Método: %s con distancia %s\n", CLUSTERING_METHOD, DISTANCE_METHOD))

  d_rows <- safe_dist(qs_by_category, method = DISTANCE_METHOD, transpose = FALSE)
  row_clust <- safe_hclust(d_rows, method = CLUSTERING_METHOD)

  d_cols_cat <- safe_dist(qs_by_category, method = DISTANCE_METHOD, transpose = TRUE)
  col_clust_cat <- safe_hclust(d_cols_cat, method = CLUSTERING_METHOD)

  d_cols_name <- safe_dist(qs_name_filt, method = DISTANCE_METHOD, transpose = TRUE)
  col_clust_name <- safe_hclust(d_cols_name, method = CLUSTERING_METHOD)

  d_rows_genus <- safe_dist(qs_genus_category, method = DISTANCE_METHOD, transpose = FALSE)
  row_clust_genus <- safe_hclust(d_rows_genus, method = CLUSTERING_METHOD)

  d_cols_genus_cat <- safe_dist(qs_genus_category, method = DISTANCE_METHOD, transpose = TRUE)
  col_clust_genus_cat <- safe_hclust(d_cols_genus_cat, method = CLUSTERING_METHOD)

  d_cols_genus_genes <- safe_dist(qs_genus_genes_filt, method = DISTANCE_METHOD, transpose = TRUE)
  col_clust_genus_genes <- safe_hclust(d_cols_genus_genes, method = CLUSTERING_METHOD)
} else {
  cat("\nClustering desactivado (orden manual).\n")
}

original_mag_ids <- gsub(" \\(o__.*", "", rownames(qs_by_category))
row_info <- nodes %>%
  filter(MAG_ID %in% original_mag_ids) %>%
  arrange(match(MAG_ID, original_mag_ids))

col_fun_cat  <- make_col_fun(qs_by_category, scheme = COLOR_SCHEME)
col_fun_name <- make_col_fun(qs_name_filt,   scheme = COLOR_SCHEME)
col_fun_genus_cat <- make_col_fun(qs_genus_category, scheme = COLOR_SCHEME)
col_fun_genus_genes <- make_col_fun(qs_genus_genes_filt, scheme = COLOR_SCHEME)

row_ha <- rowAnnotation(
  Group     = row_info$Group,
  `QS genes`= row_info$n_qs_genes,
  Degree    = row_info$degree,
  col = list(
    Group      = group_colors,
    `QS genes` = colorRamp2(c(min(row_info$n_qs_genes, na.rm=TRUE),
                              median(row_info$n_qs_genes, na.rm=TRUE),
                              max(row_info$n_qs_genes, na.rm=TRUE)),
                            c("#F7F7F7", "#969696", "#252525")),
    Degree     = colorRamp2(c(min(row_info$degree, na.rm=TRUE),
                              median(row_info$degree, na.rm=TRUE),
                              max(row_info$degree, na.rm=TRUE)),
                            c("#F7FBFF", "#6BAED6", "#08519C"))
  ),
  annotation_width = unit(c(3, 3, 3), "mm"),
  annotation_name_side = "bottom",
  annotation_name_gp = gpar(fontsize = 7),
  simple_anno_size_adjust = TRUE,
  gap = unit(0.5, "mm")
)

genus_rows <- rownames(qs_genus_category)
genus_dist <- genus_group_info[genus_rows]
genus_dist[is.na(genus_dist)] <- "unknown"

row_ha_genus <- rowAnnotation(
  Distribution = genus_dist,
  col = list(Distribution = group_colors),
  annotation_width = unit(3, "mm"),
  annotation_name_side = "bottom",
  annotation_name_gp = gpar(fontsize = 7)
)

col_ha_cat <- HeatmapAnnotation(
  `Total genes` = anno_barplot(
    colSums(qs_by_category),
    height = unit(10, "mm"),
    gp = gpar(fill = "#525252", col = NA),
    ylim = c(0, max(colSums(qs_by_category)) * 1.1)
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 7)
)

col_ha_genus_cat <- HeatmapAnnotation(
  `Total by genus` = anno_barplot(
    colSums(qs_genus_category),
    height = unit(10, "mm"),
    gp = gpar(fill = "#525252", col = NA),
    ylim = c(0, max(colSums(qs_genus_category)) * 1.1)
  ),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 7)
)

cell_fun_empty <- function(j, i, x, y, width, height, fill) {

}

ht_cat <- Heatmap(
  qs_by_category,
  name = "Count_Cat",
  col = col_fun_cat,
  cluster_rows    = row_clust,
  cluster_columns = col_clust_cat,
  show_row_dend = USE_CLUSTERING,
  show_column_dend = USE_CLUSTERING,
  row_dend_width = unit(10, "mm"),
  column_dend_height = unit(10, "mm"),
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_names_rot = 45,
  left_annotation = row_ha,
  top_annotation  = col_ha_cat,
  row_title = sprintf("MAGs (n = %d)", nrow(qs_by_category)),
  column_title = "Functional Categories",
  row_title_gp = gpar(fontsize = 9),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  cell_fun = cell_fun_empty,
  border = TRUE,
  rect_gp = gpar(col = "gray90", lwd = 0.25),
  width = unit(40, "mm"),
  height = unit(120, "mm")
)

ht_name <- Heatmap(
  qs_name_filt,
  name = "Count_Genes",
  col = col_fun_name,
  cluster_rows    = row_clust,
  cluster_columns = col_clust_name,
  show_row_dend   = FALSE,
  show_column_dend = USE_CLUSTERING,
  column_dend_height = unit(10, "mm"),
  row_names_gp = gpar(fontsize = 4),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 45,
  row_title = NULL,
  column_title = sprintf("Individual QS Genes (n = %d)", ncol(qs_name_filt)),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  cell_fun = cell_fun_empty,
  border = TRUE,
  rect_gp = gpar(col = "gray90", lwd = 0.25),
  show_heatmap_legend = FALSE,
  width = unit(120, "mm"),
  height = unit(120, "mm")
)

ht_genus_cat <- Heatmap(
  qs_genus_category,
  name = "Count_Genus_Cat",
  col = col_fun_genus_cat,
  cluster_rows = row_clust_genus,
  cluster_columns = col_clust_genus_cat,
  show_row_dend = USE_CLUSTERING,
  show_column_dend = USE_CLUSTERING,
  row_dend_width = unit(10, "mm"),
  column_dend_height = unit(10, "mm"),
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_names_rot = 45,
  left_annotation = row_ha_genus,
  top_annotation = col_ha_genus_cat,
  row_title = sprintf("Genera (n = %d)", nrow(qs_genus_category)),
  column_title = "QS Categories by Genus",
  row_title_gp = gpar(fontsize = 9),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  cell_fun = cell_fun_empty,
  border = TRUE,
  rect_gp = gpar(col = "gray90", lwd = 0.25),
  width = unit(50, "mm"),
  height = unit(80, "mm")
)

ht_genus_genes <- Heatmap(
  qs_genus_genes_filt,
  name = "Count_Genus_Genes",
  col = col_fun_genus_genes,
  cluster_rows = row_clust_genus,
  cluster_columns = col_clust_genus_genes,
  show_row_dend = FALSE,
  show_column_dend = USE_CLUSTERING,
  column_dend_height = unit(10, "mm"),
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 6),
  column_names_rot = 45,
  row_title = NULL,
  column_title = sprintf("QS Genes by Genus (n = %d)", ncol(qs_genus_genes_filt)),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  cell_fun = cell_fun_empty,
  border = TRUE,
  rect_gp = gpar(col = "gray90", lwd = 0.25),
  show_heatmap_legend = FALSE,
  width = unit(140, "mm"),
  height = unit(80, "mm")
)

pdf("Figure_QS_Combined_MAGs_WithTaxonomy_mBio.pdf", width = (FIG_WIDTH+20)/25.4, height = FIG_HEIGHT/25.4)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(2, 18), "cm"))))

pushViewport(viewport(layout.pos.row = 1))
grid.text("Quorum Sensing Gene Distribution in Kelp-Associated MAGs",
          gp = gpar(fontsize = 14, fontface = "bold"))
grid.text("\n(MAG names include taxonomic annotation: o__Order; g__Genus)",
          gp = gpar(fontsize = 10, fontface = "italic"))
if (USE_CLUSTERING) {
  grid.text(sprintf("\n\n(Clustered by %s distance, %s linkage)",
                    DISTANCE_METHOD, CLUSTERING_METHOD),
            gp = gpar(fontsize = 9, fontface = "italic"))
}
popViewport()

pushViewport(viewport(layout.pos.row = 2))
draw(ht_cat + ht_name, ht_gap = unit(5, "mm"),
     padding = unit(c(5, 5, 5, 5), "mm"))

ComplexHeatmap::decorate_heatmap_body("Count_Cat", {
  grid.text("A", x = unit(-3, "mm"), y = unit(1, "npc") + unit(4, "mm"),
            gp = gpar(fontsize = 12, fontface = "bold"))
})
ComplexHeatmap::decorate_heatmap_body("Count_Genes", {
  grid.text("B", x = unit(-3, "mm"), y = unit(1, "npc") + unit(4, "mm"),
            gp = gpar(fontsize = 12, fontface = "bold"))
})

popViewport()
dev.off()

pdf("Figure_QS_Combined_Genera_mBio.pdf", width = FIG_WIDTH/25.4, height = (FIG_HEIGHT-30)/25.4)

grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(2, 16), "cm"))))

pushViewport(viewport(layout.pos.row = 1))
grid.text("Quorum Sensing Gene Distribution by Genus",
          gp = gpar(fontsize = 14, fontface = "bold"))
if (USE_CLUSTERING) {
  grid.text(sprintf("\n(Clustered by %s distance, %s linkage)",
                    DISTANCE_METHOD, CLUSTERING_METHOD),
            gp = gpar(fontsize = 10, fontface = "italic"))
}
popViewport()

pushViewport(viewport(layout.pos.row = 2))
draw(ht_genus_cat + ht_genus_genes, ht_gap = unit(5, "mm"),
     padding = unit(c(5, 5, 5, 5), "mm"))

ComplexHeatmap::decorate_heatmap_body("Count_Genus_Cat", {
  grid.text("C", x = unit(-3, "mm"), y = unit(1, "npc") + unit(4, "mm"),
            gp = gpar(fontsize = 12, fontface = "bold"))
})
ComplexHeatmap::decorate_heatmap_body("Count_Genus_Genes", {
  grid.text("D", x = unit(-3, "mm"), y = unit(1, "npc") + unit(4, "mm"),
            gp = gpar(fontsize = 12, fontface = "bold"))
})

popViewport()
dev.off()

pdf("Figure_QS_Categories_Individual_WithTaxonomy.pdf", width = 10, height = 12)
draw(ht_cat, padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()

pdf("Figure_QS_Genes_Individual_WithTaxonomy.pdf", width = 16, height = 12)
draw(ht_name, padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()

pdf("Figure_QS_Genera_Categories.pdf", width = 8, height = 8)
draw(ht_genus_cat, padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()

pdf("Figure_QS_Genera_Genes.pdf", width = 14, height = 8)
draw(ht_genus_genes, padding = unit(c(10, 10, 10, 10), "mm"))
dev.off()

cat("\n=== RESUMEN FINAL ===\n")
cat("✓ Figura combinada MAGs (con taxonomía): Figure_QS_Combined_MAGs_WithTaxonomy_mBio.pdf\n")
cat("✓ Figura combinada Géneros: Figure_QS_Combined_Genera_mBio.pdf\n")
cat("✓ Figuras individuales generadas (con taxonomía donde corresponde)\n")

cat("\n=== DISTRIBUCIÓN POR CATEGORÍAS (MAGs) ===\n")
cat_summary_mags <- qs_by_category %>%
  as.data.frame() %>%
  summarise(across(everything(), list(
    present = ~sum(. > 0),
    mean    = ~mean(.),
    max     = ~max(.)
  )))
print(t(cat_summary_mags))

cat("\n=== DISTRIBUCIÓN POR CATEGORÍAS (Géneros) ===\n")
cat_summary_genus <- qs_genus_category %>%
  as.data.frame() %>%
  summarise(across(everything(), list(
    present = ~sum(. > 0),
    mean    = ~mean(.),
    max     = ~max(.)
  )))
print(t(cat_summary_genus))

supp_table_mags <- qs_by_category %>%
  as.data.frame() %>%
  rownames_to_column("MAG_ID_with_taxonomy") %>%
  mutate(
    MAG_ID = gsub(" \\(o__.*", "", MAG_ID_with_taxonomy),
    Taxonomy = gsub(".*( \\(o__.*\\))", "\\1", MAG_ID_with_taxonomy)
  ) %>%
  left_join(nodes %>% select(MAG_ID, Group, Site, Order, Genus), by = "MAG_ID") %>%
  select(MAG_ID_with_taxonomy, MAG_ID, Taxonomy, Group, Site, Order, Genus, everything()) %>%
  select(-MAG_ID)
write.csv(supp_table_mags, "Table_S_QS_Categories_MAGs_WithTaxonomy.csv", row.names = FALSE)

supp_table_genus <- qs_genus_category %>%
  as.data.frame() %>%
  rownames_to_column("Genus") %>%
  left_join(genus_distribution %>% select(Genus_clean, Group) %>%
              rename(Genus = Genus_clean, Distribution = Group), by = "Genus") %>%
  select(Genus, Distribution, everything())
write.csv(supp_table_genus, "Table_S_QS_Categories_Genera.csv", row.names = FALSE)

cat("\n=== ESTADÍSTICAS COMPARATIVAS ===\n")
cat(sprintf("MAGs totales analizados: %d\n", nrow(qs_by_category)))
cat(sprintf("Géneros únicos: %d\n", nrow(qs_genus_category)))
cat(sprintf("Genes QS filtrados (prevalencia >= %d): %d\n", MIN_PREVALENCE, ncol(qs_name_filt)))
cat(sprintf("Categorías funcionales: %d\n", ncol(qs_by_category)))

wide_mags <- sum(row_info$Group == "wide", na.rm = TRUE)
rare_mags <- sum(row_info$Group == "rare", na.rm = TRUE)
wide_genera <- sum(genus_distribution$Group == "wide", na.rm = TRUE)
rare_genera <- sum(genus_distribution$Group == "rare", na.rm = TRUE)

cat(sprintf("\nDistribución MAGs - Wide: %d, Rare: %d\n", wide_mags, rare_mags))
cat(sprintf("Distribución Géneros - Wide: %d, Rare: %d\n", wide_genera, rare_genera))

gene_prevalence <- colSums(qs_name_filt > 0)
top_genes <- sort(gene_prevalence, decreasing = TRUE)
cat("\nTop 10 genes QS más prevalentes:\n")
print(head(top_genes, 10))

cat("\n=== ANÁLISIS DE DIVERSIDAD TAXONÓMICA ===\n")
tax_diversity <- nodes %>%
  filter(MAG_ID %in% original_mag_ids) %>%
  summarise(
    unique_orders = n_distinct(Order, na.rm = TRUE),
    unique_genera = n_distinct(Genus, na.rm = TRUE),
    mags_with_genus = sum(!is.na(Genus) & Genus != "" & Genus != "Unknown"),
    mags_with_order = sum(!is.na(Order) & Order != "" & Order != "Unknown")
  )

cat(sprintf("Órdenes únicos representados: %d\n", tax_diversity$unique_orders))
cat(sprintf("Géneros únicos representados: %d\n", tax_diversity$unique_genera))
cat(sprintf("MAGs con género identificado: %d (%.1f%%)\n",
            tax_diversity$mags_with_genus,
            100 * tax_diversity$mags_with_genus / length(original_mag_ids)))
cat(sprintf("MAGs con orden identificado: %d (%.1f%%)\n",
            tax_diversity$mags_with_order,
            100 * tax_diversity$mags_with_order / length(original_mag_ids)))

cat("\n✓ Tablas suplementarias generadas:\n")
cat("  - Table_S_QS_Categories_MAGs_WithTaxonomy.csv (ACTUALIZADA con nombres taxonómicos)\n")
cat("  - Table_S_QS_Categories_Genera.csv\n")

cat("\n=== PROCESO COMPLETADO ===\n")
cat("Mejoras implementadas:\n")
cat("• Nombres de MAGs incluyen taxonomía: MAG_ID (o__Order; g__Genus)\n")
cat("• Tamaño de fuente ajustado para acomodar nombres más largos\n")
cat("• Ancho de figuras ajustado para mejor visualización\n")
cat("• Eliminados números de celdas (SHOW_NUMBERS = FALSE)\n")
cat("• Colores consistentes con esquema AMR\n")
cat("• Análisis agregado por género\n")
cat("• Heatmaps adicionales por género\n")
cat("• Tablas suplementarias actualizadas con información taxonómica\n")
cat("• Análisis de diversidad taxonómica agregado\n")

cat("\n=== EJEMPLOS DE NOMBRES CON TAXONOMÍA ===\n")
cat("Primeros 5 nombres de MAGs con taxonomía:\n")
example_names <- head(rownames(qs_by_category), 5)
for(i in seq_along(example_names)) {
  cat(sprintf("  %d. %s\n", i, example_names[i]))
}
