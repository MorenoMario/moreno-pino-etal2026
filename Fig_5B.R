# Fig_5B.R
# Description: Creates alluvial (Sankey) diagrams showing the flow of CRISPR-Cas genes across
#              functional categories and taxonomic levels, colored by biogeographic group.
# Purpose: Visualizes CRISPR gene diversity by mapping each gene to its CRISPR functional category
#          (Adaptation, Processing, Class I/II Effectors, Accessory) and linking them to a selected
#          taxonomic level. Produces vector-format outputs (PDF, SVG, EPS) optimized for Adobe Illustrator.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(ggalluvial)
  library(readr)
  library(RColorBrewer)
})

separate_axis3_lines <- TRUE

create_separate_plots <- FALSE

export_svg <- TRUE
export_eps <- FALSE

simplify_for_illustrator <- TRUE

line_width_multiplier <- 1.5

crispr_categories <- tribble(
  ~Protein_Name, ~Category,
  "CRISPR-associated endoribonuclease Cas6 [EC:3.1.-.-]", "Processing",
  "CRISPR-associated endonuclease Csy4 [EC:3.1.-.-]", "Processing",
  "CRISPR-associated protein Cas1", "Adaptation",
  "CRISPR-associated protein Cas2", "Adaptation",
  "CRISPR-associated exonuclease Cas4 [EC:3.1.12.1]", "Adaptation",
  "CRISPR-associated endonuclease/helicase Cas3 [EC:3.1.-.- 5.6.2.4]", "Class_I_Effectors",
  "CRISPR-associated protein Cas5d", "Class_I_Effectors",
  "CRISPR-associated protein Csd1", "Class_I_Effectors",
  "CRISPR-associated protein Csc1", "Class_I_Effectors",
  "CRISPR-associated protein Csm3", "Class_I_Effectors",
  "CRISPR-associated protein Csm5", "Class_I_Effectors",
  "CRISPR-associated protein Csc2", "Class_I_Effectors",
  "CRISPR-associated protein Csd2", "Class_I_Effectors",
  "CRISPR-associated protein Csy2", "Class_I_Effectors",
  "CRISPR-associated protein Csy3", "Class_I_Effectors",
  "CRISPR-associated protein Cmr1", "Class_I_Effectors",
  "CRISPR-associated protein Cmr2", "Class_I_Effectors",
  "CRISPR-associated protein Cmr3", "Class_I_Effectors",
  "CRISPR-associated protein Cmr4", "Class_I_Effectors",
  "CRISPR-associated protein Cmr5", "Class_I_Effectors",
  "CRISPR-associated protein Cmr6", "Class_I_Effectors",
  "CRISPR-associated protein Csa3", "Class_I_Effectors",
  "CRISPR-associated protein Csm1", "Class_I_Effectors",
  "CRISPR-associated protein Csm2", "Class_I_Effectors",
  "CRISPR-associated protein Csm4", "Class_I_Effectors",
  "CRISPR-associated protein Csy1", "Class_I_Effectors",
  "CRISPR-associated endonuclease Csn1 [EC:3.1.-.-]", "Class_II_Effectors",
  "CRISPR-associated protein Csx17", "Class_II_Effectors",
  "CRISPR-associated protein Cst2", "Accessory",
  "CRISPR-associated protein Csx10", "Accessory",
  "CRISPR-associated protein Csx3", "Accessory",
  "CRISPR-associated protein Csb1", "Accessory",
  "CRISPR-associated protein Csb2", "Accessory",
  "CRISPR-associated protein Csx16", "Accessory",
  "CRISPR-associated protein Csx14", "Accessory",
  "CRISPR-associated protein Csh1", "Accessory"
)

tax_level <- "class"

gene_order_mode <- "by_tsv_order"
custom_gene_order <- c()
custom_category_order <- c("Adaptation", "Processing", "Class_I_Effectors")
custom_taxon_order <- c()
custom_group_order <- c()

args <- commandArgs(trailingOnly = TRUE)
input_file <- if (length(args) >= 1) args[1] else NULL
categories_file <- if (length(args) >= 2) args[2] else NULL

if (interactive() && is.null(input_file)) {
  input_file <- "crispr_table.tsv"
}
if (is.null(input_file)) stop("Debes especificar la ruta del archivo TSV como primer argumento.")

if (!is.null(categories_file) && file.exists(categories_file)) {
  cat("Cargando categorías desde archivo:", categories_file, "\n")
  crispr_categories <- read_delim(categories_file, delim = "\t", col_names = TRUE, show_col_types = FALSE)
  if (!all(c("Protein_Name", "Category") %in% names(crispr_categories))) {
    stop("El archivo de categorías debe contener las columnas 'Protein_Name' y 'Category'")
  }
}

normalize_colnames <- function(x) {
  x %>%
    gsub("\\[.*?\\]", "", ., perl = TRUE) %>%
    gsub("[^A-Za-z0-9]+", "_", ., perl = TRUE) %>%
    gsub("_+", "_", ., perl = TRUE) %>%
    gsub("^_|_$", "", ., perl = TRUE)
}

extract_taxonomy <- function(df) {
  df %>%
    select(name, taxon) %>%
    separate(taxon, into = c("domain","phylum","class","order","family","genus","species"),
             sep = ";", fill = "right", remove = FALSE) %>%
    mutate(
      phylum = sub("^p__", "", phylum),
      class  = sub("^c__", "", class),
      order  = sub("^o__", "", order),
      family = sub("^f__", "", family),
      genus  = sub("^g__", "", genus),
      species= sub("^s__", "", species),
      genus  = if_else(is.na(genus)  | genus  == "", "Unknown_genus", genus),
      family = if_else(is.na(family) | family == "", "Unknown_family", family),
      order  = if_else(is.na(order)  | order == "", "Unknown_order", order),
      class  = if_else(is.na(class)  | class == "", "Unknown_class", class),
      phylum = if_else(is.na(phylum) | phylum == "", "Unknown_phylum", phylum)
    ) %>%
    select(name, taxon, phylum, class, order, family, genus, species)
}

cat("Leyendo archivo: ", input_file, "\n", sep = "")
raw <- read_delim(input_file, delim = "\t", col_names = TRUE, guess_max = 1e5, show_col_types = FALSE)
orig_names <- names(raw)
original_gene_cols_order <- setdiff(orig_names, c("name","taxon","total","biogeographic_group"))

original_to_normalized <- setNames(normalize_colnames(original_gene_cols_order), original_gene_cols_order)
normalized_to_original <- setNames(original_gene_cols_order, normalize_colnames(original_gene_cols_order))

names(raw) <- normalize_colnames(names(raw))

required_cols <- c("name","taxon","total","biogeographic_group")
missing <- setdiff(required_cols, names(raw))
if (length(missing)) stop("Faltan columnas requeridas en el TSV: ", paste(missing, collapse = ", "), call. = FALSE)

gene_cols <- setdiff(names(raw), c("name","taxon","total","biogeographic_group"))
normalized_gene_order_from_tsv <- original_to_normalized[original_gene_cols_order]
normalized_gene_order_from_tsv <- normalized_gene_order_from_tsv[normalized_gene_order_from_tsv %in% gene_cols]

raw[gene_cols] <- lapply(gene_cols, \(cn) suppressWarnings(as.numeric(raw[[cn]])))
crispr <- raw

tax_info <- extract_taxonomy(crispr)
tax_info$biogeographic_group <- factor(crispr$biogeographic_group)

valid_levels <- c("genus","family","phylum","class","order","species")
if (!tax_level %in% valid_levels) stop("tax_level inválido: ", tax_level)

if (length(custom_group_order) > 0) {
  tax_info$biogeographic_group <- factor(tax_info$biogeographic_group, levels = custom_group_order)
}
present_groups <- levels(droplevels(tax_info$biogeographic_group))

crispr_categories_normalized <- crispr_categories %>%
  mutate(
    Protein_Name_Normalized = normalize_colnames(Protein_Name)
  )

gene_to_category <- setNames(crispr_categories_normalized$Category,
                             crispr_categories_normalized$Protein_Name_Normalized)

genes_with_categories <- intersect(names(gene_to_category), gene_cols)
genes_without_categories <- setdiff(gene_cols, names(gene_to_category))

if (length(genes_without_categories) > 0) {
  cat("⚠️ Genes sin categoría asignada:", paste(genes_without_categories, collapse = ", "), "\n")

  for (gene in genes_without_categories) {
    gene_to_category[gene] <- "Unknown"
  }
}

cat("✅ Genes con categorías:", length(genes_with_categories), "\n")
cat("📋 Categorías encontradas:", paste(unique(gene_to_category), collapse = ", "), "\n")

long_genes <- crispr %>%
  select(name, all_of(gene_cols)) %>%
  pivot_longer(cols = all_of(gene_cols), names_to = "Gene", values_to = "Value") %>%
  mutate(Value = suppressWarnings(as.numeric(Value))) %>%
  replace_na(list(Value = 0)) %>%
  left_join(tax_info %>% select(name, biogeographic_group, all_of(tax_level)), by = "name") %>%
  rename(Taxon = !!tax_level) %>%

  mutate(
    CRISPR_Category = gene_to_category[Gene],
    CRISPR_Category = if_else(is.na(CRISPR_Category), "Unknown", CRISPR_Category)
  )

if (gene_order_mode == "by_tsv_order") {
  gene_order <- intersect(normalized_gene_order_from_tsv, unique(long_genes$Gene))
  remaining_genes <- setdiff(unique(long_genes$Gene), gene_order)
  if (length(remaining_genes) > 0) {
    gene_order <- c(gene_order, sort(remaining_genes))
  }
} else {
  gene_order <- sort(unique(long_genes$Gene))
}

long_genes$Gene <- factor(long_genes$Gene, levels = gene_order)

category_order <- if(length(custom_category_order) > 0) {
  c(custom_category_order, setdiff(unique(long_genes$CRISPR_Category), custom_category_order))
} else {
  sort(unique(long_genes$CRISPR_Category))
}
long_genes$CRISPR_Category <- factor(long_genes$CRISPR_Category, levels = category_order)

taxon_order <- if(length(custom_taxon_order) > 0) {
  c(custom_taxon_order, setdiff(unique(long_genes$Taxon), custom_taxon_order))
} else {
  sort(unique(long_genes$Taxon))
}
long_genes$Taxon <- factor(long_genes$Taxon, levels = taxon_order)

alluvial_df <- long_genes %>%
  group_by(Gene, CRISPR_Category, Taxon, biogeographic_group) %>%
  summarise(Total = sum(Value, na.rm = TRUE), .groups = "drop") %>%
  filter(Total > 0)

default_group_colors <- c(Widespread="#1F77B4", Restricted="#FF7F0E", Intermediate="#6BAED6", Unknown="#808080")
pg <- levels(factor(alluvial_df$biogeographic_group))
group_cols2 <- default_group_colors[pg]

category_colors <- c(
  "Adaptation" = "#E31A1C",
  "Processing" = "#1F78B4",
  "Class_I_Effectors" = "#33A02C",
  "Unknown" = "#999999"
)

create_illustrator_alluvial <- function(data, separate_lines = TRUE, simplify = TRUE) {

  base_plot <- ggplot(data, aes(axis1 = Gene, axis2 = CRISPR_Category, axis3 = Taxon, y = Total))

  if (separate_lines) {

    base_plot <- base_plot +
      geom_alluvium(aes(fill = biogeographic_group, group = interaction(Gene, CRISPR_Category, Taxon)),
                    alpha = 0.8,
                    width = 1/4,
                    knot.pos = 0.2) +
      geom_stratum(alpha = 0.6, color = "white", width = 0.25, size = 0.3)
  } else {

    base_plot <- base_plot +
      geom_alluvium(aes(fill = biogeographic_group), alpha = 0.7) +
      geom_stratum(alpha = 0.4, color = "white", width = 0.3)
  }

  if (simplify) {

    base_plot <- base_plot +
      geom_text(stat = "stratum", aes(label = after_stat(stratum)),
                size = 3, fontface = "bold", color = "black")
  } else {
    base_plot <- base_plot +
      geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2)
  }

  base_plot +
    scale_fill_manual(values = group_cols2, drop = FALSE) +
    scale_x_discrete(limits = c("Gene", "CRISPR_Category", "Taxon"),
                     labels = c("Gene", "CRISPR Category", tax_level),
                     expand = c(0.1, 0.1)) +
    labs(title = paste0("CRISPR gene flow: Gene → Category → ", tax_level),
         subtitle = "Optimizado para Adobe Illustrator",
         y = "Gene count", x = NULL, fill = "Biogeography") +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray50"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold", size = 12),
      panel.background = element_rect(fill = "white", color = "white")
    )
}

p_illustrator <- create_illustrator_alluvial(alluvial_df,
                                             separate_lines = separate_axis3_lines,
                                             simplify = simplify_for_illustrator)

if (create_separate_plots) {
  taxon_levels <- unique(alluvial_df$Taxon)

  for (level in taxon_levels) {
    subset_data <- alluvial_df %>% filter(Taxon == level)

    if (nrow(subset_data) > 0) {
      p_subset <- create_illustrator_alluvial(subset_data,
                                              separate_lines = separate_axis3_lines,
                                              simplify = simplify_for_illustrator) +
        labs(title = paste0("CRISPR genes → Category → ", level),
             subtitle = paste0("Subset: ", tax_level, " = ", level))

      filename_base <- paste0("crispr_subset_", gsub("[^A-Za-z0-9]", "_", level))

      if (export_svg) {
        ggsave(paste0(filename_base, ".svg"), p_subset,
               width = 16, height = 10, dpi = 300, device = "svg")
      }
      if (export_eps) {
        ggsave(paste0(filename_base, ".eps"), p_subset,
               width = 16, height = 10, dpi = 300, device = "eps")
      }
      ggsave(paste0(filename_base, ".pdf"), p_subset,
             width = 16, height = 10, dpi = 300)
    }
  }
}

message("Guardando archivos optimizados para Ilustrator...")

if (export_svg) {
  tryCatch({
    if (requireNamespace("svglite", quietly = TRUE)) {
      ggsave("crispr_categories_illustrator.svg", p_illustrator,
             width = 20, height = 12, dpi = 300, device = "svg")
      message("✓ Archivo SVG creado: crispr_categories_illustrator.svg")
    } else {
      message("⚠ svglite no disponible. Instalando...")
      install.packages("svglite")
      ggsave("crispr_categories_illustrator.svg", p_illustrator,
             width = 20, height = 12, dpi = 300, device = "svg")
      message("✓ Archivo SVG creado: crispr_categories_illustrator.svg")
    }
  }, error = function(e) {
    message("✗ Error exportando SVG: ", e$message)
    message("→ Usando PDF como alternativa vectorial")
    export_svg <<- FALSE
  })
}

if (export_eps) {
  tryCatch({
    ggsave("crispr_categories_illustrator2.eps", p_illustrator,
           width = 20, height = 12, dpi = 300, device = "eps")
    message("✓ Archivo EPS creado: crispr_categories_illustrator.eps")
  }, error = function(e) {
    message("✗ Error exportando EPS: ", e$message)
    message("→ EPS puede requerir ghostscript instalado en el sistema")
  })
}

ggsave("crispr_categories_illustrator2.pdf", p_illustrator,
       width = 20, height = 12, dpi = 300)

write_csv(alluvial_df, "crispr_categories_data.csv")

write_csv(crispr_categories_normalized, "crispr_categories_mapping.csv")

cat("\n=== RESUMEN DE DATOS ===\n")
cat("📊 Total de observaciones:", nrow(alluvial_df), "\n")
cat("🧬 Genes únicos:", length(unique(alluvial_df$Gene)), "\n")
cat("🏷️  Categorías CRISPR:", length(unique(alluvial_df$CRISPR_Category)), "\n")
cat("   -", paste(unique(alluvial_df$CRISPR_Category), collapse = ", "), "\n")
cat("🦠 Taxa únicos (", tax_level, "):", length(unique(alluvial_df$Taxon)), "\n")
cat("🌍 Grupos biogeográficos:", length(unique(alluvial_df$biogeographic_group)), "\n")
cat("   -", paste(unique(alluvial_df$biogeographic_group), collapse = ", "), "\n")

print(p_illustrator)

cat("\n=== ARCHIVOS CREADOS PARA ILUSTRATOR ===\n")
cat("📄 Gráfico principal (VECTORIAL): crispr_categories_illustrator.pdf\n")
if (export_svg && requireNamespace("svglite", quietly = TRUE)) {
  cat("🎨 Gráfico SVG: crispr_categories_illustrator.svg\n")
}
if (export_eps) {
  cat("📐 Gráfico EPS: crispr_categories_illustrator.eps\n")
}
cat("📊 Datos CSV: crispr_categories_data.csv\n")
cat("🏷️  Mapeo de categorías: crispr_categories_mapping.csv\n")

if (create_separate_plots) {
  cat("📁 Gráficos separados por", tax_level, "creados\n")
}

cat("\n=== USO DEL SCRIPT ===\n")
cat("📝 USO BÁSICO:\n")
cat("   Rscript script.R crispr_table.tsv\n\n")
cat("📝 CON ARCHIVO DE CATEGORÍAS EXTERNO:\n")
cat("   Rscript script.R crispr_table.tsv categorias.tsv\n\n")
cat("🔧 CONFIGURACIÓN:\n")
cat("   • tax_level:", tax_level, "(línea 55)\n")
cat("   • custom_category_order: personalizable (línea 58)\n")
cat("   • separate_axis3_lines:", separate_axis3_lines, "(línea 16)\n")

cat("\n=== RECOMENDACIONES PARA ILLUSTRATOR ===\n")
cat("✅ ARCHIVO PRINCIPAL: crispr_categories_illustrator.pdf\n")
cat("   → Flujo: Gene → Categoría CRISPR →", tax_level, "\n")
cat("   → Coloreado por grupo biogeográfico\n")
cat("   → Completamente vectorial y editable\n")
