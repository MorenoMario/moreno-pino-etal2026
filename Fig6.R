# Fig6.R
# Description: Builds and visualizes a global co-occurrence network of quorum sensing (QS) genes across MAGs
#              using the CoNet methodology (Pearson + Spearman correlation intersection).
# Purpose: Constructs a QS gene co-occurrence network, computes network metrics (degree, betweenness,
#          modularity), and generates multiple visualization layouts (force-directed, circular, sectoral).
#          Also produces statistical comparisons between widespread and restricted MAG groups,
#          normalized QS gene counts per genome by site, and supplementary TSV tables.

suppressPackageStartupMessages({
  library(tidyr); library(dplyr); library(stringr)
  library(ggplot2); library(ggrepel); library(patchwork)
  library(igraph)
  library(vegan)
  library(reshape2)
})

output_file <- "QS_CoNet_global_network_analysis2.pdf"

min_prevalence <- 0
correlation_threshold <- 0.7
p_value_threshold <- 0.01
max_double_zeros <- 0.95

normalize_BC_by_genome_count <- TRUE

node_filter <- "all"
node_filter_method <- "n_qs_genes"

label_mode <- "all"
label_topdegree_quantile <- 0.0
label_text_field <- "tax"
label_size <- 2.5
label_segment_color <- "grey50"

layout_scale_factor <- 0.8
layout_target_radius <- 0.5
layout_iterations <- 100

node_size_range <- c(6, 10)
node_size_transform <- "linear"

edge_linewidth_range <- c(0.1, 1)
edge_alpha <- 0.9
edge_color <- "grey30"

taxa_colors <- c(
  "wide" = "#1F77B4",
  "rare" = "#F28A1C",
  "unknown" = "#808080"
)

site_colors <- c(
  "Algarrobo" = "#E41A1C",
  "Los_Chonos" = "#377EB8",
  "Las_Docas" = "#4DAF4A",
  "Ilque" = "#984EA3",
  "Pargua" = "#FF7F00",
  "San_Antonio" = "#FFFF33",
  "Topocalma" = "#A65628",
  "Navidad" = "#F781BF"
)

calculate_correlations <- function(data_matrix, method = "pearson") {
  n_items <- ncol(data_matrix)
  cor_matrix <- matrix(0, n_items, n_items)
  p_matrix <- matrix(1, n_items, n_items)
  colnames(cor_matrix) <- colnames(data_matrix)
  rownames(cor_matrix) <- colnames(data_matrix)
  colnames(p_matrix) <- colnames(data_matrix)
  rownames(p_matrix) <- colnames(data_matrix)

  for (i in 1:(n_items-1)) {
    for (j in (i+1):n_items) {
      both_zero <- sum(data_matrix[,i] == 0 & data_matrix[,j] == 0)
      double_zero_prop <- both_zero / nrow(data_matrix)
      if (double_zero_prop <= max_double_zeros) {
        if (sd(data_matrix[,i]) > 0 && sd(data_matrix[,j]) > 0) {
          tryCatch({
            test_result <- cor.test(data_matrix[,i], data_matrix[,j],
                                    method = method, exact = FALSE)
            cor_matrix[i,j] <- cor_matrix[j,i] <- as.numeric(test_result$estimate)
            p_matrix[i,j] <- p_matrix[j,i] <- test_result$p.value
          }, error = function(e) {})
        }
      }
    }
  }
  return(list(correlation = cor_matrix, p_value = p_matrix))
}

apply_fdr_correction <- function(p_matrix) {
  upper_tri <- upper.tri(p_matrix)
  p_values <- p_matrix[upper_tri]
  adjusted_p <- p.adjust(p_values, method = "BH")
  p_adjusted_matrix <- matrix(1, nrow(p_matrix), ncol(p_matrix))
  colnames(p_adjusted_matrix) <- colnames(p_matrix)
  rownames(p_adjusted_matrix) <- rownames(p_matrix)
  p_adjusted_matrix[upper_tri] <- adjusted_p
  p_adjusted_matrix[lower.tri(p_adjusted_matrix)] <- t(p_adjusted_matrix)[lower.tri(p_adjusted_matrix)]
  diag(p_adjusted_matrix) <- 0
  return(p_adjusted_matrix)
}

intersect_networks <- function(pearson_edges, spearman_edges) {
  if (nrow(pearson_edges) == 0 || nrow(spearman_edges) == 0) return(data.frame())
  pearson_edges$edge_id <- paste(pmin(pearson_edges$from, pearson_edges$to),
                                 pmax(pearson_edges$from, pearson_edges$to), sep = "_")
  spearman_edges$edge_id <- paste(pmin(spearman_edges$from, spearman_edges$to),
                                  pmax(spearman_edges$from, spearman_edges$to), sep = "_")
  common_edges <- intersect(pearson_edges$edge_id, spearman_edges$edge_id)
  if (length(common_edges) == 0) return(data.frame())
  intersected <- pearson_edges[pearson_edges$edge_id %in% common_edges, ]
  spearman_matches <- match(intersected$edge_id, spearman_edges$edge_id)
  intersected$spearman_cor <- spearman_edges$correlation[spearman_matches]
  intersected$avg_correlation <- (intersected$correlation + intersected$spearman_cor) / 2
  intersected$edge_id <- NULL
  return(intersected)
}

extract_site <- function(mag_id) {
  dplyr::case_when(
    str_detect(mag_id, "^al\\d+") ~ "Algarrobo",
    str_detect(mag_id, "^chon\\d+") ~ "Los_Chonos",
    str_detect(mag_id, "^doc\\d+") ~ "Las_Docas",
    str_detect(mag_id, "^ilq\\d+") ~ "Ilque",
    str_detect(mag_id, "^par\\d+") ~ "Pargua",
    str_detect(mag_id, "^sant\\d+") ~ "San_Antonio",
    str_detect(mag_id, "^top\\d+") ~ "Topocalma",
    str_detect(mag_id, "^nav\\d+") ~ "Navidad",
    TRUE ~ "Unknown"
  )
}

extract_taxonomy <- function(tax_string) {
  tax_levels <- unlist(strsplit(tax_string, ";"))
  order_match <- grep("^o__", tax_levels, value = TRUE)
  order <- if(length(order_match) > 0) gsub("^o__", "", order_match[1]) else "Unknown_order"
  genus_match <- grep("^g__", tax_levels, value = TRUE)
  genus <- if(length(genus_match) > 0) {
    genus_clean <- gsub("^g__", "", genus_match[1])
    if(genus_clean == "" || genus_clean == " ") "Unknown_genus" else genus_clean
  } else "Unknown_genus"
  if(order %in% c("", " ")) {
    class_match <- grep("^c__", tax_levels, value = TRUE)
    if(length(class_match) > 0) {
      order <- gsub("^c__", "", class_match[1])
      if(order %in% c("", " ")) order <- "Unknown_order"
    }
  }
  if(genus == "Unknown_genus") {
    family_match <- grep("^f__", tax_levels, value = TRUE)
    if(length(family_match) > 0) {
      genus <- gsub("^f__", "", family_match[1])
      if(genus %in% c("", " ")) genus <- "Unknown_genus"
    }
  }
  return(c(order = order, genus = genus))
}

create_tax_label <- function(mag_id, order, genus) {
  if(order == "Unknown_order" && genus == "Unknown_genus") {
    return(mag_id)
  } else if(order == "Unknown_order") {
    return(paste0(mag_id, " (", genus, ")"))
  } else if(genus == "Unknown_genus") {
    return(paste0(mag_id, " (", order, ")"))
  } else {
    return(paste0(mag_id, " (", order, "; ", genus, ")"))
  }
}

filter_nodes <- function(graph, node_metrics, filter_type, filter_method, n_nodes) {
  if(filter_type == "all") return(list(graph = graph, node_metrics = node_metrics))
  if(filter_method == "degree") {
    top_nodes <- node_metrics %>% arrange(desc(degree)) %>% head(n_nodes) %>% pull(MAG_ID)
  } else if(filter_method == "betweenness") {
    top_nodes <- node_metrics %>% arrange(desc(betweenness)) %>% head(n_nodes) %>% pull(MAG_ID)
  } else if(filter_method == "n_qs_genes") {
    top_nodes <- node_metrics %>% arrange(desc(n_qs_genes)) %>% head(n_nodes) %>% pull(MAG_ID)
  } else stop("filter_method debe ser 'degree', 'betweenness' o 'n_qs_genes'")
  filtered_graph <- induced_subgraph(graph, top_nodes)
  filtered_metrics <- node_metrics[node_metrics$MAG_ID %in% top_nodes, ]
  if(vcount(filtered_graph) > 0) {
    filtered_metrics$degree <- degree(filtered_graph)
    filtered_metrics$betweenness <- betweenness(filtered_graph, normalized = TRUE)
    filtered_metrics$closeness <- closeness(filtered_graph, normalized = TRUE)
    filtered_metrics$eigenvector <- eigen_centrality(filtered_graph)$vector
    threshold_hub <- mean(filtered_metrics$degree) + 2*sd(filtered_metrics$degree)
    filtered_metrics$is_hub <- filtered_metrics$degree > threshold_hub
  }
  return(list(graph = filtered_graph, node_metrics = filtered_metrics))
}

create_improved_layout <- function(graph, scale_factor = 0.8, iterations = 100, start_temp = NULL) {
  if(is.null(start_temp)) {

    start_temp <- sqrt(vcount(graph))
  }

  layout_coords_raw <- layout_with_fr(
    graph,
    weights = E(graph)$weight,
    niter = iterations,
    start.temp = start_temp,
    grid = "nogrid"
  )

  layout_coords <- layout_coords_raw * scale_factor

  return(layout_coords)
}

normalize_layout <- function(layout, target_radius = 0.8) {
  distances <- sqrt(layout[,1]^2 + layout[,2]^2)
  max_dist <- max(distances)

  if(max_dist > 0) {
    layout_norm <- layout * (target_radius / max_dist)
  } else {
    layout_norm <- layout
  }
  return(layout_norm)
}

create_linewidth_scale <- function(correlations, min_width = 0.6, max_width = 2.0) {
  abs_cor <- abs(correlations)
  if(length(unique(abs_cor)) == 1) {
    return(rep(mean(c(min_width, max_width)), length(correlations)))
  }

  scaled_widths <- (abs_cor - min(abs_cor)) / (max(abs_cor) - min(abs_cor))
  final_widths <- min_width + scaled_widths * (max_width - min_width)

  return(final_widths)
}

cat("=== CARGANDO DATOS ===\n")

qs_raw <- read.delim("../Matriz_qs_short.tsv", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

qs_matrix <- as.data.frame(t(qs_raw[,-1]))
colnames(qs_matrix) <- qs_raw$accession
qs_matrix$MAG_ID <- rownames(qs_matrix)
qs_matrix$Site <- extract_site(qs_matrix$MAG_ID)

qs_genes <- setdiff(colnames(qs_matrix), c("MAG_ID", "Site"))
qs_matrix$n_qs_genes <- rowSums(qs_matrix[, qs_genes] > 0)

meta <- read.delim("./metadata2.tsv", header = TRUE, stringsAsFactors = FALSE) %>%
  rename(MAG_ID = name, Group = group_per_genus, Taxonomy = tax) %>%
  filter(Group %in% c("wide","rare"))

cat("Procesando información taxonómica...\n")
taxonomy_info <- meta %>%
  rowwise() %>%
  do({
    tax_info <- extract_taxonomy(.$Taxonomy)
    data.frame(
      MAG_ID = .$MAG_ID,
      Group = .$Group,
      Taxonomy = .$Taxonomy,
      Order = tax_info["order"],
      Genus = tax_info["genus"],
      stringsAsFactors = FALSE
    )
  }) %>% ungroup()

taxonomy_info$Tax_Label <- mapply(create_tax_label,
                                  taxonomy_info$MAG_ID,
                                  taxonomy_info$Order,
                                  taxonomy_info$Genus)
cat(sprintf("Información taxonómica procesada para %d MAGs\n", nrow(taxonomy_info)))

qs_with_meta <- merge(qs_matrix, taxonomy_info[,c("MAG_ID", "Group", "Order", "Genus", "Tax_Label")],
                      by = "MAG_ID", all.x = TRUE)

cat("\n=== ANÁLISIS DE RED GLOBAL ===\n")

qs_genes <- setdiff(colnames(qs_with_meta), c("MAG_ID", "Site", "Group", "n_qs_genes", "Order", "Genus", "Tax_Label"))
abundance_matrix <- as.matrix(qs_with_meta[, qs_genes])
rownames(abundance_matrix) <- qs_with_meta$MAG_ID

cat(sprintf("MAGs totales: %d\n", nrow(abundance_matrix)))
cat(sprintf("Genes QS totales: %d\n", length(qs_genes)))
cat(sprintf("Sparsity (proporción de ceros): %.3f\n", sum(abundance_matrix == 0) / length(abundance_matrix)))

mag_variance <- apply(abundance_matrix, 1, function(x) sum(x > 0))
valid_mags <- names(mag_variance[mag_variance > 0])
filtered_matrix <- abundance_matrix[valid_mags, , drop = FALSE]
cat(sprintf("MAGs válidos después de filtrado: %d\n", length(valid_mags)))

mag_matrix <- t(filtered_matrix)

cat("\nCalculando correlaciones Pearson...\n")
pearson_results <- calculate_correlations(mag_matrix, method = "pearson")
cat("Calculando correlaciones Spearman...\n")
spearman_results <- calculate_correlations(mag_matrix, method = "spearman")

cat("Aplicando corrección Benjamini-Hochberg...\n")
pearson_p_adjusted <- apply_fdr_correction(pearson_results$p_value)
spearman_p_adjusted <- apply_fdr_correction(spearman_results$p_value)

extract_edges <- function(cor_matrix, p_matrix, method_name) {
  edges <- data.frame()
  for (i in 1:(ncol(cor_matrix)-1)) {
    for (j in (i+1):ncol(cor_matrix)) {
      cor_val <- cor_matrix[i,j]
      p_val <- p_matrix[i,j]
      if (!is.na(cor_val) && !is.na(p_val) &&
          abs(cor_val) >= correlation_threshold && p_val <= p_value_threshold) {
        edges <- rbind(edges, data.frame(
          from = colnames(cor_matrix)[i],
          to = colnames(cor_matrix)[j],
          correlation = cor_val,
          p_value = p_val,
          method = method_name,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  return(edges)
}

pearson_edges <- extract_edges(pearson_results$correlation, pearson_p_adjusted, "pearson")
spearman_edges <- extract_edges(spearman_results$correlation, spearman_p_adjusted, "spearman")

cat(sprintf("\nAristas Pearson significativas: %d\n", nrow(pearson_edges)))
cat(sprintf("Aristas Spearman significativas: %d\n", nrow(spearman_edges)))

if (nrow(pearson_edges) > 0 && nrow(spearman_edges) > 0) {
  final_edges <- intersect_networks(pearson_edges, spearman_edges)
  cat(sprintf("Aristas en la intersección (Pearson ∩ Spearman): %d\n", nrow(final_edges)))
} else {
  final_edges <- data.frame()
}

if (nrow(final_edges) > 0) {
  g <- graph_from_data_frame(final_edges[,c("from","to")], directed = FALSE)
  E(g)$weight <- abs(final_edges$avg_correlation)
  E(g)$correlation <- final_edges$avg_correlation

  assign_taxonomy_to_graph <- function(graph, qs_data) {
    node_names <- V(graph)$name
    n_nodes <- length(node_names)
    groups <- character(n_nodes); sites <- character(n_nodes)
    n_qs_genes_vec <- numeric(n_nodes); orders <- character(n_nodes)
    genera <- character(n_nodes); tax_labels <- character(n_nodes)
    for(i in seq_len(n_nodes)) {
      node_name <- node_names[i]
      node_row <- qs_data[qs_data$MAG_ID == node_name, ]
      if(nrow(node_row) > 0) {
        groups[i] <- ifelse(is.na(node_row$Group[1]), "unknown", as.character(node_row$Group[1]))
        sites[i] <- ifelse(is.na(node_row$Site[1]), "Unknown", as.character(node_row$Site[1]))
        n_qs_genes_vec[i] <- ifelse(is.na(node_row$n_qs_genes[1]), 0, node_row$n_qs_genes[1])
        orders[i] <- ifelse(is.na(node_row$Order[1]), "Unknown_order", as.character(node_row$Order[1]))
        genera[i] <- ifelse(is.na(node_row$Genus[1]), "Unknown_genus", as.character(node_row$Genus[1]))
        if(!is.na(node_row$Tax_Label[1])) {
          tax_labels[i] <- as.character(node_row$Tax_Label[1])
        } else {
          if(orders[i] != "Unknown_order" && genera[i] != "Unknown_genus") {
            tax_labels[i] <- paste0(node_name, " (", orders[i], "; ", genera[i], ")")
          } else if(orders[i] != "Unknown_order") {
            tax_labels[i] <- paste0(node_name, " (", orders[i], ")")
          } else if(genera[i] != "Unknown_genus") {
            tax_labels[i] <- paste0(node_name, " (", genera[i], ")")
          } else {
            tax_labels[i] <- node_name
          }
        }
      } else {
        groups[i] <- "unknown"; sites[i] <- "Unknown"; n_qs_genes_vec[i] <- 0
        orders[i] <- "Unknown_order"; genera[i] <- "Unknown_genus"; tax_labels[i] <- node_name
      }
    }
    V(graph)$Group <- groups; V(graph)$Site <- sites
    V(graph)$n_qs_genes <- n_qs_genes_vec
    V(graph)$Order <- orders; V(graph)$Genus <- genera; V(graph)$Tax_Label <- tax_labels
    return(graph)
  }

  cat("\nAsignando información taxonómica a los nodos del grafo...\n")
  g <- assign_taxonomy_to_graph(g, qs_with_meta)

  cat("\n=== MÉTRICAS DE RED GLOBAL (ANTES DEL FILTRADO) ===\n")
  network_metrics <- data.frame(
    n_nodes = vcount(g),
    n_edges = ecount(g),
    density = edge_density(g),
    mean_degree = mean(degree(g)),
    clustering = transitivity(g, type = "global"),
    diameter = ifelse(is.connected(g), diameter(g), NA),
    avg_path_length = ifelse(is.connected(g), mean_distance(g), NA),
    modularity = modularity(cluster_louvain(g))
  )
  print(network_metrics)

  node_metrics <- data.frame(
    MAG_ID = V(g)$name,
    Group = V(g)$Group,
    Site = V(g)$Site,
    n_qs_genes = V(g)$n_qs_genes,
    Order = V(g)$Order,
    Genus = V(g)$Genus,
    Tax_Label = V(g)$Tax_Label,
    degree = degree(g),
    betweenness = betweenness(g, normalized = TRUE),
    closeness = closeness(g, normalized = TRUE),
    eigenvector = eigen_centrality(g)$vector,
    stringsAsFactors = FALSE
  )
  threshold_hub <- mean(node_metrics$degree) + 2*sd(node_metrics$degree)
  node_metrics$is_hub <- node_metrics$degree > threshold_hub
  cat(sprintf("\nNodos hub (degree > %.1f): %d\n", threshold_hub, sum(node_metrics$is_hub)))

  if(node_filter != "all") {
    cat(sprintf("\n=== APLICANDO FILTRADO DE NODOS (Top %d por %s) ===\n", node_filter, node_filter_method))
    filtered_result <- filter_nodes(g, node_metrics, node_filter, node_filter_method, node_filter)
    g_viz <- filtered_result$graph
    node_metrics_viz <- filtered_result$node_metrics
    network_metrics_viz <- data.frame(
      n_nodes = vcount(g_viz),
      n_edges = ecount(g_viz),
      density = edge_density(g_viz),
      mean_degree = mean(degree(g_viz)),
      clustering = transitivity(g_viz, type = "global"),
      diameter = ifelse(is.connected(g_viz), diameter(g_viz), NA),
      avg_path_length = ifelse(is.connected(g_viz), mean_distance(g_viz), NA),
      modularity = ifelse(vcount(g_viz) > 1, modularity(cluster_louvain(g_viz)), NA)
    )
    print(network_metrics_viz)
  } else {
    g_viz <- g
    node_metrics_viz <- node_metrics
    network_metrics_viz <- network_metrics
  }

} else {
  cat("\n!!! No se encontraron aristas significativas en la intersección !!!\n")
  g <- NULL; g_viz <- NULL
  network_metrics <- NULL; network_metrics_viz <- NULL
  node_metrics <- NULL; node_metrics_viz <- NULL
}

cat("\n=== GENERANDO VISUALIZACIONES ===\n")

if (!is.null(g_viz) && !is.null(node_metrics_viz)) {

  panel_A <- final_edges %>%
    mutate(cor_type = ifelse(avg_correlation > 0, "Positive", "Negative")) %>%
    ggplot(aes(x = avg_correlation, fill = cor_type)) +
    geom_histogram(bins = 30, alpha = 0.8, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Positive" = "#2166AC", "Negative" = "#B2182B")) +
    labs(
      title = "A) Distribution of Global Network Correlations",
      subtitle = sprintf("Total edges: %d (Pearson ∩ Spearman, |r| ≥ %.1f, FDR < %.2f)",
                         nrow(final_edges), correlation_threshold, p_value_threshold),
      x = "Average Correlation Coefficient",
      y = "Count"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10, face = "italic"),
      legend.position = c(0.85, 0.85)
    )

  if (normalize_BC_by_genome_count) {
    df_B <- node_metrics_viz %>%
      filter(Group %in% c("wide", "rare")) %>%
      group_by(Group) %>%
      summarise(
        n_genomes = n(),
        total_qs_genes = sum(n_qs_genes, na.rm = TRUE),
        mean_qs_genes_per_genome = total_qs_genes / n_genomes,
        .groups = "drop"
      )
    panel_B <- ggplot(df_B, aes(x = Group, y = mean_qs_genes_per_genome, fill = Group)) +
      geom_col(alpha = 0.85) +
      geom_text(aes(label = sprintf("n=%d", n_genomes)), vjust = -0.6, size = 3) +
      scale_fill_manual(values = taxa_colors) +
      labs(
        title = "B) QS genes per genome (normalized)",
        subtitle = "Promedio de genes QS por genoma (por grupo taxonómico)",
        x = "Taxonomic Group",
        y = "Mean QS genes / genome"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10, face = "italic"),
        legend.position = "none"
      )
  } else {
    panel_B <- node_metrics_viz %>%
      filter(Group %in% c("wide", "rare")) %>%
      ggplot(aes(x = Group, y = n_qs_genes, fill = Group)) +
      geom_boxplot(alpha = 0.7, outlier.shape = 21, outlier.fill = "white") +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) +
      scale_fill_manual(values = taxa_colors) +
      labs(
        title = "B) Number of QS Genes by Taxonomic Group",
        subtitle = ifelse(node_filter != "all",
                          sprintf("Top %d nodes by %s", node_filter, node_filter_method),
                          "Per-MAG distribution"),
        x = "Taxonomic Group",
        y = "Number of QS Genes Detected"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10, face = "italic"),
        legend.position = "none"
      )
  }

  if (normalize_BC_by_genome_count) {
    df_C <- node_metrics_viz %>%
      filter(Group %in% c("wide", "rare")) %>%
      group_by(Site, Group) %>%
      summarise(
        n_genomes = n(),
        total_qs_genes = sum(n_qs_genes, na.rm = TRUE),
        mean_qs_genes_per_genome = ifelse(n_genomes > 0, total_qs_genes / n_genomes, NA_real_),
        .groups = "drop"
      )
    panel_C <- ggplot(df_C, aes(x = Site, y = mean_qs_genes_per_genome, fill = Group)) +
      geom_col(position = position_dodge(width = 0.8), alpha = 0.9) +
      geom_text(aes(label = sprintf("n=%d", n_genomes)),
                position = position_dodge(width = 0.8), vjust = -0.6, size = 3) +
      scale_fill_manual(values = taxa_colors,
                        labels = c("wide" = "Widespread", "rare" = "Restricted", "unknown" = "Unknown")) +
      labs(
        title = "C) QS genes per genome by Site and Group (normalized)",
        subtitle = ifelse(node_filter != "all",
                          sprintf("Filtered network (%d nodes)", vcount(g_viz)),
                          "Promedio de genes QS por genoma en cada sitio"),
        x = "Site",
        y = "Mean QS genes / genome"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10, face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank()
      )
  } else {
    panel_C <- node_metrics_viz %>%
      filter(!is.na(Group)) %>%
      ggplot(aes(x = Site, y = n_qs_genes, fill = Group)) +
      geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
      scale_fill_manual(values = taxa_colors,
                        labels = c("wide" = "Widespread", "rare" = "Restricted", "unknown" = "Unknown")) +
      labs(
        title = "C) QS Gene Distribution by Site and Taxonomic Group",
        subtitle = ifelse(node_filter != "all",
                          sprintf("Filtered network (%d nodes)", vcount(g_viz)),
                          "Per-MAG distribution"),
        x = "Site",
        y = "Number of QS Genes"
      ) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10, face = "italic"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        legend.title = element_blank()
      )
  }

  set.seed(123)

  layout_coords <- create_improved_layout(
    g_viz,
    scale_factor = layout_scale_factor,
    iterations = layout_iterations,
    start_temp = NULL
  )

  node_viz_data <- dplyr::tibble(
    x = layout_coords[,1],
    y = layout_coords[,2],
    name = V(g_viz)$name,
    tax_label = V(g_viz)$Tax_Label,
    group = factor(V(g_viz)$Group, levels = c("wide", "rare", "unknown")),
    site = V(g_viz)$Site,
    degree = degree(g_viz),
    n_qs_genes = V(g_viz)$n_qs_genes,
    is_hub = ifelse(node_metrics_viz$is_hub[match(V(g_viz)$name, node_metrics_viz$MAG_ID)], TRUE, FALSE)
  ) %>% mutate(is_hub_chr = ifelse(is_hub, "TRUE", "FALSE"))

  edge_viz_data <- as.data.frame(get.edgelist(g_viz))
  names(edge_viz_data) <- c("from", "to")
  edge_viz_data$correlation <- E(g_viz)$correlation
  edge_viz_data$abs_correlation <- abs(edge_viz_data$correlation)
  edge_viz_data$linewidth <- create_linewidth_scale(
    edge_viz_data$correlation,
    min_width = edge_linewidth_range[1],
    max_width = edge_linewidth_range[2]
  )
  edge_viz_data$x1 <- node_viz_data$x[match(edge_viz_data$from, node_viz_data$name)]
  edge_viz_data$y1 <- node_viz_data$y[match(edge_viz_data$from, node_viz_data$name)]
  edge_viz_data$x2 <- node_viz_data$x[match(edge_viz_data$to, node_viz_data$name)]
  edge_viz_data$y2 <- node_viz_data$y[match(edge_viz_data$to, node_viz_data$name)]

  deg_vec <- node_viz_data$degree
  is_hub_vec <- node_viz_data$is_hub
  is_topdeg_vec <- deg_vec >= stats::quantile(deg_vec, probs = label_topdegree_quantile, na.rm = TRUE)
  n_nodes_viz <- nrow(node_viz_data)
  if (label_mode == "none") {
    label_mask <- rep(FALSE, n_nodes_viz)
  } else if (label_mode == "all") {
    label_mask <- rep(TRUE, n_nodes_viz)
  } else if (label_mode == "hubs") {
    label_mask <- is_hub_vec
  } else if (label_mode == "topdegree") {
    label_mask <- is_topdeg_vec
  } else {
    if (n_nodes_viz <= 30)      label_mask <- rep(TRUE, n_nodes_viz)
    else if (n_nodes_viz <= 100) label_mask <- (is_hub_vec | is_topdeg_vec)
    else                         label_mask <- is_hub_vec
  }

  label_text <- dplyr::case_when(
    label_text_field == "mag_id" ~ node_viz_data$name,
    label_text_field == "genus"  ~ dplyr::if_else(is.na(V(g_viz)$Genus) | V(g_viz)$Genus == "Unknown_genus",
                                                  node_viz_data$name, V(g_viz)$Genus),
    label_text_field == "order"  ~ dplyr::if_else(is.na(V(g_viz)$Order) | V(g_viz)$Order == "Unknown_order",
                                                  node_viz_data$name, V(g_viz)$Order),
    TRUE                         ~ node_viz_data$tax_label
  )

  labels_df <- node_viz_data %>%
    dplyr::transmute(x, y, name, label = dplyr::if_else(label_mask, label_text, "")) %>%
    dplyr::filter(label != "")

  panel_D <- ggplot() +

    geom_segment(
      data = edge_viz_data,
      aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = abs_correlation),
      color = edge_color, alpha = edge_alpha
    ) +

    scale_linewidth_continuous(
      name = "Correlation\n(absolute)",
      range = edge_linewidth_range,
      breaks = pretty(edge_viz_data$abs_correlation, n = 4),
      labels = function(x) sprintf("%.2f", x)
    ) +
    geom_point(
      data = node_viz_data,
      aes(x = x, y = y, fill = group, size = n_qs_genes, shape = is_hub_chr),
      alpha = 0.8, stroke = 0.5
    ) +
    scale_fill_manual(values = taxa_colors,
                      name = "Group",
                      labels = c("wide" = "Widespread", "rare" = "Restricted", "unknown" = "Unknown")) +
    scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 24),
                       name = "Hub", labels = c("FALSE" = "No", "TRUE" = "Yes")) +
    scale_size_continuous(name = "QS Genes",
                          range = node_size_range,
                          breaks = pretty(node_viz_data$n_qs_genes, n = 4)) +
    { if (nrow(labels_df) > 0)
      ggrepel::geom_text_repel(
        data = labels_df,
        aes(x = x, y = y, label = label),
        size = label_size,
        segment.color = label_segment_color,
        max.overlaps = Inf,
        force = 2,
        box.padding = 0.5
      )
      else NULL } +
    labs(
      title = "D) Global QS Network - Line Thickness = Correlation Strength",
      subtitle = sprintf("%d nodes, %d edges (Line thickness = |correlation|, Node size = QS genes)",
                         vcount(g_viz), ecount(g_viz))
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10, face = "italic"),
      legend.position = "right",
      legend.box = "vertical",
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_equal()

  combined_plot <- (panel_A | panel_B) / (panel_C | panel_D) +
    plot_layout(heights = c(1, 1))

  pdf(output_file, width = 8, height = 6, useDingbats = FALSE)
  print(combined_plot)
  dev.off()
  cat(sprintf("\n✓ Figura guardada en: %s\n", output_file))

} else {
  cat("\n!!! No se pudo generar la visualización - No hay datos de red !!!\n")
}

cat("\n=== EXPORTANDO RESULTADOS ===\n")

if (!is.null(network_metrics_viz)) {
  network_metrics_export <- network_metrics_viz
  network_metrics_export$filtered <- ifelse(node_filter == "all", "no", "yes")
  network_metrics_export$filter_method <- ifelse(node_filter == "all", NA, node_filter_method)
  network_metrics_export$filter_n_nodes <- ifelse(node_filter == "all", NA, node_filter)
  write.table(network_metrics_export, "global_network_metrics.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  cat("✓ global_network_metrics.tsv\n")
}

if (!is.null(final_edges) && nrow(final_edges) > 0) {
  write.table(final_edges, "global_network_edges.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  cat("✓ global_network_edges.tsv\n")
}

if (!is.null(node_metrics_viz)) {
  write.table(node_metrics_viz, "global_node_metrics.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  cat("✓ global_node_metrics.tsv\n")
}

if (!is.null(node_metrics_viz) && normalize_BC_by_genome_count) {
  df_B_out <- node_metrics_viz %>%
    filter(Group %in% c("wide", "rare")) %>%
    group_by(Group) %>%
    summarise(
      n_genomes = n(),
      total_qs_genes = sum(n_qs_genes, na.rm = TRUE),
      mean_qs_genes_per_genome = total_qs_genes / n_genomes,
      .groups = "drop"
    )
  write.table(df_B_out, "group_mean_qs_genes_per_genome.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  cat("✓ group_mean_qs_genes_per_genome.tsv\n")

  df_C_out <- node_metrics_viz %>%
    filter(Group %in% c("wide", "rare")) %>%
    group_by(Site, Group) %>%
    summarise(
      n_genomes = n(),
      total_qs_genes = sum(n_qs_genes, na.rm = TRUE),
      mean_qs_genes_per_genome = ifelse(n_genomes > 0, total_qs_genes / n_genomes, NA_real_),
      .groups = "drop"
    )
  write.table(df_C_out, "site_group_mean_qs_genes_per_genome.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
  cat("✓ site_group_mean_qs_genes_per_genome.tsv\n")
}

if (!is.null(node_metrics_viz)) {
  cat("\n=== ANÁLISIS ESTADÍSTICO ===\n")
  wide_data <- filter(node_metrics_viz, Group == "wide")
  rare_data <- filter(node_metrics_viz, Group == "rare")
  if (nrow(wide_data) > 2 && nrow(rare_data) > 2) {
    cat("\n--- Tests de Wilcoxon (wide vs rare, por MAG) ---\n")
    test_qs_genes <- wilcox.test(wide_data$n_qs_genes, rare_data$n_qs_genes, exact = FALSE)
    cat(sprintf("N° genes QS por MAG: p-value = %.4f %s\n",
                test_qs_genes$p.value, ifelse(test_qs_genes$p.value < 0.05, "***", "")))
    test_degree <- wilcox.test(wide_data$degree, rare_data$degree, exact = FALSE)
    cat(sprintf("Degree por MAG: p-value = %.4f %s\n",
                test_degree$p.value, ifelse(test_degree$p.value < 0.05, "***", "")))
  }
}

if (!is.null(g_viz) && !is.null(node_metrics_viz)) {
  cat("\n=== GENERANDO VISUALIZACIÓN CIRCULAR CON GROSOR DE LÍNEA ===\n")

  node_order <- node_metrics_viz %>% arrange(Group, desc(n_qs_genes)) %>% pull(MAG_ID)
  ordered_vertices <- match(node_order, V(g_viz)$name); ordered_vertices <- ordered_vertices[!is.na(ordered_vertices)]

  layout_circular_raw <- layout_in_circle(g_viz, order = ordered_vertices)
  layout_circular <- normalize_layout(layout_circular_raw, target_radius = layout_target_radius)

  node_viz_circular <- dplyr::tibble(
    x = layout_circular[,1],
    y = layout_circular[,2],
    name = V(g_viz)$name,
    tax_label = V(g_viz)$Tax_Label,
    group = factor(V(g_viz)$Group, levels = c("wide", "rare", "unknown")),
    site = V(g_viz)$Site,
    degree = degree(g_viz),
    n_qs_genes = V(g_viz)$n_qs_genes
  )

  edge_viz_circular <- as.data.frame(get.edgelist(g_viz))
  names(edge_viz_circular) <- c("from", "to")
  edge_viz_circular$correlation <- E(g_viz)$correlation
  edge_viz_circular$abs_correlation <- abs(edge_viz_circular$correlation)
  edge_viz_circular$linewidth <- create_linewidth_scale(
    edge_viz_circular$correlation,
    min_width = edge_linewidth_range[1],
    max_width = edge_linewidth_range[2]
  )
  edge_viz_circular$x1 <- node_viz_circular$x[match(edge_viz_circular$from, node_viz_circular$name)]
  edge_viz_circular$y1 <- node_viz_circular$y[match(edge_viz_circular$from, node_viz_circular$name)]
  edge_viz_circular$x2 <- node_viz_circular$x[match(edge_viz_circular$to, node_viz_circular$name)]
  edge_viz_circular$y2 <- node_viz_circular$y[match(edge_viz_circular$to, node_viz_circular$name)]

  deg_vec <- node_viz_circular$degree
  is_hub_vec <- deg_vec > (mean(deg_vec) + 2*sd(deg_vec))
  is_topdeg_vec <- deg_vec >= stats::quantile(deg_vec, probs = label_topdegree_quantile, na.rm = TRUE)
  n_nodes_circ <- nrow(node_viz_circular)
  if (label_mode == "none") {
    label_mask_circ <- rep(FALSE, n_nodes_circ)
  } else if (label_mode == "all") {
    label_mask_circ <- rep(TRUE, n_nodes_circ)
  } else if (label_mode == "hubs") {
    label_mask_circ <- is_hub_vec
  } else if (label_mode == "topdegree") {
    label_mask_circ <- is_topdeg_vec
  } else {
    if (n_nodes_circ <= 50)      label_mask_circ <- rep(TRUE, n_nodes_circ)
    else if (n_nodes_circ <= 150) label_mask_circ <- (is_hub_vec | is_topdeg_vec)
    else                          label_mask_circ <- is_hub_vec
  }
  label_text_circ <- dplyr::case_when(
    label_text_field == "mag_id" ~ node_viz_circular$name,
    label_text_field == "genus"  ~ dplyr::if_else(is.na(V(g_viz)$Genus) | V(g_viz)$Genus == "Unknown_genus",
                                                  node_viz_circular$name, V(g_viz)$Genus),
    label_text_field == "order"  ~ dplyr::if_else(is.na(V(g_viz)$Order) | V(g_viz)$Order == "Unknown_order",
                                                  node_viz_circular$name, V(g_viz)$Order),
    TRUE                         ~ node_viz_circular$tax_label
  )

  labels_df_circ <- node_viz_circular %>%
    dplyr::transmute(x, y, name, label = dplyr::if_else(label_mask_circ, label_text_circ, "")) %>%
    dplyr::filter(label != "")

  circular_plot <- ggplot() +

    geom_segment(
      data = edge_viz_circular,
      aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = abs_correlation),
      color = edge_color, alpha = edge_alpha * 0.6
    ) +
    scale_linewidth_continuous(
      name = "Correlation\n(absolute)",
      range = edge_linewidth_range,
      breaks = pretty(edge_viz_circular$abs_correlation, n = 4),
      labels = function(x) sprintf("%.2f", x)
    ) +
    geom_point(
      data = node_viz_circular,
      aes(x = x, y = y, fill = group, size = n_qs_genes),
      shape = 21, alpha = 0.9, stroke = 1, color = "black"
    ) +
    scale_fill_manual(values = c("wide" = "#2E86AB", "rare" = "#F24236", "unknown" = "#95A99C"),
                      name = "Taxonomic\nGroup",
                      labels = c("wide" = "Widespread", "rare" = "Restricted", "unknown" = "Unknown")) +
    scale_size_continuous(name = "QS Genes",
                          range = node_size_range,
                          breaks = pretty(node_viz_circular$n_qs_genes, n = 4),
                          guide = guide_legend(override.aes = list(shape = 21,
                                                                   fill = "grey70",
                                                                   stroke = 1))) +
    { if (nrow(labels_df_circ) > 0)
      ggrepel::geom_text_repel(
        data = labels_df_circ,
        aes(x = x, y = y, label = label),
        size = ifelse(nrow(node_viz_circular) <= 20, 3.5, ifelse(nrow(node_viz_circular) <= 50, 3, 2.5)),
        segment.color = "grey40",
        segment.size = 0.2,
        segment.alpha = 0.5,
        max.overlaps = Inf,
        force = 2,
        force_pull = 0.5,
        box.padding = 0.5,
        point.padding = 0.3
      )
      else NULL } +
    labs(
      title = "Global QS Network - Circular Layout with Line Thickness = Correlation",
      subtitle = sprintf("Network with %d MAGs and %d connections | Line thickness = |correlation| | Node size = QS genes",
                         vcount(g_viz), ecount(g_viz))
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.box = "vertical",
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    coord_equal(clip = "off")

  circular_output <- "QS_global_network_circular_linewidth.pdf"
  pdf(circular_output, width = 16, height = 16, useDingbats = FALSE)
  print(circular_plot); dev.off()
  cat(sprintf("✓ Visualización circular con grosor de línea guardada en: %s\n", circular_output))

  png_output <- "QS_global_network_circular_linewidth.png"
  png(png_output, width = 3200, height = 3200, res = 200)
  print(circular_plot); dev.off()
  cat(sprintf("✓ Versión PNG guardada en: %s\n", png_output))
}

cat("\n=== ANÁLISIS COMPLETADO ===\n")
if(node_filter != "all") {
  cat(sprintf("\n✓ Red filtrada a %d nodos usando criterio: %s\n", node_filter, node_filter_method))
} else {
  cat("\n✓ Se mostró la red completa\n")
}
cat(ifelse(normalize_BC_by_genome_count,
           "✓ Paneles B y C normalizados por número de genomas (promedio genes/genoma)\n",
           "✓ Paneles B y C muestran distribución por MAG (no normalizado)\n"))
cat("✅ Grosor de línea representa correlación (en lugar de color)\n")
cat("✅ Leyenda de grosor de línea y tamaño de burbujas incluidas\n")
cat(sprintf("✅ Layout con proximidad controlada (scale: %.1f, target_radius: %.1f)\n",
            layout_scale_factor, layout_target_radius))
cat("✅ ERROR g_viz corregido - temperatura inicial calculada dinámicamente\n")

cat("\n=== GENERANDO ESTADÍSTICAS PARA MANUSCRITO ===\n")

if(!is.null(abundance_matrix)) {

  total_possible_pairs <- choose(ncol(abundance_matrix), 2)
  total_mags <- ncol(abundance_matrix)
  total_samples <- nrow(abundance_matrix)

  cat(sprintf("Dataset characteristics:\n"))
  cat(sprintf("- Total MAGs: %d\n", total_mags))
  cat(sprintf("- Total samples: %d\n", total_samples))
  cat(sprintf("- Total possible MAG pairs: %d\n", total_possible_pairs))
  cat(sprintf("- Current max_double_zeros threshold: %.2f\n", max_double_zeros))

  eligible_pairs <- 0
  double_zero_props <- c()

  cat("\nCalculating double-zero statistics...\n")

  for (i in 1:(ncol(abundance_matrix)-1)) {
    for (j in (i+1):ncol(abundance_matrix)) {
      both_zero <- sum(abundance_matrix[,i] == 0 & abundance_matrix[,j] == 0)
      double_zero_prop <- both_zero / nrow(abundance_matrix)
      double_zero_props <- c(double_zero_props, double_zero_prop)

      if (double_zero_prop <= max_double_zeros) {
        eligible_pairs <- eligible_pairs + 1
      }
    }
  }

  cat("\n--- DOUBLE-ZERO STATISTICS ---\n")
  cat(sprintf("Mean proportion of double-zeros: %.3f\n", mean(double_zero_props)))
  cat(sprintf("Median proportion of double-zeros: %.3f\n", median(double_zero_props)))
  cat(sprintf("Min proportion of double-zeros: %.3f\n", min(double_zero_props)))
  cat(sprintf("Max proportion of double-zeros: %.3f\n", max(double_zero_props)))

  cat("\n--- FILTERING STATISTICS ---\n")
  cat(sprintf("Pairs after double-zero filter (≤%.2f): %d (%.1f%% retained)\n",
              max_double_zeros, eligible_pairs, (eligible_pairs/total_possible_pairs)*100))
  cat(sprintf("Pairs excluded by double-zero filter: %d (%.1f%% excluded)\n",
              total_possible_pairs - eligible_pairs,
              ((total_possible_pairs - eligible_pairs)/total_possible_pairs)*100))

  if(!is.null(final_edges) && nrow(final_edges) > 0) {
    cat("\n--- FINAL NETWORK STATISTICS ---\n")
    cat(sprintf("Final network edges: %d\n", nrow(final_edges)))
    cat(sprintf("Percentage of eligible pairs that became edges: %.1f%%\n",
                (nrow(final_edges)/eligible_pairs)*100))
    cat(sprintf("Percentage of total possible pairs that became edges: %.1f%%\n",
                (nrow(final_edges)/total_possible_pairs)*100))

    cat(sprintf("Mean correlation strength: %.3f\n", mean(abs(final_edges$avg_correlation))))
    cat(sprintf("Range of correlations: %.3f to %.3f\n",
                min(final_edges$avg_correlation), max(final_edges$avg_correlation)))

    pos_corr <- sum(final_edges$avg_correlation > 0)
    neg_corr <- sum(final_edges$avg_correlation < 0)
    cat(sprintf("Positive correlations: %d (%.1f%%)\n", pos_corr, (pos_corr/nrow(final_edges))*100))
    cat(sprintf("Negative correlations: %d (%.1f%%)\n", neg_corr, (neg_corr/nrow(final_edges))*100))
  } else {
    cat("\n--- FINAL NETWORK STATISTICS ---\n")
    cat("No significant edges found in final network\n")
  }

  cat("\n--- THRESHOLD COMPARISON ---\n")
  thresholds <- c(0.50, 0.70, 0.80, 0.90, 0.95)

  for(threshold in thresholds) {
    eligible_at_threshold <- sum(double_zero_props <= threshold)
    cat(sprintf("At threshold %.2f: %d pairs eligible (%.1f%% of total)\n",
                threshold, eligible_at_threshold,
                (eligible_at_threshold/total_possible_pairs)*100))
  }

  cat("\n--- TABLE FOR MANUSCRIPT ---\n")
  cat("| Filtering Step | MAG Pairs | Percentage Retained |\n")
  cat("|---------------|-----------|-------------------|\n")
  cat(sprintf("| Total possible pairs | %d | 100.0%% |\n", total_possible_pairs))
  cat(sprintf("| After double-zero filter (≤%.2f) | %d | %.1f%% |\n",
              max_double_zeros, eligible_pairs, (eligible_pairs/total_possible_pairs)*100))

  if(!is.null(final_edges) && nrow(final_edges) > 0) {
    cat(sprintf("| After correlation threshold (|r| ≥ %.1f) | %d | %.1f%% |\n",
                correlation_threshold, nrow(final_edges), (nrow(final_edges)/total_possible_pairs)*100))
    cat(sprintf("| Final network (Pearson ∩ Spearman) | %d | %.1f%% |\n",
                nrow(final_edges), (nrow(final_edges)/total_possible_pairs)*100))
  }

  cat("\n--- TEXT FOR METHODS SECTION ---\n")
  cat(sprintf("\"MAG pairs exhibiting >%.0f%% co-absence across all samples ", (max_double_zeros*100)))
  cat(sprintf("(max_double_zeros = %.2f) were excluded from correlation analysis. ", max_double_zeros))
  cat(sprintf("This resulted in %d of %d possible MAG pairs (%.1f%%) being retained ",
              eligible_pairs, total_possible_pairs, (eligible_pairs/total_possible_pairs)*100))
  cat("for correlation analysis.\"\n")

} else {
  cat("abundance_matrix not found. Please run the data loading section first.\n")
}

cat("\n=== ESTADÍSTICAS COMPLETADAS ===\n")

if (!is.null(g_viz) && !is.null(node_metrics_viz)) {
  cat("\n=== GENERANDO VISUALIZACIÓN ESTILO MICROBIOME JOURNAL ===\n")

  create_microbiome_layout <- function(graph, node_metrics) {

    group_info <- node_metrics %>%
      group_by(Group) %>%
      summarise(
        n_nodes = n(),
        avg_degree = mean(degree),
        total_qs = sum(n_qs_genes),
        .groups = "drop"
      ) %>%
      arrange(desc(avg_degree))

    layout_coords_raw <- layout_with_graphopt(
      graph,
      niter = 500,
      charge = 0.001,
      mass = 30,
      spring.length = 30,
      spring.constant = 1,
      max.sa.movement = 5
    )

    node_data_temp <- data.frame(
      name = V(graph)$name,
      group = V(graph)$Group,
      x = layout_coords_raw[,1],
      y = layout_coords_raw[,2],
      stringsAsFactors = FALSE
    )

    group_centroids <- node_data_temp %>%
      group_by(group) %>%
      summarise(
        cx = mean(x),
        cy = mean(y),
        .groups = "drop"
      )

    target_positions <- data.frame(
      group = c("wide", "rare", "unknown"),
      target_x = c(-2, 2, 0),
      target_y = c(-1, -1, 2),
      stringsAsFactors = FALSE
    )

    layout_final <- layout_coords_raw

    for(i in 1:nrow(group_centroids)) {
      group_name <- group_centroids$group[i]
      target_info <- target_positions[target_positions$group == group_name, ]

      if(nrow(target_info) > 0) {

        group_nodes <- which(V(graph)$Group == group_name)

        if(length(group_nodes) > 0) {

          current_cx <- group_centroids$cx[i]
          current_cy <- group_centroids$cy[i]
          offset_x <- target_info$target_x - current_cx
          offset_y <- target_info$target_y - current_cy

          layout_final[group_nodes, 1] <- layout_final[group_nodes, 1] + offset_x
          layout_final[group_nodes, 2] <- layout_final[group_nodes, 2] + offset_y
        }
      }
    }

    return(layout_final)
  }

  create_separated_circular_layout <- function(graph, node_metrics) {

    wide_nodes <- node_metrics[node_metrics$Group == "wide", "MAG_ID"]
    rare_nodes <- node_metrics[node_metrics$Group == "rare", "MAG_ID"]
    unknown_nodes <- node_metrics[node_metrics$Group == "unknown", "MAG_ID"]

    layout_coords <- matrix(0, nrow = vcount(graph), ncol = 2)
    rownames(layout_coords) <- V(graph)$name

    if(length(wide_nodes) > 0) {
      wide_indices <- match(wide_nodes, V(graph)$name)
      wide_indices <- wide_indices[!is.na(wide_indices)]

      angles_wide <- seq(pi/2, 3*pi/2, length.out = length(wide_indices) + 1)[1:length(wide_indices)]
      radius_wide <- 3

      layout_coords[wide_indices, 1] <- radius_wide * cos(angles_wide)
      layout_coords[wide_indices, 2] <- radius_wide * sin(angles_wide)
    }

    if(length(rare_nodes) > 0) {
      rare_indices <- match(rare_nodes, V(graph)$name)
      rare_indices <- rare_indices[!is.na(rare_indices)]

      angles_rare <- seq(-pi/2, pi/2, length.out = length(rare_indices) + 1)[1:length(rare_indices)]
      radius_rare <- 3

      layout_coords[rare_indices, 1] <- radius_rare * cos(angles_rare)
      layout_coords[rare_indices, 2] <- radius_rare * sin(angles_rare)
    }

    if(length(unknown_nodes) > 0) {
      unknown_indices <- match(unknown_nodes, V(graph)$name)
      unknown_indices <- unknown_indices[!is.na(unknown_indices)]

      angles_unknown <- seq(0, pi, length.out = length(unknown_indices) + 1)[1:length(unknown_indices)]
      radius_unknown <- 2

      layout_coords[unknown_indices, 1] <- radius_unknown * cos(angles_unknown)
      layout_coords[unknown_indices, 2] <- radius_unknown * sin(angles_unknown)
    }

    return(layout_coords)
  }

  set.seed(123)

  layout_method1 <- create_microbiome_layout(g_viz, node_metrics_viz)

  layout_method2 <- create_separated_circular_layout(g_viz, node_metrics_viz)

  layout_coords <- layout_method2

  node_viz_microbiome <- dplyr::tibble(
    x = layout_coords[,1],
    y = layout_coords[,2],
    name = V(g_viz)$name,
    tax_label = V(g_viz)$Tax_Label,
    group = factor(V(g_viz)$Group, levels = c("wide", "rare", "unknown")),
    site = V(g_viz)$Site,
    degree = degree(g_viz),
    n_qs_genes = V(g_viz)$n_qs_genes,
    betweenness = betweenness(g_viz, normalized = TRUE)
  )

  edge_viz_microbiome <- as.data.frame(get.edgelist(g_viz))
  names(edge_viz_microbiome) <- c("from", "to")
  edge_viz_microbiome$correlation <- E(g_viz)$correlation
  edge_viz_microbiome$abs_correlation <- abs(edge_viz_microbiome$correlation)
  edge_viz_microbiome$correlation_type <- ifelse(edge_viz_microbiome$correlation > 0, "Positive", "Negative")

  edge_viz_microbiome$linewidth <- create_linewidth_scale(
    edge_viz_microbiome$correlation,
    min_width = 0.3,
    max_width = 2.0
  )

  edge_viz_microbiome$x1 <- node_viz_microbiome$x[match(edge_viz_microbiome$from, node_viz_microbiome$name)]
  edge_viz_microbiome$y1 <- node_viz_microbiome$y[match(edge_viz_microbiome$from, node_viz_microbiome$name)]
  edge_viz_microbiome$x2 <- node_viz_microbiome$x[match(edge_viz_microbiome$to, node_viz_microbiome$name)]
  edge_viz_microbiome$y2 <- node_viz_microbiome$y[match(edge_viz_microbiome$to, node_viz_microbiome$name)]

  from_group <- node_viz_microbiome$group[match(edge_viz_microbiome$from, node_viz_microbiome$name)]
  to_group <- node_viz_microbiome$group[match(edge_viz_microbiome$to, node_viz_microbiome$name)]

  edge_viz_microbiome$connection_type <- case_when(
    from_group == to_group ~ "Intra-group",
    from_group != to_group ~ "Inter-group"
  )

  deg_vec <- node_viz_microbiome$degree

  labels_df_microbiome <- node_viz_microbiome %>%
    dplyr::mutate(
      order_clean = ifelse(is.na(V(g_viz)$Order) | V(g_viz)$Order == "Unknown_order" | V(g_viz)$Order == "",
                           "Unknown", V(g_viz)$Order),
      genus_clean = ifelse(is.na(V(g_viz)$Genus) | V(g_viz)$Genus == "Unknown_genus" | V(g_viz)$Genus == "",
                           "Unknown", V(g_viz)$Genus),

      label = paste0(name, " (", order_clean, "; ", genus_clean, ")"),

      show_label = TRUE
    ) %>%
    dplyr::filter(show_label) %>%
    dplyr::select(x, y, name, label)

  microbiome_colors <- c(
    "wide" = "#3498db",
    "rare" = "#e74c3c",
    "unknown" = "#95a5a6"
  )

  connection_colors <- c(
    "Intra-group" = "#34495e",
    "Inter-group" = "#e67e22"
  )

  microbiome_plot <- ggplot() +

    geom_segment(
      data = edge_viz_microbiome %>% filter(connection_type == "Intra-group"),
      aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = abs_correlation),
      color = connection_colors["Intra-group"], alpha = 0.6
    ) +
    geom_segment(
      data = edge_viz_microbiome %>% filter(connection_type == "Inter-group"),
      aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = abs_correlation),
      color = connection_colors["Inter-group"], alpha = 0.8
    ) +

    scale_linewidth_continuous(
      name = "Correlation\n(absolute)",
      range = c(0.3, 2.0),
      breaks = pretty(edge_viz_microbiome$abs_correlation, n = 4),
      labels = function(x) sprintf("%.2f", x)
    ) +

    geom_point(
      data = node_viz_microbiome,
      aes(x = x, y = y, fill = group, size = n_qs_genes),
      shape = 21, alpha = 0.9, stroke = 0.5, color = "black"
    ) +

    scale_fill_manual(
      values = microbiome_colors,
      name = "Taxonomic Group",
      labels = c("wide" = "Widespread", "rare" = "Restricted", "unknown" = "Unknown")
    ) +
    scale_size_continuous(
      name = "QS Genes",
      range = c(2, 8),
      breaks = pretty(node_viz_microbiome$n_qs_genes, n = 4)
    ) +

    { if (nrow(labels_df_microbiome) > 0)
      ggrepel::geom_text_repel(
        data = labels_df_microbiome,
        aes(x = x, y = y, label = label),
        size = 2.2,
        segment.color = "grey40",
        segment.size = 0.2,
        max.overlaps = 50,
        force = 2,
        force_pull = 0.2,
        box.padding = 0.35,
        point.padding = 0.25,
        min.segment.length = 0.1,
        xlim = c(NA, NA),
        ylim = c(NA, NA),
        direction = "both"
      )
      else NULL } +

    labs(
      title = "Global QS Network - Microbiome Journal Style",
      subtitle = sprintf("Network with %d MAGs and %d connections | Labels: MAG_ID (Order; Genus)",
                         vcount(g_viz), ecount(g_viz)),
      caption = "Dark lines = Intra-group | Orange lines = Inter-group | Node size = QS genes"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
      plot.caption = element_text(size = 10, hjust = 0.5, color = "grey40"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 9),
      legend.box = "vertical",
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    coord_equal(clip = "off")

  set.seed(456)
  layout_compact <- layout_with_fr(g_viz, niter = 300) * 0.8

  node_viz_compact <- node_viz_microbiome
  node_viz_compact$x <- layout_compact[,1]
  node_viz_compact$y <- layout_compact[,2]

  edge_viz_compact <- edge_viz_microbiome
  edge_viz_compact$x1 <- node_viz_compact$x[match(edge_viz_compact$from, node_viz_compact$name)]
  edge_viz_compact$y1 <- node_viz_compact$y[match(edge_viz_compact$from, node_viz_compact$name)]
  edge_viz_compact$x2 <- node_viz_compact$x[match(edge_viz_compact$to, node_viz_compact$name)]
  edge_viz_compact$y2 <- node_viz_compact$y[match(edge_viz_compact$to, node_viz_compact$name)]

  labels_df_compact <- node_viz_compact %>%
    dplyr::mutate(
      order_clean = ifelse(is.na(V(g_viz)$Order) | V(g_viz)$Order == "Unknown_order" | V(g_viz)$Order == "",
                           "Unknown", V(g_viz)$Order),
      genus_clean = ifelse(is.na(V(g_viz)$Genus) | V(g_viz)$Genus == "Unknown_genus" | V(g_viz)$Genus == "",
                           "Unknown", V(g_viz)$Genus),

      label = paste0(name, " (", order_clean, "; ", genus_clean, ")"),

      show_label = TRUE
    ) %>%
    dplyr::filter(show_label) %>%
    dplyr::select(x, y, name, label)

  compact_plot <- ggplot() +

    geom_segment(
      data = edge_viz_compact,
      aes(x = x1, y = y1, xend = x2, yend = y2, linewidth = abs_correlation),
      color = "grey50", alpha = 0.5
    ) +
    scale_linewidth_continuous(
      name = "Correlation\n(absolute)",
      range = c(0.2, 1.5),
      guide = guide_legend(override.aes = list(color = "black", alpha = 1))
    ) +

    geom_point(
      data = node_viz_compact,
      aes(x = x, y = y, fill = group, size = n_qs_genes),
      shape = 21, alpha = 0.8, stroke = 0.3, color = "white"
    ) +
    scale_fill_manual(
      values = microbiome_colors,
      name = "Group"
    ) +
    scale_size_continuous(
      name = "QS Genes",
      range = c(3, 10)
    ) +

    { if (nrow(labels_df_compact) > 0)
      ggrepel::geom_text_repel(
        data = labels_df_compact,
        aes(x = x, y = y, label = label),
        size = 2.0,
        segment.color = "grey40",
        segment.size = 0.2,
        max.overlaps = 50,
        force = 2.5,
        force_pull = 0.1,
        box.padding = 0.4,
        point.padding = 0.3,
        min.segment.length = 0.1,
        direction = "both"
      )
      else NULL } +
    labs(
      title = "Global QS Network - Compact Layout",
      subtitle = sprintf("%d MAGs, %d connections | Labels: MAG_ID (Order; Genus)",
                         vcount(g_viz), ecount(g_viz))
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_equal()

  microbiome_output <- "QS_global_network_microbiome_style.pdf"
  pdf(microbiome_output, width = 12, height = 10, useDingbats = FALSE)
  print(microbiome_plot)
  dev.off()
  cat(sprintf("✓ Visualización estilo Microbiome Journal guardada en: %s\n", microbiome_output))

  compact_output <- "QS_global_network_compact_style.pdf"
  pdf(compact_output, width = 10, height = 8, useDingbats = FALSE)
  print(compact_plot)
  dev.off()
  cat(sprintf("✓ Visualización compacta guardada en: %s\n", compact_output))

  png(gsub(".pdf", ".png", microbiome_output), width = 2400, height = 2000, res = 200)
  print(microbiome_plot)
  dev.off()

  png(gsub(".pdf", ".png", compact_output), width = 2000, height = 1600, res = 200)
  print(compact_plot)
  dev.off()

  cat("\n=== VISUALIZACIÓN ESTILO MICROBIOME JOURNAL COMPLETADA ===\n")
  cat("✅ Layout sectorial con grupos separados\n")
  cat("✅ Conexiones intra-grupo vs inter-grupo diferenciadas\n")
  cat("✅ Versión compacta adicional incluida\n")
  cat("✅ Colores inspirados en papers de microbioma\n")
  cat("✅ Etiquetas muestran: MAG_ID (Orden; Genus)\n")
  cat("✅ Todas las etiquetas visibles (ajustar 'show_label' si es necesario)\n")

} else {
  cat("\n!!! No se pudo generar la visualización - No hay datos de red !!!\n")
}
