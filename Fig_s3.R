# Fig_s3.R
# Description: Compares polysaccharide lyases (PLs) and vitamin biosynthesis gene counts between
#              widespread and restricted MAG groups using bar plots, dot plots, and statistical tests.
# Purpose: Performs Wilcoxon rank-sum tests and generates four visualization variants for panels A-D:
#          mean bar plots with SE (A, B) and per-MAG dot plots (C, D) in clean, numbered, shortened-name,
#          and horizontal orientations. Exports PDF, TIFF, and CSV files for mBio supplementary figures.

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(cowplot)

if (!dir.exists("figures_output")) {
  dir.create("figures_output")
}

font_sizes <- list(
  panel_A = list(base = 11, axis_text = 10, axis_title = 11, legend = 9),
  panel_B = list(base = 11, axis_text = 10, axis_title = 11, legend = 9),
  panel_C = list(base = 10, axis_text_y = 6, axis_text_x = 9, axis_title = 10, legend = 8),
  panel_D = list(base = 10, axis_text_y = 6, axis_text_x = 9, axis_title = 10, legend = 8)
)

vitamin_pl_data <- read.delim("pls_vitamins_test2.txt",
                              sep = "\t",
                              header = TRUE,
                              stringsAsFactors = FALSE)

names(vitamin_pl_data) <- c("MAG_ID", "Group", "Vitamins", "PLs")

vitamin_pl_data$Group <- trimws(vitamin_pl_data$Group)

colors_publication <- c("wide" = "#125DA1",
                        "rare" = "#F28A1A")

stats_by_group <- vitamin_pl_data %>%
  group_by(Group) %>%
  summarise(
    N_MAGs = n(),
    Total_PLs = sum(PLs, na.rm = TRUE),
    Total_Vitamins = sum(Vitamins, na.rm = TRUE),
    Mean_PLs = round(mean(PLs, na.rm = TRUE), 2),
    Mean_Vitamins = round(mean(Vitamins, na.rm = TRUE), 2),
    SE_PLs = round(sd(PLs, na.rm = TRUE) / sqrt(n()), 2),
    SE_Vitamins = round(sd(Vitamins, na.rm = TRUE) / sqrt(n()), 2),
    .groups = 'drop'
  ) %>%
  mutate(
    Group_label = case_when(
      Group == "wide" ~ paste0("Widespread\n(n=", N_MAGs, ")"),
      Group == "rare" ~ paste0("Rare\n(n=", N_MAGs, ")"),
      TRUE ~ Group
    )
  )

print("Estadísticas normalizadas por grupo:")
print(stats_by_group)

vitamin_pl_data <- read.delim("pls_vitamins_test2.txt",
                              sep = "\t",
                              header = TRUE,
                              stringsAsFactors = FALSE)

names(vitamin_pl_data) <- c("MAG_ID", "Group", "Vitamins", "PLs")

vitamin_pl_data$Group <- trimws(vitamin_pl_data$Group)

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("=== REALIZANDO ANÁLISIS ESTADÍSTICO ===\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

cat("Estructura de los datos:\n")
print(str(vitamin_pl_data))
cat("Primeras filas:\n")
print(head(vitamin_pl_data))

cat("\nNúmero de MAGs por grupo:\n")
print(table(vitamin_pl_data$Group))

cat("\n--- TEST PARA POLYSACCHARIDE LYASES ---\n")
pl_test <- wilcox.test(PLs ~ Group, data = vitamin_pl_data,
                       alternative = "two.sided", exact = FALSE)

cat("Wilcoxon rank-sum test para PLs:\n")
cat("W =", pl_test$statistic, "\n")
cat("p-value =", format(pl_test$p.value, digits = 4), "\n")
cat("Significativo (p < 0.05):", pl_test$p.value < 0.05, "\n")

cat("\n--- TEST PARA VITAMIN BIOSYNTHESIS ---\n")
vitamin_test <- wilcox.test(Vitamins ~ Group, data = vitamin_pl_data,
                            alternative = "two.sided", exact = FALSE)

cat("Wilcoxon rank-sum test para Vitamins:\n")
cat("W =", vitamin_test$statistic, "\n")
cat("p-value =", format(vitamin_test$p.value, digits = 4), "\n")
cat("Significativo (p < 0.05):", vitamin_test$p.value < 0.05, "\n")

cat("\n--- ESTADÍSTICAS DESCRIPTIVAS ---\n")
descriptive_stats <- vitamin_pl_data %>%
  group_by(Group) %>%
  summarise(
    n = n(),
    PLs_mean = round(mean(PLs, na.rm = TRUE), 2),
    PLs_median = median(PLs, na.rm = TRUE),
    PLs_sd = round(sd(PLs, na.rm = TRUE), 2),
    PLs_se = round(sd(PLs, na.rm = TRUE) / sqrt(n()), 2),
    Vitamins_mean = round(mean(Vitamins, na.rm = TRUE), 2),
    Vitamins_median = median(Vitamins, na.rm = TRUE),
    Vitamins_sd = round(sd(Vitamins, na.rm = TRUE), 2),
    Vitamins_se = round(sd(Vitamins, na.rm = TRUE) / sqrt(n()), 2),
    .groups = 'drop'
  )

print(descriptive_stats)

statistical_summary <- data.frame(
  Variable = c("Polysaccharide Lyases", "Vitamin Biosynthesis"),
  Group_Rare_Mean = c(descriptive_stats$PLs_mean[descriptive_stats$Group == "rare"],
                      descriptive_stats$Vitamins_mean[descriptive_stats$Group == "rare"]),
  Group_Wide_Mean = c(descriptive_stats$PLs_mean[descriptive_stats$Group == "wide"],
                      descriptive_stats$Vitamins_mean[descriptive_stats$Group == "wide"]),
  W_statistic = c(pl_test$statistic, vitamin_test$statistic),
  p_value = c(pl_test$p.value, vitamin_test$p.value),
  Significant = c(pl_test$p.value < 0.05, vitamin_test$p.value < 0.05)
)

cat("\n--- RESUMEN DE TESTS ESTADÍSTICOS ---\n")
print(statistical_summary)

write.csv(statistical_summary, "figures_output/statistical_results.csv", row.names = FALSE)
write.csv(descriptive_stats, "figures_output/descriptive_statistics.csv", row.names = FALSE)

cat("\nArchivos guardados:\n")
cat("- statistical_results.csv\n")
cat("- descriptive_statistics.csv\n")

cat(paste(rep("=", 60), collapse = ""), "\n")
cat("=== ANÁLISIS ESTADÍSTICO COMPLETADO ===\n")
cat(paste(rep("=", 60), collapse = ""), "\n")

panel_A_normalized <- ggplot(stats_by_group, aes(x = Group_label, y = Mean_PLs, fill = Group)) +
  geom_col(width = 0.7, alpha = 0.8, color = "black", size = 0.5) +
  geom_errorbar(aes(ymin = Mean_PLs - SE_PLs, ymax = Mean_PLs + SE_PLs),
                width = 0.25, size = 0.7, color = "black") +
  geom_text(aes(label = Mean_PLs), vjust = -1.2, size = 4, fontface = "bold", color = "black") +
  scale_fill_manual(values = colors_publication, name = "MAG Group") +
  scale_y_continuous(
    name = "Mean polysaccharide lyases per genome",
    breaks = pretty_breaks(n = 6),
    limits = c(0, max(stats_by_group$Mean_PLs + stats_by_group$SE_PLs) * 1.2),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(name = "MAG Group", expand = expansion(add = c(0.5, 0.5))) +
  theme_classic(base_size = font_sizes$panel_A$base) +
  theme(
    axis.text.x = element_text(size = font_sizes$panel_A$axis_text, color = "black", face = "bold"),
    axis.text.y = element_text(size = font_sizes$panel_A$axis_text, color = "black"),
    axis.title = element_text(size = font_sizes$panel_A$axis_title, color = "black", face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 30, b = 15, l = 15, unit = "pt")  )

panel_B_normalized <- ggplot(stats_by_group, aes(x = Group_label, y = Mean_Vitamins, fill = Group)) +
  geom_col(width = 0.7, alpha = 0.8, color = "black", size = 0.5) +
  geom_errorbar(aes(ymin = Mean_Vitamins - SE_Vitamins, ymax = Mean_Vitamins + SE_Vitamins),
                width = 0.25, size = 0.7, color = "black") +
  geom_text(aes(label = Mean_Vitamins), vjust = -1.2, size = 4, fontface = "bold", color = "black") +
  scale_fill_manual(values = colors_publication, name = "MAG Group") +
  scale_y_continuous(
    name = "Mean vitamin biosynthesis genes per genome",
    breaks = pretty_breaks(n = 6),
    limits = c(0, max(stats_by_group$Mean_Vitamins + stats_by_group$SE_Vitamins) * 1.2),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(name = "MAG Group", expand = expansion(add = c(0.5, 0.5))) +
  theme_classic(base_size = font_sizes$panel_B$base) +
  theme(
    axis.text.x = element_text(size = font_sizes$panel_B$axis_text, color = "black", face = "bold"),
    axis.text.y = element_text(size = font_sizes$panel_B$axis_text, color = "black"),
    axis.title = element_text(size = font_sizes$panel_B$axis_title, color = "black", face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 30, b = 15, l = 15, unit = "pt")
    )

vitamin_pl_data_no_zeros_PLs <- vitamin_pl_data %>% filter(PLs > 0)
vitamin_pl_data_no_zeros_Vitamins <- vitamin_pl_data %>% filter(Vitamins > 0)

mags_with_data <- vitamin_pl_data %>%
  filter(PLs > 0 | Vitamins > 0) %>%
  mutate(Total_Combined = PLs + Vitamins) %>%
  arrange(desc(Total_Combined)) %>%
  pull(MAG_ID)

cat("MAGs con datos > 0 para paneles C y D:", length(mags_with_data), "de", nrow(vitamin_pl_data), "totales\n")

mag_mapping <- data.frame(
  MAG_ID = mags_with_data,
  MAG_Number = 1:length(mags_with_data)
)

vitamin_pl_data_numbered_PLs <- vitamin_pl_data_no_zeros_PLs %>%
  left_join(mag_mapping, by = "MAG_ID")

vitamin_pl_data_numbered_Vitamins <- vitamin_pl_data_no_zeros_Vitamins %>%
  left_join(mag_mapping, by = "MAG_ID")

shorten_mag_names <- function(mag_ids, max_length = 6) {
  shortened <- ifelse(
    nchar(mag_ids) > max_length,
    paste0(substr(mag_ids, 1, max_length), "..."),
    mag_ids
  )
  return(shortened)
}

panel_C_clean <- ggplot(vitamin_pl_data_no_zeros_PLs, aes(
  x = PLs,
  y = factor(MAG_ID, levels = rev(mags_with_data)),
  color = Group
)) +
  geom_point(size = 2.5, alpha = 0.7, stroke = 0) +
  scale_color_manual(
    values = colors_publication,
    name = "MAG Group",
    labels = c("Rare", "Widespread")
  ) +
  scale_x_continuous(
    name = "Number of polysaccharide lyases",
    breaks = pretty_breaks(n = 6),
    limits = c(0.5, NA),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_discrete(
    name = "Individual MAGs (ordered by total metabolic activity)",
    expand = expansion(add = c(0.5, 0.5)),
    labels = NULL
  ) +
  theme_classic(base_size = font_sizes$panel_C$base) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = font_sizes$panel_C$axis_text_x, color = "black"),
    axis.title = element_text(size = font_sizes$panel_C$axis_title, color = "black", face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = font_sizes$panel_C$legend, face = "bold"),
    legend.text = element_text(size = font_sizes$panel_C$legend),
    legend.key.size = unit(0.8, "cm"),
    legend.margin = margin(t = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 30, b = 15, l = 15, unit = "pt")  ) +
  guides(
    color = guide_legend(
      title = "MAG Group",
      override.aes = list(size = 3.5, alpha = 1),
      nrow = 1
    )
  )

panel_D_clean <- ggplot(vitamin_pl_data_no_zeros_Vitamins, aes(
  x = Vitamins,
  y = factor(MAG_ID, levels = rev(mags_with_data)),
  color = Group
)) +
  geom_point(size = 2.5, alpha = 0.7, stroke = 0) +
  scale_color_manual(
    values = colors_publication,
    name = "MAG Group",
    labels = c("Rare", "Widespread")
  ) +
  scale_x_continuous(
    name = "Number of vitamin biosynthesis genes",
    breaks = pretty_breaks(n = 6),
    limits = c(0.5, NA),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_discrete(
    name = "Individual MAGs (ordered by total metabolic activity)",
    expand = expansion(add = c(0.5, 0.5)),
    labels = NULL
  ) +
  theme_classic(base_size = font_sizes$panel_D$base) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = font_sizes$panel_D$axis_text_x, color = "black"),
    axis.title = element_text(size = font_sizes$panel_D$axis_title, color = "black", face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = font_sizes$panel_D$legend, face = "bold"),
    legend.text = element_text(size = font_sizes$panel_D$legend),
    legend.key.size = unit(0.8, "cm"),
    legend.margin = margin(t = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 30, b = 15, l = 15, unit = "pt")  ) +
  guides(
    color = guide_legend(
      title = "MAG Group",
      override.aes = list(size = 3.5, alpha = 1),
      nrow = 1
    )
  )

panel_C_numbered <- ggplot(vitamin_pl_data_numbered_PLs, aes(
  x = PLs,
  y = MAG_Number,
  color = Group
)) +
  geom_point(size = 2.5, alpha = 0.7, stroke = 0) +
  scale_color_manual(
    values = colors_publication,
    name = "MAG Group",
    labels = c("Rare", "Widespread")
  ) +
  scale_x_continuous(
    name = "Number of polysaccharide lyases",
    breaks = pretty_breaks(n = 6),
    limits = c(0.5, NA),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name = "MAG rank (by total metabolic activity)",
    breaks = pretty_breaks(n = 8),
    expand = expansion(add = c(0.5, 0.5))
  ) +
  theme_classic(base_size = font_sizes$panel_C$base) +
  theme(
    axis.text.y = element_text(size = font_sizes$panel_C$axis_text_x, color = "black"),
    axis.text.x = element_text(size = font_sizes$panel_C$axis_text_x, color = "black"),
    axis.title = element_text(size = font_sizes$panel_C$axis_title, color = "black", face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = font_sizes$panel_C$legend, face = "bold"),
    legend.text = element_text(size = font_sizes$panel_C$legend),
    legend.key.size = unit(0.8, "cm"),
    legend.margin = margin(t = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 30, b = 15, l = 15, unit = "pt")
    ) +
  guides(
    color = guide_legend(
      title = "MAG Group",
      override.aes = list(size = 3.5, alpha = 1),
      nrow = 1
    )
  )

panel_D_numbered <- ggplot(vitamin_pl_data_numbered_Vitamins, aes(
  x = Vitamins,
  y = MAG_Number,
  color = Group
)) +
  geom_point(size = 2.5, alpha = 0.7, stroke = 0) +
  scale_color_manual(
    values = colors_publication,
    name = "MAG Group",
    labels = c("Rare", "Widespread")
  ) +
  scale_x_continuous(
    name = "Number of vitamin biosynthesis genes",
    breaks = pretty_breaks(n = 6),
    limits = c(0.5, NA),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_continuous(
    name = "MAG rank (by total metabolic activity)",
    breaks = pretty_breaks(n = 8),
    expand = expansion(add = c(0.5, 0.5))
  ) +
  theme_classic(base_size = font_sizes$panel_D$base) +
  theme(
    axis.text.y = element_text(size = font_sizes$panel_D$axis_text_x, color = "black"),
    axis.text.x = element_text(size = font_sizes$panel_D$axis_text_x, color = "black"),
    axis.title = element_text(size = font_sizes$panel_D$axis_title, color = "black", face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = font_sizes$panel_D$legend, face = "bold"),
    legend.text = element_text(size = font_sizes$panel_D$legend),
    legend.key.size = unit(0.8, "cm"),
    legend.margin = margin(t = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 30, b = 15, l = 15, unit = "pt")
    ) +
  guides(
    color = guide_legend(
      title = "MAG Group",
      override.aes = list(size = 3.5, alpha = 1),
      nrow = 1
    )
  )

panel_C_shortened <- ggplot(vitamin_pl_data_no_zeros_PLs, aes(
  x = PLs,
  y = factor(MAG_ID, levels = rev(mags_with_data)),
  color = Group
)) +
  geom_point(size = 2.5, alpha = 0.7, stroke = 0) +
  scale_color_manual(
    values = colors_publication,
    name = "MAG Group",
    labels = c("Rare", "Widespread")
  ) +
  scale_x_continuous(
    name = "Number of polysaccharide lyases",
    breaks = pretty_breaks(n = 6),
    limits = c(0.5, NA),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_discrete(
    name = "MAG ID",
    expand = expansion(add = c(0.5, 0.5)),
    labels = shorten_mag_names(rev(mags_with_data), max_length = 6)
  ) +
  theme_classic(base_size = font_sizes$panel_C$base) +
  theme(
    axis.text.y = element_text(size = 6, color = "black", family = "mono"),
    axis.text.x = element_text(size = font_sizes$panel_C$axis_text_x, color = "black"),
    axis.title = element_text(size = font_sizes$panel_C$axis_title, color = "black", face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = font_sizes$panel_C$legend, face = "bold"),
    legend.text = element_text(size = font_sizes$panel_C$legend),
    legend.key.size = unit(0.8, "cm"),
    legend.margin = margin(t = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 30, b = 15, l = 15, unit = "pt")  ) +
  guides(
    color = guide_legend(
      title = "MAG Group",
      override.aes = list(size = 3.5, alpha = 1),
      nrow = 1
    )
  )

panel_D_shortened <- ggplot(vitamin_pl_data_no_zeros_Vitamins, aes(
  x = Vitamins,
  y = factor(MAG_ID, levels = rev(mags_with_data)),
  color = Group
)) +
  geom_point(size = 2.5, alpha = 0.7, stroke = 0) +
  scale_color_manual(
    values = colors_publication,
    name = "MAG Group",
    labels = c("Rare", "Widespread")
  ) +
  scale_x_continuous(
    name = "Number of vitamin biosynthesis genes",
    breaks = pretty_breaks(n = 6),
    limits = c(0.5, NA),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_y_discrete(
    name = "MAG ID",
    expand = expansion(add = c(0.5, 0.5)),
    labels = shorten_mag_names(rev(mags_with_data), max_length = 6)
  ) +
  theme_classic(base_size = font_sizes$panel_D$base) +
  theme(
    axis.text.y = element_text(size = 6, color = "black", family = "mono"),
    axis.text.x = element_text(size = font_sizes$panel_D$axis_text_x, color = "black"),
    axis.title = element_text(size = font_sizes$panel_D$axis_title, color = "black", face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black", size = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = font_sizes$panel_D$legend, face = "bold"),
    legend.text = element_text(size = font_sizes$panel_D$legend),
    legend.key.size = unit(0.8, "cm"),
    legend.margin = margin(t = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 30, b = 15, l = 15, unit = "pt")  ) +
  guides(
    color = guide_legend(
      title = "MAG Group",
      override.aes = list(size = 3.5, alpha = 1),
      nrow = 1
    )
  )

panel_C_horizontal <- ggplot(vitamin_pl_data_no_zeros_PLs, aes(
  y = PLs,
  x = factor(MAG_ID, levels = mags_with_data),
  color = Group
)) +
  geom_point(size = 2.5, alpha = 0.7, stroke = 0) +
  scale_color_manual(
    values = colors_publication,
    name = "MAG Group",
    labels = c("Rare", "Widespread")
  ) +
  scale_y_continuous(
    name = "Number of polysaccharide lyases",
    breaks = pretty_breaks(n = 6),
    limits = c(0.5, NA),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_discrete(
    name = "Individual MAGs (ordered by total metabolic activity)",
    expand = expansion(add = c(0.5, 0.5)),
    labels = NULL
  ) +
  theme_classic(base_size = font_sizes$panel_C$base) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = font_sizes$panel_C$axis_text_x, color = "black"),
    axis.title = element_text(size = font_sizes$panel_C$axis_title, color = "black", face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = font_sizes$panel_C$legend, face = "bold"),
    legend.text = element_text(size = font_sizes$panel_C$legend),
    legend.key.size = unit(0.8, "cm"),
    legend.margin = margin(t = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 30, b = 15, l = 15, unit = "pt")  ) +
  guides(
    color = guide_legend(
      title = "MAG Group",
      override.aes = list(size = 3.5, alpha = 1),
      nrow = 1
    )
  )

panel_D_horizontal <- ggplot(vitamin_pl_data_no_zeros_Vitamins, aes(
  y = Vitamins,
  x = factor(MAG_ID, levels = mags_with_data),
  color = Group
)) +
  geom_point(size = 2.5, alpha = 0.7, stroke = 0) +
  scale_color_manual(
    values = colors_publication,
    name = "MAG Group",
    labels = c("Rare", "Widespread")
  ) +
  scale_y_continuous(
    name = "Number of vitamin biosynthesis genes",
    breaks = pretty_breaks(n = 6),
    limits = c(0.5, NA),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_x_discrete(
    name = "Individual MAGs (ordered by total metabolic activity)",
    expand = expansion(add = c(0.5, 0.5)),
    labels = NULL
  ) +
  theme_classic(base_size = font_sizes$panel_D$base) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = font_sizes$panel_D$axis_text_x, color = "black"),
    axis.title = element_text(size = font_sizes$panel_D$axis_title, color = "black", face = "bold"),
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_line(color = "black", size = 0.5),
    legend.position = "bottom",
    legend.title = element_text(size = font_sizes$panel_D$legend, face = "bold"),
    legend.text = element_text(size = font_sizes$panel_D$legend),
    legend.key.size = unit(0.8, "cm"),
    legend.margin = margin(t = 10),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 30, b = 15, l = 15, unit = "pt")  ) +
  guides(
    color = guide_legend(
      title = "MAG Group",
      override.aes = list(size = 3.5, alpha = 1),
      nrow = 1
    )
  )

if (!dir.exists("figures_output")) {
  dir.create("figures_output")
}

combined_figure_v1_clean <- plot_grid(
  panel_A_normalized, panel_B_normalized,
  panel_C_clean, panel_D_clean,
  labels = c("A", "B", "C", "D"),
  label_size = 14,
  label_fontface = "bold",
  ncol = 2,
  nrow = 2,
  rel_widths = c(1, 1),
  align = "hv",
  axis = "tb"
)

combined_figure_v2_numbered <- plot_grid(
  panel_A_normalized, panel_B_normalized,
  panel_C_numbered, panel_D_numbered,
  labels = c("A", "B", "C", "D"),
  label_size = 14,
  label_fontface = "bold",
  ncol = 2,
  nrow = 2,
  rel_widths = c(1, 1),
  align = "hv",
  axis = "tb"
)

combined_figure_v3_shortened <- plot_grid(
  panel_A_normalized, panel_B_normalized,
  panel_C_shortened, panel_D_shortened,
  labels = c("A", "B", "C", "D"),
  label_size = 14,
  label_fontface = "bold",
  ncol = 2,
  nrow = 2,
  rel_widths = c(1, 1),
  align = "hv",
  axis = "tb"
)

combined_figure_v4_horizontal <- plot_grid(
  panel_A_normalized, panel_B_normalized,
  panel_C_horizontal, panel_D_horizontal,
  labels = c("A", "B", "C", "D"),
  label_size = 14,
  label_fontface = "bold",
  ncol = 2,
  nrow = 2,
  rel_widths = c(1, 1),
  align = "hv",
  axis = "tb"
)

combined_CD_v1_clean <- plot_grid(
  panel_C_clean, panel_D_clean,
  labels = c("C", "D"),
  label_size = 14,
  label_fontface = "bold",
  ncol = 2,
  rel_widths = c(1, 1),
  align = "h",
  axis = "tb"
)

combined_CD_v2_numbered <- plot_grid(
  panel_C_numbered, panel_D_numbered,
  labels = c("C", "D"),
  label_size = 14,
  label_fontface = "bold",
  ncol = 2,
  rel_widths = c(1, 1),
  align = "h",
  axis = "tb"
)

combined_CD_v3_shortened <- plot_grid(
  panel_C_shortened, panel_D_shortened,
  labels = c("C", "D"),
  label_size = 14,
  label_fontface = "bold",
  ncol = 2,
  rel_widths = c(1, 1),
  align = "h",
  axis = "tb"
)

combined_CD_v4_horizontal <- plot_grid(
  panel_C_horizontal, panel_D_horizontal,
  labels = c("C", "D"),
  label_size = 14,
  label_fontface = "bold",
  ncol = 2,
  rel_widths = c(1, 1),
  align = "h",
  axis = "tb"
)

combined_figure_AB <- plot_grid(
  panel_A_normalized, panel_B_normalized,
  labels = c("A", "B"),
  label_size = 14,
  label_fontface = "bold",
  ncol = 2,
  rel_widths = c(1, 1),
  align = "h",
  axis = "tb"
)

cat("=== GUARDANDO PANELES INDIVIDUALES (A y B) ===\n")

ggsave(
  filename = "figures_output/Figure_PanelA_Mean_PLs_by_Group.pdf",
  plot = panel_A_normalized,
  width = 6, height = 5, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_PanelB_Mean_Vitamins_by_Group.pdf",
  plot = panel_B_normalized,
  width = 6, height = 5, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_Combined_AB_Normalized.pdf",
  plot = combined_figure_AB,
  width = 12, height = 5, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

cat("=== GUARDANDO VARIACIÓN 1: SIN NOMBRES (RECOMENDADA) ===\n")

ggsave(
  filename = "figures_output/Figure_V1_PanelC_PLs_Clean.pdf",
  plot = panel_C_clean,
  width = 6, height = 8, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_V1_PanelD_Vitamins_Clean.pdf",
  plot = panel_D_clean,
  width = 6, height = 8, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_V1_Combined_CD_Clean.pdf",
  plot = combined_CD_v1_clean,
  width = 12, height = 8, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_V1_Complete_ABCD_Clean.pdf",
  plot = combined_figure_v1_clean,
  width = 12, height = 12, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

cat("=== GUARDANDO VARIACIÓN 2: CON NÚMEROS DE RANKING ===\n")

ggsave(
  filename = "figures_output/Figure_V2_PanelC_PLs_Numbered.pdf",
  plot = panel_C_numbered,
  width = 6, height = 8, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_V2_PanelD_Vitamins_Numbered.pdf",
  plot = panel_D_numbered,
  width = 6, height = 8, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_V2_Combined_CD_Numbered.pdf",
  plot = combined_CD_v2_numbered,
  width = 12, height = 10, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_V2_Complete_ABCD_Numbered.pdf",
  plot = combined_figure_v2_numbered,
  width = 12, height = 10, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

cat("=== GUARDANDO VARIACIÓN 3: CON NOMBRES ACORTADOS ===\n")

ggsave(
  filename = "figures_output/Figure_V3_PanelC_PLs_Shortened.pdf",
  plot = panel_C_shortened,
  width = 6, height = 8, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_V3_PanelD_Vitamins_Shortened.pdf",
  plot = panel_D_shortened,
  width = 6, height = 8, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_V3_Combined_CD_Shortened.pdf",
  plot = combined_CD_v3_shortened,
  width = 12, height = 8, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_V3_Complete_ABCD_Shortened.pdf",
  plot = combined_figure_v3_shortened,
  width = 12, height = 12, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

cat("=== GUARDANDO VARIACIÓN 4: ORIENTACIÓN HORIZONTAL ===\n")

ggsave(
  filename = "figures_output/Figure_V4_PanelC_PLs_Horizontal.pdf",
  plot = panel_C_horizontal,
  width = 8, height = 6, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_V4_PanelD_Vitamins_Horizontal.pdf",
  plot = panel_D_horizontal,
  width = 8, height = 6, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_V4_Combined_CD_Horizontal.pdf",
  plot = combined_CD_v4_horizontal,
  width = 12, height = 8, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

ggsave(
  filename = "figures_output/Figure_V4_Complete_ABCD_Horizontal.pdf",
  plot = combined_figure_v4_horizontal,
  width = 12, height = 12, units = "in", dpi = 300, device = "pdf", useDingbats = FALSE
)

cat("=== GUARDANDO VERSIONES TIFF PARA LAS MEJORES OPCIONES ===\n")

ggsave(
  filename = "figures_output/Figure_V1_Combined_AB.tiff",
  plot = combined_figure_AB,
  width = 10, height = 4, units = "in", dpi = 300, compression = "lzw"
)

ggsave(
  filename = "figures_output/Figure_V1_Combined_CD_Clean.tiff",
  plot = combined_CD_v1_clean,
  width = 10, height = 6, units = "in", dpi = 300, compression = "lzw"
)

ggsave(
  filename = "figures_output/Figure_V1_Complete_ABCD_Clean.tiff",
  plot = combined_figure_v1_clean,
  width = 10, height = 10, units = "in", dpi = 300, compression = "lzw"
)

ggsave(
  filename = "figures_output/Figure_V2_Complete_ABCD_Numbered.tiff",
  plot = combined_figure_v2_numbered,
  width = 10, height = 10, units = "in", dpi = 300, compression = "lzw"
)

cat("\n=== MOSTRANDO FIGURAS PARA COMPARACIÓN ===\n")

print("VARIACIÓN 1 - Sin nombres (MÁS RECOMENDADA):")
print(combined_figure_v1_clean)

print("VARIACIÓN 2 - Con números de ranking:")
print(combined_figure_v2_numbered)

print("VARIACIÓN 3 - Nombres acortados:")
print(combined_figure_v3_shortened)

print("VARIACIÓN 4 - Orientación horizontal:")
print(combined_figure_v4_horizontal)

additional_stats <- vitamin_pl_data %>%
  group_by(Group) %>%
  summarise(
    N_MAGs = n(),
    Mean_PLs_per_MAG = round(mean(PLs, na.rm = TRUE), 2),
    Median_PLs_per_MAG = median(PLs, na.rm = TRUE),
    Mean_Vitamins_per_MAG = round(mean(Vitamins, na.rm = TRUE), 2),
    Median_Vitamins_per_MAG = median(Vitamins, na.rm = TRUE),
    Total_PLs = sum(PLs, na.rm = TRUE),
    Total_Vitamins = sum(Vitamins, na.rm = TRUE),
    .groups = 'drop'
  )

print("\nEstadísticas detalladas por grupo:")
print(additional_stats)

total_pls_all <- sum(vitamin_pl_data$PLs, na.rm = TRUE)
total_vitamins_all <- sum(vitamin_pl_data$Vitamins, na.rm = TRUE)

proportions <- stats_by_group %>%
  mutate(
    Prop_PLs = round(Total_PLs / total_pls_all * 100, 1),
    Prop_Vitamins = round(Total_Vitamins / total_vitamins_all * 100, 1)
  )

print("\nProporciones de totales por grupo (%):")
print(proportions[c("Group_label", "Prop_PLs", "Prop_Vitamins")])

cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("=== RESUMEN COMPLETO DE FIGURAS GENERADAS ===\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

cat("\nPANELES INDIVIDUALES (A y B - iguales para todas las variaciones):\n")
cat("- Figure_PanelA_Mean_PLs_by_Group.pdf (6x5 pulgadas)\n")
cat("- Figure_PanelB_Mean_Vitamins_by_Group.pdf (6x5 pulgadas)\n")
cat("- Figure_Combined_AB_Normalized.pdf (12x5 pulgadas)\n")

cat("\nVARIACIÓN 1 - SIN NOMBRES (RECOMENDADA):\n")
cat("INDIVIDUAL:\n")
cat("- Figure_V1_PanelC_PLs_Clean.pdf (6x8 pulgadas)\n")
cat("- Figure_V1_PanelD_Vitamins_Clean.pdf (6x8 pulgadas)\n")
cat("COMBINADAS:\n")
cat("- Figure_V1_Combined_CD_Clean.pdf (12x8 pulgadas)\n")
cat("- Figure_V1_Complete_ABCD_Clean.pdf (12x12 pulgadas)\n")

cat("\nVARIACIÓN 2 - CON NÚMEROS DE RANKING:\n")
cat("INDIVIDUAL:\n")
cat("- Figure_V2_PanelC_PLs_Numbered.pdf (6x8 pulgadas)\n")
cat("- Figure_V2_PanelD_Vitamins_Numbered.pdf (6x8 pulgadas)\n")
cat("COMBINADAS:\n")
cat("- Figure_V2_Combined_CD_Numbered.pdf (12x8 pulgadas)\n")
cat("- Figure_V2_Complete_ABCD_Numbered.pdf (12x12 pulgadas)\n")

cat("\nVARIACIÓN 3 - CON NOMBRES ACORTADOS:\n")
cat("INDIVIDUAL:\n")
cat("- Figure_V3_PanelC_PLs_Shortened.pdf (6x8 pulgadas)\n")
cat("- Figure_V3_PanelD_Vitamins_Shortened.pdf (6x8 pulgadas)\n")
cat("COMBINADAS:\n")
cat("- Figure_V3_Combined_CD_Shortened.pdf (12x8 pulgadas)\n")
cat("- Figure_V3_Complete_ABCD_Shortened.pdf (12x12 pulgadas)\n")

cat("\nVARIACIÓN 4 - ORIENTACIÓN HORIZONTAL:\n")
cat("INDIVIDUAL:\n")
cat("- Figure_V4_PanelC_PLs_Horizontal.pdf (8x6 pulgadas)\n")
cat("- Figure_V4_PanelD_Vitamins_Horizontal.pdf (8x6 pulgadas)\n")
cat("COMBINADAS:\n")
cat("- Figure_V4_Combined_CD_Horizontal.pdf (12x8 pulgadas)\n")
cat("- Figure_V4_Complete_ABCD_Horizontal.pdf (12x12 pulgadas)\n")

cat("\nVERSIONES TIFF PARA ENVÍO (mejores opciones):\n")
cat("- Figure_V1_Combined_AB.tiff (10x4 pulgadas)\n")
cat("- Figure_V1_Combined_CD_Clean.tiff (10x6 pulgadas)\n")
cat("- Figure_V1_Complete_ABCD_Clean.tiff (10x10 pulgadas)\n")
cat("- Figure_V2_Complete_ABCD_Numbered.tiff (10x10 pulgadas)\n")

cat("\n" + "="*80 + "\n")
cat("RECOMENDACIONES:\n")
cat("="*80 + "\n")
cat("1. VARIACIÓN 1 (Sin nombres): Más limpia y profesional\n")
cat("2. VARIACIÓN 2 (Con números): Buena si necesitas referenciar MAGs específicos\n")
cat("3. VARIACIÓN 3 (Nombres acortados): Solo si es esencial mostrar IDs\n")
cat("4. VARIACIÓN 4 (Horizontal): Diferente pero puede ser menos legible\n")
cat("\nTodas las figuras están optimizadas para estándares de publicación.\n")
cat("Archivos guardados en el directorio 'figures_output/'\n")

cat("\nTotal de MAGs analizados:", nrow(vitamin_pl_data), "\n")
cat("MAGs con PLs > 0:", nrow(vitamin_pl_data_no_zeros_PLs), "\n")
cat("MAGs con Vitaminas > 0:", nrow(vitamin_pl_data_no_zeros_Vitamins), "\n")
cat("MAGs con datos > 0 (PLs O Vitaminas):", length(mags_with_data), "\n")
