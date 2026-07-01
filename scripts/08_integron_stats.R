#!/usr/bin/env Rscript
# ==============================================================
# Script: 08_integron_stats.R
# Purpose: Deep Integron Visualization & Statistical Analysis
#
# Input:
#   --class_dist   : {pop}/class_distribution.csv
#   --cassette     : {pop}/cassette_matrix.csv
#   --specificity  : {pop}/region_specificity.csv
#   --amr_link     : {pop}/amr_integron_link.csv
#   --output_dir   : {pop}/ (folder output)
#   --pop          : nama populasi
#
# Output:
#   figures/Fig_Int1_class_distribution.pdf  - Stacked bar: integron class per country
#   figures/Fig_Int2_cassette_heatmap.pdf    - Heatmap top cassettes x country
#   figures/Fig_Int3_amr_drugclass.pdf       - Bar: AMR drug class di integron
#   tables/integron_fisher_test.csv          - Fisher's Exact Test: class distribution
#   tables/top_cassettes.csv                 - Top 20 gene cassettes by frequency
# ==============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(RColorBrewer)
  library(scales)
})

# --- Parse arguments ---
args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) > 0 && idx + 1 <= length(args)) return(args[idx + 1])
  return(default)
}

class_dist_file  <- parse_arg("--class_dist")
cassette_file    <- parse_arg("--cassette")
specificity_file <- parse_arg("--specificity")
amr_link_file    <- parse_arg("--amr_link")
output_dir       <- parse_arg("--output_dir")
pop_name         <- parse_arg("--pop", "population")

# --- Setup output directories ---
fig_dir   <- file.path(output_dir, "figures")
table_dir <- file.path(output_dir, "tables")
dir.create(fig_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

cat("========================================\n")
cat(sprintf("[08_integron_stats.R] Population: %s\n", pop_name))
cat("========================================\n")

# ============================================================
# PLOT 1: Integron Class Distribution per Country (Stacked Bar)
# ============================================================
cat("\n[1] Membuat Stacked Bar — Integron Class Distribution...\n")

tryCatch({
  class_df <- read_csv(class_dist_file, show_col_types = FALSE)

  if (nrow(class_df) > 0 && "integron_class" %in% colnames(class_df)) {
    # Pilih kolom yang relevan
    plot_df <- class_df %>%
      filter(!is.na(integron_class), integron_class != "") %>%
      group_by(country_name, region, integron_class) %>%
      summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
      group_by(country_name) %>%
      mutate(proportion = count / sum(count)) %>%
      ungroup()

    # Warna untuk class
    class_colors <- c(
      "complete" = "#2166AC",
      "calin"    = "#F4A582",
      "in0"      = "#BABABA",
      "unknown"  = "#D9D9D9"
    )

    p1 <- ggplot(plot_df, aes(x = country_name, y = proportion, fill = integron_class)) +
      geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
      scale_fill_manual(
        values = class_colors,
        name   = "Integron Class",
        labels = c(
          "complete" = "Complete (Class 1/2)",
          "calin"    = "CALIN (Historical HGT)",
          "in0"      = "In0",
          "unknown"  = "Unknown"
        )
      ) +
      scale_y_continuous(labels = percent_format()) +
      labs(
        title    = sprintf("Integron Class Distribution — %s", pop_name),
        subtitle = "Proportion of integron types per country (IntegronFinder)",
        x        = "Country",
        y        = "Proportion (%)"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title       = element_text(face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(hjust = 0.5, color = "grey50"),
        axis.text.x      = element_text(angle = 30, hjust = 1),
        legend.position  = "right",
        panel.grid.minor = element_blank()
      )

    ggsave(
      filename = file.path(fig_dir, "Fig_Int1_class_distribution.pdf"),
      plot     = p1,
      width    = 8, height = 5, device = "pdf"
    )
    cat(sprintf("   -> Tersimpan: Fig_Int1_class_distribution.pdf\n"))
  } else {
    cat("   [SKIP] Data class distribution kosong atau kolom tidak lengkap.\n")
  }
}, error = function(e) {
  cat(sprintf("   [WARN] Gagal membuat Plot 1: %s\n", conditionMessage(e)))
})


# ============================================================
# PLOT 2: Gene Cassette Heatmap (Top 20 cassettes x Country)
# ============================================================
cat("\n[2] Membuat Gene Cassette Heatmap...\n")

tryCatch({
  cassette_df <- read_csv(cassette_file, show_col_types = FALSE)

  if (nrow(cassette_df) > 0 && "annotation" %in% colnames(cassette_df)) {
    # Pilih top 20 cassette berdasarkan total kemunculan
    top_cassettes <- cassette_df %>%
      group_by(annotation) %>%
      summarise(total = sum(count, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(total)) %>%
      slice_head(n = 20) %>%
      pull(annotation)

    heat_df <- cassette_df %>%
      filter(annotation %in% top_cassettes) %>%
      group_by(annotation, country_name) %>%
      summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = country_name, values_from = count, values_fill = 0) %>%
      pivot_longer(-annotation, names_to = "country_name", values_to = "count")

    # Simpan tabel top cassettes
    top_cassette_table <- cassette_df %>%
      filter(annotation %in% top_cassettes) %>%
      group_by(annotation, region) %>%
      summarise(
        total_count = sum(count, na.rm = TRUE),
        n_samples   = sum(n_samples, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(total_count))

    write_csv(top_cassette_table, file.path(table_dir, "top_cassettes.csv"))
    cat(sprintf("   -> Top cassettes table tersimpan: tables/top_cassettes.csv\n"))

    p2 <- ggplot(heat_df, aes(x = country_name, y = annotation, fill = count)) +
      geom_tile(color = "white", linewidth = 0.4) +
      scale_fill_gradient(
        low  = "#F7FBFF",
        high = "#08306B",
        name = "Count"
      ) +
      labs(
        title    = sprintf("Top Gene Cassette Distribution — %s", pop_name),
        subtitle = "Top 20 gene cassettes from complete integrons per country",
        x        = "Country",
        y        = "Gene Cassette"
      ) +
      theme_bw(base_size = 11) +
      theme(
        plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
        axis.text.x   = element_text(angle = 30, hjust = 1),
        axis.text.y   = element_text(size = 8),
        legend.position = "right"
      )

    ggsave(
      filename = file.path(fig_dir, "Fig_Int2_cassette_heatmap.pdf"),
      plot     = p2,
      width    = 9, height = 7, device = "pdf"
    )
    cat(sprintf("   -> Tersimpan: Fig_Int2_cassette_heatmap.pdf\n"))
  } else {
    cat("   [SKIP] Data cassette matrix kosong.\n")
    write_csv(data.frame(), file.path(table_dir, "top_cassettes.csv"))
  }
}, error = function(e) {
  cat(sprintf("   [WARN] Gagal membuat Plot 2: %s\n", conditionMessage(e)))
  write_csv(data.frame(), file.path(table_dir, "top_cassettes.csv"))
})


# ============================================================
# PLOT 3: AMR Drug Class on Integrons (Bar Chart)
# ============================================================
cat("\n[3] Membuat Bar Chart — AMR Drug Class di Integrons...\n")

tryCatch({
  amr_df <- read_csv(amr_link_file, show_col_types = FALSE)

  if (nrow(amr_df) > 0 && "Drug Class" %in% colnames(amr_df)) {
    region_col <- if ("Region" %in% colnames(amr_df)) "Region" else "Country"

    amr_summary <- amr_df %>%
      rename(drug_class = `Drug Class`) %>%
      group_by(drug_class, !!sym(region_col)) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(drug_class) %>%
      mutate(total = sum(count)) %>%
      ungroup() %>%
      arrange(desc(total))

    # Top 15 drug classes
    top_classes <- amr_summary %>%
      distinct(drug_class, total) %>%
      slice_max(total, n = 15) %>%
      pull(drug_class)

    plot_amr <- amr_summary %>%
      filter(drug_class %in% top_classes) %>%
      mutate(drug_class = factor(drug_class, levels = rev(top_classes)))

    # Color palette per region
    n_regions <- length(unique(plot_amr[[region_col]]))
    region_colors <- brewer.pal(max(3, min(n_regions, 8)), "Set2")[seq_len(n_regions)]
    names(region_colors) <- unique(plot_amr[[region_col]])

    p3 <- ggplot(plot_amr, aes(x = drug_class, y = count, fill = !!sym(region_col))) +
      geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
      coord_flip() +
      scale_fill_manual(values = region_colors, name = "Region") +
      labs(
        title    = sprintf("AMR Drug Classes Carried by Integrons — %s", pop_name),
        subtitle = "Number of AMR genes co-localized with complete integrons",
        x        = "Drug Class",
        y        = "Count"
      ) +
      theme_bw(base_size = 12) +
      theme(
        plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey50"),
        legend.position = "right",
        panel.grid.minor = element_blank()
      )

    ggsave(
      filename = file.path(fig_dir, "Fig_Int3_amr_drugclass.pdf"),
      plot     = p3,
      width    = 9, height = 6, device = "pdf"
    )
    cat(sprintf("   -> Tersimpan: Fig_Int3_amr_drugclass.pdf\n"))
  } else {
    cat("   [SKIP] Data AMR-integron link kosong.\n")
  }
}, error = function(e) {
  cat(sprintf("   [WARN] Gagal membuat Plot 3: %s\n", conditionMessage(e)))
})


# ============================================================
# STATISTIK: Fisher's Exact Test — Class Distribution
# ============================================================
cat("\n[4] Menjalankan Fisher's Exact Test — Integron Class Distribution...\n")

tryCatch({
  class_df <- read_csv(class_dist_file, show_col_types = FALSE)

  if (nrow(class_df) > 0 && nrow(class_df %>% filter(integron_class == "complete")) > 0) {
    # Buat contingency table: complete vs bukan complete per country
    contingency <- class_df %>%
      mutate(is_complete = ifelse(integron_class == "complete", "complete", "non_complete")) %>%
      group_by(country_name, is_complete) %>%
      summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = is_complete, values_from = count, values_fill = 0)

    if (nrow(contingency) >= 2 && all(c("complete", "non_complete") %in% colnames(contingency))) {
      mat <- as.matrix(contingency[, c("complete", "non_complete")])
      rownames(mat) <- contingency$country_name

      fisher_result <- fisher.test(mat, simulate.p.value = TRUE)

      fisher_df <- data.frame(
        test       = "Fisher's Exact Test (simulated)",
        comparison = "complete vs non-complete integron across countries",
        p_value    = fisher_result$p.value,
        note       = ifelse(fisher_result$p.value < 0.05,
                            "Significant: integron class distribution differs across countries",
                            "Not significant: similar integron class distribution across countries")
      )

      write_csv(fisher_df, file.path(table_dir, "integron_fisher_test.csv"))
      cat(sprintf("   -> Fisher p-value: %.4f\n", fisher_result$p.value))
      cat(sprintf("   -> Tersimpan: tables/integron_fisher_test.csv\n"))
    } else {
      cat("   [SKIP] Tabel tidak memiliki cukup kolom untuk Fisher test.\n")
      write_csv(data.frame(), file.path(table_dir, "integron_fisher_test.csv"))
    }
  } else {
    cat("   [SKIP] Data tidak cukup untuk Fisher test.\n")
    write_csv(data.frame(), file.path(table_dir, "integron_fisher_test.csv"))
  }
}, error = function(e) {
  cat(sprintf("   [WARN] Fisher test gagal: %s\n", conditionMessage(e)))
  write_csv(data.frame(), file.path(table_dir, "integron_fisher_test.csv"))
})


cat("\n========================================\n")
cat("[SELESAI] 08_integron_stats.R complete.\n")
cat(sprintf("Output: %s\n", output_dir))
cat("========================================\n")
