#!/usr/bin/env Rscript
# Script: 05_network_analysis.R
# Purpose: Cross-population AMR-MGE network comparison.
#          Membangun bipartite network per populasi (node = AMR gene | MGE type,
#          edge = co-occurrence) dan membandingkan strukturnya antar DNK, CHN, IND.
#
# Ini adalah NOVELTY FIGURE utama paper:
#   "Apakah ada hub AMR gene yang universal, atau ada yang region-specific?"
#
# Input:
#   - results/colocalization_summary.csv
#
# Output:
#   - results/figures/Fig3_Network_DNK.pdf
#   - results/figures/Fig3_Network_CHN.pdf
#   - results/figures/Fig3_Network_IND.pdf
#   - results/figures/Fig3_Network_Combined.pdf
#   - results/tables/Table_Network_Metrics.csv
#   - results/tables/Table_Jaccard_Similarity.csv

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(igraph)
  library(ggraph)
  library(RColorBrewer)
})

cat("============================================================\n")
cat("[Phase 5] Cross-Population AMR-MGE Network Comparison\n")
cat("============================================================\n\n")

# --- Base R Argument Parser ---
args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  params <- list()
  i <- 1
  while(i <= length(args)) {
    if(startsWith(args[i], "--")) {
      key <- sub("^--", "", args[i])
      val <- args[i+1]
      params[[key]] <- val
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  return(params)
}
params <- parse_args(args)

# Fallback defaults jika dijalankan manual tanpa arguments
input_file      <- ifelse(!is.null(params$input), params$input, "results/amr_mge_association_matrix.csv")
output_dir      <- ifelse(!is.null(params$output_dir), params$output_dir, "results")
pop_name        <- ifelse(!is.null(params$pop), params$pop, "combined")
MIN_COOCCURRENCE <- as.numeric(ifelse(!is.null(params$min_cooccurrence), params$min_cooccurrence, "2"))
LAYOUT_METHOD   <- ifelse(!is.null(params$layout), params$layout, "stress")

# --- Setup direktori output ---
fig_dir <- file.path(output_dir, "figures")
tab_dir <- file.path(output_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load data ---
cat("[1] Membaca data input...\n")
agg_dir <- dirname(input_file)

# File name berdasarkan prefix populasi
coloc_path <- file.path(agg_dir, paste0(pop_name, "_all_coloc.csv"))

if (!file.exists(coloc_path)) {
  cat("[WARN] File input tidak ditemukan. Mungkin populasi ini kosong.\n")
  cat("       Mencoba membaca fallback dari folder results/ jika ada...\n")
  coloc_path <- "results/colocalization_summary.csv"
}

# Membaca data secara aman
df <- tryCatch(read.csv(coloc_path), error = function(e) data.frame())

if (nrow(df) == 0) {
  cat("[WARN] Data kosong. Menghentikan analisis secara aman.\n")
  quit(status = 0)
}

# Hanya ambil gen AMR yang berasosiasi dengan MGE (exclude Chromosomal)
df_mobile <- df %>%
  filter(MGE_Type != "Chromosomal") %>%
  rename(Drug_Class = Drug.Class, AMR_Gene = Best_Hit_ARO)

cat("    Total gen AMR mobile:", nrow(df_mobile), "\n")
cat("    Populasi            :", paste(unique(df_mobile$Country), collapse = ", "), "\n\n")

# ============================================================
# Fungsi: Bangun graph per populasi
# ============================================================
build_network <- function(data, country_id) {
  sub <- data %>% filter(Country == country_id)
  if (nrow(sub) == 0) return(NULL)
  
  # Edge list: AMR_Gene — MGE_Type (bobot = jumlah co-occurrence)
  edges <- sub %>%
    count(AMR_Gene, MGE_Type, name = "weight") %>%
    rename(from = AMR_Gene, to = MGE_Type)
  
  if (nrow(edges) == 0) return(NULL)
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  # Tandai tipe node (AMR atau MGE)
  V(g)$node_type <- ifelse(
    V(g)$name %in% unique(sub$AMR_Gene), "AMR Gene", "MGE Type"
  )
  V(g)$degree <- degree(g)
  E(g)$weight <- edges$weight
  
  return(g)
}

# ============================================================
# Fungsi: Plot network dengan ggraph
# ============================================================
plot_network <- function(g, title, color_amr = "#1565C0", color_mge = "#B71C1C") {
  if (is.null(g)) {
    cat("    [SKIP] Graph kosong untuk:", title, "\n")
    return(NULL)
  }
  
  ggraph(g, layout = LAYOUT_METHOD) +
    geom_edge_link(aes(width = weight), alpha = 0.3, color = "gray40") +
    geom_node_point(aes(color = node_type, size = degree)) +
    geom_node_label(aes(label = name, color = node_type),
                    size = 2.5, repel = TRUE, max.overlaps = 20) +
    scale_color_manual(values = c("AMR Gene" = color_amr, "MGE Type" = color_mge)) +
    scale_edge_width(range = c(0.5, 3)) +
    scale_size(range = c(3, 10)) +
    labs(title = title,
         subtitle = "Edge = co-occurrence; Node size = degree centrality",
         color = "Node Type", size = "Degree") +
    theme_graph(base_family = "sans", base_size = 11) +
    theme(legend.position = "right")
}

# ============================================================
# Build & Plot per populasi
# ============================================================
countries <- unique(df_mobile$Country)
graphs <- list()
network_metrics <- list()

for (ctr in countries) {
  cat(paste0("[2.", which(countries == ctr), "] Processing network: ", ctr, "...\n"))
  
  g <- build_network(df_mobile, ctr)
  graphs[[ctr]] <- g
  
  if (!is.null(g)) {
    p <- plot_network(g, title = paste("AMR-MGE Network:", ctr))
    out_path <- file.path(fig_dir, paste0("Fig3_Network_", ctr, ".pdf"))
    ggsave(out_path, p, width = 10, height = 8)
    cat("    -> Saved:", out_path, "\n")
    
    # Hitung network metrics
    network_metrics[[ctr]] <- data.frame(
      Country           = ctr,
      N_nodes           = vcount(g),
      N_edges           = ecount(g),
      Density           = edge_density(g),
      Avg_Degree        = mean(degree(g)),
      Max_Degree        = max(degree(g)),
      Hub_Node          = V(g)$name[which.max(degree(g))],
      stringsAsFactors  = FALSE
    )
    cat(paste0("    Nodes=", vcount(g), " | Edges=", ecount(g),
               " | Density=", round(edge_density(g), 3),
               " | Hub=", V(g)$name[which.max(degree(g))], "\n"))
  }
}

# ============================================================
# Jaccard Similarity antar network
# ============================================================
cat("\n[3] Menghitung Jaccard Similarity antar jaringan populasi...\n")

get_edge_set <- function(g) {
  if (is.null(g)) return(character(0))
  apply(get.edgelist(g), 1, function(x) paste(sort(x), collapse = "--"))
}

jaccard <- function(a, b) {
  if (length(union(a, b)) == 0) return(NA)
  length(intersect(a, b)) / length(union(a, b))
}

jaccard_results <- list()
for (i in seq_along(countries)) {
  for (j in seq_along(countries)) {
    if (i < j) {
      ea <- get_edge_set(graphs[[countries[i]]])
      eb <- get_edge_set(graphs[[countries[j]]])
      jaccard_results[[length(jaccard_results) + 1]] <- data.frame(
        Pop1 = countries[i], Pop2 = countries[j],
        Jaccard = round(jaccard(ea, eb), 3),
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(jaccard_results) > 0) {
  jac_df <- do.call(rbind, jaccard_results)
  cat("    Jaccard Similarity antar network:\n")
  print(jac_df)
  write.csv(jac_df, file.path(tab_dir, "Table_Jaccard_Similarity.csv"), row.names = FALSE)
  cat("    -> Saved: Table_Jaccard_Similarity.csv\n")
}

# ============================================================
# Tabel Network Metrics
# ============================================================
if (length(network_metrics) > 0) {
  metrics_df <- do.call(rbind, network_metrics)
  write.csv(metrics_df, file.path(tab_dir, "Table_Network_Metrics.csv"), row.names = FALSE)
  cat("\n    Network Metrics Summary:\n")
  print(metrics_df)
  cat("    -> Saved: Table_Network_Metrics.csv\n")
}

cat("\n============================================================\n")
cat("SELESAI - Phase 5 Network Analysis\n")
cat("  figures/Fig3_Network_*.pdf  : Network per populasi\n")
cat("  tables/Table_Network_*.csv  : Metrics dan Jaccard similarity\n")
cat("============================================================\n")
