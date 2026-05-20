#!/usr/bin/env Rscript
# Script: 04_run_stats.R
# Purpose: Analisis statistik komprehensif dan visualisasi untuk manuskrip.
#          Setiap blok analisis dipetakan eksplisit ke sub-pertanyaan riset.
#
# Sub-question 1: Apakah diversitas dan abundansi AMR genes berbeda antar populasi?
# Sub-question 2: MGE tipe apa yang dominan per populasi?
# Sub-question 3: Apakah AMR-MGE association berbeda antar populasi (jalur HGT)?
#
# Input:
#   - results/amr_abundance_matrix.csv
#   - results/mge_distribution_matrix.csv
#   - results/amr_mge_association_matrix.csv
#   - results/colocalization_summary.csv
#
# Output:
#   - results/figures/  (semua plot)
#   - results/tables/   (tabel statistik)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(vegan)
  library(RColorBrewer)
  library(writexl)
})

# Coba load package opsional
tryCatch({
  suppressPackageStartupMessages(library(ComplexHeatmap))
  HAVE_HEATMAP <- TRUE
}, error = function(e) { HAVE_HEATMAP <<- FALSE })

cat("============================================================\n")
cat("[Phase 4] Analisis Statistik Komprehensif - AMR Gene Mobility\n")
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
input_file   <- ifelse(!is.null(params$input), params$input, "results/amr_mge_association_matrix.csv")
output_dir   <- ifelse(!is.null(params$output_dir), params$output_dir, "results")
pop_name     <- ifelse(!is.null(params$pop), params$pop, "combined")
p_adjust     <- ifelse(!is.null(params$p_adjust_method), params$p_adjust_method, "BH")
alpha_val    <- as.numeric(ifelse(!is.null(params$alpha), params$alpha, "0.05"))
n_permutations <- as.numeric(ifelse(!is.null(params$permanova_permutations), params$permanova_permutations, "999"))

# --- Setup direktori output ---
fig_dir <- file.path(output_dir, "figures")
tab_dir <- file.path(output_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load data ---
cat("[1] Membaca data input...\n")
agg_dir <- dirname(input_file)

# File names berdasarkan prefix populasi
coloc_path <- file.path(agg_dir, paste0(pop_name, "_all_coloc.csv"))
abund_path <- file.path(agg_dir, paste0(pop_name, "_amr_abundance.csv"))
mge_path   <- file.path(agg_dir, paste0(pop_name, "_mge_distribution.csv"))
assoc_path <- input_file

if (!file.exists(coloc_path) || !file.exists(abund_path) || !file.exists(mge_path)) {
  cat("[WARN] Salah satu file input tidak ditemukan. Mungkin populasi ini kosong.\n")
  cat("       Mencoba membaca fallback dari folder results/ jika ada...\n")
  coloc_path <- "results/colocalization_summary.csv"
  abund_path <- "results/amr_abundance_matrix.csv"
  mge_path   <- "results/mge_distribution_matrix.csv"
  assoc_path <- "results/amr_mge_association_matrix.csv"
}

# Membaca data secara aman
coloc_df  <- tryCatch(read.csv(coloc_path), error = function(e) data.frame())
abund_df  <- tryCatch(read.csv(abund_path), error = function(e) data.frame())
mge_df    <- tryCatch(read.csv(mge_path), error = function(e) data.frame())
assoc_df  <- tryCatch(read.csv(assoc_path), error = function(e) data.frame())

if (nrow(coloc_df) == 0 || nrow(abund_df) == 0) {
  cat("[WARN] Data kosong untuk populasi ini. Menghentikan analisis secara aman.\n")
  # Buat dummy file agar Snakemake tidak error
  writeLines("Analisis statistik di-skip karena data kosong untuk populasi ini.", file.path(output_dir, "summary_stats.txt"))
  quit(status = 0)
}

cat("    colocalization_summary: ", nrow(coloc_df), "baris\n")
cat("    amr_abundance_matrix  : ", nrow(abund_df), "sampel\n\n")

# Warna konsisten per populasi (mendukung negara baru secara dinamis)
base_colors <- c("DNK" = "#2196F3", "CHN" = "#F44336", "IND" = "#4CAF50", "CMR" = "#FF9800", "MDG" = "#9C27B0")
unique_countries <- unique(c(abund_df$Country, coloc_df$Country))
pop_colors <- base_colors
missing_countries <- setdiff(unique_countries, names(pop_colors))
if (length(missing_countries) > 0) {
  # Gunakan RColorBrewer untuk warna tambahan jika ada negara baru
  extra_colors <- colorRampPalette(suppressWarnings(brewer.pal(8, "Set2")))(length(missing_countries))
  names(extra_colors) <- missing_countries
  pop_colors <- c(pop_colors, extra_colors)
}

#=============================================================
# 4A: DIVERSITY ANALYSIS -> Sub-question 1
#=============================================================
cat("[4A] Diversity Analysis (Menjawab Sub-question 1)...\n")

# Buat matrix numerik: baris = sampel, kolom = drug class
meta_cols <- c("Sample_ID", "Country", "Region", "Age", "Gender", "Study")
num_cols <- setdiff(colnames(abund_df), meta_cols)

if (length(num_cols) > 0 && nrow(abund_df) >= 3) {
  species_matrix <- abund_df[, num_cols]
  rownames(species_matrix) <- abund_df$Sample_ID

  # Alpha diversity: Shannon Index
  abund_df$Shannon <- diversity(species_matrix, index = "shannon")
  abund_df$Richness <- rowSums(species_matrix > 0)

  cat("    Alpha diversity (Shannon) per sampel:\n")
  # Print conditionally to avoid errors if Age/Gender are not loaded
  print_cols <- intersect(c("Sample_ID", "Country", "Shannon", "Richness", "Age", "Gender"), colnames(abund_df))
  print(abund_df %>% select(all_of(print_cols)))

  # Plot alpha diversity
  p_alpha <- ggplot(abund_df, aes(x = Country, y = Shannon, fill = Country)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.1, size = 2) +
    scale_fill_manual(values = pop_colors) +
    labs(title = "Alpha Diversity AMR Genes (Shannon Index)",
         subtitle = "Sub-question 1: Apakah diversitas AMR berbeda antar populasi?",
         x = "Populasi", y = "Shannon Index") +
    theme_bw(base_size = 12) +
    theme(legend.position = "none")
  ggsave(file.path(fig_dir, "Fig_S1_alpha_diversity.pdf"), p_alpha, width = 6, height = 5)
  cat("    -> Saved: Fig_S1_alpha_diversity.pdf\n")

  # Beta diversity: Bray-Curtis + PERMANOVA + PERMDISP
  if (nrow(abund_df) >= 3) {
    bc_dist <- vegdist(species_matrix, method = "bray")
    
    # Simple PERMANOVA
    permanova_result <- adonis2(bc_dist ~ Country, data = abund_df, permutations = n_permutations)
    cat("\n    PERMANOVA (Beta diversity - Bray-Curtis):\n")
    print(permanova_result)
    write.csv(as.data.frame(permanova_result), file.path(tab_dir, "Table_PERMANOVA.csv"))
    cat("    -> Saved: Table_PERMANOVA.csv\n")
    
    # Homogeneity of Dispersion (PERMDISP)
    cat("\n    PERMDISP (Multivariate Homogeneity of Groups Dispersions):\n")
    disp_result <- tryCatch({
      betadisper(bc_dist, abund_df$Country)
    }, error = function(e) {
      cat("    [WARN] betadisper failed:", e$message, "\n")
      NULL
    })
    
    if (!is.null(disp_result)) {
      permutest_disp <- tryCatch({
        permutest(disp_result, permutations = n_permutations)
      }, error = function(e) {
        cat("    [WARN] permutest on betadisper failed:", e$message, "\n")
        NULL
      })
      
      if (!is.null(permutest_disp)) {
        print(permutest_disp)
        write.csv(as.data.frame(permutest_disp$tab), file.path(tab_dir, "Table_PERMDISP.csv"))
        cat("    -> Saved: Table_PERMDISP.csv\n")
      }
      
      # Plot PERMDISP
      pdf(file.path(fig_dir, "Fig_S2_permdisp_plot.pdf"), width = 6, height = 5)
      plot(disp_result, main = "Multivariate Dispersions (PERMDISP)", sub = "")
      dev.off()
      cat("    -> Saved: Fig_S2_permdisp_plot.pdf\n")
    }

    # Pairwise PERMANOVA (Post-hoc comparisons between countries)
    unique_countries <- unique(abund_df$Country)
    if (length(unique_countries) > 2) {
      cat("\n    Pairwise PERMANOVA (Post-hoc comparisons between countries):\n")
      pairwise_results <- list()
      pairs <- combn(unique_countries, 2, simplify = FALSE)
      
      for (pair in pairs) {
        pair_df <- abund_df %>% filter(Country %in% pair)
        if (nrow(pair_df) >= 3) {
          pair_matrix <- pair_df[, num_cols]
          pair_bc <- vegdist(pair_matrix, method = "bray")
          
          pair_perm <- tryCatch({
            adonis2(pair_bc ~ Country, data = pair_df, permutations = n_permutations)
          }, error = function(e) NULL)
          
          if (!is.null(pair_perm)) {
            p_val <- pair_perm$`Pr(>F)`[1]
            r2 <- pair_perm$R2[1]
            f_stat <- pair_perm$F[1]
            
            pairwise_results[[paste(pair, collapse = " vs ")]] <- data.frame(
              Comparison = paste(pair, collapse = " vs "),
              F_value = f_stat,
              R2 = r2,
              p_value = p_val,
              stringsAsFactors = FALSE
            )
          }
        }
      }
      
      if (length(pairwise_results) > 0) {
        pairwise_df <- do.call(rbind, pairwise_results)
        pairwise_df$p_adjusted <- p.adjust(pairwise_df$p_value, method = p_adjust)
        print(pairwise_df)
        write.csv(pairwise_df, file.path(tab_dir, "Table_Pairwise_PERMANOVA.csv"), row.names = FALSE)
        cat("    -> Saved: Table_Pairwise_PERMANOVA.csv\n")
      }
    }

    # Multi-factor PERMANOVA (controlling for covariates Age, Gender)
    factors_to_try <- c("Gender", "Age")
    valid_factors <- c()
    for (f in factors_to_try) {
      if (f %in% colnames(abund_df)) {
        val_vector <- abund_df[[f]]
        non_na_count <- sum(!is.na(val_vector) & val_vector != "Unknown" & val_vector != "")
        unique_count <- length(unique(val_vector[!is.na(val_vector) & val_vector != "Unknown" & val_vector != ""]))
        if (non_na_count >= 3 && unique_count >= 2) {
          valid_factors <- c(valid_factors, f)
        }
      }
    }

    if (length(valid_factors) > 0) {
      clean_df <- abund_df
      for (f in valid_factors) {
        clean_df <- clean_df %>% filter(!is.na(.data[[f]]) & .data[[f]] != "Unknown" & .data[[f]] != "")
      }
      
      if (nrow(clean_df) >= 4) {
        clean_matrix <- clean_df[, num_cols]
        clean_bc <- vegdist(clean_matrix, method = "bray")
        
        multi_formula <- as.formula(paste("clean_bc ~ Country +", paste(valid_factors, collapse = " + ")))
        cat("\n    Multi-factor PERMANOVA (controlling for covariates):\n")
        cat("    Formula:", deparse(multi_formula), "\n")
        
        multi_permanova <- tryCatch({
          adonis2(multi_formula, data = clean_df, permutations = n_permutations)
        }, error = function(e) {
          cat("    [WARN] Multi-factor PERMANOVA failed:", e$message, "\n")
          NULL
        })
        
        if (!is.null(multi_permanova)) {
          print(multi_permanova)
          write.csv(as.data.frame(multi_permanova), file.path(tab_dir, "Table_MultiFactor_PERMANOVA.csv"))
          cat("    -> Saved: Table_MultiFactor_PERMANOVA.csv\n")
        }
      }
    }

    # --- Mantel Test & Procrustes: Correlation between AMR and Taxonomy ---
    tax_file <- "data/metadata/taxonomy_abundance.csv"
    if (file.exists(tax_file)) {
      cat("\n    [INFO] Membaca data taksonomi komunitas (MetaPhlAn) untuk analisis Mantel...\n")
      tax_df <- tryCatch(read.csv(tax_file, check.names = FALSE), error = function(e) data.frame())
      if (nrow(tax_df) > 0) {
        # Filter tax_df untuk mencocokkan Sample_ID yang ada di abund_df saat ini
        common_samples <- intersect(abund_df$Sample_ID, tax_df$Sample_ID)
        if (length(common_samples) >= 3) {
          # Subset dan urutkan kedua data frame agar urutan sampel sama persis
          amr_sub <- abund_df %>% filter(Sample_ID %in% common_samples) %>% arrange(Sample_ID)
          tax_sub <- tax_df %>% filter(Sample_ID %in% common_samples) %>% arrange(Sample_ID)
          
          # Ekstrak matriks numerik
          amr_mat_sub <- amr_sub[, num_cols, drop = FALSE]
          tax_mat_sub <- tax_sub[, setdiff(colnames(tax_sub), "Sample_ID"), drop = FALSE]
          
          # Hitung matriks jarak Bray-Curtis
          amr_bc_sub <- vegdist(amr_mat_sub, method = "bray")
          tax_bc_sub <- vegdist(tax_mat_sub, method = "bray")
          
          # Jalankan Mantel Test
          cat("    Menjalankan Mantel Test (AMR vs Profil Taksonomi Species)...\n")
          mantel_res <- tryCatch({
            mantel(amr_bc_sub, tax_bc_sub, method = "spearman", permutations = n_permutations)
          }, error = function(e) {
            cat("    [WARN] Mantel test failed:", e$message, "\n")
            NULL
          })
          
          if (!is.null(mantel_res)) {
            print(mantel_res)
            mantel_df <- data.frame(
              Analysis = "Mantel Test (AMR vs Taxonomic Abundance)",
              Correlation_r = mantel_res$statistic,
              p_value = mantel_res$signif,
              Permutations = mantel_res$permutations,
              stringsAsFactors = FALSE
            )
            write.csv(mantel_df, file.path(tab_dir, "Table_Mantel_Taxonomy.csv"), row.names = FALSE)
            cat("    -> Saved: Table_Mantel_Taxonomy.csv\n")
            
            # --- Dual-Mantel Test (Decoupling Analysis) ---
            cat("    Menjalankan Dual-Mantel Test (Decoupling Analysis: Mobile vs Non-Mobile AMR vs Taxonomy)...\n")
            
            # List of all unique samples
            all_samples <- coloc_df %>% 
              select(Sample_ID, Country) %>% 
              distinct()

            # Mobile AMR Matrix (all samples included, missing filled with 0)
            mobile_counts <- coloc_df %>%
              filter(MGE_Type != "Chromosomal") %>%
              count(Sample_ID, `Drug.Class`, name = "Count")

            if (nrow(mobile_counts) > 0) {
              mobile_pivot <- mobile_counts %>%
                pivot_wider(names_from = `Drug.Class`, values_from = Count, values_fill = 0)
              mobile_abund <- all_samples %>% 
                left_join(mobile_pivot, by = "Sample_ID") %>%
                mutate(across(-c(Sample_ID, Country), ~ replace_na(.x, 0)))
            } else {
              mobile_abund <- all_samples
            }

            # Non-Mobile AMR Matrix (all samples included, missing filled with 0)
            nonmobile_counts <- coloc_df %>%
              filter(MGE_Type == "Chromosomal") %>%
              count(Sample_ID, `Drug.Class`, name = "Count")

            if (nrow(nonmobile_counts) > 0) {
              nonmobile_pivot <- nonmobile_counts %>%
                pivot_wider(names_from = `Drug.Class`, values_from = Count, values_fill = 0)
              nonmobile_abund <- all_samples %>% 
                left_join(nonmobile_pivot, by = "Sample_ID") %>%
                mutate(across(-c(Sample_ID, Country), ~ replace_na(.x, 0)))
            } else {
              nonmobile_abund <- all_samples
            }

            # Filter for common samples
            mobile_sub <- mobile_abund %>% filter(Sample_ID %in% common_samples) %>% arrange(Sample_ID)
            nonmobile_sub <- nonmobile_abund %>% filter(Sample_ID %in% common_samples) %>% arrange(Sample_ID)
            
            # Extract numeric matrices
            mobile_mat_sub <- mobile_sub[, setdiff(colnames(mobile_sub), c("Sample_ID", "Country")), drop = FALSE]
            nonmobile_mat_sub <- nonmobile_sub[, setdiff(colnames(nonmobile_sub), c("Sample_ID", "Country")), drop = FALSE]
            
            # Run Mantel for Non-Mobile AMR vs Taxonomy
            nonmobile_bc <- vegdist(nonmobile_mat_sub, method = "bray")
            mantel_nonmobile <- tryCatch({
              mantel(nonmobile_bc, tax_bc_sub, method = "spearman", permutations = n_permutations)
            }, error = function(e) {
              cat("    [WARN] Non-mobile Mantel test failed:", e$message, "\n")
              NULL
            })
            
            # Run Mantel for Mobile AMR vs Taxonomy
            mantel_mobile <- NULL
            if (ncol(mobile_mat_sub) > 0 && sum(mobile_mat_sub) > 0) {
              mobile_bc <- vegdist(mobile_mat_sub, method = "bray")
              mantel_mobile <- tryCatch({
                mantel(mobile_bc, tax_bc_sub, method = "spearman", permutations = n_permutations)
              }, error = function(e) {
                cat("    [WARN] Mobile Mantel test failed:", e$message, "\n")
                NULL
              })
            }
            
            # Save results
            dual_mantel_results <- data.frame(
              Analysis = c("Non-Mobile AMR vs Taxonomy (Vertical)", "Mobile AMR vs Taxonomy (Horizontal/HGT)"),
              Correlation_r = c(
                ifelse(!is.null(mantel_nonmobile), mantel_nonmobile$statistic, NA),
                ifelse(!is.null(mantel_mobile), mantel_mobile$statistic, NA)
              ),
              p_value = c(
                ifelse(!is.null(mantel_nonmobile), mantel_nonmobile$signif, NA),
                ifelse(!is.null(mantel_mobile), mantel_mobile$signif, NA)
              ),
              Permutations = c(
                ifelse(!is.null(mantel_nonmobile), mantel_nonmobile$permutations, NA),
                ifelse(!is.null(mantel_mobile), mantel_mobile$permutations, NA)
              ),
              stringsAsFactors = FALSE
            )
            write.csv(dual_mantel_results, file.path(tab_dir, "Table_Dual_Mantel_Decoupling.csv"), row.names = FALSE)
            cat("    -> Saved: Table_Dual_Mantel_Decoupling.csv\n")
            
            # Procrustes Analysis
            cat("\n    Menjalankan Procrustes Analysis...\n")
            amr_pcoa <- cmdscale(amr_bc_sub, k = 2)
            tax_pcoa <- cmdscale(tax_bc_sub, k = 2)
            pro_res <- tryCatch({
              procrustes(amr_pcoa, tax_pcoa, symmetric = TRUE)
            }, error = function(e) {
              cat("    [WARN] Procrustes failed:", e$message, "\n")
              NULL
            })
            
            if (!is.null(pro_res)) {
              pro_test <- protest(amr_pcoa, tax_pcoa, permutations = n_permutations)
              print(pro_test)
              
              procrustes_df <- data.frame(
                Analysis = "Procrustes (AMR PCoA vs Taxonomy PCoA)",
                Correlation_R = pro_test$t0,
                p_value = pro_test$signif,
                Permutations = pro_test$permutations,
                stringsAsFactors = FALSE
              )
              write.csv(procrustes_df, file.path(tab_dir, "Table_Procrustes_Taxonomy.csv"), row.names = FALSE)
              cat("    -> Saved: Table_Procrustes_Taxonomy.csv\n")
              
              # Plot Procrustes
              pdf(file.path(fig_dir, "Fig_S3_procrustes_plot.pdf"), width = 6, height = 5)
              plot(pro_res, main = "Procrustes Superimposition (AMR vs Taxonomy)", sub = paste("R =", round(pro_test$t0, 3), ", p =", pro_test$signif))
              points(pro_res, display = "target", col = "red", pch = 16)
              points(pro_res, display = "rotated", col = "blue", pch = 17)
              dev.off()
              cat("    -> Saved: Fig_S3_procrustes_plot.pdf\n")
            }
          }
        } else {
          cat("    [WARN] Jumlah sampel iris untuk Mantel terlalu sedikit (butuh >=3 sampel)\n")
        }
      }
    } else {
      cat("    [INFO] data/metadata/taxonomy_abundance.csv tidak ditemukan. Skip analisis Mantel.\n")
    }
  }
} else {
  cat("    [WARN] Data sampel terlalu sedikit untuk diversity analysis (butuh >=3 sampel)\n")
}

#=============================================================
# 4B: AMR CLASS COMPARISON -> Sub-question 1
#=============================================================
cat("\n[4B] AMR Class Comparison - Kruskal-Wallis (Sub-question 1)...\n")

amr_long <- coloc_df %>%
  count(Sample_ID, Country, Region, `Drug.Class`, name = "Count") %>%
  rename(Drug_Class = `Drug.Class`)

# Kruskal-Wallis per kelas AMR
kw_results <- amr_long %>%
  group_by(Drug_Class) %>%
  summarise(
    n_samples = n(),
    kw_p_value = tryCatch({
      if (n_distinct(Country) >= 2 && n() >= 3)
        kruskal.test(Count ~ Country)$p.value
      else NA_real_
    }, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(p_adjusted = p.adjust(kw_p_value, method = p_adjust)) %>%
  arrange(kw_p_value)

cat("    Top AMR Classes by Kruskal-Wallis significance:\n")
print(head(kw_results, 10))
write.csv(kw_results, file.path(tab_dir, "Table1_KruskalWallis_AMR_Class.csv"), row.names = FALSE)
cat("    -> Saved: Table1_KruskalWallis_AMR_Class.csv\n")

#=============================================================
# 4C: MGE TYPE DISTRIBUTION -> Sub-question 2
#=============================================================
cat("\n[4C] MGE Type Distribution (Sub-question 2)...\n")

mge_colors <- c(
  "Plasmid"      = "#9C27B0",
  "Integron"     = "#FF9800",
  "IS/Transposon"= "#00BCD4",
  "Chromosomal"  = "#9E9E9E"
)

p_mge <- ggplot(mge_df, aes(x = Country, y = Proportion, fill = MGE_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = mge_colors) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Distribusi Tipe MGE yang Membawa Gen AMR",
    subtitle = "Sub-question 2: MGE tipe apa yang dominan per populasi?",
    x = "Populasi", y = "Proporsi (%)", fill = "Tipe MGE"
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(fig_dir, "Fig2_MGE_distribution.pdf"), p_mge, width = 7, height = 5)
cat("    -> Saved: Fig2_MGE_distribution.pdf\n")

# Chi-square test: apakah distribusi MGE berbeda antar populasi?
mge_contingency <- mge_df %>%
  select(Country, MGE_Type, Count) %>%
  pivot_wider(names_from = MGE_Type, values_from = Count, values_fill = 0) %>%
  column_to_rownames("Country")

if (nrow(mge_contingency) >= 2 && ncol(mge_contingency) >= 2) {
  chi_result <- chisq.test(as.matrix(mge_contingency))
  cat("    Chi-square test MGE distribution:\n")
  cat("      Chi2 =", round(chi_result$statistic, 3),
      "| df =", chi_result$parameter,
      "| p =", round(chi_result$p.value, 4), "\n")
  write.csv(data.frame(
    Chi2 = chi_result$statistic,
    df   = chi_result$parameter,
    p_value = chi_result$p.value
  ), file.path(tab_dir, "Table_ChiSquare_MGE.csv"), row.names = FALSE)
  cat("    -> Saved: Table_ChiSquare_MGE.csv\n")
}

#=============================================================
# 4D: AMR-MGE ASSOCIATION TEST -> Sub-question 3 (INTI NOVELTY)
#=============================================================
cat("\n[4D] AMR-MGE Association Test - Fisher's Exact (Sub-question 3)...\n")

# Fisher's exact untuk setiap pasangan AMR class × MGE type (mobile vs chromosomal)
coloc_df <- coloc_df %>%
  mutate(Is_Mobile = ifelse(MGE_Type != "Chromosomal", "Mobile", "Chromosomal"))

fisher_results <- list()

for (country_name in unique(coloc_df$Country)) {
  sub <- coloc_df %>% filter(Country == country_name)
  drug_classes <- unique(sub$Drug.Class)
  
  for (drug in drug_classes) {
    drug_sub <- sub %>%
      mutate(Is_Drug = ifelse(Drug.Class == drug, "Yes", "No"))
    
    ctab <- table(drug_sub$Is_Drug, drug_sub$Is_Mobile)
    if (all(dim(ctab) >= 2)) {
      ft <- tryCatch(fisher.test(ctab, simulate.p.value = TRUE), error = function(e) NULL)
      if (!is.null(ft)) {
        fisher_results[[length(fisher_results) + 1]] <- data.frame(
          Country   = country_name,
          Drug_Class = drug,
          OddsRatio = ifelse(is.numeric(ft$estimate), ft$estimate, NA),
          p_value   = ft$p.value,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

if (length(fisher_results) > 0) {
  fisher_df <- do.call(rbind, fisher_results) %>%
    mutate(p_adjusted = p.adjust(p_value, method = p_adjust)) %>%
    arrange(p_adjusted)
  
  cat("    Top significant AMR-MGE associations (FDR < 0.05):\n")
  sig <- fisher_df %>% filter(p_adjusted < 0.05)
  if (nrow(sig) > 0) print(head(sig, 10)) else cat("    (Tidak ada yang signifikan - wajar untuk pilot kecil)\n")
  
  write.csv(fisher_df, file.path(tab_dir, "Table2_Fisher_AMR_MGE.csv"), row.names = FALSE)
  cat("    -> Saved: Table2_Fisher_AMR_MGE.csv\n")
}

#=============================================================
# 4F: RESISTOME MOBILITY INDEX ANALYSIS (Enhancement for Novelty)
#=============================================================
cat("\n[4F] Resistome Mobility Index Analysis...\n")

mobility_df <- coloc_df %>%
  group_by(Sample_ID, Country) %>%
  summarise(
    Total_AMR = n(),
    Mobile_AMR = sum(MGE_Type != "Chromosomal"),
    Mobility_Index = ifelse(Total_AMR > 0, Mobile_AMR / Total_AMR, 0),
    .groups = "drop"
  )

write.csv(mobility_df, file.path(tab_dir, "Table_Mobility_Index.csv"), row.names = FALSE)
cat("    -> Saved: Table_Mobility_Index.csv\n")

# Run Kruskal-Wallis test on Resistome Mobility Index across countries
if (n_distinct(mobility_df$Country) >= 2 && nrow(mobility_df) >= 3) {
  mobility_kw <- tryCatch({
    kruskal.test(Mobility_Index ~ Country, data = mobility_df)
  }, error = function(e) NULL)
  
  if (!is.null(mobility_kw)) {
    cat("    Kruskal-Wallis test on Resistome Mobility Index:\n")
    print(mobility_kw)
    
    # Save Kruskal-Wallis result
    mobility_kw_df <- data.frame(
      Metric = "Resistome Mobility Index",
      Chi2 = mobility_kw$statistic,
      df = mobility_kw$parameter,
      p_value = mobility_kw$p.value,
      stringsAsFactors = FALSE
    )
    write.csv(mobility_kw_df, file.path(tab_dir, "Table_Mobility_KruskalWallis.csv"), row.names = FALSE)
    cat("    -> Saved: Table_Mobility_KruskalWallis.csv\n")
    
    # Post-hoc pairwise Wilcoxon test if significant (built-in R, no package dependency)
    if (mobility_kw$p.value < alpha_val) {
      cat("    Running Post-hoc Pairwise Wilcoxon Test on Mobility Index...\n")
      wt <- tryCatch({
        pairwise.wilcox.test(mobility_df$Mobility_Index, mobility_df$Country, p.adjust.method = p_adjust)
      }, error = function(e) NULL)
      
      if (!is.null(wt) && !is.null(wt$p.value)) {
        p_mat <- wt$p.value
        dunn_df <- as.data.frame(as.table(p_mat)) %>%
          filter(!is.na(Freq)) %>%
          rename(Pop1 = Var1, Pop2 = Var2, p_adjusted = Freq) %>%
          mutate(Comparison = paste0(Pop1, " vs ", Pop2)) %>%
          select(Comparison, p_adjusted)
        
        write.csv(dunn_df, file.path(tab_dir, "Table_Mobility_DunnPostHoc.csv"), row.names = FALSE)
        cat("    -> Saved: Table_Mobility_DunnPostHoc.csv\n")
      }
    }
  }
}

# Plot Mobility Index per Country
p_mobility <- ggplot(mobility_df, aes(x = Country, y = Mobility_Index, fill = Country)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  scale_fill_manual(values = pop_colors) +
  labs(
    title = "Resistome Mobility Index across Populations",
    subtitle = "Mobility Index = Mobile AMR Genes / Total AMR Genes",
    x = "Country", y = "Mobility Index"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "Fig_S4_mobility_index.pdf"), p_mobility, width = 6, height = 5)
cat("    -> Saved: Fig_S4_mobility_index.pdf\n")

#=============================================================
# 4E: HEATMAP UTAMA (Figure 1 Manuskrip)
#=============================================================
cat("\n[4E] Membuat Heatmap Utama (Figure 1)...\n")

heatmap_data <- coloc_df %>%
  count(Sample_ID, Country, `Drug.Class`, name = "Count") %>%
  pivot_wider(names_from = `Drug.Class`, values_from = Count, values_fill = 0)

heatmap_mat <- heatmap_data %>%
  select(-Sample_ID, -Country) %>%
  as.matrix()
# Hindari duplikasi awalan Country_ jika Sample_ID sudah memilikinya
rownames(heatmap_mat) <- sapply(1:nrow(heatmap_data), function(i) {
  sid <- heatmap_data$Sample_ID[i]
  ctr <- heatmap_data$Country[i]
  if (startsWith(sid, paste0(ctr, "_"))) {
    return(sid)
  } else {
    return(paste0(ctr, "_", sid))
  }
})

if (HAVE_HEATMAP && nrow(heatmap_mat) >= 2 && ncol(heatmap_mat) >= 2) {
  pdf(file.path(fig_dir, "Fig1_Heatmap_AMR.pdf"), width = 14, height = 8)
  ht <- Heatmap(
    log10(heatmap_mat + 1),
    name = "log10(count+1)",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 7),
    column_title = "Drug Class",
    row_title = "Samples (by Country)",
    col = colorRampPalette(c("white", "#1565C0"))(100)
  )
  draw(ht)
  dev.off()
  cat("    -> Saved: Fig1_Heatmap_AMR.pdf\n")
} else {
  cat("    [SKIP] ComplexHeatmap tidak tersedia atau data terlalu kecil.\n")
}

#=============================================================
# EKSPOR EXCEL TERINTEGRASI
#=============================================================
cat("\n[5] Mengekspor semua tabel ke Excel terintegrasi...\n")

sheets <- list(
  "Kruskal-Wallis AMR Class" = kw_results,
  "MGE Distribution"         = mge_df
)
if (exists("fisher_df")) sheets[["Fisher AMR-MGE"]] <- fisher_df
if (exists("permanova_result")) sheets[["PERMANOVA Beta Diversity"]] <- as.data.frame(permanova_result)
if (exists("multi_permanova") && !is.null(multi_permanova)) sheets[["Multi-Factor PERMANOVA"]] <- as.data.frame(multi_permanova)
if (exists("permutest_disp") && !is.null(permutest_disp)) sheets[["PERMDISP Homogeneity"]] <- as.data.frame(permutest_disp$tab)
if (exists("pairwise_df")) sheets[["Pairwise PERMANOVA"]] <- pairwise_df
if (exists("mantel_df")) sheets[["Mantel Taxonomy Correlation"]] <- mantel_df
if (exists("procrustes_df")) sheets[["Procrustes Taxonomy"]] <- procrustes_df
if (exists("mobility_df")) sheets[["Mobility Index per Sample"]] <- mobility_df
if (exists("mobility_kw_df")) sheets[["Mobility Kruskal-Wallis"]] <- mobility_kw_df
if (exists("dunn_df")) sheets[["Mobility Dunn PostHoc"]] <- dunn_df
if (exists("dual_mantel_results")) sheets[["Dual-Mantel Decoupling"]] <- dual_mantel_results

write_xlsx(sheets, file.path(tab_dir, "Supplementary_Statistics.xlsx"))
cat("    -> Saved: Supplementary_Statistics.xlsx\n")

cat("\n============================================================\n")
cat("SELESAI. Semua output tersimpan di:\n")
cat("  figures/  : Grafik untuk manuskrip\n")
cat("  tables/   : Tabel statistik\n")
cat("============================================================\n")

# Buat file summary_stats.txt untuk Snakemake
writeLines("Analisis statistik selesai successfully.", file.path(output_dir, "summary_stats.txt"))
