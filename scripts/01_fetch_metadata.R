#!/usr/bin/env Rscript
# Script: 01_fetch_metadata.R
# Purpose: Kurasi metadata sampel metagenom gut dari curatedMetagenomicData.
#          Semua parameter seleksi dibaca dari config.yaml — tidak ada hardcoded filter.
#          Output: full_samples.csv + sample_map.csv + QC report

suppressPackageStartupMessages(library(curatedMetagenomicData))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(yaml))

cat("========================================================\n")
cat("[Fase 1] Kurasi Metadata - curatedMetagenomicData\n")
cat("========================================================\n\n")

# =============================================================================
# STEP 0: Load config.yaml
# =============================================================================
config_path <- "config.yaml"
if (!file.exists(config_path)) {
  stop(sprintf("[ERROR] config.yaml tidak ditemukan di '%s'. Pastikan kamu menjalankan script dari root direktori proyek.", config_path))
}

cfg <- yaml::read_yaml(config_path)
ss  <- cfg$sample_selection   # shorthand untuk sample_selection

# Ekstrak semua parameter dari config
n_per_country       <- ss$n_per_country
random_seed         <- ss$random_seed
body_site_filter    <- ss$body_site
countries_map       <- ss$countries          # list: {DNK: "Europe", CHN: "East_Asia", ...}
country_names_map   <- ss$country_names      # list: {DNK: "Denmark", ...}
min_reads           <- ss$min_reads
disease_status      <- ss$disease_status
antibiotics_filter  <- ss$antibiotics_current_use
platforms           <- ss$sequencing_platforms  # vector/list
study_cond          <- ss$study_condition
infant_age_thresh   <- ss$infant_age_threshold
max_infant          <- ss$max_infant_per_country

# Derive vectors dari config
target_countries <- names(countries_map)

cat("--- Konfigurasi yang digunakan (dari config.yaml) ---\n")
cat(sprintf("  Negara target      : %s\n", paste(target_countries, collapse=", ")))
cat(sprintf("  n per negara       : %d\n", n_per_country))
cat(sprintf("  Random seed        : %d\n", random_seed))
cat(sprintf("  body_site          : %s\n", body_site_filter))
cat(sprintf("  min_reads          : %s\n", format(min_reads, big.mark=",")))
cat(sprintf("  disease            : %s\n", disease_status))
cat(sprintf("  antibiotics_use    : %s\n", antibiotics_filter))
cat(sprintf("  platforms          : %s\n", paste(platforms, collapse=", ")))
cat(sprintf("  study_condition    : %s\n", study_cond))
cat("-----------------------------------------------------\n\n")

# =============================================================================
# STEP 1: Load metadata
# =============================================================================
cat("1. Me-load database metadata (sampleMetadata)...\n")
meta <- sampleMetadata
cat(sprintf("   Total sampel di database: %s\n\n", format(nrow(meta), big.mark=",")))

# =============================================================================
# STEP 2: Terapkan Hard Filters dari config
# =============================================================================
cat("2. Menerapkan 9 Hard Filters dari config.yaml...\n")

target_data <- meta %>%
  # [1] body_site
  filter(body_site == body_site_filter) %>%

  # [2] Negara target
  filter(country %in% target_countries) %>%

  # [3] Harus punya NCBI accession
  filter(!is.na(NCBI_accession)) %>%

  # [4] Harus punya data jumlah reads
  filter(!is.na(number_reads)) %>%

  # [5] Minimum reads
  filter(number_reads >= min_reads) %>%

  # [6] Status kesehatan
  filter(disease == disease_status) %>%

  # [7] Tidak sedang konsumsi antibiotik
  filter(antibiotics_current_use %in% antibiotics_filter) %>%

  # [8] Platform Illumina only
  filter(sequencing_platform %in% platforms) %>%

  # [9] Kondisi studi = control
  filter(study_condition %in% study_cond) %>%

  select(country, NCBI_accession, study_name, sample_id, subject_id,
         number_reads, disease, age, gender, sequencing_platform,
         antibiotics_current_use, study_condition, everything())

cat("\n--- Distribusi sampel SETELAH hard filter ---\n")
eligible_summary <- target_data %>% count(country, name = "n_eligible")
print(eligible_summary)
cat("\n")

# Tampilkan peringatan per negara yang kekurangan sampel
for (cc in target_countries) {
  n_avail <- eligible_summary %>% filter(country == cc) %>% pull(n_eligible)
  if (length(n_avail) == 0) n_avail <- 0
  if (n_avail == 0) {
    cat(sprintf("  ⛔ ERROR : %s tidak punya sampel yang memenuhi kriteria!\n", cc))
  } else if (n_avail < n_per_country) {
    cat(sprintf("  ⚠ WARNING: %s hanya punya %d sampel (target %d). Semua akan digunakan.\n",
                cc, n_avail, n_per_country))
  } else {
    cat(sprintf("  ✓ OK: %s memiliki %d sampel eligible (akan diambil %d).\n",
                cc, n_avail, n_per_country))
  }
}
cat("\n")

# =============================================================================
# STEP 3: Stratified Random Sampling
# =============================================================================
cat(sprintf("3. Stratified Random Sampling (seed=%d, maks %d per negara)...\n",
            random_seed, n_per_country))
set.seed(random_seed)

final_samples <- target_data %>%
  group_by(country) %>%
  slice_sample(n = min(n_per_country, n())) %>%
  ungroup()

cat("\n--- SAMPEL FINAL TERPILIH ---\n")
print(final_samples %>% count(country, name = "n_selected"))
cat(sprintf("Total keseluruhan: %d sampel\n\n", nrow(final_samples)))

# =============================================================================
# STEP 4: Post-sampling QC Report
# =============================================================================
cat("4. POST-SAMPLING QC REPORT\n")
cat("----------------------------------\n")

# 4a. Distribusi usia
cat("[QC] Distribusi usia per negara:\n")
age_summary <- final_samples %>%
  group_by(country) %>%
  summarise(
    n          = n(),
    n_missing  = sum(is.na(age)),
    age_median = round(median(age, na.rm = TRUE), 1),
    age_min    = suppressWarnings(min(age, na.rm = TRUE)),
    age_max    = suppressWarnings(max(age, na.rm = TRUE)),
    n_infant   = sum(age < infant_age_thresh, na.rm = TRUE),
    .groups    = "drop"
  )
print(age_summary)
infant_warning <- age_summary %>% filter(n_infant > max_infant)
if (nrow(infant_warning) > 0) {
  cat(sprintf("  ⚠ WARNING: Negara berikut punya >%d infant (usia < %d tahun): %s\n",
              max_infant, infant_age_thresh,
              paste(infant_warning$country, collapse=", ")))
  cat("    Pertimbangkan filter age > 3 atau stratifikasi usia di analisis.\n")
}
cat("\n")

# 4b. Keseimbangan gender
cat("[QC] Distribusi gender per negara:\n")
gender_summary <- final_samples %>%
  group_by(country, gender) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = gender, values_from = n, values_fill = 0)
print(gender_summary)
cat("\n")

# 4c. Study diversity (waspadai study effect)
cat("[QC] Jumlah studi yang berkontribusi per negara:\n")
study_summary <- final_samples %>%
  group_by(country) %>%
  summarise(
    n_studies = n_distinct(study_name),
    studies   = paste(unique(study_name), collapse = " | "),
    .groups   = "drop"
  )
print(study_summary)
single_study <- study_summary %>% filter(n_studies == 1)
if (nrow(single_study) > 0) {
  cat(sprintf("  ⚠ WARNING: Negara berikut hanya punya 1 studi: %s\n",
              paste(single_study$country, collapse=", ")))
  cat("    Cantumkan ini sebagai 'Limitation of Study' di naskah.\n")
}
cat("\n")

# =============================================================================
# STEP 5: Simpan Output
# =============================================================================
cat("5. Menyimpan output...\n")
out_dir <- "data/metadata"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Full metadata
write.csv(final_samples,
          file.path(out_dir, "full_samples.csv"),
          row.names = FALSE, quote = FALSE)
cat("   ✓ data/metadata/full_samples.csv\n")

# Bangun country_map dari config (bukan hardcoded!)
country_map <- data.frame(
  country      = names(countries_map),
  region       = unlist(countries_map),
  country_name = unlist(country_names_map[names(countries_map)]),
  stringsAsFactors = FALSE
)

# Sample map untuk Snakemake pipeline
sample_map <- final_samples %>%
  select(country, accession = NCBI_accession) %>%
  left_join(country_map, by = "country") %>%
  mutate(sample_id = paste0(country, "_", accession)) %>%
  select(sample_id, country, country_name, region, accession)

write.csv(sample_map,
          file.path(out_dir, "sample_map.csv"),
          row.names = FALSE, quote = FALSE)
cat("   ✓ data/metadata/sample_map.csv\n\n")

# =============================================================================
# RINGKASAN AKHIR
# =============================================================================
cat("========================================================\n")
cat("RINGKASAN FINAL\n")
cat("========================================================\n")
cat(sprintf("Total sampel terpilih  : %d\n", nrow(final_samples)))
cat(sprintf("Total negara           : %d\n", n_distinct(final_samples$country)))
cat(sprintf("Total wilayah (region) : %d\n", n_distinct(country_map$region)))
cat("\nDistribusi per wilayah:\n")
print(
  final_samples %>%
    left_join(country_map, by = "country") %>%
    count(region, country, country_name, name = "n_samples")
)
cat("\nOutput files:\n")
cat("  - data/metadata/full_samples.csv  : metadata lengkap\n")
cat("  - data/metadata/sample_map.csv    : mapping untuk Snakemake pipeline\n\n")
cat("LANGKAH SELANJUTNYA:\n")
cat("  1. Review QC report di atas (usia, gender, study effect)\n")
cat("  2. Download FASTQ menggunakan kolom 'accession' di SRA Toolkit\n")
cat("  3. Assembly dengan MEGAHIT di HPC BRIN / Galaxy\n")
cat("  4. Rename contig FASTA: format '{sample_id}.fasta' → letakkan di data/contigs/\n")
