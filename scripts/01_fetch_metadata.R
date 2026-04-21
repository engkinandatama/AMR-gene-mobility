#!/usr/bin/env Rscript
# Script: 01_fetch_metadata.R
# Purpose: Mengambil metadata sampel metagenom gut dari curatedMetagenomicData,
#          fokus pada negara Denmark (Eropa), China (Asia Timur), dan India (Asia Selatan).
#          Output: pilot_samples.csv + sample_map.csv (country/region mapping untuk pipeline)

suppressPackageStartupMessages(library(curatedMetagenomicData))
suppressPackageStartupMessages(library(dplyr))

cat("========================================================\n")
cat("[Fase 1] Menarik metadata dari curatedMetagenomicData...\n")
cat("========================================================\n\n")

# 1. Menarik kerangka metadata seluruh database
cat("1. Me-load database metadata...\n")
meta <- sampleMetadata

# 2. Filter berdasarkan region & tipe sampel
cat("2. Memfilter sampel (stool) dari Denmark, China, dan India...\n")
target_data <- meta %>%
  filter(body_site == "stool") %>%
  filter(country %in% c("DNK", "CHN", "IND")) %>%
  filter(!is.na(NCBI_accession)) %>%
  select(country, NCBI_accession, study_name, sample_id, subject_id, everything())

cat("   Total sampel yang memenuhi kriteria awal:", nrow(target_data), "sampel.\n\n")

# 3. Pengambilan Acak (Random Sampling) untuk Pilot
cat("3. Melakukan seleksi representasi untuk PILOT RUN (n=3)...\n")
set.seed(123) # Mengunci seed agar sampel yang terpilih fix

pilot_samples <- target_data %>%
  group_by(country) %>%
  sample_n(1) %>%
  ungroup()

cat("----------------------------------\n")
cat("SAMPEL UJI COBA (PILOT RUN) TERPILIH:\n")
cat("----------------------------------\n")
print(pilot_samples %>% select(country, NCBI_accession, study_name))
cat("\n")

# 4. Menyimpan pilot_samples.csv
cat("4. Menyimpan daftar ke 'data/metadata/pilot_samples.csv' ...\n")
dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)
write.csv(pilot_samples, "data/metadata/pilot_samples.csv", row.names = FALSE, quote = FALSE)

# 5. Auto-generate sample_map.csv
#    (Single Source of Truth: sample_id -> country -> region)
cat("5. Membuat 'data/metadata/sample_map.csv' ...\n")

country_map <- data.frame(
  country      = c("DNK", "CHN", "IND"),
  country_name = c("Denmark", "China", "India"),
  region       = c("Europe", "East_Asia", "South_Asia"),
  stringsAsFactors = FALSE
)

sample_map <- pilot_samples %>%
  select(country, accession = NCBI_accession) %>%
  left_join(country_map, by = "country") %>%
  mutate(
    # sample_id = prefix negara + accession => digunakan sebagai nama file .fasta di pipeline
    sample_id = paste0(country, "_", accession)
  ) %>%
  select(sample_id, country, country_name, region, accession)

cat("----------------------------------\n")
cat("SAMPLE MAP (digunakan oleh Snakemake pipeline):\n")
cat("----------------------------------\n")
print(sample_map)
cat("\n")

write.csv(sample_map, "data/metadata/sample_map.csv", row.names = FALSE, quote = FALSE)

cat("=== SELESAI ===\n")
cat("Output:\n")
cat("  - data/metadata/pilot_samples.csv : metadata lengkap curatedMetagenomicData\n")
cat("  - data/metadata/sample_map.csv    : country/region mapping untuk Snakemake\n\n")
cat("LANGKAH SELANJUTNYA:\n")
cat("  Rename file .fasta dari Galaxy sesuai kolom sample_id di sample_map.csv,\n")
cat("  lalu letakkan di folder data/contigs/\n")
