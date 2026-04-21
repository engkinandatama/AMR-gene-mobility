#!/usr/bin/env Rscript
# Script: 01_fetch_metadata.R
# Purpose: Mengambil metadata sampel metagenom gut dari curatedMetagenomicData,
#          fokus pada negara Denmark (Eropa), China (Asia Timur), dan India (Asia Selatan).
#          Mengekstrak 3 sampel acak sebagai Pilot Run untuk dikerjakan di UseGalaxy.

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

# 3. Pengambilan Serampangan / Acak (Random Sampling) untuk Pilot
cat("3. Melakukan seleksi representasi untuk PILOT RUN (n=3)...\n")
set.seed(123) # Mengunci seed agar sampel yang terpilih fix dan tidak berubah-ubah saat di rerun.

pilot_samples <- target_data %>%
  group_by(country) %>%
  sample_n(1) %>%
  ungroup()

cat("----------------------------------\n")
cat("SAMPEL UJI COBA (PILOT RUN) TERPILIH:\n")
cat("----------------------------------\n")
print(pilot_samples %>% select(country, NCBI_accession, study_name))
cat("\n")

# 4. Menyimpan output
cat("4. Menyimpan daftar ke 'data/metadata/pilot_samples.csv' ...\n")
dir.create("data/metadata", recursive = TRUE, showWarnings = FALSE)
output_file <- "data/metadata/pilot_samples.csv"
write.csv(pilot_samples, output_file, row.names = FALSE, quote = FALSE)

cat("SELESAI. File metada berhasil ditulis!\n")
