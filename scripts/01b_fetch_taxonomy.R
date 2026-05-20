#!/usr/bin/env Rscript
# Script: 01b_fetch_taxonomy.R
# Purpose: Fetch species-level relative abundance (MetaPhlAn) from curatedMetagenomicData
#          for the selected samples and save as a clean abundance matrix.

suppressPackageStartupMessages({
  library(curatedMetagenomicData)
  library(dplyr)
  library(tibble)
  library(SummarizedExperiment)
})

cat("========================================================\n")
cat("[Fase 1b] Fetch Taxonomic Abundance (MetaPhlAn)\n")
cat("========================================================\n\n")

# Load full selected samples metadata
full_samples_path <- "data/metadata/full_samples.csv"
sample_map_path   <- "data/metadata/sample_map.csv"

if (!file.exists(full_samples_path) || !file.exists(sample_map_path)) {
  stop("[ERROR] full_samples.csv atau sample_map.csv tidak ditemukan di data/metadata/. Jalankan 01_fetch_metadata.R terlebih dahulu.")
}

full_samples <- read.csv(full_samples_path, stringsAsFactors = FALSE)
sample_map   <- read.csv(sample_map_path, stringsAsFactors = FALSE)

# curatedMetagenomicData query using returnSamples
cat("1. Menarik data kelimpahan taksonomi menggunakan returnSamples...\n")
tse <- tryCatch({
  returnSamples(full_samples, dataType = "relative_abundance", counts = FALSE, rownames = "long")
}, error = function(e) {
  stop("[ERROR] Gagal menarik data dari curatedMetagenomicData: ", e$message)
})

cat("   ✓ Berhasil mendownload data taksonomi.\n\n")

# Ekstrak matriks kelimpahan relatif
tax_mat <- assay(tse)
cat(sprintf("   Jumlah clades taksonomi terdeteksi: %d\n", nrow(tax_mat)))
cat(sprintf("   Jumlah sampel terdeteksi          : %d\n\n", ncol(tax_mat)))

# Filter hanya tingkat spesies (s__ tetapi tidak ada t__)
cat("2. Memfilter data ke tingkat spesies (Species level only)...\n")
species_idx <- grepl("s__[^|]+$", rownames(tax_mat))
species_mat <- tax_mat[species_idx, , drop = FALSE]

# Bersihkan nama spesies: ambil teks setelah 's__'
clean_species_names <- sub(".*s__", "", rownames(species_mat))
rownames(species_mat) <- clean_species_names

# Transpose matriks agar Baris = Sampel, Kolom = Spesies
tax_df <- as.data.frame(t(species_mat)) %>%
  rownames_to_column(var = "original_sample_id")

# Buat mapping antara original_sample_id dan customized sample_id kita
# original_sample_id di tse dicocokkan dengan sample_id di full_samples
# customized sample_id di sample_map dipetakan lewat accession (SRA accession)
mapping_df <- full_samples %>%
  select(original_sample_id = sample_id, accession = NCBI_accession) %>%
  mutate(accession_clean = sub(";.*", "", accession)) %>%
  filter(!is.na(accession_clean) & accession_clean != "") %>%
  left_join(
    sample_map %>% 
      select(sample_id, accession) %>% 
      mutate(accession_clean = sub(";.*", "", accession)) %>%
      filter(!is.na(accession_clean) & accession_clean != ""),
    by = "accession_clean"
  ) %>%
  select(original_sample_id, final_sample_id = sample_id)

# Gabungkan mapping ke tax_df
final_tax_df <- mapping_df %>%
  inner_join(tax_df, by = "original_sample_id") %>%
  select(-original_sample_id) %>%
  rename(Sample_ID = final_sample_id)

# Simpan hasil matriks kelimpahan taksonomi
out_dir <- "data/metadata"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, "taxonomy_abundance.csv")
write.csv(final_tax_df, out_file, row.names = FALSE)
cat(sprintf("✓ Sukses menyimpan matriks kelimpahan taksonomi ke: %s\n", out_file))
