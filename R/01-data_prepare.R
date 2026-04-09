## ---------------------------- SETUP ------------------------------------------
# Initialize reproducibility state for any downstream stochastic steps.
set.seed(16748991)
seed_before <- .Random.seed

# Load required packages.
# Package installation is handled by renv
packages <- c("tidyverse", "data.table", "DBI", "duckdb", "phiper")

missing_packages <- packages[!packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  stop("Missing packages: ", paste(missing_packages, collapse = ", "))
}

invisible(lapply(packages, library, character.only = TRUE))
rm(packages, missing_packages)

## ---------------------------- PATHS ------------------------------------------
# Input archive containing PRIMM PhIP-Seq measurements and sample metadata.
zip_path <- "/lisc/data/work/ccr/SHARED_RESOURCES/external_data/PRIMM/FinalB_2.zip"

# Output directory for processed tables used by downstream analysis.
out_dir <- file.path("data", "processed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Data input: for Carlos' automated report generator
exist_csv_path    <- file.path(out_dir, "exist.csv")
samples_meta_path <- file.path(out_dir, "samples_meta.csv")
# Metadata samples: for Rmd report
metadata_all_path <- file.path(out_dir, "metadata_all.csv")
metadata_path     <- file.path(out_dir, "metadata.csv")
# Metadata peptides: for phiper analysis
peplib_path       <- file.path(out_dir, "peptide_library.csv")
# Parquet file: for phiper analysis
parquet_path      <- file.path(out_dir, "primm_full.parquet")

## ---------------------------- LOAD PHIP-SEQ DATA -----------------------------
# Read wide-format peptide fold-change matrix from the ZIP archive.
primm_wide <- read.csv(unz(zip_path, "FinalB/Fig1_table1/phipseq_data_fixed_all.csv"))

# Convert to long format expected by downstream processing.
# Add `exist` indicator for peptide presence based on non-zero fold change.
primm <- primm_wide %>%
  pivot_longer(-name, names_to = "peptide_id", values_to = "fold_change") %>%
  rename(sample_id = name) %>%
  mutate(
    fold_change = as.double(fold_change),
    exist = as.integer(!is.na(fold_change) & fold_change != 0)
  ) %>%
  relocate(sample_id, peptide_id, fold_change, exist)

## ---------------------------- LOAD METADATA ----------------------------------
# Read sample-level metadata from the ZIP archive.
metadata <- read.csv(unz(zip_path, "FinalB/Fig1_table1/metadata_fixed_all.csv"))
# The read_csv() function is better at handling booleans
peplib <- read_csv(unz(zip_path, "FinalB/Fig1_table1/oligos_metadata.csv"))
# Extract the peptide library from phiper (downloaded from GitHub)
peplib_vogl <- readRDS("data/raw/combined_library_15.01.26.rds")

# Create one-hot encoded group indicators for response and clinical covariates
# stratified by timepoint. Missing and empty values are excluded so that no
# columns such as antibiotics_NA_T1 are created.
group_vars <- c("response", "ORR", "toxicity", "colitis", "combiIO", "antibiotics", "ppi")

metadata <- metadata %>%
  mutate(response = if_else(PFS12 == "yes", "R", "NR")) %>%
  pivot_longer(
    cols = all_of(group_vars),
    names_to = "group_type",
    values_to = "group_value",
    values_drop_na = TRUE
  ) %>%
  filter(group_value != "") %>%
  mutate(
    group = if_else(
      group_type == "response",
      paste(group_value, timepoint, sep = "_"),
      paste(group_type, group_value, timepoint, sep = "_")
    ),
    in_group = 1L
  ) %>%
  select(-group_type, -group_value) %>%
  pivot_wider(
    names_from = group,
    values_from = in_group,
    values_fill = 0L,
    values_fn = max,
    names_sort = TRUE
  ) %>%
  mutate(response = if_else(PFS12 == "yes", "R", "NR")) %>%
  relocate(response, .before = "timepoint") %>%
  relocate(paste0("R_T", 1:4), .after = "age") %>%
  relocate(paste0("NR_T", 1:4), .after = "R_T4")
  
# Preserve the unfiltered encoded metadata for manual review.
metadata_all <- metadata

# This sample appears twice with completely different metadata
# Need to double-check with collaborators
# MVCC011 T3 doesn't correlate with T1+T2, so I assume it's the other patient
metadata <- metadata %>%
  filter(!(LV_code == "LV1016273557" & patientid == "MVCC011"))

## ---------------------------- JOIN DATA --------------------------------------
# Retain only samples present in both the assay data and metadata.
primm <- primm %>%
  inner_join(metadata, by = c("sample_id" = "LV_code"))

## ---------------------------- EXPORT CSV FILES -------------------------------
# Wide binary peptide-by-sample matrix for Carlos' automated report generator.
primm_exist <- primm %>%
  pivot_wider(
    id_cols = peptide_id,
    names_from = sample_id,
    values_from = exist
  ) %>%
  column_to_rownames("peptide_id")

# Sample metadata table with standardized column names expected downstream.
samples_meta <- metadata %>%
  rename(
    SampleName = LV_code,
    Sex = sex,
    Age = age
  ) %>%
  relocate(Sex, Age, .after = last_col())

# Make peptide metadata compatible with downstream analysis
peplib <- peplib %>%
  rename(peptide_id = oligo) %>%
  mutate(
    is_bac_flagella = case_when(
      is_bac_flagella == FALSE ~ "no",
      is_bac_flagella == TRUE ~ "yes",
      TRUE ~ NA
    )
  ) %>%
  left_join(
    peplib_vogl %>%
      select(peptide_id, phylum, class, order, family, genus, species),
    by = "peptide_id"
  )

# Data input for Carlos' automated report generator
write.csv(primm_exist, exist_csv_path)
write.csv(samples_meta, samples_meta_path, row.names = FALSE)
# Sample metadata for Rmd report
write.csv(metadata_all, metadata_all_path, row.names = FALSE)
write.csv(metadata, metadata_path, row.names = FALSE)
# Peptide metadata
write.csv(peplib, peplib_path, row.names = FALSE)

## ---------------------------- EXPORT PARQUET ---------------------------------
# Export the full joined long-format table to Parquet for phiper-compatible use.
primm_export <- primm

# Fail explicitly if the primary sample-peptide key is not unique.
primm_export <- primm_export %>%
  distinct(sample_id, peptide_id, .keep_all = TRUE)

if (nrow(primm_export) != nrow(primm)) {
  stop("Duplicate sample_id/peptide_id combinations found in primm_export")
}

# Open (or create) a DuckDB database
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:", read_only = FALSE)

# Expose the data frame to DuckDB
duckdb::duckdb_register(con, "df_tbl", primm_export)

# Persist it to Parquet
dbExecute(
  con,
  paste0("COPY df_tbl TO '", parquet_path, "' (FORMAT PARQUET);")
)

# Close connection
dbDisconnect(con, shutdown = TRUE)

## ---------------------------- CLEANUP ----------------------------------------
# Remove large intermediate objects no longer needed in the session.
rm(
  peplib,
  peplib_vogl,
  primm_wide,
  primm,
  primm_exist,
  primm_export,
  con
)
gc()
