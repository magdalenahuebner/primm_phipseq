## ----------------------------READING PACKAGES---------------------------------
# Setting random generator seed
set.seed(16748991)
seed_before <- .Random.seed

# Creating vector of necessary packages
packages <- c(
  "tidyverse",
  "data.table",
  "DBI",
  "duckdb"
)

# Load development version of phiper
library(phiper)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  pak::pkg_install(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

# Removing unnecessary variables
rm(list = c('installed_packages', 'packages'))

## ---------------------LOADING DATA--------------------------------------------
# Read in PhIPSeq data
zip_path <- "L:/lisc/data/work/ccr/SHARED_RESOURCES/external_data/PRIMM/FinalB_2.zip"
primm_wide <- read.csv(unz(zip_path, "FinalB/Fig1_table1/phipseq_data_fixed_all.csv"))

# Convert to long format
primm <- primm_wide %>%
  pivot_longer(-name, names_to = "peptide_id", values_to = "fold_change") %>%
  rename(sample_id = name) %>%
  mutate(fold_change = as.double(fold_change)) %>%
  mutate(exist = if_else(fold_change == 0, 0, 1)) %>%
  relocate(sample_id, peptide_id, fold_change, exist)

# Read the metadata
metadata <- read.csv(unz(zip_path, "FinalB/Fig1_table1/metadata_fixed_all.csv"))

# Add one-hot encoded groups
metadata <- metadata %>%
  mutate(
    response = if_else(PFS12 == "yes", "R", "NR"),
    group = paste0(response, "_", timepoint),
    in_group = 1L
  ) %>%
  pivot_wider(
    names_from = group, 
    values_from = in_group, 
    names_sort = TRUE,
    values_fill = 0L
  ) %>%
  # This sample appears twice with completely different metadata
  # Need to double-check with collaborators
  filter(LV_code != "LV1016273557")

# Join to the PRIMM on sample_id
primm <- primm %>% right_join(metadata, by = c("sample_id" = "LV_code"))

## ---------------------------------- DATA EXPORT ------------------------------
# Export for Carlos' report generator
primm_exist <- primm %>%
  pivot_wider(
    id_cols = "peptide_id", 
    names_from = "sample_id", 
    values_from = "exist"
  ) %>%
  column_to_rownames("peptide_id")

samples_meta <- metadata %>%
  rename(
    SampleName = LV_code,
    Sex = sex,
    Age = age
  ) %>%
  relocate(Sex, Age, .after = last_col())

write.csv(primm_exist, "data/processed/exist.csv")
write.csv(samples_meta, "data/processed/samples_meta.csv", row.names = FALSE)
rm(list = c("primm_exist", "samples_meta"))

# Exporting the data to a .parquet file --> for phiper use
primm_export <- primm

# Safety-check
primm_export <- primm_export %>%
  distinct(sample_id, peptide_id, .keep_all = TRUE)

# 1. Open (or create) a DuckDB database ----------------------------------------
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:", read_only = FALSE)

# 2. Expose the data frame to DuckDB -------------------------------------------
duckdb::duckdb_register(con, "df_tbl", primm_export)

# -- df is now a virtual table called "df_tbl" inside the DB
# 3. Persist it to Parquet -----------------------------------------------------
dbExecute(
  con,
  "COPY df_tbl TO 'data/processed/primm_full.parquet' (FORMAT PARQUET);"
)

# 4. Clean up ------------------------------------------------------------------
dbDisconnect(con, shutdown = TRUE)
rm(list = c("primm_wide", "primm", "primm_export", "metadata", "con"))
gc()
