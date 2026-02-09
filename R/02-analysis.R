# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
# Load required R packages
library(phiper)
library(rlang)
library(ggplot2)
library(Cairo)
library(openxlsx)
library(dplyr)
library(purrr)
library(locfdr)

set.seed(632961)

# ------------------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------------------
# Parse CLI inputs and optional parameters (e.g., paths, filters, flags),
# validate values, and populate defaults used throughout the script
# ------------------------------------------------------------------------------
N_CORES <- 6
LOG <- TRUE
LOG_FILE <- NULL
MAX_GB <- 10
FORCE <- TRUE
args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
  if (grepl("=", arg)) {
    parts <- strsplit(arg, "=")[[1]]
    key <- parts[1]; value <- parts[2]
    if (key == "N_CORES") {
      val_num <- suppressWarnings(as.numeric(value))
      if (!is.na(val_num)) N_CORES <- val_num
    } else if (key == "MAX_GB") {
      val_num <- suppressWarnings(as.numeric(value))
      if (!is.na(val_num)) MAX_GB <- val_num
    } else if (key == "LOG") {
      val_log <- tolower(value)
      if (val_log %in% c("true", "t", "1")) {
        LOG <- TRUE
      } else if (val_log %in% c("false", "f", "0")) {
        LOG <- FALSE
      }
    } else if (key == "LOG_FILE") {
      LOG_FILE <- sub("^['\\\"]|['\\\"]$", "", value)
    } else if (key == "FORCE") {
      val_log <- tolower(value)
      if (val_log %in% c("true", "t", "1")) {
        FORCE <- TRUE
      } else if (val_log %in% c("false", "f", "0")) {
        FORCE <- FALSE
      }
    }
  }
}

# ------------------------------------------------------------------------------
# Build PHIP data object (reproducible)
# ------------------------------------------------------------------------------
# Create the PHIP data object from the input files and set a fixed seed to
# ensure reproducible sampling, random splits, and any stochastic steps
# ------------------------------------------------------------------------------
withr::with_preserve_seed({
  ps <- phip_convert(
    data_long_path    = "data/processed/primm_full.parquet",
    sample_id         = "sample_id",
    peptide_id        = "peptide_id",
    exist             = "exist",
    fold_change       = "fold_change",
    counts_input      = NULL,
    counts_hit        = NULL,
    peptide_library   = TRUE,
    materialise_table = TRUE,
    auto_expand       = TRUE,
    n_cores           = 10
  )
})

# ------------------------------------------------------------------------------
# Results directory + peptide library snapshot
# ------------------------------------------------------------------------------
# Create a base results folder and persist the peptide library used in this run
# for provenance and reproducibility
dir.create("results", recursive = TRUE, showWarnings = FALSE)

peptide_library <- ps %>%
  get_peptide_library() %>%
  collect() %>%
  as.data.frame()

saveRDS(peptide_library, file = file.path("results", "peptide_library.rds"))

# ------------------------------------------------------------------------------
# Analysis setup
# ------------------------------------------------------------------------------
# List of comparisons (each entry is a pair of group labels)
comparisons <- list(
  c("R_T1", "NR_T1")
  # c("R_T2", "NR_T2"),
  # c("R_T3", "NR_T3"),
  # c("R_T4", "NR_T4"),
  # c("R_T1", "R_T2"),
  # c("NR_T1", "NR_T2")
)

# Columns that are always retained when exporting data
base_cols <- c("sample_id", "peptide_id", "group_char", "exist")

# ------------------------------------------------------------------------------
# I/O helpers
# ------------------------------------------------------------------------------
# Save an R object to an RDS file, creating parent directories if needed
save_rds_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file = path)
}

# Select requested columns (base + extras), optionally drop rows with NA in
# selected columns, collect to R, and save to RDS
make_and_save <- function(data, out_path, extra_vars = NULL, drop_na = TRUE) {
  cols_to_select <- unique(c(base_cols, extra_vars %||% character(0)))

  avail_cols <- names(data) # MH: THIS IS ALWAYS NULL
  missing_cols <- setdiff(cols_to_select, avail_cols)
  if (length(missing_cols) > 0) {
    message("Skipping missing columns: ", paste(missing_cols, collapse = ", "))
  }

  df <- data %>%
    dplyr::select(dplyr::any_of(cols_to_select))

  if (isTRUE(drop_na)) {
    cols_present <- intersect(cols_to_select, names(df))
    df <- df %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(cols_present), ~ !is.na(.x)))
  }

  df <- df %>% dplyr::collect()
  save_rds_safe(df, out_path)
  invisible(df)
}

# ------------------------------------------------------------------------------
# Parallel backend (DELTA)
# ------------------------------------------------------------------------------
# Configure BLAS/OpenMP to single-threaded mode to avoid CPU oversubscription
# when running multiple parallel R workers. Then set up a `future` plan that
# works across platforms (multisession on Windows; multicore/sequential on Unix)
Sys.setenv(
  OMP_NUM_THREADS      = "1",
  MKL_NUM_THREADS      = "1",
  OPENBLAS_NUM_THREADS = "1"
)

options(
  future.globals.maxSize = MAX_GB * 1024^3,
  future.scheduling      = 1
)

# Store the current plan so it can be restored later.
original_plan <- future::plan()

if (.Platform$OS.type == "windows") {
  future::plan(future::multisession, workers = N_CORES)
} else if (N_CORES > 1L) {
  future::plan(future::multicore, workers = N_CORES)
} else {
  future::plan(future::sequential)
}

# ------------------------------------------------------------------------------
# Comparisons loop
# ------------------------------------------------------------------------------
# Run the same pipeline for each (var1 vs var2) comparison:
#   1) Filter + export comparison subset
#   2) Enrichment counts
#   3) Alpha diversity
#   4) Beta diversity (distances, ordinations, PERMANOVA, dispersion, t-SNE)
#   5) POP framework (prevalence comparison + plots)
#   6) DELTA framework (permutation-based differential prevalence) MH: SKIP
for (cmp in comparisons) {
  var1 <- cmp[[1]]
  var2 <- cmp[[2]]

  message("Running comparison: ", var1, " vs ", var2)

  label_dir <- paste(var1, "vs", var2, sep = "_")
  out_dir   <- file.path("results", label_dir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ----------------------------------------------------------------------------
  # Filter to the two groups and add analysis-friendly grouping columns
  # ----------------------------------------------------------------------------
  ps_cmp <- ps %>%
    dplyr::filter((.data[[var1]] == 1L) | (.data[[var2]] == 1L)) %>%
    dplyr::mutate(
      group_char  = dplyr::if_else(.data[[var1]] == 1L, var1, var2),
      group_dummy = dplyr::if_else(.data[[var1]] == 1L, 1L, 0L)
    )

  # Persist the filtered dataset used for downstream frameworks.
  make_and_save(
    data       = ps_cmp,
    out_path   = file.path(out_dir, paste0(label_dir, "_data.rds")),
    extra_vars = "fold_change"
  )

  # ----------------------------------------------------------------------------
  # Enrichment counts
  # ----------------------------------------------------------------------------
  CairoSVG(file.path(out_dir, "enrichment_counts.svg"), dpi = 300,
           height = 30, width = 40, unit = "cm", bg = "white")
  p_enrich <- plot_enrichment_counts(ps_cmp, group_cols = "group_char") +
    # theme(text = element_text(family = "DejaVu Sans")) +
    facet_wrap(
      vars(forcats::fct_relevel(as.factor(Cohort), "R_T1", "NR_T1")),
      ncol = 2,
      scales = "free_x",
      drop = FALSE
    )
  print(p_enrich)
  dev.off()

  # ----------------------------------------------------------------------------
  # Alpha diversity
  # ----------------------------------------------------------------------------
  alpha_div <- compute_alpha_diversity(ps_cmp,
                                       group_cols = "group_char",
                                       carry_cols = c("sex", "age"))

  dir.create(file.path(out_dir, "alpha_diversity"),
             recursive = TRUE, showWarnings = FALSE)

  write.xlsx(alpha_div, file.path(out_dir, "alpha_diversity", "table.xlsx"))

  # Richness
  CairoSVG(file.path(out_dir, "alpha_diversity", "plot.svg"), dpi = 300,
           height = 30, width = 40, unit = "cm", bg = "white")

  p_alpha_r <- plot_alpha_diversity(alpha_div,
                                    metric = "richness",
                                    group_col = "group_char",
                                    x_order = cmp)
  print(p_alpha_r)
  dev.off()
  
  # Shannon diversity
  CairoSVG(file.path(out_dir, "alpha_diversity", "plot_shannon.svg"), dpi = 300,
           height = 30, width = 40, unit = "cm", bg = "white")
  
  p_alpha_sh <- plot_alpha_diversity(alpha_div,
                                     metric = "shannon_diversity",
                                     group_col = "group_char",
                                     x_order = cmp)
  print(p_alpha_sh)
  dev.off()

  # ----------------------------------------------------------------------------
  # Beta diversity
  # ----------------------------------------------------------------------------
  dist_bc <- phiper:::compute_distance(ps_cmp, value_col = "exist",
                                       method_normalization = "hellinger",
                                       distance = "bray", n_threads = 10)

  dir.create(file.path(out_dir, "beta_diversity"),
             recursive = TRUE,
             showWarnings = FALSE)

  dist_mat <- as.matrix(dist_bc)
  openxlsx::write.xlsx(dist_mat,
                       file = file.path(out_dir, "beta_diversity",
                                        "distance_matrix.xlsx"),
                       rowNames = TRUE)

  pcoa_res <- phiper:::compute_pcoa(dist_bc,
                                    neg_correction = "none",
                                    n_axes = 109)
  saveRDS(pcoa_res, file.path(out_dir, "beta_diversity", "pcoa_results.rds"))

  cap_res <- phiper:::compute_capscale(dist_bc,
                                       ps = ps_cmp,
                                       formula = ~ group_char)
  saveRDS(cap_res, file.path(out_dir, "beta_diversity", "capscale_results.rds"))

  permanova_res <- phiper:::compute_permanova(dist_bc,
                                              ps = ps_cmp,
                                              group_col = "group_char")
  saveRDS(permanova_res, file.path(out_dir, "beta_diversity",
                                   "permanova_results.rds"))

  disp_res <- phiper:::compute_dispersion(dist_bc,
                                          ps = ps_cmp,
                                          group_col = "group_char")
  saveRDS(disp_res, file.path(out_dir, "beta_diversity",
                              "dispersion_results.rds"))
  print(disp_res)

  tsne_res <- phiper:::compute_tsne(ps = ps_cmp,
                                    dist_obj = dist_bc,
                                    dims = 2L,
                                    perplexity = 15,
                                    meta_cols = c("group_char"))
  openxlsx::write.xlsx(tsne_res, file = file.path(out_dir,
                                                  "beta_diversity",
                                                  "tsne2d_results.xlsx"),
                       rowNames = TRUE)

  CairoSVG(file.path(out_dir, "beta_diversity", "tsne2d_plot.svg"), dpi = 300,
           height = 30, width = 30, unit = "cm", bg = "white")
  p_tsne2d <- phiper:::plot_tsne(
    tsne_res %>% mutate(group_char = factor(group_char, levels = cmp)),
    view = "2d",
    colour = "group_char",
    palette = phip_palette[1:2]
  )
  print(p_tsne2d)
  dev.off()

  tsne_res <- phiper:::compute_tsne(ps = ps_cmp,
                                    dist_obj = dist_bc,
                                    dims = 3L,
                                    perplexity = 20,
                                    meta_cols = c("group_char"))
  openxlsx::write.xlsx(tsne_res, file = file.path(out_dir,
                                                  "beta_diversity",
                                                  "tsne3d_results.xlsx"),
                       rowNames = TRUE)

  p3d <- phiper:::plot_tsne(
    tsne_res %>% mutate(group_char = factor(group_char, levels = cmp)),
    view = "3d",
    colour = "group_char",
    palette = phip_palette[1:2]
  )
  htmlwidgets::saveWidget(p3d, file = file.path(out_dir, "beta_diversity",
                                                "tsne3d_plot.html"),
                          selfcontained = TRUE)

  # Add group information to PCoA sample coordinates
  pcoa_res$sample_coords <- pcoa_res$sample_coords %>%
    dplyr::left_join(
      ps_cmp$data_long %>%
        dplyr::select(sample_id, group_char) %>%
        dplyr::distinct(),
      by = "sample_id",
      copy = TRUE
    )
  pcoa_res$sample_coords <- pcoa_res$sample_coords %>%
    mutate(group_char = factor(group_char, levels = cmp))

  # PCoA plot with group centroids and ellipses
  CairoSVG(file.path(out_dir, "beta_diversity", "pcoa_plot.svg"), dpi = 300,
           height = 30, width = 30, unit = "cm", bg = "white")
  p_pcoa <- phiper:::plot_pcoa(
    pcoa_res,
    axes = c(1, 2),
    group_col = "group_char",
    ellipse_by = "group",
    show_centroids = TRUE
  ) +
    scale_colour_manual(values = setNames(phip_palette[1:2], cmp))
  print(p_pcoa)
  dev.off()

  # Scree plot for first 15 axes of PCoA
  CairoSVG(file.path(out_dir, "beta_diversity", "scree_plot.svg"), dpi = 300,
           height = 30, width = 30, unit = "cm", bg = "white")
  p_scree <- phiper:::plot_scree(pcoa_res, n_axes = 15, type = "line") +
    theme()
  print(p_scree)
  dev.off()

  # Determine which contrast label actually exists in the dispersion object
  available_contrasts <- unique(disp_res$distances$contrast)

  # Preferred pairwise label: var1 vs var2 (e.g. "control vs MCI")
  pair_contrast <- paste(var1, "vs", var2)

  if (pair_contrast %in% available_contrasts) {
    contrast_to_use <- pair_contrast
  } else if ("<global>" %in% available_contrasts) {
    # fallback: use global dispersion if pairwise distances are not stored
    contrast_to_use <- "<global>"
  } else {
    # last resort: just take the first available contrast and warn
    contrast_to_use <- available_contrasts[1]
    message(
      "Warning: requested contrast '", pair_contrast,
      "' not found in disp_res$distances$contrast. Using '",
      contrast_to_use, "' instead."
    )
  }

  disp_res$distances <- disp_res$distances %>%
    mutate(level = factor(level, levels = cmp))
  
  CairoSVG(file.path(out_dir, "beta_diversity", "dispersion_plot.svg"),
           dpi = 300, height = 30, width = 30, unit = "cm", bg = "white")
  p_disp <- phiper:::plot_dispersion(
    disp_res,
    scope        = "group",
    contrast     = contrast_to_use,
    show_violin  = TRUE,
    show_box     = TRUE,
    show_points  = TRUE
  ) +
    scale_colour_manual(values = setNames(phip_palette[1:2], cmp)) +
    scale_fill_manual(values = setNames(phip_palette[1:2], cmp))
  
  print(p_disp)
  dev.off()

  # ----------------------------------------------------------------------------
  # POP framework
  # ----------------------------------------------------------------------------
  data_frameworks <- readRDS(file.path(out_dir, paste0(label_dir, "_data.rds")))
  data_frameworks$group_char <- factor(data_frameworks$group_char, levels = cmp)
  peplib <- readRDS(file.path("results", "peptide_library.rds"))

  extract_tbl <- function(obj) {
    if (is.data.frame(obj)) {
      return(tibble::as_tibble(obj))
    }
    for (nm in c("data", "table", "tbl", "df", "result", "results")) {
      if (!is.null(obj[[nm]])) {
        return(tibble::as_tibble(obj[[nm]]))
      }
    }
    out <- try(tibble::as_tibble(obj), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }
    out <- try(as.data.frame(obj), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(tibble::as_tibble(out))
    }
    stop("Cannot extract a data table from the POP result object.")
  }

  dir.create(file.path(out_dir, "POP_framework"), recursive = TRUE,
             showWarnings = FALSE)
  
  pep_tbl_file <- file.path(out_dir, "POP_framework", "single_peptide.csv")
  if(!file.exists(pep_tbl_file) | FORCE) {
    prev_res_pep <- phiper::ph_prevalence_compare(
      x                 = data_frameworks,
      group_cols        = "group_char",
      rank_cols         = "peptide_id",
      compute_ratios_db = TRUE,
      parallel          = TRUE,
      collect           = TRUE
    )
    pep_tbl <- extract_tbl(prev_res_pep)
    write.csv(pep_tbl, pep_tbl_file)
  } else {
    pep_tbl <- read.csv(pep_tbl_file)
  }
  
  ranks_tax <- c("phylum", "class", "order", "family", "genus", "species")
  
  rank_tbl_file <- file.path(out_dir, "POP_framework", "taxa_ranks.csv")
  if(!file.exists(pep_tbl_file) | FORCE) {
    prev_res_rank <- phiper::ph_prevalence_compare(
      x                 = data_frameworks,
      group_cols        = "group_char",
      rank_cols         = ranks_tax,
      compute_ratios_db = FALSE,
      parallel          = TRUE,
      peptide_library   = peplib,
      collect           = TRUE
    )
    rank_tbl <- extract_tbl(prev_res_rank)
    write.csv(rank_tbl, rank_tbl_file)
  } else {
    rank_tbl <- read.csv(rank_tbl_file)
  }

  ranks_combined <- c(ranks_tax, "peptide_id")
  plots_dir <- file.path(out_dir, "POP_framework", "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

  for (rank_name in ranks_combined) {
    rank_chr <- as.character(rank_name)
    out_name <- file.path(plots_dir, rank_chr)
    df_rank <- if (rank_chr == "peptide_id") {
      pep_tbl
    } else {
      rank_tbl %>% filter(rank == rank_chr)
    }
    p_static <- scatter_static(
      df   = df_rank,
      rank = rank_chr,
      xlab = df_rank$group1[1],
      ylab = df_rank$group2[1],
      color_by = "is_flagellum",
      color_title = "Flagellins",
      point_size       = 2,
      jitter_width_pp  = 0.15,
      jitter_height_pp = 0.15,
      point_alpha      = 0.85,
      font_size        = 12
    ) +
      ggplot2::coord_cartesian(xlim = c(-2, 102), ylim = c(-2, 102),
                               expand = TRUE) +
      ggplot2::theme(
        plot.margin = grid::unit(c(12, 12, 12, 12), "pt"),
        text        = ggplot2::element_text(family = "Montserrat")
      )
    ggsave(paste0(out_name, "_static.svg"), p_static, dpi = 300,
           height = 30, width = 30, unit = "cm", bg = "white")

    p_inter <- scatter_interactive(
      df   = df_rank,
      rank = rank_chr,
      xlab = df_rank$group1[1],
      ylab = df_rank$group2[1],
      peplib = peplib,
      point_size = 10,
      jitter_width_pp  = 0.25,
      jitter_height_pp = 0.25,
      point_alpha = 0.85,
      font_size = 12
    )

    p_inter <- plotly::layout(
      p_inter,
      autosize = TRUE,
      margin   = list(l = 70, r = 30, t = 10, b = 70),
      xaxis    = list(range = c(-2, 102), automargin = TRUE),
      yaxis    = list(range = c(-2, 102), automargin = TRUE)
    )
    htmlwidgets::saveWidget(
      p_inter,
      file = paste0(out_name, "_interactive.html"),
      selfcontained = TRUE
    )
  }
}

# ------------------------------------------------------------------------------
# Restore original future plan
# ------------------------------------------------------------------------------
future::plan(original_plan)
