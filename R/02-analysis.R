# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------
# Load required R packages
library(phiper)
library(tidyverse)
library(Cairo)
library(ggpubr)
library(locfdr)
library(mgcv)
library(openxlsx)
library(patchwork)
library(rlang)
library(showtext)
library(vegan)

set.seed(632961)

source(file.path("R", "helpers.R"))

# Add font + phiper use font in all plots
font_add_google("Montserrat", "monte")
phip_use_montserrat()
showtext_auto()

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
    peptide_library   = FALSE,  # using library from Weizmann Institute
    materialise_table = TRUE,
    auto_expand       = TRUE,
    n_cores           = 10
  )
})

# Adding custom peptide library
peplib <- read.csv("data/processed/peptide_library.csv")

con <- ps$data_long %>% dbplyr::remote_con()
ps$peptide_library <- copy_to(
  dest = con,
  df = peplib,
  name = "peptide_library",
  temporary = FALSE,
  overwrite = TRUE
)
rm(con)

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
comp_main <- list(
  c("R_T0", "NR_T0"),
  c("R_T1", "NR_T1"),
  c("R_T2", "NR_T2"),
  c("R_T3", "NR_T3"),
  c("R_T0", "R_T1"),
  c("NR_T0", "NR_T1")
)

# Specify timepoint comparisons
timepoints <- paste0("T", 0:3)
tp_pairs <- combn(timepoints, 2, simplify = FALSE)
comp_tps <- c(
  lapply(tp_pairs, \(p) paste0("R_",  p)),
  lapply(tp_pairs, \(p) paste0("NR_", p)),
  list(
    c("R_T0", "NR_T1"),
    c("NR_T0", "R_T1")
  )
)

# Add more comparisons:
# Vector of group variables and timepoints
group_vars <- c("ORR", "toxicity", "colitis", "combiIO", "antibiotics", "ppi")
timepoints <- c("T0", "T1","T2","T3")
levels2 <- c("yes","no")

# All group labels that to consider
labs <- expand_grid(
  group_type = group_vars,
  group_value = levels2,
  timepoint = timepoints
) %>%
  mutate(lab = str_c(group_type, group_value, timepoint, sep = "_"))

# 1) within-timepoint: yes vs. no
cmp_within <- labs %>%
  select(group_type, group_value, timepoint, lab) %>%
  pivot_wider(names_from = group_value, values_from = lab) %>%
  transmute(cmp = map2(yes, no, c))

# Leave for now, maybe add back later for supplementary figures?
# # 2) across-timepoint within level: yes_Ti vs. yes_Tj AND no_Ti vs. no_Tj
# cmp_across <- labs %>%
#   group_by(group_type, group_value) %>%
#   summarise(cmp = list(combn(lab, 2, simplify = FALSE)), .groups = "drop") %>%
#   unnest(cmp)
# comp_add <- c(cmp_within$cmp, cmp_across$cmp)
# comparisons <- c(comp_main, comp_add)

comparisons <- unique(c(comp_main, comp_tps, cmp_within$cmp))

# Columns that are always retained when exporting data
base_cols <- c("sample_id", "peptide_id", "group_char", "exist")

# Colours scheme for plotting
labels <- unique(unlist(comparisons))
group_palette <- setNames(
  ifelse(
    grepl("(^|_)NR(_|$)|(^|_)no(_|$)", labels),
    phip_palette[9],
    ifelse(
      grepl("(^|_)R(_|$)|(^|_)yes(_|$)", labels),
      phip_palette[1],
      NA_character_
    )
  ),
  labels
)
group_palette[["NR"]] <- phip_palette[9]
group_palette[["R"]] <- phip_palette[1]

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
# Run the same pipeline for each (var1 vs. var2) comparison:
#   1) Filter + export comparison subset
#   2) Enrichment counts
#   3) Alpha diversity
#   4) Beta diversity (distances, ordinations, PERMANOVA, dispersion, t-SNE)
#   5) POP framework (prevalence comparison + plots)
#   6) DELTA framework (permutation-based differential prevalence)
for (cmp in comparisons) {
  var1 <- cmp[[1]]
  var2 <- cmp[[2]]

  message("Running comparison: ", var1, " vs. ", var2)

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
  # MH: Formatting of y axis, log10, bold -> phiper issue
  # MH: Facet wrap in order of vars -> phiper issue
  p_enrich <- plot_enrichment_counts(
    ps_cmp, 
    group_cols = "group_char",
    custom_colors = group_palette
  ) +
    facet_wrap(
      vars(forcats::fct_relevel(as.factor(Cohort), var1, var2)),
      ncol = 2,
      scales = "free_x"
    ) +
    theme(text = element_text(size = 12)) +
    labs(title = NULL)
  
  ggsave(
    filename = file.path(out_dir, "enrichment_counts.pdf"),
    plot = p_enrich,
    width = 18, height = 9, units = "cm",
    device = cairo_pdf, bg = "white"
  )

  # ----------------------------------------------------------------------------
  # Alpha diversity
  # ----------------------------------------------------------------------------
  alpha_div <- compute_alpha_diversity(
    ps_cmp,
    group_cols = "group_char",
    carry_cols = c("sex", "age")
  )
  alpha_div <- alpha_div$group_char

  dir.create(file.path(out_dir, "alpha_diversity"), recursive = T, showWarnings = F)
  write.xlsx(alpha_div, file.path(out_dir, "alpha_diversity", "table.xlsx"))

  # Alpha diversity plot (Richness + Shannon)
  metrics <- c("richness", "shannon_diversity")

  p_alpha <- lapply(metrics, function(m) {
    # Calculate and format p-values first
    pv <- compare_means(
      formula = reformulate("group_char", response = m),
      data = alpha_div,
      method = "wilcox.test",
      comparisons = list(cmp)
    ) %>%
      mutate(
        p = format_pval(p),
        y.position = max(alpha_div[[m]], na.rm = TRUE) * 1.1
      )
    
    p_alpha <- phiper::plot_alpha_diversity(
      alpha_div,
      metric = m,
      group_col = "group_char",
      x_order = cmp
    ) +
      stat_pvalue_manual(
        pv,
        label = "p",
        y.position = "y.position",
        vjust = -0.3,  # increase space between bracket and text
        tip.length = 0.02,
        label.size = 4.23,
        family = "Montserrat"
      ) +
      scale_fill_manual(values = group_palette) +
      coord_cartesian(clip = "off") + 
      theme(
        text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "plain"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")
      )
    
    if (m == "richness") {
      p_alpha <- p_alpha + labs(y = "# of significantly bound peptides")
    }
    
    p_alpha
  })
  p_alpha <- wrap_plots(p_alpha)

  ggsave(
    filename = file.path(out_dir, "alpha_diversity", "plot_alpha_diversity.pdf"),
    plot = p_alpha,
    width = 18, height = 9, units = "cm",
    device = cairo_pdf, bg = "white"
  )
  
  # ----------------------------------------------------------------------------
  # Beta diversity
  # ----------------------------------------------------------------------------
  dist_bc <- phiper:::compute_distance(
    ps_cmp, value_col = "exist",
    method_normalization = "auto",
    distance = "bray", n_threads = 10
  )

  dir.create(file.path(out_dir, "beta_diversity"), recursive = T, showWarnings = F)
  dist_mat <- as.matrix(dist_bc)
  openxlsx::write.xlsx(
    dist_mat,
    file = file.path(out_dir, "beta_diversity", "distance_matrix.xlsx"),
    rowNames = TRUE
  )

  pcoa_res <- phiper:::compute_pcoa(
    dist_bc,
    neg_correction = "none",
    n_axes = 109
  )
  saveRDS(pcoa_res, file.path(out_dir, "beta_diversity", "pcoa_results.rds"))

  permanova_res <- phiper:::compute_permanova(
    dist_bc,
    ps = ps_cmp,
    group_col = "group_char"
  )
  saveRDS(permanova_res, file.path(out_dir, "beta_diversity", "permanova_results.rds"))

  disp_res <- phiper:::compute_dispersion(
    dist_bc,
    ps = ps_cmp,
    group_col = "group_char"
  )
  saveRDS(disp_res, file.path(out_dir, "beta_diversity", "dispersion_results.rds"))

  # Add group information to PCoA sample coordinates
  pcoa_res$sample_coords <- pcoa_res$sample_coords %>%
    dplyr::left_join(
      ps_cmp$data_long %>%
        dplyr::select(sample_id, group_char) %>%
        dplyr::distinct(),
      by = "sample_id",
      copy = TRUE
    ) %>%
    mutate(group_char = factor(group_char, levels = cmp))

  lab_perm <- paste0("PERMANOVA p = ", format_pval(permanova_res$p_adjust))
  lab_disp <- paste0( "Dispersion   p = ", format_pval(disp_res$tests$p_adjust))
  
  # PCoA plot with group centroids and ellipses
  p_pcoa <- phiper:::plot_pcoa(
    pcoa_res,
    axes = c(1, 2),
    group_col = "group_char",
    show_ellipses = FALSE,
    show_centroids = TRUE,
    point_size = 2
  ) +
    stat_ellipse(
      aes(group = group_char),
      colour = "grey70",
      linetype = "dashed"
    ) +
    scale_colour_manual(values = group_palette) +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      axis.title = element_text(face = "plain"),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")
    ) +
    annotate(
      "text",
      x = Inf,
      y = -Inf,
      label = paste(lab_perm, lab_disp, sep = "\n"),
      hjust = 1,  # nudge left a bit
      vjust = -0.2,  # nudge up a bit
      size = 4.23
    )
  
  p_pcoa_noleg <- p_pcoa + theme(legend.position = "none")
  
  ggsave(
    filename = file.path(out_dir, "beta_diversity", "pcoa_plot.pdf"),
    plot = p_pcoa,
    width = 9, height = 9, units = "cm",
    device = cairo_pdf, bg = "white"
  )
  ggsave(
    filename = file.path(out_dir, "beta_diversity", "pcoa_plot_noleg.pdf"),
    plot = p_pcoa_noleg,
    width = 9, height = 9, units = "cm",
    device = cairo_pdf, bg = "white"
  )
  
  # Scree plot for first 15 axes of PCoA
  p_scree <- phiper:::plot_scree(pcoa_res, n_axes = 15, type = "line")
  
  ggsave(
    filename = file.path(out_dir, "beta_diversity", "scree_plot.pdf"),
    plot = p_scree,
    width = 9, height = 9, units = "cm",
    device = cairo_pdf, bg = "white"
  )

  # Pearson correlation of log(fold change)

  if (any(vapply(comp_tps, identical, logical(1), cmp))) {
    # Calculate correlation of fold changes between different timepoints
    resp1 <- ifelse(str_starts(var1, "NR"), "NR", "R")
    resp2 <- ifelse(str_starts(var2, "NR"), "NR", "R")
    tp1 <- str_extract(var1, "T[0-9]+")
    tp2 <- str_extract(var2, "T[0-9]+")

    ps_1 <- ps %>%
      filter(response == resp1, timepoint == tp1) %>%
      select(patientid, peptide_id, fold_change, exist)
    ps_2 <- ps %>%
      filter(response == resp2, timepoint == tp2) %>%
      select(patientid, peptide_id, fold_change, exist)

    # Fold change matrices (Pearson)
    mat1_fc <- ps_1 %>%
      distinct() %>%
      arrange(patientid) %>%
      collect() %>%
      pivot_wider(
        id_cols = peptide_id, 
        names_from = patientid, 
        values_from = fold_change
      ) %>%
      column_to_rownames("peptide_id") %>%
      as.matrix()
    mat2_fc <- ps_2 %>%
      distinct() %>%
      arrange(patientid) %>%
      collect() %>%
      pivot_wider(
        id_cols = peptide_id, 
        names_from = patientid, 
        values_from = fold_change
      ) %>%
      column_to_rownames("peptide_id") %>%
      as.matrix()

    # Align peptides across comparisons
    common_peptides <- intersect(rownames(mat1_fc), rownames(mat2_fc))
    mat1_fc <- mat1_fc[common_peptides, ]
    mat2_fc <- mat2_fc[common_peptides, ]

    # Compute matrices
    pearson <- cor(mat1_fc, mat2_fc, use = "pairwise.complete.obs")

    # Convert to long table
    cor_mat <- pearson %>%
      as.data.frame() %>%
      rownames_to_column("pid_t1") %>%
      pivot_longer(-pid_t1, names_to = "pid_t2", values_to = "pearson")

    cor_mat_file <- file.path(out_dir, "beta_diversity", "cor_matrix_pearson.csv")
    write.csv(cor_mat, file = cor_mat_file, row.names = FALSE)
  }

  # Correlation coefficients (Histogram)

  if (any(vapply(comp_tps[1:12], identical, logical(1), cmp))) {
    cor_mat <- read.csv(cor_mat_file)

    # Keep only patients occurring in both timepoints
    common_patients <- intersect(unique(cor_mat$pid_t1), unique(cor_mat$pid_t2))

    # Extract correlations of matched/unmatched samples
    matched <- cor_mat %>%
      filter(pid_t1 == pid_t2) %>%
      pull(pearson)
    
    unmatched <- cor_mat %>%
      filter(pid_t1 %in% common_patients, pid_t2 %in% common_patients) %>%
      filter(pid_t1 != pid_t2) %>%
      pull(pearson)

    scale_fct <- 10
    bins <- 40
    group <- if_else(resp1 == "R", "Responders", "Non-responders")
    col_pal = c(um = "grey70", m = group_palette[[var1]])
    
    df_um <- data.frame(r = unmatched, type = "um")
    df_m <- data.frame(r = matched, type = "m")
    
    p_corr_hist <- ggplot() +
      geom_histogram(
        data = df_um,
        aes(x = r, y = after_stat(count), fill = type),
        bins = bins,
        alpha = 0.8
      ) +
      geom_histogram(
        data = df_m,
        aes(x = r, y = after_stat(count) * scale_fct, fill = type),
        bins = bins,
        alpha = 0.8
      ) +
      scale_y_continuous(
        name = "Unmatched individuals correlation",
        sec.axis = sec_axis(
          ~ . / scale_fct, 
          name = "Matched individuals correlation"
        )
      ) +
      scale_fill_manual(values = col_pal) +
      labs(
        subtitle = group,
        x = paste0("Pearson correlation of \nlog(fold change) ", tp1, " vs. ", tp2)
      ) +
      theme_phip(12) +
      theme(
        text = element_text(size = 12, family = "Montserrat"),
        plot.subtitle = element_text(size = 12, hjust = 0),
        axis.title = element_text(face = "plain"),
        axis.title.y.right = element_text(margin = margin(l = 8), angle = 90, vjust = 0.5),
        axis.title.y.left  = element_text(margin = margin(r = 8)),
        panel.grid = element_blank(),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt"),
        legend.position = "none"
      )
    
    ggsave(
      filename = file.path(out_dir, "beta_diversity", "corr_histogram_plot.pdf"),
      plot = p_corr_hist,
      width = 9, height = 9, units = "cm",
      device = cairo_pdf, bg = "white"
    )
  }

  # ----------------------------------------------------------------------------
  # POP framework
  # ----------------------------------------------------------------------------
  data_frameworks <- readRDS(file.path(out_dir, paste0(label_dir, "_data.rds")))
  data_frameworks$group_char <- factor(data_frameworks$group_char, levels = cmp)
  peptide_library <- readRDS(file = file.path("results", "peptide_library.rds"))

  dir.create(file.path(out_dir, "POP_framework"), recursive = T, showWarnings = F)
  
  pep_tbl_rds <- file.path(out_dir, "POP_framework", "single_peptide.rds")
  if(!file.exists(pep_tbl_rds) | FORCE) {
    pep_tbl <- phiper::ph_prevalence_compare(
      x                 = data_frameworks,
      group_cols        = "group_char",
      rank_cols         = "peptide_id",
      compute_ratios_db = TRUE,
      parallel          = TRUE,
      peptide_library   = peptide_library,
      collect           = TRUE
    )
    saveRDS(pep_tbl, pep_tbl_rds)
  } else {
    pep_tbl <- readRDS(pep_tbl_rds)
  }
  
  ranks_tax <- c("phylum", "class", "order", "family", "genus", "species")
  
  rank_tbl_file <- file.path(out_dir, "POP_framework", "taxa_ranks.csv")
  if(!file.exists(rank_tbl_file) | FORCE) {
    prev_res_rank <- phiper::ph_prevalence_compare(
      x                 = data_frameworks,
      group_cols        = "group_char",
      rank_cols         = ranks_tax,
      compute_ratios_db = FALSE,
      parallel          = TRUE,
      peptide_library   = peptide_library,
      collect           = TRUE
    )
    rank_tbl <- extract_tbl(prev_res_rank)
    write.csv(rank_tbl, rank_tbl_file)
  } else {
    rank_tbl <- read.csv(rank_tbl_file)
  }
  
  plots_dir <- file.path(out_dir, "POP_framework", "plots")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  
  out_name <- file.path(plots_dir, "peptide_id")
  
  # ----------------------------------------------------------------------------
  # MW test (all categories)
  # ----------------------------------------------------------------------------
  # Apply filter
  pep_tbl_prev10 <- pep_tbl %>% filter((prop1 + prop2 >= 0.1))
  
  # Add peptide library
  pep_tbl_peplib <- pep_tbl_prev10 %>%
    left_join(
      peptide_library,
      by = c("feature" = "peptide_id")
    ) %>%
    mutate(is_library = 1)

  N1 <- unique(pep_tbl_peplib$N1)
  N2 <- unique(pep_tbl_peplib$N2)
  
  # Prepare long format for boxplots by annotation
  family_cols <- c(
    "is_library"       = "Complete library",
    "is_MPA"           = "Microbiome",
    "is_patho"         = "Pathogenic strains",
    "is_IgA"           = "IgA",
    "is_probio"        = "Probiotic",
    "is_bac_flagella"  = "Flagellins",
    "is_infect"        = "VFDB",
    "is_phage"         = "Phages",
    "is_allergens"     = "Allergens",
    "is_IEDB_or_cntrl" = "IEDB controls"
  )
  
  pep_tbl_peplib_long <- pep_tbl_peplib %>%
    select(feature, ratio, all_of(names(family_cols))) %>%
    pivot_longer(
      cols = all_of(names(family_cols)),
      names_to = "family",
      values_to = "present"
    ) %>%
    filter(present == 1) %>%
    mutate(family_lab = factor(family_cols[family], levels = family_cols))
  
  # Set colours
  peptide_palette <- c(
    "Complete library"   = "#9467bd",
    "Microbiome"         = "#ff9896",
    "Pathogenic strains" = "#d62728",
    "IgA"                = "#8c564b",
    "Probiotic"          = "#98df8a",
    "Flagellins"         = "#2ca02c",
    "VFDB"               = "#e377c2",
    "Phages"             = "#bcbd22",
    "Allergens"          = "#17becf",
    "IEDB controls"      = "#7f7f7f"
  )
  
  # MW test (all categories)
  pval_mw_df <- compare_means(
    ratio ~ family_lab,
    data = pep_tbl_peplib_long,
    method = "wilcox.test",
    p.adjust.method = "bonferroni"
  ) %>%
  rstatix::add_significance("p.adj") %>%
  select(-p.format)
  
  pval_mw_file <- file.path(out_dir, "POP_framework", "pval_mw_df.csv")
  write.csv(pval_mw_df, file = pval_mw_file, row.names = FALSE)
  
  ## Flagellins scatter plot
  # Extract adjusted p-value
  pval_mw <- pval_mw_df %>% 
    filter(group1 == "Complete library", group2 == "Flagellins") %>%
    pull(p.adj)
  
  # Removing datapoints where is_flagellum is NA, so NAs won't be plotted
  # Alternatively could set NAs to FALSE in peptide_library to plot NAs too
  keep_flag <- peptide_library %>% filter(!is.na(is_bac_flagella))
  pep_tbl_flag <- pep_tbl_prev10 %>% filter(feature %in% keep_flag$peptide_id)
  
  legend_lab <- paste0("Flagellins (p = ", format_pval(pval_mw), ")")
  
  p_static <- scatter_static(
    df = pep_tbl_flag,
    rank = "peptide_id",
    xlab = var1,
    ylab = var2,
    color_by = "is_bac_flagella",
    color_title = NULL,
    point_size = 2,
    point_alpha = 0,
    font_size = 12
  ) +
    coord_cartesian(xlim = c(-2, 102), ylim = c(-2, 102), expand = TRUE) +
    theme(
      text = element_text(size = 12, family = "Montserrat"),
      legend.position = c(0, 1), 
      legend.justification = c("left", "top"), 
      legend.title = element_blank(),
      axis.title = element_text(face = "plain"),
      legend.text = element_text(size = 12),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")
    ) +
    labs(
      x = paste0(
        "% of ", var1, " with a peptide \nsignificantly bound ",
        "(n = ", N1, ")"
      ),
      y = paste0(
        "% of ", var2, " with a peptide \nsignificantly bound ",
        "(n = ", N2, ")"
      )
    )
  
  p_static$data <- p_static$data %>%
    mutate(
      pt_size = if_else(is_bac_flagella == 1, 2, 1),
      legend_group = if_else(is_bac_flagella == 1, legend_lab, NA_character_)
    )

  p_static <- p_static +
    geom_point(
      data = p_static$data %>% filter(is_bac_flagella != 1),
      aes(x = percent1, y = percent2),
      colour = "grey70",
      size = 1.2,
      alpha = 0.5,
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +
    geom_point(
      data = p_static$data %>% filter(is_bac_flagella == 1),
      aes(x = percent1, y = percent2, colour = legend_group),
      size = 2,
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    scale_colour_manual(
      values = c(setNames(peptide_palette[["Flagellins"]], legend_lab))
    ) +
    scale_size_identity() +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 0.5)))
  
  ggsave(
    filename = file.path(paste0(out_name, "_static.pdf")),
    plot = p_static,
    width = 9, height = 9, units = "cm",
    device = cairo_pdf, bg = "white"
  )
  
  ## Boxplot - peptide categories (each dot = peptide)
  # Remove outliers for plotting (nicer plots)
  pep_tbl_peplib_long <- pep_tbl_peplib_long %>%
    group_by(family_lab) %>%
    filter(!ratio %in% boxplot.stats(ratio)$out)
  
  # Stacked y positions so brackets don't overlap
  max_y <- max(pep_tbl_peplib_long$ratio, na.rm = TRUE)

  # Keep factor order exactly as used in the plot
  x_levels <- levels(factor(pep_tbl_peplib_long$family_lab))
  flag_pos <- match("Flagellins", x_levels)

  # Build p-value table for only comparisons involving Flagellins
  pdat <- pval_mw_df %>%
    filter(group1 == "Flagellins" | group2 == "Flagellins") %>%
    mutate(
      other_group = if_else(group1 == "Flagellins", group2, group1),
      other_pos   = match(other_group, x_levels),
      side        = if_else(other_pos < flag_pos, "left", "right"),
      dist        = abs(other_pos - flag_pos)
    ) %>%
    # Make sure xmin is always the left category and xmax the right category
    mutate(
      xmin = if_else(other_pos < flag_pos, other_group, "Flagellins"),
      xmax = if_else(other_pos < flag_pos, "Flagellins", other_group)
    ) %>%
    # Alternate left/right from inside to outside
    arrange(dist, factor(side, levels = c("left", "right"))) %>%
    mutate(
      y.position = max_y * (1 + 0.11 * row_number())
    )

  # Boxplot (each dot = peptide)
  p_pep_box <- ggplot(
    pep_tbl_peplib_long, 
    aes(x = family_lab, y = ratio, fill = family_lab)
  ) +
    geom_boxplot(outlier.shape = NA) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    scale_fill_manual(values = peptide_palette) +
    stat_pvalue_manual(
      pdat,
      label = "p.adj.signif",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.02
    ) +
    coord_cartesian(ylim = c(0, 6)) +
    labs(
      x = NULL, 
      y = paste0(
        "Ratios of Ig response in \n", label_dir,
        " (", N1, ", ", N2, ")"
      )
    ) +
    theme_phip(12) +
    theme(
      text = element_text(size = 12, family = "Montserrat"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.title = element_text(face = "plain"),
      plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")
    )
  
  ggsave(
    filename = file.path(paste0(out_name, "_pep_box.pdf")),
    plot = p_pep_box,
    width = 18, height = 9, units = "cm",
    device = cairo_pdf, bg = "white"
  )
  
  if (any(vapply(comp_main, identical, logical(1), cmp))) {
    # ----------------------------------------------------------------------------
    # DELTA framework
    # ----------------------------------------------------------------------------
    dir.create(file.path(out_dir, "DELTA_framework"), recursive = T, showWarnings = F)
    
    data_frameworks$subject_id <- data_frameworks$sample_id
    peptide_library[] <- lapply(peptide_library, as.character)
    log_file_current <- if (is.null(LOG_FILE)) {
      file.path(out_dir, "DELTA_framework", "log.txt")
    } else {
      LOG_FILE
    }
    
    delta_file <- file.path(out_dir, "DELTA_framework", "delta_table.csv")
    if(!file.exists(delta_file) | FORCE) {
      res <- phiper::compute_delta(
        x = data_frameworks,
        exist_col = "exist",
        rank_cols = c(
          "phylum", "class", "order", "family", "genus", "species",
          "is_auto", "is_infect", "is_EBV", "is_toxin", "is_PNP", "is_EM",
          "is_MPA", "is_patho", "is_probio", "is_IgA", "is_bac_flagella",
          "is_allergens"
        ),
        group_cols          = "group_char",
        peptide_library     = peptide_library,
        B_permutations      = 150000L,
        smooth_eps_num      = 0.5,
        smooth_eps_den_mult = 2.0,
        min_max_prev        = 0.0,
        weight_mode         = "n_eff_sqrt",
        stat_mode           = "asin",
        prev_strat          = "none",
        winsor_z            = Inf,
        rank_feature_keep   = list(
          phylum  = NULL, class = NULL, order = NULL, family = NULL, genus = NULL,
          species = NULL,
          is_auto = "1", is_infect = "1", is_EBV = "1", is_toxin = "1",
          is_PNP = "1", is_EM = "1", is_MPA  = "1", is_patho = "1", 
          is_probio = "1", is_IgA = "1", is_bac_flagella = "1", 
          is_allergens = "1"
        ),
        log                 = LOG,
        log_file            = log_file_current,
        fold_change         = "sum",
        cross_prev          = "mean"
      )
      res <- as.data.frame(res)
      write.csv(res, file = delta_file)
    } else {
      res <- read.csv(delta_file)
    }
  }
}

# ------------------------------------------------------------------------------
# Restore original future plan
# ------------------------------------------------------------------------------
future::plan(original_plan)

# ------------------------------------------------------------------------------
# Add some additional cross-comparison plots
# ------------------------------------------------------------------------------

# Heatmap: Pearson (T0/T1) - R and NR combined, matched patients only

cor_mat <- list()
for(cmp in c(
  "R_T0_vs_R_T1", 
  "R_T0_vs_NR_T1", 
  "NR_T0_vs_R_T1", 
  "NR_T0_vs_NR_T1"
)) {
  cor_mat[[cmp]] <- read.csv(
    file.path("results", cmp, "beta_diversity", "cor_matrix_pearson.csv")
  )
}

# Four blocks
df_rr <- cor_mat[["R_T0_vs_R_T1"]]
df_rn <- cor_mat[["R_T0_vs_NR_T1"]]
df_nr <- cor_mat[["NR_T0_vs_R_T1"]]
df_nn <- cor_mat[["NR_T0_vs_NR_T1"]]

# Keep only patients occurring in both timepoints within each response group
common_r <- sort(intersect(unique(df_rr$pid_t1), unique(df_rr$pid_t2)))
common_nr <- sort(intersect(unique(df_nn$pid_t1), unique(df_nn$pid_t2)))

# Create sequential labels
r_ids  <- setNames(as.character(seq_along(common_r)),  common_r)
nr_ids <- setNames(as.character(seq_along(common_nr)), common_nr)

# Filter all four blocks consistently
df_rr_f <- df_rr %>%
  filter(pid_t1 %in% common_r, pid_t2 %in% common_r) %>%
  mutate(
    pid_t1_num = factor(r_ids[pid_t1], levels = as.character(seq_along(common_r))),
    pid_t2_num = factor(r_ids[pid_t2], levels = as.character(seq_along(common_r)))
  )

df_rn_f <- df_rn %>%
  filter(pid_t1 %in% common_r, pid_t2 %in% common_nr) %>%
  mutate(
    pid_t1_num = factor(r_ids[pid_t1], levels = as.character(seq_along(common_r))),
    pid_t2_num = factor(nr_ids[pid_t2], levels = as.character(seq_along(common_nr)))
  )

df_nr_f <- df_nr %>%
  filter(pid_t1 %in% common_nr, pid_t2 %in% common_r) %>%
  mutate(
    pid_t1_num = factor(nr_ids[pid_t1], levels = as.character(seq_along(common_nr))),
    pid_t2_num = factor(r_ids[pid_t2], levels = as.character(seq_along(common_r)))
  )

df_nn_f <- df_nn %>%
  filter(pid_t1 %in% common_nr, pid_t2 %in% common_nr) %>%
  mutate(
    pid_t1_num = factor(nr_ids[pid_t1], levels = as.character(seq_along(common_nr))),
    pid_t2_num = factor(nr_ids[pid_t2], levels = as.character(seq_along(common_nr)))
  )

# Axis ordering
r_x_ids  <- levels(df_rr_f$pid_t1_num)
nr_x_ids <- levels(df_nn_f$pid_t1_num)
r_y_ids  <- levels(df_rr_f$pid_t2_num)
nr_y_ids <- levels(df_nn_f$pid_t2_num)

x_levels <- c(paste0("R_", r_x_ids), paste0("NR_", nr_x_ids))
y_levels <- c(paste0("R_", r_y_ids), paste0("NR_", nr_y_ids))

x_labels <- c(
  setNames(ifelse(seq_along(r_x_ids) %in% c(1, 10, 20, 30, 40), seq_along(r_x_ids), ""),
           paste0("R_", r_x_ids)),
  setNames(ifelse(seq_along(nr_x_ids) %in% c(1, 10, 20, 30, 40), seq_along(nr_x_ids), ""),
           paste0("NR_", nr_x_ids))
)

y_labels <- c(
  setNames(ifelse(seq_along(r_y_ids) %in% c(1, 10, 20, 30, 40), seq_along(r_y_ids), ""),
           paste0("R_", r_y_ids)),
  setNames(ifelse(seq_along(r_y_ids) %in% c(1, 10, 20, 30, 40), seq_along(nr_y_ids), ""),
           paste0("NR_", nr_y_ids))
)

# Map each filtered block onto the combined axes
df_comb <- bind_rows(
  df_rr_f %>% mutate(
    x_comb = paste0("R_", pid_t1_num),
    y_comb = paste0("R_", pid_t2_num)
  ),
  df_rn_f %>% mutate(
    x_comb = paste0("R_", pid_t1_num),
    y_comb = paste0("NR_", pid_t2_num)
  ),
  df_nr_f %>% mutate(
    x_comb = paste0("NR_", pid_t1_num),
    y_comb = paste0("R_", pid_t2_num)
  ),
  df_nn_f %>% mutate(
    x_comb = paste0("NR_", pid_t1_num),
    y_comb = paste0("NR_", pid_t2_num)
  )
) %>%
  mutate(
    x_comb = factor(x_comb, levels = x_levels),
    y_comb = factor(y_comb, levels = y_levels)
  )

p_corr_heat <- ggplot(df_comb, aes(x = x_comb, y = y_comb, fill = pearson)) +
  geom_tile() +
  scale_fill_viridis_c(
    option = "plasma",
    breaks = seq(0, 1, 0.1),
    guide = guide_colorbar(
      barheight = unit(0.75, "npc"),
      barwidth  = unit(0.3, "cm"),
      frame.colour = "black",
      frame.linewidth = 0.2
    )
  ) +
  scale_x_discrete(
    breaks = c(
      paste0("R_", r_x_ids[c(1, seq(5, length(r_x_ids), 5))]),
      paste0("NR_", nr_x_ids[c(1, seq(5, length(nr_x_ids), 5))])
    ),
    labels = c(
      x_labels[paste0("R_", r_x_ids[c(1, seq(5, length(r_x_ids), 5))])],
      x_labels[paste0("NR_", nr_x_ids[c(1, seq(5, length(nr_x_ids), 5))])]
    )
  ) +
  scale_y_discrete(
    breaks = c(
      paste0("R_", r_y_ids[c(1, seq(5, length(r_y_ids), 5))]),
      paste0("NR_", nr_y_ids[c(1, seq(5, length(nr_y_ids), 5))])
    ),
    labels = c(
      y_labels[paste0("R_", r_y_ids[c(1, seq(5, length(r_y_ids), 5))])],
      y_labels[paste0("NR_", nr_y_ids[c(1, seq(5, length(nr_y_ids), 5))])]
    )
  ) +
  coord_fixed() +
  geom_vline(
    xintercept = length(common_r) + 0.5,
    linewidth = 0.35, 
    colour = "white", 
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = length(common_r) + 0.5,
    linewidth = 0.35, 
    colour = "white", 
    linetype = "dashed"
  ) +
  theme_bw(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(size = 12),
    plot.subtitle = element_text(size = 10),
    axis.title = element_blank(),
    axis.text = element_text(size = 10, colour = "black"),
    axis.ticks = element_line(colour = "black", linewidth = 0.3),
    axis.ticks.length = unit(1.5, "mm"),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    legend.margin = margin(0.5, 0.5, 0.5, 0.5),
    legend.box.margin = margin(0, 0, 0, -5),
    panel.grid = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")
  ) +
  labs(
    title = "T0 vs. T1 Ig epitope repertoires",
    subtitle = "Pearson correlation of log(fold change)"
  )

ggsave(
  filename = file.path("results", "corr_heatmap_plot.pdf"),
  plot = p_corr_heat,
  width = 8.8, height = 8.8, units = "cm",
  device = cairo_pdf, bg = "white"
)

# Correlation coefficients (Barplots)

cor_stats <- list()
for(cmp in comp_tps[1:12]) {
  cmp_lab <- paste(cmp[[1]], "vs", cmp[[2]], sep = "_")
  cor_mat <- read.csv(
    file.path("results", cmp_lab, "beta_diversity", "cor_matrix_pearson.csv")
  )
  common_patients <- intersect(unique(cor_mat$pid_t1), unique(cor_mat$pid_t2))
  matched <- cor_mat %>%
    filter(pid_t1 == pid_t2) %>%
    pull(pearson)
  unmatched <- cor_mat %>%
    filter(pid_t1 %in% common_patients, pid_t2 %in% common_patients) %>%
    filter(pid_t1 != pid_t2) %>%
    pull(pearson)
  cor_stats[[cmp_lab]] <- list(matched = matched, unmatched = unmatched)
}

# Extract just the correlations of the matched samples
matched_df <- bind_rows(
  lapply(names(cor_stats), function(nm) {
    vals <- cor_stats[[nm]]$matched
    data.frame(contrast = nm, value = vals)
  })
) %>%
  mutate(
    group = str_extract(contrast, "(R|NR)"),
    time_comp = str_match(contrast, "T(\\d+)_vs_(?:NR|R)_T(\\d+)") %>% 
      (\(m) paste0("T", m[,2], "/T", m[,3]))()  
  )

# Compute significance between R/NR at each timepoint comparison
stats_tests <- compare_means(
  value ~ group,
  data  = matched_df,
  group.by = "time_comp",
  method = "wilcox.test",
  p.adjust.method = "bonferroni"
) %>%
  rstatix::add_significance("p.adj")

# Put everything together
matched_df <- matched_df %>%
  group_by(time_comp, group) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(stats_tests, by = "time_comp") %>%
  mutate(group = factor(group, levels = c("R", "NR")))

p_corr_bars <- ggplot(matched_df, aes(x = time_comp, y = mean, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  geom_text(
    aes(x = time_comp, y = 1.1, label = p.adj.signif),
    inherit.aes = FALSE,
    vjust = 0,
    size = 3.51
  ) +
  scale_fill_manual(values = group_palette) +
  labs(
    x = NULL,
    y = "Mean Pearson correlation \n(matched samples)"
  ) +
  theme_phip(12) +
  theme(
    text = element_text(size = 12, family = "Arial"),
    axis.title = element_text(face = "plain"),
    axis.text = element_text(size = 10),
    panel.grid = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path("results", "corr_bars_plot.pdf"),
  plot = p_corr_bars,
  width = 9, height = 9, units = "cm",
  device = cairo_pdf, bg = "white"
)

# Ratios of response to flagellins by class at different timepoints

pep_tbl_all <- list()
pval_mw_df_all <- list()
prev10 <- c()
for(cmp_pair in c(comp_main[1:4], cmp_within$cmp)) {
  cmp <- paste0(cmp_pair[1], "_vs_", cmp_pair[2])
  pep_tbl_all[[cmp]] <- readRDS(file.path("results", cmp, "POP_framework", "single_peptide.rds"))
  pval_mw_df_all[[cmp]] <- read.csv(file.path("results", cmp, "POP_framework", "pval_mw_df.csv"))
  
  pep_tbl_prev10 <- pep_tbl_all[[cmp]] %>% filter((prop1 + prop2 >= 0.1))
  prev10 <- union(prev10, pep_tbl_prev10$feature)
}

flagellins <- peptide_library %>% filter(is_bac_flagella == 1) %>% pull(peptide_id)
class_lvls <- c("response (R/NR)", unique(labs$group_type))

class_tp_df <- bind_rows(pep_tbl_all, .id = "contrast") %>%
  filter(feature %in% flagellins) %>%
  filter(feature %in% prev10) %>%
  mutate(
    tp = str_extract(group1, "T\\d+$"),
    grp = str_remove(group1, "_T\\d+$"),
    class = case_when(
      grp %in% c("R", "NR") ~ "response (R/NR)",
      str_detect(grp, "_") ~ str_extract(grp, "^[^_]+"),
      TRUE ~ grp
    ),
    class = factor(class, levels = class_lvls)
  )

# Remove outliers for plotting (nicer plots)
class_tp_df <- class_tp_df %>%
  group_by(class, tp) %>%
  filter(!ratio %in% boxplot.stats(ratio)$out)

# Calculate the top of the whisker
ypos_df <- class_tp_df %>%
  group_by(class, tp) %>%
  summarise(
    y.position = boxplot.stats(ratio)$stats[5] + 0.5, 
    .groups = "drop"
  )

# Build p-value table for flagellins vs. complete library only
pdat <- bind_rows(pval_mw_df_all, .id = "contrast") %>%
  filter(
    (group1 == "Flagellins" & group2 == "Complete library") | 
    (group1 == "Complete library" & group2 == "Flagellins")
  ) %>%
  mutate(p.adj.signif = if_else(p.adj.signif == "ns", "", p.adj.signif)) %>%
  select(contrast, p.adj.signif) %>%
  left_join(
    class_tp_df %>% distinct(contrast, group1, group2, class, tp),
    by = "contrast"
  ) %>%
  left_join(ypos_df, by = c("class", "tp"))

p_class_tp_box <- class_tp_df %>%
  ggplot(aes(x = class, y = ratio, fill = tp)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  labs(
    x = NULL,
    y = "Ratios of Ig response to \nflagellins by class (yes/no)", 
    fill = "Timepoint"
  ) +
  scale_fill_manual(values = c(
    "T0" = "grey70",
    "T1" = "#6BBF58",
    "T2" = "#59A14F",
    "T3" = "#3D702D"
  )) +
  geom_text(
    data = pdat,
    aes(x = class, y = y.position, label = p.adj.signif, group = tp),
    position = position_dodge(width = 0.8),
    inherit.aes = FALSE,
    size = 3.51,
    angle = 90,
    hjust = 0,
    vjust = 0.7
  ) +
  theme_phip(12) +
  theme(
    text = element_text(size = 12, family = "Montserrat"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(face = "plain"),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.title = element_text(face = "plain")
  )

ggsave(
  filename = file.path("results", "class_tp_box_plot.pdf"),
  plot = p_class_tp_box,
  width = 18, height = 9, units = "cm",
  device = cairo_pdf, bg = "white"
)

# Time interval between treatment start and sample collection
# Sample metadata
meta <- read.csv(file.path("data/processed/", "metadata.csv")) %>%
  mutate(group_char = paste(response, timepoint, sep = "_"))

p_time_diff <- meta %>%
  mutate(
    baseline_date = as.Date(start_io, "%d/%m/%Y"),
    timepoint_date = as.Date(serum_collection_date, "%d/%m/%Y"),
    diff_tp_base_week = difftime(timepoint_date, baseline_date, units = "weeks")
  ) %>%
  ggplot(aes(x = timepoint, y = diff_tp_base_week)) +
  geom_boxplot(outliers = TRUE) +
  labs(
    x = NULL,
    y = "Time interval between treatment \nstart and sample collection (weeks)"
  ) +
  theme_phip(12) + 
  theme(
    text = element_text(size = 12, family = "Montserrat"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "plain"),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, unit = "pt")
  )

ggsave(
  filename = file.path("results", "time_diff_plot.pdf"),
  plot = p_time_diff,
  width = 9, height = 9, units = "cm",
  device = cairo_pdf, bg = "white"
)
