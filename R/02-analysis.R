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

# 1) within-timepoint: yes vs no
cmp_within <- labs %>%
  select(group_type, group_value, timepoint, lab) %>%
  pivot_wider(names_from = group_value, values_from = lab) %>%
  transmute(cmp = map2(yes, no, c))

# Leave for now, maybe add back later for supplementary figures?
# # 2) across-timepoint within level: yes_Ti vs yes_Tj AND no_Ti vs no_Tj
# cmp_across <- labs %>%
#   group_by(group_type, group_value) %>%
#   summarise(cmp = list(combn(lab, 2, simplify = FALSE)), .groups = "drop") %>%
#   unnest(cmp)
# comp_add <- c(cmp_within$cmp, cmp_across$cmp)
# comparisons <- c(comp_main, comp_add)

comparisons <- c(comp_main, cmp_within)

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
#   6) DELTA framework (permutation-based differential prevalence)
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
    mutate(is_library = TRUE)

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
  )
  pval_mw_file <- file.path(out_dir, "POP_framework", "pval_mw_df.csv")
  write.csv(pval_mw_df, file = pval_mw_file)
  
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
    xlab = pep_tbl$group1[1],
    ylab = pep_tbl$group2[1],
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
  p_pep_box <- ggplot(pep_tbl_peplib_long, aes(x = family_lab, y = ratio, fill = family_lab)) +
    geom_boxplot(outlier.shape = NA) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    scale_fill_manual(values = peptide_palette) +
    stat_pvalue_manual(
      pdat,
      label = "p.signif",
      xmin = "xmin",
      xmax = "xmax",
      y.position = "y.position",
      tip.length = 0.02
    ) +
    coord_cartesian(ylim = c(0, 6)) +
    labs(
      x = NULL, 
      y = paste0(
        "Antibody response ratios \nin ", label_dir,
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
          is_auto = "TRUE", is_infect = "TRUE", is_EBV = "TRUE", is_toxin = "TRUE",
          is_PNP = "TRUE", is_EM = "TRUE", is_MPA  = "TRUE", is_patho = "TRUE", 
          is_probio = "TRUE", is_IgA = "TRUE", is_bac_flagella = "yes", 
          is_allergens = "TRUE"
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
