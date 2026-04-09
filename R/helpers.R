# Function to calculate Kulczynski similarity ----------------------------------
kulczynski_mat <- function(A, B) {
  a <- t(A) %*% B   # (nA x nB) shared presences
  nA <- colSums(A)  # length nA
  nB <- colSums(B)  # length nB
  termA <- sweep(a, 1, nA, "/")  # divide each row by nA
  termB <- sweep(a, 2, nB, "/")  # divide each col by nB
  kulc <- 0.5 * (termA + termB)
  
  # If a sample has nA==0 or nB==0, define similarity as 0
  kulc[!is.finite(kulc)] <- 0
  kulc
}

# Function to format p-values nicely -------------------------------------------
format_pval <- function(p, alpha = 0.05) {
  
  # helper to drop trailing zeros (e.g. "1.00" → "1", "0.500" → "0.5")
  drop_zeros <- function(x) sub("\\.?0+$", "", x)
  
  if (is.na(p)) return("NA")
  
  # -------------------------------------------------------------------------
  # 1. Non-significant (p > alpha)
  # -------------------------------------------------------------------------
  if (p > alpha) {
    raw <- formatC(p, digits = 2, format = "f")
    raw <- drop_zeros(raw)
    return(paste0("ns [", raw, "]"))
  }
  
  # -------------------------------------------------------------------------
  # 2. Normal fixed-decimal formatting (0.001 ≤ p ≤ alpha)
  # -------------------------------------------------------------------------
  if (p >= 0.001) {
    raw <- formatC(p, digits = 3, format = "f")
    raw <- drop_zeros(raw)
    return(raw)
  }
  
  # -------------------------------------------------------------------------
  # 3. Scientific notation (< 0.001)
  # -------------------------------------------------------------------------
  raw <- formatC(p, digits = 2, format = "e")
  
  # remove unnecessary zeros: "1.00e-05" → "1e-05"
  raw <- sub("([0-9]+)\\.0+e", "\\1e", raw)
  
  # remove trailing zeros inside "1.10e-04" → "1.1e-04"
  raw <- sub("([0-9]+\\.[0-9]*[1-9])0+e", "\\1e", raw)
  
  return(raw)
}


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
  
  avail_cols <- names(data) # MH: THIS IS ALWAYS NULL -> phiper issue
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

# Extract peptide table from POP result
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
