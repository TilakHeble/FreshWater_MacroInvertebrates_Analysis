# ================================================================
# PHASE 0 — Loading Libraries and Datasets
# ================================================================

# ---------------------------------------------------------------
# 0) Package setup — install if missing, then load quietly
# (All packages are CRAN; install once per machine.)
# ---------------------------------------------------------------
load_if_needed <- function(pkgs) {
  missing <- setdiff(pkgs, rownames(installed.packages()))
  if (length(missing)) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, dependencies = TRUE)
  }
  invisible(lapply(pkgs, function(p)
    suppressPackageStartupMessages(library(p, character.only = TRUE))))
}

req_pkgs <- c(
  # --- core tidy/data handling ---
  "arrow",           # Fast Parquet/Feather IO + lazy datasets
  "dplyr",           # Data wrangling (filter/mutate/summarise/join)
  "magrittr",        # Pipes (%>%, %<>%) & helpers
  "lubridate",       # Dates/times (ymd(), year(), decimal_date())
  "tidyr",           # Reshaping (pivot_longer/wider, separate, unite)
  "forcats",         # Factor tools (fct_* helpers)
  "purrr",           # Functional tools (map, pmap, safely)

  # --- modelling & visualisation (used later in the pipeline) ---
  "ggplot2",         # Plotting
  "scales",          # Axis/label formatting
  "mgcv",            # GAM/BAM fitting
  "gratia",          # GAM diagnostics & partial-effects
  "readr",           # Fast read/write for CSV
  "MASS",            # Misc. stats (e.g., mvrnorm)
  "patchwork",       # Compose multiple ggplots

  # --- light mapping (used in later EDA) ---
  "sf",
  "rnaturalearth",
  "rnaturalearthdata"
)

load_if_needed(req_pkgs)

# ---------------------------------------------------------------
# Reproducibility settings
# ---------------------------------------------------------------
options(stringsAsFactors = FALSE)
set.seed(20250816)

# ================================================================
# 1) Data location — point to the Parquet files
# *Adjust DATA_DIR if your files live elsewhere.*
# ================================================================
DATA_DIR <- getwd()  # e.g. change to "D:/EA_open_data" if needed

paths <- list(
  sites   = file.path(DATA_DIR, "INV_OPEN_DATA_SITE.parquet"),
  metrics = file.path(DATA_DIR, "INV_OPEN_DATA_METRICS.parquet"),
  whpt    = file.path(DATA_DIR, "R_INV_WHPT_METRICS_B.parquet"),
  taxa    = file.path(DATA_DIR, "INV_OPEN_DATA_TAXA.parquet")
)

# ---------------------------------------------------------------
# Helper: fail fast if files are missing
# ---------------------------------------------------------------
assert_parquet <- function(path) {
  if (!file.exists(path)) {
    stop(
      sprintf(
        "Parquet not found:\n  %s\n\nFix: set DATA_DIR correctly or download the file.",
        normalizePath(path, winslash = "/", mustWork = FALSE)
      ),
      call. = FALSE
    )
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

# ================================================================
# 2) Open raw datasets (lazy loading with Arrow)
# Nothing is loaded into memory until `collect()` is used
# ================================================================
sites_raw   <- arrow::open_dataset(assert_parquet(paths$sites))
metrics_raw <- arrow::open_dataset(assert_parquet(paths$metrics))
whpt_raw    <- arrow::open_dataset(assert_parquet(paths$whpt))
taxa_raw    <- arrow::open_dataset(assert_parquet(paths$taxa))

# ---------------------------------------------------------------
# Optional sanity checks (uncomment if needed)
# ---------------------------------------------------------------
# dplyr::glimpse(dplyr::collect(head(sites_raw,   3)))
# dplyr::glimpse(dplyr::collect(head(metrics_raw, 3)))
# dplyr::glimpse(dplyr::collect(head(whpt_raw,    3)))
# dplyr::glimpse(dplyr::collect(head(taxa_raw,    3)))

# ---------------------------------------------------------------
# Note:
# Keep operations lazy as long as possible.
# Use `collect()` only when necessary for performance.
# ---------------------------------------------------------------
