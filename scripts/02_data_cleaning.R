# ============================================================================
# PHASE 1 — Configure & Clean
# ============================================================================

# Purpose:
# Build the three core cleaned tables used throughout the analysis:
#   - metrics_clean: one row per (SITE_ID, SAMPLE_DATE)
#   - whpt_clean   : harmonised family counts
#   - sites_clean  : site metadata

# Notes:
#   • Keep only ANAA / ANLA / ANLE methods
#   • Default to S3PO sampling
#   • Drop pre-1990 data (standardisation point)
#   • Use Arrow for efficient data handling

# ----------------------------------------------------------------------------
# 0) Project Options
# ----------------------------------------------------------------------------
USE_S3PO   <- TRUE
DATE_FLOOR <- as.Date("1990-01-01")

FAMILIES <- c(
  "Aphelocheiridae",
  "Brachycentridae",
  "Odontoceridae",
  "Cordulegastridae"
)

# ----------------------------------------------------------------------------
# 1) Clean METRICS table
# ----------------------------------------------------------------------------
metrics_clean <- metrics_raw %>%
  dplyr::filter(ANALYSIS_METHOD %in% c("ANAA","ANLA","ANLE")) %>%
  dplyr::mutate(
    SITE_ID     = arrow::cast(SITE_ID, arrow::int32()),
    SAMPLE_ID   = arrow::cast(SAMPLE_ID, arrow::int32()),
    SAMPLE_DATE = arrow::cast(SAMPLE_DATE, arrow::date32())
  ) %>%
  dplyr::select(
    SITE_ID, SAMPLE_ID, SAMPLE_DATE,
    SAMPLE_TYPE, SAMPLE_METHOD,
    ANALYSIS_TYPE, ANALYSIS_METHOD
  ) %>%
  dplyr::collect() %>%
  dplyr::mutate(SAMPLE_DATE = as.Date(SAMPLE_DATE)) %>%
  dplyr::arrange(SITE_ID, SAMPLE_DATE, SAMPLE_ID) %>%
  { if (USE_S3PO) dplyr::filter(., SAMPLE_METHOD == "S3PO") else . } %>%
  dplyr::filter(SAMPLE_DATE >= DATE_FLOOR) %>%
  dplyr::group_by(SITE_ID, SAMPLE_DATE) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup()

# ----------------------------------------------------------------------------
# 2) Clean WHPT table
# ----------------------------------------------------------------------------
whpt_clean <- whpt_raw %>%
  dplyr::mutate(
    SITE_ID      = arrow::cast(SITE_ID, arrow::int32()),
    SAMPLE_ID    = arrow::cast(SAMPLE_ID, arrow::int32()),
    TOTAL_NUMBER = arrow::cast(TOTAL_NUMBER, arrow::int32())
  ) %>%
  dplyr::select(SITE_ID, SAMPLE_ID, EQ_TAXON_UNIT, TOTAL_NUMBER) %>%
  dplyr::collect()

# ----------------------------------------------------------------------------
# 3) Clean SITES table
# ----------------------------------------------------------------------------
sites_clean <- sites_raw %>%
  dplyr::mutate(SITE_ID = arrow::cast(SITE_ID, arrow::int32())) %>%
  dplyr::select(
    SITE_ID, REPORTING_AREA, CATCHMENT, WATER_BODY,
    FULL_EASTING, FULL_NORTHING
  ) %>%
  dplyr::collect()
