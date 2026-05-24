# ============================================================================
# PHASE 2 — Classify BIN vs COUNT
# ============================================================================

# ----------------------------------------------------------------------------
# 0) Load TAXA
# ----------------------------------------------------------------------------
taxa_clean <- taxa_raw %>%
  dplyr::mutate(
    SAMPLE_ID       = arrow::cast(SAMPLE_ID, arrow::int32()),
    TOTAL_ABUNDANCE = arrow::cast(TOTAL_ABUNDANCE, arrow::int32())
  ) %>%
  dplyr::select(SAMPLE_ID, TOTAL_ABUNDANCE) %>%
  dplyr::collect()

# ----------------------------------------------------------------------------
# 1) Classification Logic
# ----------------------------------------------------------------------------
bin_placeholders <- c(1L, 3L, 33L, 333L, 3333L, 33333L)

taxa_nonzero <- taxa_clean %>%
  dplyr::filter(!is.na(TOTAL_ABUNDANCE), TOTAL_ABUNDANCE != 0L)

taxa_flags_initial <- taxa_nonzero %>%
  dplyr::group_by(SAMPLE_ID) %>%
  dplyr::summarise(
    all_placeholders = all(TOTAL_ABUNDANCE %in% bin_placeholders),
    any_placeholder  = any(TOTAL_ABUNDANCE %in% bin_placeholders),
    any_numeric      = any(!(TOTAL_ABUNDANCE %in% bin_placeholders)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    data_type_taxa = dplyr::case_when(
      all_placeholders ~ "bin",
      any_placeholder & any_numeric ~ "bin",
      TRUE ~ "count"
    )
  )

sample_flags_taxa <- metrics_clean %>%
  dplyr::select(SAMPLE_ID, ANALYSIS_METHOD) %>%
  dplyr::left_join(taxa_flags_initial, by = "SAMPLE_ID")

resolve_data_type <- function(analysis_method, from_taxa) {
  if (is.na(analysis_method)) return(NA_character_)
  if (analysis_method == "ANLE") return("bin")
  if (analysis_method == "ANAA") return("count")
  if (analysis_method == "ANLA") {
    if (!is.na(from_taxa)) return(from_taxa)
    return("bin")
  }
  "bin"
}

sample_flags_final <- sample_flags_taxa %>%
  dplyr::mutate(
    data_type = purrr::pmap_chr(
      list(ANALYSIS_METHOD, data_type_taxa),
      ~ resolve_data_type(..1, ..2)
    )
  ) %>%
  dplyr::select(SAMPLE_ID, data_type)

roster_flagged <- metrics_clean %>%
  dplyr::left_join(sample_flags_final, by = "SAMPLE_ID") %>%
  dplyr::mutate(data_type = tidyr::replace_na(data_type, "bin"))

# ----------------------------------------------------------------------------
# 2) Season Feature
# ----------------------------------------------------------------------------
KEEP_OFFSEASON <- FALSE

roster_flagged <- roster_flagged %>%
  mutate(
    month  = month(SAMPLE_DATE),
    season = case_when(
      month %in% 3:5  ~ "spring",
      month %in% 9:11 ~ "autumn",
      TRUE ~ NA_character_
    ),
    season_f = factor(season, levels = c("spring","autumn"))
  ) %>%
  { if (!KEEP_OFFSEASON) dplyr::filter(., !is.na(season_f)) else . } %>%
  dplyr::select(
    SITE_ID, SAMPLE_ID, SAMPLE_DATE,
    SAMPLE_METHOD, ANALYSIS_METHOD,
    data_type, season_f, month
  )

# ============================================================================
# PHASE 3 — Build Modelling Tables
# ============================================================================

families <- c(
  "Aphelocheiridae",
  "Brachycentridae",
  "Odontoceridae",
  "Cordulegastridae"
)

make_family_tables <- function(family_name, roster_tbl, whpt_tbl) {

  sites_with_fam <- whpt_tbl %>%
    dplyr::filter(EQ_TAXON_UNIT == family_name) %>%
    dplyr::distinct(SITE_ID) %>%
    dplyr::pull(SITE_ID)

  roster_family <- roster_tbl %>%
    dplyr::filter(SITE_ID %in% sites_with_fam)

  family_counts <- whpt_tbl %>%
    dplyr::filter(EQ_TAXON_UNIT == family_name) %>%
    dplyr::select(SITE_ID, SAMPLE_ID, TOTAL_NUMBER)

  fam_data <- roster_family %>%
    dplyr::left_join(family_counts, by = c("SITE_ID","SAMPLE_ID")) %>%
    dplyr::mutate(
      TOTAL_NUMBER  = tidyr::replace_na(TOTAL_NUMBER, 0L),
      taxon_present = TOTAL_NUMBER > 0,
      decimal_date  = lubridate::decimal_date(SAMPLE_DATE),
      SITE_ID.F     = factor(SITE_ID)
    )

  pa_tbl <- fam_data

  ord_tbl <- fam_data %>%
    dplyr::filter(data_type == "bin") %>%
    dplyr::mutate(
      ord_class = dplyr::case_when(
        TOTAL_NUMBER == 0L ~ "AB0",
        TOTAL_NUMBER == 1L ~ "AB1",
        TOTAL_NUMBER == 3L ~ "AB2",
        TOTAL_NUMBER == 33L ~ "AB3",
        TOTAL_NUMBER >= 333L ~ "AB4",
        TRUE ~ NA_character_
      ),
      ord_class = forcats::fct_collapse(ord_class, AB3p = c("AB3","AB4"))
    ) %>%
    dplyr::filter(!is.na(ord_class)) %>%
    dplyr::mutate(
      ord_class = factor(ord_class,
                         levels = c("AB0","AB1","AB2","AB3p"),
                         ordered = TRUE)
    )

  cpois_tbl <- fam_data %>%
    dplyr::filter(data_type == "count") %>%
    dplyr::mutate(
      lower = dplyr::case_when(
        TOTAL_NUMBER <= 9L ~ as.numeric(TOTAL_NUMBER) - 0.5,
        TRUE ~ floor(TOTAL_NUMBER / 10) * 10 - 0.5
      ),
      upper = dplyr::case_when(
        TOTAL_NUMBER <= 9L ~ as.numeric(TOTAL_NUMBER) + 0.5,
        TRUE ~ floor(TOTAL_NUMBER / 10) * 10 + 9.5
      )
    )

  cnorm_tbl <- cpois_tbl

  list(
    pa      = pa_tbl,
    ordinal = ord_tbl,
    cpois   = cpois_tbl,
    cnorm   = cnorm_tbl
  )
}

models_data <- purrr::map(
  families,
  make_family_tables,
  roster_tbl = roster_flagged,
  whpt_tbl   = whpt_clean
) %>%
  rlang::set_names(families)
