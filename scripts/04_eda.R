# -------------------------------------------------------------------
# Polished EDA plots for the EA macroinvertebrate project
# Expects: metrics_clean, whpt_clean, roster_flagged already created
# Uses: families vector (defined here if missing)
# Purpose: Produce clear, publication-ready figures directly tied to the
#          project questions (methods over time, family frequencies,
#          categorical vs numeric mix, abundance distributions, seasonality,
#          and presence rates).
# -------------------------------------------------------------------

# If families wasn't defined earlier
if (!exists("families")) {
  families <- c("Aphelocheiridae","Brachycentridae","Odontoceridae","Cordulegastridae")
}

# ---------- Styling ----------
# Minimal, consistent house style for all EDA figures
theme_ea <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = element_blank(),                 # de-clutter
      panel.grid.major.x = element_line(size = 0.25, colour = "grey85"),
      panel.grid.major.y = element_line(size = 0.25, colour = "grey90"),
      plot.title        = element_text(face = "bold", size = base_size + 2),
      plot.subtitle     = element_text(colour = "grey30"),
      plot.caption      = element_text(colour = "grey40", size = base_size - 2),
      strip.text        = element_text(face = "bold")
    )
}
# Compact numeric labels (e.g., 1.2k, 3.4M) — avoids crowded axes
short_num  <- scales::label_number(scale_cut = scales::cut_short_scale())
# Colour palettes used consistently across plots
method_cols <- c(ANAA = "#E86E58", ANLA = "#2A9D8F", ANLE = "#457B9D")
dtype_cols  <- c(bin = "#F4A261", count = "#2A9D8F")

# ---------- 1) Sampling effort over time by analysis method ----------
# What it shows: volume of sampling each year split by EA analysis method.
# Why it matters: motivates date filtering (>=1990) and shows the method mix over time.
p_effort <- metrics_clean %>%
  dplyr::mutate(year = lubridate::year(SAMPLE_DATE)) %>%
  dplyr::group_by(year, ANALYSIS_METHOD) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  ggplot2::ggplot(ggplot2::aes(year, n, fill = ANALYSIS_METHOD)) +
  ggplot2::geom_col(width = 0.9) +
  ggplot2::scale_fill_manual(values = method_cols, name = "Analysis Method") +
  ggplot2::scale_x_continuous(breaks = seq(1990, lubridate::year(Sys.Date()) + 1, 5)) +
  ggplot2::scale_y_continuous(labels = short_num) +
  ggplot2::labs(
    title = "Sampling Effort Over Time by Analysis Method",
    subtitle = "Environment Agency macroinvertebrate samples (≥1990)",
    x = "Year", y = "Number of samples",
    caption = "Note: Around 2000, some ANLA remain categorical; method labelling not a hard cut-off."
  ) +
  theme_ea()
print(p_effort)

# ---------- 2) Frequency of selected families ----------
# What it shows: overall counts of records for the four target families.
# Why it matters: supports the choice of “not too rare, not too common” families.
p_freq <- whpt_clean %>%
  dplyr::filter(EQ_TAXON_UNIT %in% families) %>%
  dplyr::count(EQ_TAXON_UNIT, name = "n") %>%
  dplyr::mutate(EQ_TAXON_UNIT = forcats::fct_reorder(EQ_TAXON_UNIT, n)) %>%
  ggplot2::ggplot(ggplot2::aes(n, EQ_TAXON_UNIT)) +
  ggplot2::geom_col(fill = "#3D5A80") +
  ggplot2::scale_x_continuous(labels = short_num) +
  ggplot2::labs(
    title = "Frequency of Selected Families in the EA Dataset",
    x = "Total records", y = "Family",
    caption = "Counts from harmonised WHPT family table."
  ) +
  theme_ea()
print(p_freq)

# ---------- 3) Categorical vs numerical proportion over time (by family) ----------
# What it shows: yearly mix of categorical vs numeric data for each family.
# Why it matters: demonstrates the ANLA/ANLE/ANAA transition and justifies
#                 modelling strategy (ordered categorical vs censored/count).
p_dtype <- roster_flagged %>%
  dplyr::left_join(whpt_clean %>% dplyr::select(SAMPLE_ID, EQ_TAXON_UNIT), by = "SAMPLE_ID") %>%
  dplyr::filter(EQ_TAXON_UNIT %in% families) %>%
  dplyr::mutate(year = lubridate::year(SAMPLE_DATE)) %>%
  dplyr::group_by(EQ_TAXON_UNIT, year, data_type) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop_last") %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(ggplot2::aes(year, prop, fill = data_type)) +
  ggplot2::geom_col(position = "fill", width = 0.95) +
  ggplot2::facet_wrap(~ EQ_TAXON_UNIT, ncol = 2, scales = "free_x") +
  ggplot2::scale_fill_manual(values = dtype_cols, name = "Data Type",
                    labels = c("Categorical (bin)", "Numerical (count)")) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_x_continuous(breaks = seq(1990, lubridate::year(Sys.Date()) + 1, 5)) +
  ggplot2::labs(
    title = "Categorical vs Numerical Proportion Over Time",
    subtitle = "Per-sample type from TAXA placeholders vs numerics with ANLE/ANAA/ANLA overrides",
    x = "Year", y = "Proportion of samples",
    caption = "Mixed samples treated as categorical (safe rule)."
  ) +
  theme_ea()
print(p_dtype)

# ---------- 4) Log–log abundance distributions by family & method (no warnings) ----------
# What it shows: distribution of non-zero abundances by method for each family on log–log axes.
# Why it matters: highlights categorical spikes (3/33/333...) vs smoother numeric counts.
abund_df <- whpt_clean %>%
  dplyr::filter(EQ_TAXON_UNIT %in% families) %>%
  dplyr::left_join(metrics_clean %>% dplyr::select(SAMPLE_ID, ANALYSIS_METHOD), by = "SAMPLE_ID") %>%
  dplyr::filter(!is.na(ANALYSIS_METHOD), !is.na(TOTAL_NUMBER), TOTAL_NUMBER > 0)

p_hist <- ggplot2::ggplot(abund_df, ggplot2::aes(x = TOTAL_NUMBER)) +
  ggplot2::geom_histogram(
    bins = 40,
    ggplot2::aes(y = ggplot2::after_stat(ifelse(count == 0, NA, count))),  # drop empty bins before log transform
    alpha = 0.9, fill = "#F4A261", colour = "white", linewidth = 0.15
  ) +
  ggplot2::scale_x_log10(labels = short_num) +
  ggplot2::scale_y_log10(labels = short_num) +
  ggplot2::facet_grid(EQ_TAXON_UNIT ~ ANALYSIS_METHOD, scales = "free_y") +
  ggplot2::labs(
    title = "Log–Log Abundance Distribution by Family & Analysis Method",
    x = "Total number (log scale)", y = "Frequency (log scale)",
    caption = "Zeros removed; empty bins omitted for log scaling. Spikes often indicate categorical (3/33/333...) bins."
  ) +
  theme_ea()
print(p_hist)

# ---------- 5) Seasonality of sampling (by method) ----------
# What it shows: months when sampling occurs, split by analysis method.
# Why it matters: validates spring/autumn focus and the season factor used in models.
p_season <- metrics_clean %>%
  dplyr::mutate(month = factor(lubridate::month(SAMPLE_DATE, label = TRUE), ordered = TRUE)) %>%
  dplyr::group_by(ANALYSIS_METHOD, month) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  ggplot2::ggplot(ggplot2::aes(month, n)) +
  ggplot2::geom_col(fill = "#2A9D8F") +
  ggplot2::facet_wrap(~ ANALYSIS_METHOD, ncol = 3, scales = "free_y") +
  ggplot2::scale_y_continuous(labels = short_num) +
  ggplot2::labs(
    title = "Seasonality of Sampling by Analysis Method",
    x = "Month", y = "Number of samples",
    caption = "Sampling concentrates in spring (Mar–May) and autumn (Sep–Nov)."
  ) +
  theme_ea()
print(p_season)


# ---------- 6) Presence rate over time (per family) ----------
# helper: presence table for one family using site roster where family ever observed
# Rationale: compute annual presence proportion only across sites that have
#            recorded the family at least once (avoids inflating absences at
#            sites where the family never occurs).
presence_rate_tbl <- function(fam) {
  sites_with <- whpt_clean %>%
    dplyr::filter(EQ_TAXON_UNIT == fam) %>%
    dplyr::distinct(SITE_ID) %>% dplyr::pull()

  roster <- metrics_clean %>%
    dplyr::filter(SITE_ID %in% sites_with) %>%           # only sites where fam seen at least once
    dplyr::select(SITE_ID, SAMPLE_ID, SAMPLE_DATE)

  counts <- whpt_clean %>%
    dplyr::filter(EQ_TAXON_UNIT == fam) %>%
    dplyr::select(SITE_ID, SAMPLE_ID, TOTAL_NUMBER)

  roster %>%
    dplyr::left_join(counts, by = c("SITE_ID","SAMPLE_ID")) %>%
    dplyr::mutate(TOTAL_NUMBER = tidyr::replace_na(TOTAL_NUMBER, 0L),
           present      = TOTAL_NUMBER > 0,
           year         = lubridate::year(SAMPLE_DATE)) %>%
    dplyr::group_by(EQ_TAXON_UNIT = fam, year) %>%
    dplyr::summarise(p_present = mean(present),
              n_samples = dplyr::n(), .groups = "drop")
}

# Bind presence-rate time series for all requested families
presence_df <- dplyr::bind_rows(lapply(families, presence_rate_tbl))

# Plot annual presence proportion with 0–100% scale per panel
p_presence <- presence_df %>%
  ggplot2::ggplot(ggplot2::aes(year, p_present)) +
  ggplot2::geom_line(linewidth = 0.9, colour = "#3D5A80") +
  ggplot2::geom_point(size = 1.6, colour = "#3D5A80") +
  ggplot2::facet_wrap(~ EQ_TAXON_UNIT, ncol = 2) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  ggplot2::scale_x_continuous(breaks = seq(1990, lubridate::year(Sys.Date()) + 1, 5)) +
  ggplot2::labs(
    title = "Presence Rate Over Time (per Family)",
    subtitle = "Share of samples with the family present, using only sites where that family has ever been observed",
    x = "Year", y = "Presence rate",
    caption = "Zeros added for non-detections within the family’s site roster."
  ) +
  ggplot2::theme_minimal(base_size = 12)
print(p_presence)

```

```{r}
# -------------------------------------------------------------------
# EDA add-ons for abundance & data-type patterns
# (Uses: whpt_clean, roster_flagged, sites_clean, families,
#        plus `theme_ea()` and `short_num` defined earlier.)
# -------------------------------------------------------------------

# ----------------------------
# Helpers
# ----------------------------

# Return the roster (metrics rows) for sites where a given family
# has been observed at least once (so we add zeros correctly).
family_roster <- function(fam) {
  sites_with <- whpt_clean %>%
    filter(EQ_TAXON_UNIT == fam) %>% 
    distinct(SITE_ID) %>% 
    pull()
  roster_flagged %>% filter(SITE_ID %in% sites_with)
}

# Join family counts to that roster and fill non-detections with 0.
# (Qualified dplyr::select to avoid conflicts with other packages.)
join_counts_zero <- function(roster_tbl, fam) {
  counts <- whpt_clean %>%
    filter(EQ_TAXON_UNIT == fam) %>%
    dplyr::select(SITE_ID, SAMPLE_ID, TOTAL_NUMBER)
  roster_tbl %>%
    left_join(counts, by = c("SITE_ID","SAMPLE_ID")) %>%
    mutate(TOTAL_NUMBER = replace_na(TOTAL_NUMBER, 0L))
}

# ----------------------------
# 1) Numeric-only abundance trend (median + IQR) per family
#    • Uses only samples classified as numeric ("count")
#    • Adds zeros (explicit absences) within the family’s site roster
# ----------------------------
num_trend <- bind_rows(lapply(families, function(fam) {
  roster <- family_roster(fam) %>% filter(data_type == "count")
  join_counts_zero(roster, fam) %>%
    mutate(year = year(SAMPLE_DATE)) %>%
    group_by(EQ_TAXON_UNIT = fam, year) %>%
    summarise(
      n   = n(),
      med = median(TOTAL_NUMBER),
      q25 = quantile(TOTAL_NUMBER, 0.25),
      q75 = quantile(TOTAL_NUMBER, 0.75),
      .groups = "drop"
    )
}))

p_num_trend <- ggplot(num_trend, aes(year, med)) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, fill = "#2A9D8F") +
  geom_line(colour = "#2A9D8F", linewidth = 0.9) +
  facet_wrap(~ EQ_TAXON_UNIT, ncol = 2, scales = "free_y") +
  scale_y_continuous(labels = short_num) +
  scale_x_continuous(breaks = seq(1990, year(Sys.Date()) + 1, 5)) +
  labs(
    title    = "Numeric-only Abundance Trend (Median with IQR)",
    subtitle = "Counts from samples classified as numeric (data_type = 'count'), zeros included",
    x = "Year", y = "Abundance (median, with IQR)"
  ) +
  theme_ea()
print(p_num_trend)

# ----------------------------
# 2) Categorical composition over time (AB0/AB1/AB2/AB3+)
#    • Uses only samples classified as categorical ("bin")
#    • Collapses high bins to AB3+ and shows proportions per year
# ----------------------------
cat_comp <- bind_rows(lapply(families, function(fam) {
  roster <- family_roster(fam) %>% filter(data_type == "bin")
  df <- join_counts_zero(roster, fam) %>%
    mutate(
      year = year(SAMPLE_DATE),
      ord  = case_when(
        TOTAL_NUMBER == 0L   ~ "AB0",
        TOTAL_NUMBER == 1L   ~ "AB1",
        TOTAL_NUMBER == 3L   ~ "AB2",
        TOTAL_NUMBER >= 33L  ~ "AB3+",
        TRUE                 ~ NA_character_
      )
    ) %>%
    filter(!is.na(ord)) %>%
    group_by(EQ_TAXON_UNIT = fam, year, ord) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(EQ_TAXON_UNIT, year) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  df
}))

p_cat_comp <- ggplot(cat_comp, aes(year, prop, fill = ord)) +
  geom_col(width = 0.95, position = "fill") +
  facet_wrap(~ EQ_TAXON_UNIT, ncol = 2) +
  scale_fill_manual(
    values = c(AB0 = "#E0E0E0", AB1 = "#F4A261", AB2 = "#E76F51", `AB3+` = "#2A9D8F"),
    name   = "Category"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(1990, year(Sys.Date()) + 1, 5)) +
  labs(
    title    = "Categorical Composition Over Time (AB0/AB1/AB2/AB3+)",
    subtitle = "Only samples classified as categorical (data_type = 'bin')",
    x = "Year", y = "Proportion"
  ) +
  theme_ea()
print(p_cat_comp)

# ----------------------------
# 3) Rounding fingerprint for numeric counts (last digit distribution)
#    • Focus on TOTAL_NUMBER ≥ 10 where 1 s.f. rounding is typical
#    • Optional pre/post-2012 era split for visual comparison
# ----------------------------
round_fp <- bind_rows(lapply(families, function(fam) {
  roster <- family_roster(fam) %>% filter(data_type == "count")
  join_counts_zero(roster, fam) %>%
    filter(TOTAL_NUMBER >= 10) %>%                           # restrict to 1 s.f. zone
    mutate(
      last_digit = TOTAL_NUMBER %% 10,
      era        = if_else(year(SAMPLE_DATE) < 2012, "<2012", "≥2012")
    ) %>%
    group_by(EQ_TAXON_UNIT = fam, era, last_digit) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(EQ_TAXON_UNIT, era) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
}))

p_round_fp <- ggplot(round_fp, aes(factor(last_digit), prop)) +
  geom_col(fill = "#3D5A80") +
  facet_grid(EQ_TAXON_UNIT ~ era) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title    = "Rounding Fingerprint for Numeric Counts",
    subtitle = "Distribution of last digit for TOTAL_NUMBER ≥ 10 (expect spikes at 0 and often 5)",
    x = "Last digit", y = "Proportion"
  ) +
  theme_ea()
print(p_round_fp)

# ----------------------------
# 4) Seasonal abundance (numeric counts, sqrt scale)
#    • Boxplots of sqrt(count) by month for each family
# ----------------------------
season_df <- bind_rows(lapply(families, function(fam) {
  roster <- family_roster(fam) %>% filter(data_type == "count")
  join_counts_zero(roster, fam) %>%
    mutate(
      month         = factor(month(SAMPLE_DATE, label = TRUE), ordered = TRUE),
      EQ_TAXON_UNIT = fam
    )
}))

p_season_abund <- ggplot(season_df, aes(month, sqrt(TOTAL_NUMBER))) +
  geom_boxplot(outlier.alpha = 0.25, fill = "#A8DADC") +
  facet_wrap(~ EQ_TAXON_UNIT, ncol = 2, scales = "free_y") +
  labs(
    title = "Seasonal Distribution of Numeric Abundance (sqrt scale)",
    x = "Month", y = "sqrt(Abundance)"
  ) +
  theme_ea()
print(p_season_abund)

# ----------------------------
# 5) Regional heterogeneity: data-type mix by EA area (top 12)
#    • Shows proportion of categorical vs numeric samples by area
#    • Uses sites_clean to fetch REPORTING_AREA for each SITE_ID
# ----------------------------
area_mix <- roster_flagged %>%
  left_join(sites_clean %>% dplyr::select(SITE_ID, REPORTING_AREA), by = "SITE_ID") %>%
  filter(!is.na(REPORTING_AREA)) %>%
  count(REPORTING_AREA, data_type, name = "n") %>%
  group_by(REPORTING_AREA) %>%
  mutate(total = sum(n), prop = n / total) %>%
  ungroup() %>%
  slice_max(order_by = total, n = 12, with_ties = FALSE) %>%
  mutate(REPORTING_AREA = fct_reorder(REPORTING_AREA, total))

p_area_mix <- ggplot(area_mix, aes(REPORTING_AREA, prop, fill = data_type)) +
  geom_col(position = "fill") +
  coord_flip() +
  scale_fill_manual(values = c(bin = "#F4A261", count = "#2A9D8F"), name = "Data Type") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Data-Type Mix by EA Reporting Area (Top 12 by sample volume)",
    x = "Reporting area", y = "Proportion of samples"
  ) +
  theme_ea()
print(p_area_mix)

```
```{r, fig.width=9, fig.height=6}
# ===================================================================
# Lightweight theme + number formatter 
#   - 1) ANLA usage shift over time
#   - 2) Ordinal sparsity (AB0/AB1/AB2/AB3+) — full + zoomed views
#   - 3) Rounding fingerprint for numeric counts
#   - 4) Site-level presence map for one family
#   - 5) Same-day replicates (pre de-dup) by year
# ===================================================================

# --- Short number axis formatter (e.g., 1k, 2.5k) ----------------------------
if (!exists("short_num")) short_num <- label_number(scale_cut = cut_short_scale())

# ===================================================================
# 1) ANLA meaning shift through time (categorical vs numeric)
#    • Restrict to ANALYSIS_METHOD == "ANLA"
#    • Compute yearly share of bin vs count
# ===================================================================
p_anla_shift <- roster_flagged %>%
  mutate(year = year(SAMPLE_DATE)) %>%
  filter(ANALYSIS_METHOD == "ANLA") %>%
  count(year, data_type, name = "n") %>%
  group_by(year) %>% mutate(prop = n / sum(n)) %>% ungroup() %>%
  ggplot(aes(year, prop, fill = data_type)) +
  geom_col(width = 0.95, position = "fill") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(1990, year(Sys.Date()) + 1, 5)) +
  scale_fill_manual(values = c(bin = "#F4A261", count = "#2A9D8F"),
                    name = "Data Type",
                    labels = c("Categorical (bin)", "Numerical (count)")) +
  labs(
    title = "How ANLA is Used Over Time",
    subtitle = "Share of ANLA-labelled samples behaving as categorical vs numeric",
    x = "Year", y = "Proportion"
  ) + theme_ea()
print(p_anla_shift)

# ===================================================================
# 2) Ordinal sparsity: AB0 / AB1 / AB2 / AB3+ by family (all years)
#    • Work on categorical samples only (data_type == 'bin')
#    • Add explicit zeros; collapse high bins to AB3+
#    Output: (A) full stacked, (B) zoomed to small categories
# ===================================================================

# -- Build once: overall category proportions per family ----------------------
cat_sparsity_df <- bind_rows(lapply(families, function(fam) {
  roster <- roster_flagged %>% filter(data_type == "bin")                 # only categorical samples
  counts <- whpt_clean %>%
    filter(EQ_TAXON_UNIT == fam) %>%
    dplyr::select(SITE_ID, SAMPLE_ID, TOTAL_NUMBER)                       # minimal join fields
  roster %>%
    left_join(counts, by = c("SITE_ID","SAMPLE_ID")) %>%                  # add counts for this family
    mutate(TOTAL_NUMBER = replace_na(TOTAL_NUMBER, 0L),                   # explicit zeros for non-detections
           ord = case_when(                                               # map raw values to ordinal bins
             TOTAL_NUMBER == 0L  ~ "AB0",
             TOTAL_NUMBER == 1L  ~ "AB1",
             TOTAL_NUMBER == 3L  ~ "AB2",
             TOTAL_NUMBER >= 33L ~ "AB3+",
             TRUE ~ NA_character_
           )) %>%
    filter(!is.na(ord)) %>%                                               # drop any unexpected values
    count(EQ_TAXON_UNIT = fam, ord, name = "n")                           # counts per bin
})) %>%
  group_by(EQ_TAXON_UNIT) %>%
  mutate(prop = n / sum(n)) %>%                                           # within-family proportions
  ungroup()

# -- (A) Full stacked view ----------------------------------------------------
p_cat_sparsity_better <- ggplot(cat_sparsity_df,
                                aes(x = fct_relevel(ord,"AB0","AB1","AB2","AB3+"),
                                    y = prop, fill = ord)) +
  geom_col(width = 0.85) +                                                # 100% stacked bars
  facet_wrap(~ EQ_TAXON_UNIT, ncol = 2) +                                 # one panel per family
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(AB0="#DADADA", AB1="#F4A261", AB2="#E76F51", `AB3+`="#2A9D8F"),
                    name = "Category") +
  labs(title = "Ordinal Category Mass by Family (full scale)",
       subtitle = "AB0 often dominates; see zoomed panel below to inspect small categories",
       x = "Ordinal class", y = "Proportion") +
  theme_ea()

# -- (B) Zoomed view (exclude AB0) -------------------------------------------
ab0_lab <- cat_sparsity_df %>%
  filter(ord == "AB0") %>%
  transmute(EQ_TAXON_UNIT, lab = paste0("AB0 = ", scales::percent(prop, 1)))  # facet annotation

p_cat_sparsity_zoom <- cat_sparsity_df %>%
  filter(ord != "AB0") %>%                                                # focus on small categories
  ggplot(aes(x = fct_relevel(ord,"AB1","AB2","AB3+"), y = prop, fill = ord)) +
  geom_col(width = 0.85) +
  geom_text(aes(label = scales::percent(prop, 1)),                         # label bars with %
            vjust = -0.25, size = 3, colour = "grey20") +
  geom_text(data = ab0_lab, aes(x = "AB3+", y = 0.095, label = lab),       # show AB0 share per facet
            inherit.aes = FALSE, hjust = 1, size = 3.2, colour = "grey30") +
  facet_wrap(~ EQ_TAXON_UNIT, ncol = 2) +
  coord_cartesian(ylim = c(0, 0.10)) +                                    # zoom to 0–10%
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = c(AB1="#F4A261", AB2="#E76F51", `AB3+`="#2A9D8F"),
                    name = "Category") +
  labs(title = "Ordinal Category Mass by Family — zoomed to small categories",
       subtitle = "Y-axis limited to 0–10% of all samples; facet text shows AB0 share",
       x = "Ordinal class", y = "Proportion (of all samples)") +
  theme_ea()

print(p_cat_sparsity_better)
print(p_cat_sparsity_zoom)

# ===================================================================
# 3) Rounding mix (numeric counts): 1–9 exact vs tens vs 5s vs other
#    • Work on numeric samples only (data_type == 'count')
#    • Categorise rounding signature and plot yearly proportions
# ===================================================================
p_rounding_mix <- bind_rows(lapply(families, function(fam) {
  roster <- roster_flagged %>% filter(data_type == "count")
  counts <- whpt_clean %>%
    filter(EQ_TAXON_UNIT == fam) %>%
    dplyr::select(SITE_ID, SAMPLE_ID, TOTAL_NUMBER)  # qualified select to avoid masking
  roster %>%
    left_join(counts, by = c("SITE_ID","SAMPLE_ID")) %>%
    mutate(TOTAL_NUMBER = replace_na(TOTAL_NUMBER, 0L),
           group = case_when(
             TOTAL_NUMBER >= 1  & TOTAL_NUMBER <= 9  ~ "1–9 exact",
             TOTAL_NUMBER >= 10 & TOTAL_NUMBER %% 10 == 0 ~ "rounded to 10s",
             TOTAL_NUMBER >= 10 & TOTAL_NUMBER %% 10 == 5 ~ "ends with 5",
             TOTAL_NUMBER >= 10                           ~ "other ≥10",
             TRUE                                          ~ NA_character_
           ),
           year = year(SAMPLE_DATE)) %>%
    filter(!is.na(group)) %>%
    count(EQ_TAXON_UNIT = fam, year, group, name = "n") %>%
    group_by(EQ_TAXON_UNIT, year) %>% mutate(prop = n/sum(n)) %>% ungroup()
})) %>%
  ggplot(aes(year, prop, fill = group)) +
  geom_col(width = 0.95, position = "fill") +
  facet_wrap(~ EQ_TAXON_UNIT, ncol = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_x_continuous(breaks = seq(1990, year(Sys.Date()) + 1, 5)) +
  scale_fill_manual(values = c("1–9 exact"="#8ecae6","rounded to 10s"="#219ebc",
                               "ends with 5"="#ffb703","other ≥10"="#fb8500"),
                    name = "Numeric rounding pattern") +
  labs(
    title = "Rounding Fingerprint in Numeric Counts",
    subtitle = "Share of numeric observations by rounding signature (per family, per year)",
    x = "Year", y = "Proportion"
  ) + theme_ea()
print(p_rounding_mix)

# ===================================================================
# 4) Site-level presence map for one family (effort & detection)
#    • Pick first family; compute visits & detection rate per site
#    • Join EA coordinates; scatter by easting/northing
# ===================================================================
fam_map <- families[1]  # pick any family to display
site_presence <- {
  sites_with <- whpt_clean %>%
    filter(EQ_TAXON_UNIT == fam_map) %>% distinct(SITE_ID) %>% pull()
  roster <- roster_flagged %>% filter(SITE_ID %in% sites_with)
  counts <- whpt_clean %>% filter(EQ_TAXON_UNIT == fam_map) %>%
    dplyr::select(SITE_ID, SAMPLE_ID, TOTAL_NUMBER)  # qualified select

  roster %>%
    left_join(counts, by = c("SITE_ID","SAMPLE_ID")) %>%
    mutate(TOTAL_NUMBER = replace_na(TOTAL_NUMBER, 0L),
           present = TOTAL_NUMBER > 0) %>%
    group_by(SITE_ID) %>%
    summarise(visits = n(), p_present = mean(present), .groups = "drop") %>%
    left_join(sites_clean %>% dplyr::select(SITE_ID, FULL_EASTING, FULL_NORTHING),  # qualified select
              by = "SITE_ID") %>%
    filter(!is.na(FULL_EASTING), !is.na(FULL_NORTHING))
}

p_map <- ggplot(site_presence, aes(FULL_EASTING, FULL_NORTHING)) +
  geom_point(aes(size = visits, colour = p_present), alpha = 0.8) +
  scale_size_continuous(range = c(0.8, 5), name = "Visits") +
  scale_colour_viridis_c(name = "Presence rate", labels = percent_format(accuracy = 1)) +
  coord_equal() +
  labs(
    title = paste0("Site Effort & Presence for ", fam_map),
    subtitle = "Size = number of visits; colour = detection rate at that site",
    x = "Easting", y = "Northing",
    caption = "British National Grid coordinates; for cartographic basemaps use {sf}."
  ) + theme_ea()
print(p_map)

# ===================================================================
# 5) Same-day replicates (pre-dedup) by year
#    • Re-derive a minimal metrics table from Arrow (as in analysis)
#    • Count days with >1 SAMPLE_ID per site/date before your de-dup
# ===================================================================
if (!exists("DATE_FLOOR")) DATE_FLOOR <- as.Date("1990-01-01")
if (!exists("USE_S3PO"))   USE_S3PO   <- TRUE

replicates_year <- {
  m_raw <- metrics_raw %>%
    dplyr::mutate(
      SITE_ID     = arrow::cast(SITE_ID,   arrow::int32()),
      SAMPLE_ID   = arrow::cast(SAMPLE_ID, arrow::int32()),
      SAMPLE_DATE = arrow::cast(SAMPLE_DATE, arrow::date32())
    ) %>%
    dplyr::select(SITE_ID, SAMPLE_ID, SAMPLE_DATE, SAMPLE_METHOD) %>%
    dplyr::collect() %>%
    dplyr::mutate(SAMPLE_DATE = as.Date(SAMPLE_DATE)) %>%
    dplyr::filter(SAMPLE_DATE >= DATE_FLOOR) %>%
    { if (USE_S3PO) dplyr::filter(., SAMPLE_METHOD == "S3PO") else . }

  m_raw %>%
    group_by(SITE_ID, SAMPLE_DATE) %>%
    summarise(n_ids = n_distinct(SAMPLE_ID), .groups = "drop") %>%
    mutate(year = year(SAMPLE_DATE),
           is_replicate = n_ids > 1) %>%
    group_by(year) %>%
    summarise(replicate_days = sum(is_replicate),
              total_days = n(),
              prop = replicate_days / total_days,
              .groups = "drop")
}

p_reps <- ggplot(replicates_year, aes(year, replicate_days)) +
  geom_col(fill = "#9b2226") +
  scale_x_continuous(breaks = seq(1990, year(Sys.Date()) + 1, 5)) +
  scale_y_continuous(labels = short_num) +
  labs(
    title = "Same-day Replicates Before De-duplication",
    subtitle = "Count of site-days with >1 SAMPLE_ID (filters matched to your analysis set)",
    x = "Year", y = "Replicate site-days"
  ) + theme_ea()
print(p_reps)

```




```{r, fig.width=9, fig.height=6}
# ================================================================
# Three EDA plots (presence/absence focus)
#   1) Presence rate over time by season
#   2) Colonisation / extinction between windows
#   3) Zero proportion over time by season
# ------------------------------------------------
# Expects in memory:
#   - roster_flagged: deduped roster with SAMPLE_DATE, SITE_ID, SAMPLE_ID, season_f
#   - whpt_clean   : harmonised WHPT with EQ_TAXON_UNIT, SITE_ID, SAMPLE_ID, TOTAL_NUMBER
#   - families     : character vector of family names
# ================================================================

# ---- Helper: build a joined table for one family ----
# Adds explicit zeros, presence flag, year & season.
fam_join <- function(fam) {
  sites_with <- whpt_clean %>%
    dplyr::filter(EQ_TAXON_UNIT == fam) %>%
    dplyr::distinct(SITE_ID) %>% dplyr::pull()

  roster <- roster_flagged %>%
    dplyr::filter(SITE_ID %in% sites_with) %>%
    dplyr::select(SITE_ID, SAMPLE_ID, SAMPLE_DATE, season_f)

  counts <- whpt_clean %>%
    dplyr::filter(EQ_TAXON_UNIT == fam) %>%
    dplyr::select(SITE_ID, SAMPLE_ID, TOTAL_NUMBER)

  roster %>%
    dplyr::left_join(counts, by = c("SITE_ID","SAMPLE_ID")) %>%
    dplyr::mutate(
      TOTAL_NUMBER  = tidyr::replace_na(TOTAL_NUMBER, 0L),
      present       = TOTAL_NUMBER > 0L,
      year          = lubridate::year(SAMPLE_DATE),
      EQ_TAXON_UNIT = fam
    )
}

# ================================================================
# 1) Presence rate over time by season
#    - Site roster restricted to where the family has ever been observed
#    - Presence = share of samples with TOTAL_NUMBER > 0
# ================================================================
presence_season_df <- dplyr::bind_rows(lapply(families, fam_join)) %>%
  dplyr::group_by(EQ_TAXON_UNIT, season_f, year) %>%
  dplyr::summarise(p_present = mean(present), n = dplyr::n(), .groups = "drop")

p_presence_season <- ggplot(presence_season_df,
                            aes(year, p_present, colour = season_f)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ EQ_TAXON_UNIT, ncol = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(1990, year(Sys.Date()) + 1, 5)) +
  scale_colour_manual(values = c(spring = "#2A9D8F", autumn = "#E76F51"),
                      name = "Season") +
  labs(title = "Presence rate over time by season",
       x = "Year", y = "Presence rate") +
  theme_ea()

print(p_presence_season)

# ================================================================
# 2) Colonisation / extinction between windows
#    - Windows: Early = 1990–2005, Recent = 2006–2024 (inclusive)
#    - Consider sites that were sampled in BOTH windows
#    - Classification per site:
#         Colonised (0 -> 1), Extinct (1 -> 0),
#         Stayed absent (0 -> 0), Stayed present (1 -> 1)
# ================================================================
early_years  <- 1990:2005
recent_years <- 2006:2024

col_ext_df <- dplyr::bind_rows(lapply(families, fam_join)) %>%
  dplyr::mutate(period = dplyr::case_when(
    year %in% early_years  ~ "early",
    year %in% recent_years ~ "recent",
    TRUE ~ NA_character_
  )) %>%
  dplyr::filter(!is.na(period)) %>%
  dplyr::group_by(EQ_TAXON_UNIT, SITE_ID, period) %>%
  dplyr::summarise(p_any = as.integer(any(present)), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = period, values_from = p_any) %>%
  # Keep sites sampled in BOTH windows
  dplyr::filter(!is.na(early), !is.na(recent)) %>%
  dplyr::mutate(class = dplyr::case_when(
    early == 0L & recent == 1L ~ "Colonised",
    early == 1L & recent == 0L ~ "Extinct",
    early == 0L & recent == 0L ~ "Stayed absent",
    early == 1L & recent == 1L ~ "Stayed present",
    TRUE ~ NA_character_
  )) %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::count(EQ_TAXON_UNIT, class, name = "n") %>%
  dplyr::group_by(EQ_TAXON_UNIT) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  dplyr::ungroup()

p_col_ext <- ggplot(col_ext_df, aes(class, prop, fill = class)) +
  geom_col(width = 0.85) +
  facet_wrap(~ EQ_TAXON_UNIT, ncol = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_fill_manual(values = c("Colonised" = "#2A9D8F",
                               "Extinct" = "#E76F51",
                               "Stayed absent" = "#BDBDBD",
                               "Stayed present" = "#345995"),
                    guide = "none") +
  labs(title = "Colonisation / extinction between windows",
       subtitle = "Sites sampled in both 1990–2005 and 2006–2024",
       x = NULL, y = "Share of sites") +
  theme_ea()

print(p_col_ext)

# ================================================================
# 3) Zero proportion over time by season
#    - Same denominator as (1): roster restricted to sites where family ever seen
#    - Metric = share of samples with TOTAL_NUMBER == 0
# ================================================================
zero_season_df <- dplyr::bind_rows(lapply(families, fam_join)) %>%
  dplyr::group_by(EQ_TAXON_UNIT, season_f, year) %>%
  dplyr::summarise(p_zero = mean(TOTAL_NUMBER == 0L), n = dplyr::n(), .groups = "drop")

p_zero_season <- ggplot(zero_season_df, aes(year, p_zero, colour = season_f)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ EQ_TAXON_UNIT, ncol = 2) +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(1990, year(Sys.Date()) + 1, 5)) +
  scale_colour_manual(values = c(spring = "#2A9D8F", autumn = "#E76F51"),
                      name = "Season") +
  labs(title = "Zero proportion over time by season",
       x = "Year", y = "Share zeros") +
  theme_ea()

print(p_zero_season)
