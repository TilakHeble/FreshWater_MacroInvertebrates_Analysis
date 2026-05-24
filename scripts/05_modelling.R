#4. MODELLING PHASE

#4.1 PRESENCE/ABSENCE BINOMIAL GAMM MODEL
# ================================================
# PHASE 4 — MODEL 1 : Presence/Absence GAMM (one family)

#   Goal: fit a binomial GAMM for presence/absence of a chosen family,
#         with season-specific temporal smooths and a site random intercept.
#   Inputs expected in memory: models_data (from Phase 3).
#   Output: pa_model object + quick diagnostics.
# ================================================

# -------- choose family --------
fam_name <- "Aphelocheiridae"          # <-- change to any of your families
pa_df    <- models_data[[fam_name]]$pa  # presence/absence table built in Phase 3

# -------- optional: restrict to S3PO (EA-recommended) --------
USE_S3PO <- TRUE
if (USE_S3PO && "SAMPLE_METHOD" %in% names(pa_df)) {
  pa_df <- dplyr::filter(pa_df, SAMPLE_METHOD == "S3PO")  # keep only S3PO samples if column exists
}

# -------- sensible k based on span of years --------
n_years <- dplyr::n_distinct(lubridate::year(pa_df$SAMPLE_DATE))  # number of distinct calendar years
k_time  <- min(20, max(8, round(0.6 * n_years)))                  # cap between 8 and 20, ~0.6*years

# -------- fit binomial GAMM --------
pa_model <- bam(
  taxon_present ~
    season_f +                                   # season-specific intercepts (spring vs autumn)
    s(decimal_date, by = season_f, k = k_time) + # two time smooths, one per season
    s(SITE_ID.F, bs = "re"),                     # random intercept for site (controls spatial heterogeneity)
  family   = binomial(link = "logit"),           # presence/absence with logit link
  data     = pa_df,                              # modelling table
  method   = "fREML",                            # stable, fast REML for smoothing selection
  discrete = TRUE,                               # speed-up for large data
  select   = TRUE,    # shrink unnecessary wiggle (adds penalties to drop unneeded basis functions)
  gamma    = 1.2      # mild extra penalty to reduce overfitting
)

# -------- minimal diagnostics / outputs --------
print(summary(pa_model))       # EDFs per smooth, parametric terms, deviance explained, etc.
gam.check(pa_model)            # residual checks + k-index (are bases big enough?)
draw(pa_model, select = 1:2)   # partial effect plots for the two season-specific time smooths


# ================================================
# PA model quick Model Diagnostics
# Requires: pa_model (bam fit), pa_df (training data)
# Purpose: run lightweight checks + a small set of
#          quantitative summaries.
# ================================================

# 1) Basis/smooth adequacy
gam.check(pa_model, rep = 0)   # rep=0 avoids expensive re-smoothing; reports k-index and residual patterns

# 2) Concurvity — robust to return shape
con_obj <- try(mgcv::concurvity(pa_model, full = TRUE), silent = TRUE)  # can return list or matrix depending on terms
if (!inherits(con_obj, "try-error")) {
  if (is.list(con_obj) && !is.null(con_obj$estimate)) {
    cat("\nConcurvity (estimate, rounded):\n")
    print(round(con_obj$estimate, 3))   # typical output: para / smooth blocks
  } else {
    cat("\nConcurvity (raw object):\n")
    print(con_obj)                      # fallback for alternate structures
  }
} else {
  message("concurvity() failed or not available; skipping.")  # don’t fail the whole diagnostic pass
}

# 3) Compact residual diagnostics
appraise(pa_model)  # produces QQ, residual vs fitted, histogram, scale-location (silent if device is off)

# 4) Residual autocorrelation (global)
res_dev <- residuals(pa_model, type = "deviance")  # deviance residuals on the response scale
acf(res_dev[is.finite(res_dev)], na.action = na.pass,
    main = "ACF of deviance residuals")           # check for temporal dependence left over

# 5) Calibration (observed vs predicted by decile)
cal <- tibble(
  phat = fitted(pa_model, type = "response"),     # predicted probabilities
  y    = as.numeric(pa_df$taxon_present)          # 0/1 outcome as numeric
) |>
  mutate(bin = cut(phat, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)) |>  # decile bins
  group_by(bin) |>
  summarise(p_hat = mean(phat), p_obs = mean(y), n = n(), .groups = "drop")        # mean predicted vs observed

ggplot(cal, aes(p_hat, p_obs, size = n)) +
  geom_point(alpha = 0.85) +                                                           # larger points for bigger bins
  geom_abline(slope = 1, intercept = 0, linetype = 2) +                                # perfect-calibration line
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_size_continuous(range = c(2, 8), guide = "none") +
  labs(title = "Calibration: observed vs predicted (deciles)",
       x = "Mean predicted probability", y = "Observed presence rate") +
  theme_minimal(base_size = 12)

# 6) Quick metrics
phat  <- fitted(pa_model, type = "response")   # vector of predicted probabilities
y     <- as.numeric(pa_df$taxon_present)       # 0/1 truth
brier <- mean((y - phat)^2)                    # Brier score (lower is better)
devex <- summary(pa_model)$dev.expl            # deviance explained (fraction)
cat(sprintf("\nBrier score: %.4f | Deviance explained: %.1f%%\n",
            brier, 100 * devex))

# 7) Season-specific time smooths (spring/autumn)
draw(pa_model, select = 1:2)   # partial effects for s(decimal_date):season_fspring & :autumn



# ================================================
# PHASE 4 — Batch PA GAMMs for 4 families
# Goal: fit the same presence/absence GAMM to each family,
#       print compact text diagnostics, and return a result list.
#       (No figures produced; gam.check plots are suppressed.)
# ================================================

# Families to run
families <- c("Aphelocheiridae","Brachycentridae","Odontoceridae","Cordulegastridae")

# Option: EA-recommended S3PO-only
USE_S3PO <- TRUE    # if TRUE, restricts to kick samples with consistent method coding

# ---- Helper: fit one family ----
fit_pa_family <- function(fam_name, md = models_data, use_s3po = TRUE) {
  # Safety: check that the models_data list has this family and a 'pa' table
  stopifnot(!is.null(md[[fam_name]]), "pa" %in% names(md[[fam_name]]))
  df <- md[[fam_name]]$pa   # training data for this family's PA model

  # Filter to S3PO if requested (only if the column exists)
  if (use_s3po && "SAMPLE_METHOD" %in% names(df)) {
    df <- dplyr::filter(df, SAMPLE_METHOD == "S3PO")
  }

  # Build 2-level season if missing (spring/autumn only)
  # Uses month(SAMPLE_DATE) to label spring = Mar–May; autumn = Sep–Nov
  if (!("season_f" %in% names(df))) {
    df <- df %>%
      mutate(m = month(SAMPLE_DATE),
             season_f = factor(if_else(m %in% 3:5, "spring", "autumn"),
                               levels = c("spring","autumn")))
  }

  # Sensible k from span of years (avoid over/under-smoothing)
  k_time <- {
    ny <- dplyr::n_distinct(lubridate::year(df$SAMPLE_DATE))
    min(20, max(8, round(0.6 * ny)))
  }

  # Fit model:
  # - season-specific intercepts (season_f)
  # - season-varying time smooth s(decimal_date, by = season_f)
  # - site random intercept s(SITE_ID.F, bs="re")
  # - binomial logit PA, fREML, discrete=TRUE for speed, select+gamma to regularise
  m <- bam(
    taxon_present ~
      season_f +
      s(decimal_date, by = season_f, k = k_time) +
      s(SITE_ID.F, bs = "re"),
    family   = binomial(link = "logit"),
    data     = df,
    method   = "fREML",
    discrete = TRUE,
    select   = TRUE,
    gamma    = 1.2
  )

  # Compact text summary to console: parametric terms, smooth EDFs, deviance explained
  cat("\n============================\n")
  cat(sprintf("Family: %s\n", fam_name))
  cat("============================\n")
  print(summary(m))

  # gam.check without plotting to screen:
  # send plots to a temp PDF device, close, then delete file
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  gc_out <- try(gam.check(m, rep = 0), silent = TRUE)  # rep=0 => no refits; k-index etc.
  dev.off()
  unlink(tmp)

  # Quick scalar metrics on training data
  phat  <- fitted(m, type = "response")   # predicted presence prob
  y     <- as.numeric(df$taxon_present)   # 0/1 vector
  brier <- mean((y - phat)^2)             # Brier score (lower is better)
  devex <- summary(m)$dev.expl            # fraction deviance explained
  aic   <- AIC(m)                         # AIC for rough comparability
  nobs_ <- nobs(m)                        # effective N used

  cat(sprintf("Brier score: %.4f | Deviance explained: %.1f%% | AIC: %.1f | N: %d\n",
              brier, 100 * devex, aic, nobs_))

  # Return a compact result object (model + meta for later use)
  list(
    family   = fam_name,
    model    = m,
    data     = df,
    k_time   = k_time,
    brier    = brier,
    dev_expl = devex,
    AIC      = aic,
    n        = nobs_
  )
}

# ---- Run for all families ----
# Produces a named list of result objects, one per family
pa_results <- lapply(families, fit_pa_family, md = models_data, use_s3po = USE_S3PO)
names(pa_results) <- families

# ---- Optional: collect a summary table ----
# Binds core metrics across families for quick comparison
pa_summary <- do.call(rbind, lapply(pa_results, function(x) {
  data.frame(
    family   = x$family,
    n        = x$n,
    k_time   = x$k_time,
    brier    = x$brier,
    dev_expl = x$dev_expl,
    AIC      = x$AIC,
    row.names = NULL
  )
}))
print(pa_summary)

# ---- Optional: save models ----
# invisible(lapply(pa_results, function(x) saveRDS(x$model, paste0("pa_model_", x$family, ".rds"))))



# PA model quick Model Diagnostics for all families
# (Produces: QQ plot, residuals vs linear predictor, residual histogram,
# response vs fitted jitter, and season-specific time smooths.)


plot_pa_diagnostics_clean <- function(res) {
  # `res` is one element of pa_results: list(model, family, data, ...)
  m   <- res$model          # mgcv::bam fitted model
  fam <- res$family         # family name string (for titles)
  df  <- res$data           # training data used to fit

  # 1) QQ plot of residuals.
  #    rep = 0 avoids costly reference simulation; pch/cex make points subtle.
  par(mfrow = c(1,1))
  mgcv::qq.gam(m, rep = 0, pch = 19, cex = 0.3, main = paste("QQ plot —", fam))

  # Prepare a tidy frame with the key diagnostics quantities:
  eta   <- predict(m, type = "link")         # linear predictor η
  resid <- residuals(m, type = "deviance")   # deviance residuals
  phat  <- fitted(m, type = "response")      # fitted probabilities on response scale
  y     <- as.numeric(df$taxon_present)      # observed 0/1
  dd    <- data.frame(eta = eta, resid = resid, phat = phat, y = y)

  # 2) Residuals vs linear predictor.
  p_res_lin <- ggplot(dd, aes(eta, resid)) +
    geom_point(alpha = 0.15, size = 0.6) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = paste("Residuals vs Linear Predictor —", fam),
         x = "Linear predictor (link scale)", y = "Deviance residuals") +
    theme_minimal(base_size = 12)
  print(p_res_lin)

  # 3) Histogram of residuals (distribution shape check).
  p_hist <- ggplot(dd, aes(resid)) +
    geom_histogram(bins = 40, fill = "#3D5A80", colour = "white", linewidth = 0.15) +
    labs(title = paste("Histogram of Residuals —", fam),
         x = "Deviance residuals", y = "Frequency") +
    theme_minimal(base_size = 12)
  print(p_hist)

  # 4) Response vs Fitted: jittered points only.
  p_resp_fit <- ggplot(dd, aes(phat, y)) +
    geom_jitter(height = 0.045, width = 0, alpha = 0.15, size = 0.5) +
    labs(title = paste("Response vs Fitted —", fam),
         x = "Fitted probability", y = "Observed (0/1)") +
    theme_minimal(base_size = 12)
  print(p_resp_fit)

  # 5) Season-specific time smooths:
  #    locate smooth terms containing "decimal_date" and draw each.
  sm_ids   <- which(grepl("decimal_date", smooths(m)))
  sm_names <- smooths(m)[sm_ids]
  for (j in seq_along(sm_ids)) {
    nm  <- sm_names[j]
    lab <- if (grepl("spring", nm, ignore.case = TRUE)) "spring"
           else if (grepl("autumn", nm, ignore.case = TRUE)) "autumn" else nm
    p_sm <- draw(m, select = sm_ids[j]) + ggtitle(paste("Time smooth (", lab, ") — ", fam, sep = ""))
    print(p_sm)
  }
}

# Run the diagnostic panels for every fitted PA model in `pa_results`.
# `invisible()` avoids printing the returned list of NULLs to the console.
invisible(lapply(pa_results, plot_pa_diagnostics_clean))



#4.2 ORDINAL(ORDERED-CATEGORICAL) GAMM MODEL
# ================================================================
# Model 2 : ORDINAL (ordered-categorical) MODEL — one family
# Pre-step: inspect class counts & sparsity; then fit ocat GAMM
# ================================================================

# ----- choose family -----
fam_name <- "Aphelocheiridae"          # set the target family (use any from your list)
pa_df    <- models_data[[fam_name]]$pa # start from PA table (has zeros + season_f + data_type)

# -------- optional: restrict to S3PO (EA-recommended) --------
USE_S3PO <- TRUE
if (USE_S3PO && "SAMPLE_METHOD" %in% names(pa_df)) {
  pa_df <- dplyr::filter(pa_df, SAMPLE_METHOD == "S3PO")  # keep S3PO-only if column exists
}

# -------- restrict to categorical samples for ordinal modelling --------
ord_base <- pa_df %>%
  filter(data_type == "bin")  # ordinal model uses the categorical (bin) subset only

# -------- build ordered classes using *intervals* (robust to odd values like 6/66/etc.) --------
# AB0 = 0, AB1 = 1–9, AB2 = 10–99, AB3 = 100–999, AB4 = 1000+
ord_df <- ord_base %>%
  mutate(
    ord_class = case_when(
      TOTAL_NUMBER == 0L                               ~ "AB0 (0)",
      TOTAL_NUMBER > 0L & TOTAL_NUMBER < 10L           ~ "AB1 (1–9)",
      TOTAL_NUMBER >= 10L & TOTAL_NUMBER < 100L        ~ "AB2 (10–99)",
      TOTAL_NUMBER >= 100L & TOTAL_NUMBER < 1000L      ~ "AB3 (100–999)",
      TOTAL_NUMBER >= 1000L                            ~ "AB4 (1000+)",
      TRUE ~ NA_character_
    ),
    ord_class = factor(
      ord_class,
      levels = c("AB0 (0)","AB1 (1–9)","AB2 (10–99)","AB3 (100–999)","AB4 (1000+)")
    )
  ) %>%
  filter(!is.na(ord_class)) %>%         # drop anything outside the intended bins
  mutate(SITE_ID.F = factor(SITE_ID))   # safety: ensure site is a factor for the RE

# --------- Plot & print class counts (sparsity check before fitting) ---------
class_counts <- ord_df %>% count(ord_class, name = "n") %>% 
  mutate(prop = n / sum(n))    # proportions to show how dominant AB0 is, etc.

print(class_counts)            # numeric table in console (n and share per class)

ggplot(class_counts, aes(x = ord_class, y = n)) +
  geom_col(fill = "#3D5A80") +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
            vjust = -0.4, size = 3.3) +
  labs(title = paste("Ordinal class counts —", fam_name),
       x = "Class", y = "Count",
       caption = "Proportions shown above bars. Consider merging AB4→AB3 if very sparse.") +
  theme_minimal(base_size = 12)

# --------- Optional merge rule (AB4→AB3 if sparse); prints a recommendation ---------
min_n    <- 50     # minimum count threshold for highest class
min_prop <- 0.02   # or <2% of data -> merge down
n_ab4    <- class_counts$n[class_counts$ord_class == "AB4 (1000+)"]
p_ab4    <- class_counts$prop[class_counts$ord_class == "AB4 (1000+)"]
if (length(n_ab4) && (is.na(n_ab4) || n_ab4 < min_n || p_ab4 < min_prop)) {
  message(sprintf("AB4 is sparse for %s (n=%s, p=%.2f%%) -> merging AB4 into AB3.",
                  fam_name, ifelse(length(n_ab4)==0, 0, n_ab4), 100*ifelse(length(p_ab4)==0, 0, p_ab4)))
  ord_df <- ord_df %>%
    mutate(ord_class = fct_collapse(ord_class, `AB3 (100–999)` = c("AB3 (100–999)","AB4 (1000+)"))) %>%
    droplevels()  # keep only used levels after merge
}

# --------- Prepare response as integers 1..K expected by mgcv::ocat ---------
ord_df <- ord_df %>%
  mutate(
    y_ord = as.integer(factor(ord_class, levels = levels(ord_class))) # ordered categories -> 1..K
  )

K <- length(levels(ord_df$ord_class))  # number of ordinal categories after any merge

# --------- Choose k for the time smooth based on span of years ---------
k_time <- {
  ny <- dplyr::n_distinct(lubridate::year(ord_df$SAMPLE_DATE))
  min(20, max(8, round(0.6 * ny)))  # pragmatic cap/floor to avoid under/overfitting
}

# --------- Fit ordered categorical GAMM (season-varying trend + site RE) ---------
# Notes:
# - season_f enters as parametric shift
# - s(decimal_date, by = season_f) allows different smooths by season
# - s(SITE_ID.F, bs = "re") is a random intercept for site
# - family = ocat(R = K) sets ordered categorical likelihood with K thresholds
ocat_model <- bam(
  y_ord ~
    season_f +
    s(decimal_date, by = season_f, k = k_time) +
    s(SITE_ID.F, bs = "re"),
  family   = ocat(R = K),     # ordered categorical family (proportional odds link by default)
  data     = ord_df,
  method   = "fREML",
  discrete = TRUE,            # speed-up for large data
  select   = TRUE,            # shrink unneeded wiggle
  gamma    = 1.2              # mild extra penalty against overfit
)

# --------- Minimal reporting (text only, no figures) ---------
summary(ocat_model)              # parametric terms, EDFs, deviance explained
gam.check(ocat_model, rep = 0)   # basis checks & residual sanity (no re-smoothing)

# Optional: store for later use (predictions/plots/report)
# saveRDS(ocat_model, paste0("ocat_model_", fam_name, ".rds"))


# =========================================================
# Ordinal (ocat) Model Diagnostics 
# Requires: ocat_model, ord_df, fam_name
# Purpose: run lightweight checks + residual visuals + smooths + thresholds
# =========================================================

inspect_ocat <- function(model, data, fam = "Family") {
  # ---- 1) Basis/smooth adequacy (text only; suppress plot device) ----
  tmp <- tempfile(fileext = ".pdf")   # open a dummy device so gam.check doesn't draw to screen
  pdf(tmp)                            # divert any plots to a temp PDF
  gc_out <- try(gam.check(model, rep = 0), silent = TRUE)  # rep=0 avoids expensive re-smoothing
  dev.off(); unlink(tmp)              # close and delete the temp device/file
  if (!inherits(gc_out, "try-error")) print(gc_out)  # print the textual part if available

  # ---- 2) QQ plot of deviance residuals (single panel, base graphics) ----
  par(mfrow = c(1,1))                 # ensure a single plot frame
  mgcv::qq.gam(model, rep = 0, pch = 19, cex = 0.3,
               main = paste("QQ plot —", fam))  # quick normality check of residuals

  # ---- Prepare one data.frame with linear predictor & residuals (used below) ----
  dd <- data.frame(
    eta   = predict(model, type = "link"),         # linear predictor (latent scale)
    resid = residuals(model, type = "deviance")    # deviance residuals for ordered logit
  )

  # ---- 3) Residuals vs linear predictor (no smoother overlay) ----
  p_res_lin <- ggplot(dd, aes(eta, resid)) +
    geom_point(alpha = 0.15, size = 0.5) +         # light scatter; large n-safe
    geom_hline(yintercept = 0, linetype = 2) +     # zero-reference line
    labs(title = paste("Residuals vs Linear Predictor —", fam),
         x = "Linear predictor (link)", y = "Deviance residuals") +
    theme_minimal(base_size = 12)
  print(p_res_lin)

  # ---- 4) Histogram of deviance residuals (shape check) ----
  p_hist <- ggplot(dd, aes(resid)) +
    geom_histogram(bins = 40, fill = "#3D5A80", colour = "white", linewidth = 0.15) +
    labs(title = paste("Histogram of Residuals —", fam),
         x = "Deviance residuals", y = "Frequency") +
    theme_minimal(base_size = 12)
  print(p_hist)

  # ---- 5) Partial effects: season-specific time smooths (spring/autumn) ----
  # Identify smooth indices that involve 'decimal_date' (the temporal smooth)
  sm_ids   <- which(grepl("decimal_date", smooths(model)))
  sm_names <- smooths(model)[sm_ids]
  for (j in seq_along(sm_ids)) {
    # Friendly facet label from smooth name
    lab <- if (grepl("spring", sm_names[j], ignore.case = TRUE)) "spring"
           else if (grepl("autumn", sm_names[j], ignore.case = TRUE)) "autumn"
           else sm_names[j]
    # Draw each selected smooth with gratia::draw()
    p_sm <- draw(model, select = sm_ids[j]) +
      ggtitle(paste("Time smooth (", lab, ") — ", fam, sep = ""))
    print(p_sm)
  }

  # ---- 6) Thresholds (cutpoints a1…aK-1) on the latent scale with 95% CI ----
  cf  <- coef(model)                  # full coefficient vector
  V   <- vcov(model)                  # covariance matrix of coefficients
  idx <- grep("^a[0-9]+$", names(cf)) # threshold parameters are named a1, a2, ...
  if (length(idx)) {
    # Assemble estimates and Wald CIs for the ordered cutpoints
    th <- data.frame(threshold = names(cf)[idx],
                     est = cf[idx],
                     se  = sqrt(diag(V)[idx]))
    th$lower <- th$est - 1.96 * th$se
    th$upper <- th$est + 1.96 * th$se
    th$threshold <- factor(th$threshold, levels = th$threshold)

    # Plot thresholds with error bars — increasing cutpoints indicates ordered categories
    p_th <- ggplot(th, aes(threshold, est)) +
      geom_point() +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
      labs(title = paste("Category thresholds (latent scale) —", fam),
           x = "Cutpoint (a1 < a2 < …)", y = "Estimate ± 95% CI") +
      theme_minimal(base_size = 12)
    print(p_th)
  }

  invisible(NULL)  # return nothing (side-effect is printed output/plots)
}

# ---- Run diagnostics on your fitted ordered categorical model ----
inspect_ocat(ocat_model, ord_df, fam_name)


# ================================================================
# PHASE 4 — Batch Ordered-Categorical (OCAT) GAMMs for 4 families
#  - Uses categorical samples only (data_type == "bin")
#  - Optional S3PO-only filter
#  - Top-down merge of sparse top class (AB4 -> AB3) when very rare
#  - Prints summary() and gam.check() output 
#  - Returns compact results + a summary table
# ================================================================

# Families to run
families <- c("Aphelocheiridae","Brachycentridae","Odontoceridae","Cordulegastridae")
# ^ Vector of target macroinvertebrate families to process.

# Project options
USE_S3PO       <- TRUE   # EA-recommended
MIN_TOP_N      <- 50     # merge AB4 into AB3 if AB4 count < MIN_TOP_N ...
MIN_TOP_PROP   <- 0.02   # ... or if AB4 share < 2%
# ^ Controls for optional collapsing of a very sparse highest category.

# ---- Helper: build ordinal classes from a PA table ----
build_ordinal_df <- function(df_pa) {
  df_pa %>%
    mutate(
      # Map TOTAL_NUMBER into ordered bins using inclusive intervals:
      # AB0: exactly 0; AB1: 1–9; AB2: 10–99; AB3: 100–999; AB4: 1000+
      ord_class = case_when(
        TOTAL_NUMBER == 0L                         ~ "AB0 (0)",
        TOTAL_NUMBER > 0L   & TOTAL_NUMBER < 10L   ~ "AB1 (1–9)",
        TOTAL_NUMBER >= 10L & TOTAL_NUMBER < 100L  ~ "AB2 (10–99)",
        TOTAL_NUMBER >= 100L & TOTAL_NUMBER < 1000L~ "AB3 (100–999)",
        TOTAL_NUMBER >= 1000L                      ~ "AB4 (1000+)",
        TRUE ~ NA_character_
      ),
      # Fix the factor level order to reflect the natural ordinal scale
      ord_class = factor(
        ord_class,
        levels = c("AB0 (0)","AB1 (1–9)","AB2 (10–99)","AB3 (100–999)","AB4 (1000+)")
      ),
      # Ensure a site-level random-effect factor exists (create if missing)
      SITE_ID.F = if (!"SITE_ID.F" %in% names(.)) factor(SITE_ID) else SITE_ID.F
    ) %>%
    filter(!is.na(ord_class))  # drop any rows that failed classification
}

# ---- Helper: fit OCAT model for one family ----
fit_ocat_family <- function(fam_name,
                            md      = models_data,
                            use_s3po = TRUE,
                            min_top_n = MIN_TOP_N,
                            min_top_prop = MIN_TOP_PROP) {
  stopifnot(!is.null(md[[fam_name]]), "pa" %in% names(md[[fam_name]]))
  df <- md[[fam_name]]$pa
  # ^ Start from the Presence/Absence table for this family (already has zeros etc.)

  # Keep categorical samples only
  df <- dplyr::filter(df, data_type == "bin")
  # ^ The ordered-categorical model is defined on categorical observations.

  # Optional: S3PO-only
  if (use_s3po && "SAMPLE_METHOD" %in% names(df)) {
    df <- dplyr::filter(df, SAMPLE_METHOD == "S3PO")
  }
  # ^ Restrict to EA’s recommended sampling method when available.

  # Season guard (should already exist if you added the preprocessing patch)
  if (!("season_f" %in% names(df))) {
    df <- df %>%
      mutate(m = month(SAMPLE_DATE),
             season_f = factor(if_else(m %in% 3:5, "spring", "autumn"),
                               levels = c("spring","autumn")))
  }
  # ^ Ensure a two-level seasonal factor (spring/autumn) is present.

  # Map to ordinal classes
  ord_df <- build_ordinal_df(df)
  # ^ Convert raw TOTAL_NUMBER into ordered categories (AB0…AB4).

  # Merge AB4 -> AB3 if very sparse
  class_counts <- ord_df %>% count(ord_class, name = "n") %>%
    mutate(prop = n / sum(n))
  n_ab4 <- class_counts$n[class_counts$ord_class == "AB4 (1000+)"]
  p_ab4 <- class_counts$prop[class_counts$ord_class == "AB4 (1000+)"]
  if (length(n_ab4) && (is.na(n_ab4) || n_ab4 < min_top_n || p_ab4 < min_top_prop)) {
    ord_df <- ord_df %>%
      mutate(ord_class = fct_collapse(ord_class,
                        `AB3 (100–999)` = c("AB3 (100–999)","AB4 (1000+)"))) %>%
      droplevels()
  }
  # ^ Optional sparsity rule: collapse a rare top class to stabilize the fit.

  # Response 1..K for ocat()
  ord_df <- ord_df %>%
    mutate(y_ord = as.integer(factor(ord_class, levels = levels(ord_class))))
  K <- nlevels(ord_df$ord_class)
  # ^ mgcv::ocat expects integer categories 1..K; keep K for family().

  # Sensible k from span of years
  k_time <- {
    ny <- dplyr::n_distinct(lubridate::year(ord_df$SAMPLE_DATE))
    min(20, max(8, round(0.6 * ny)))
  }
  # ^ Time-smooth basis dimension chosen from temporal coverage (guardrails 8..20).

  # Fit the ordered-categorical GAMM
  m <- bam(
    y_ord ~
      season_f +
      s(decimal_date, by = season_f, k = k_time) +
      s(SITE_ID.F, bs = "re"),
    family   = ocat(R = K),
    data     = ord_df,
    method   = "fREML",
    discrete = TRUE,
    select   = TRUE,
    gamma    = 1.2
  )
  # ^ Model: season-specific temporal smooth + site random intercept on the latent scale.

  # Header + summary
  cat("\n============================\n")
  cat(sprintf("OCAT Family: %s | Classes: %s\n", fam_name,
              paste(levels(ord_df$ord_class), collapse = " | ")))
  cat("============================\n")
  print(summary(m))
  # ^ Print parameteric terms, EDFs, deviance explained, etc.

  # gam.check (plots suppressed)
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  gc_out <- try(gam.check(m, rep = 0), silent = TRUE)
  dev.off(); unlink(tmp)
  if (!inherits(gc_out, "try-error")) print(gc_out)
  # ^ Report k-index and residual checks without opening the 2x2 plot grid.

  # Metrics
  devex <- summary(m)$dev.expl
  aic   <- AIC(m)
  nobs_ <- nobs(m)
  # ^ Quick scalar diagnostics: deviance explained, AIC, sample size.

  list(
    family     = fam_name,          # family name
    model      = m,                 # fitted mgcv::bam object
    data       = ord_df,            # modelling data (ordinal-coded)
    classes    = levels(ord_df$ord_class),  # human-readable class labels
    K          = K,                 # number of ordinal levels
    k_time     = k_time,            # time-smooth basis size used
    dev_expl   = devex,             # deviance explained (fraction)
    AIC        = aic,               # AIC
    n          = nobs_              # number of observations
  )
}

# ---- Run for all families ----
ocat_results <- lapply(families, fit_ocat_family,
                       md = models_data, use_s3po = USE_S3PO)
names(ocat_results) <- families
# ^ Fit OCAT models family-by-family and name the results list for convenience.

# ---- Summary table across families ----
ocat_summary <- do.call(rbind, lapply(ocat_results, function(x) {
  data.frame(
    family   = x$family,                         # family name
    n        = x$n,                              # sample size used
    K        = x$K,                              # number of ordinal classes
    k_time   = x$k_time,                         # basis size for time smooth
    dev_expl = x$dev_expl,                       # deviance explained
    AIC      = x$AIC,                            # AIC
    classes  = paste(x$classes, collapse = " | "),  # label set used
    row.names = NULL
  )
}))
print(ocat_summary)
# ^ One tidy frame summarizing all fitted families.

# ---- Optional: save models ----
# invisible(lapply(ocat_results, function(x)


      
# =========================================================
# OCAT Model Diagnostics and Plots
# Expects: ocat_results (each item has $model, $family, $data)
# Purpose:
#   For each fitted ordered-categorical GAMM, produce a compact set
#   of diagnostic visuals:
#     1) QQ plot of residuals (base graphics)
#     2) Residuals vs linear predictor (no smoother)
#     3) Histogram of residuals
#     4) Season-specific partial-effect smooths over time
# =========================================================

plot_ocat_diagnostics_clean <- function(res) {
  m   <- res$model   # fitted mgcv::bam OCAT model
  fam <- res$family  # family name (for titles)
  df  <- res$data    # modelling data (not directly plotted here)

  # --- 1) QQ plot (single, base) ---
  #     Uses mgcv::qq.gam to compare residuals against theoretical quantiles.
  #     rep=0 avoids expensive re-smoothing; small points to de-clutter.
  par(mfrow = c(1,1))
  mgcv::qq.gam(m, rep = 0, pch = 19, cex = 0.3,
               main = paste("QQ plot —", fam))

  # Prepare diagnostics frame once
  # - eta: linear predictor on the latent (link) scale for OCAT
  # - resid: deviance residuals from the model
  eta   <- as.numeric(predict(m, type = "link"))
  resid <- as.numeric(residuals(m, type = "deviance"))
  dd    <- data.frame(eta = eta, resid = resid)

  # --- 2) Residuals vs linear predictor (no smooth) ---
  #     Checks for structure in residuals as a function of the linear predictor.
  #     Only a horizontal reference line at 0 is shown (no smoother).
  p_res_lin <- ggplot(dd, aes(eta, resid)) +
    geom_point(alpha = 0.15, size = 0.6) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(title = paste("Residuals vs Linear Predictor —", fam),
         x = "Linear predictor (latent scale)", y = "Deviance residuals") +
    theme_minimal(base_size = 12)
  print(p_res_lin)

  # --- 3) Histogram of residuals ---
  #     Quick look at residual distribution; large bins for a smooth view.
  p_hist <- ggplot(dd, aes(resid)) +
    geom_histogram(bins = 40, fill = "#3D5A80", colour = "white", linewidth = 0.15) +
    labs(title = paste("Histogram of Residuals —", fam),
         x = "Deviance residuals", y = "Frequency") +
    theme_minimal(base_size = 12)
  print(p_hist)

  # --- 4) Season-specific time smooths (spring/autumn), printed separately ---
  #     Identify the smooths that involve decimal_date (time), then draw them one by one.
  sm_ids   <- which(grepl("decimal_date", smooths(m)))
  sm_names <- smooths(m)[sm_ids]
  for (j in seq_along(sm_ids)) {
    nm  <- sm_names[j]
    # Give each panel a readable season label if present in the smooth name
    lab <- if (grepl("spring", nm, ignore.case = TRUE)) "spring"
           else if (grepl("autumn", nm, ignore.case = TRUE)) "autumn" else nm
    p_sm <- draw(m, select = sm_ids[j]) +
      ggtitle(paste("Time smooth (", lab, ") — ", fam, sep = ""))
    print(p_sm)
  }
}

# Run for all fitted OCAT models (prints plots for each family)
# - invisible() to avoid printing the lapply return value;
#   each plot is printed as a side effect within the function.
invisible(lapply(ocat_results, plot_ocat_diagnostics_clean))




#4.3 CENSORED POISSON GAMM MODEL
# ================================================
# CENSORED POISSON (cpois) — one family
# Goal:
#   Fit an interval-censored Poisson GAMM to numeric (count) data
#   using 1-significant-figure censoring bounds and site random effects.
# Inputs expected:
#   - models_data (from Phase 3), containing per-family tables
#   - For the chosen family: models_data[[fam_name]]$pa with columns:
#       SITE_ID, SAMPLE_ID, SAMPLE_DATE, TOTAL_NUMBER, data_type,
#       season_f (spring/autumn; if absent it is built below),
#       and SAMPLE_METHOD (to optionally filter to S3PO)
# Notes:
#   - Code *assumes* numeric samples are flagged as data_type == "count".
#   - Bounds follow 1-s.f. rounding guidance (non-integer limits).
# ================================================

# ---- pick a family ----
# Choose one of your focus families present in models_data.
fam_name <- "Aphelocheiridae"

# ---- helper: half-up rounding & 1-sf interval bounds ----
# round_half_up(): base "half up" rounding (5 → up) at given digit place.
round_half_up <- function(x, digits = 0) {
  s <- sign(x); x <- abs(x) * 10^digits
  s * floor(x + 0.5) / 10^digits
}

# bounds_1sf(): produce lower/upper censoring bounds implied by 1-s.f. rounding.
#   - For 0        → [-0.5, 0.5]
#   - For 1..9     → [x-0.5, x+0.5]
#   - For >= 10    → centre on the 1-s.f. rounded value with half-step width,
#                    minus 0.5 to align Poisson mass at integers (e.g. 20 → [14.5, 24.5]).
bounds_1sf <- function(v) {
  # vectorised container
  out <- matrix(NA_real_, nrow = length(v), ncol = 2,
                dimnames = list(NULL, c("lower","upper")))
  # exact zeros
  z <- v == 0L
  out[z,] <- cbind(-0.5, 0.5)
  # exact 1..9
  s <- v >= 1L & v <= 9L
  out[s,] <- cbind(v[s] - 0.5, v[s] + 0.5)
  # 10+
  b <- v >= 10L
  if (any(b)) {
    k     <- floor(log10(v[b]))                     # order of magnitude
    step  <- 10^k                                   # 1-s.f. step size (…, 1, 10, 100, …)
    s1sf  <- round_half_up(v[b], digits = -k)       # 1-s.f. rounded centre (5s round up)
    lower <- s1sf - step/2 - 0.5                    # lower non-integer bound
    upper <- s1sf + step/2 - 0.5                    # upper non-integer bound
    out[b,] <- cbind(lower, upper)
  }
  out
}

# ---- build modelling data from your Phase 3 objects ----
# Start from the presence/absence table but keep only numeric-count samples,
# optionally restrict to S3PO (EA recommendation), and ensure factors exist.
cp_df <- models_data[[fam_name]]$pa %>%
  filter(data_type == "count") %>%                 # only numeric-count samples
  filter(SAMPLE_METHOD == "S3PO") %>%              # EA-recommended (toggle here if needed)
  mutate(
    SITE_ID.F = if (!"SITE_ID.F" %in% names(.)) factor(SITE_ID) else SITE_ID.F,  # random-effect factor
    # safety: build season_f if absent (should already exist from preprocessing)
    season_f  = if (!"season_f" %in% names(.))
                  factor(if_else(month(SAMPLE_DATE) %in% 3:5, "spring", "autumn"),
                         levels = c("spring","autumn"))
                else season_f
  )

# Compute non-integer censoring bounds for each observation using the 1-s.f. rule.
B <- bounds_1sf(cp_df$TOTAL_NUMBER)
cp_df$lower <- B[, "lower"]
cp_df$upper <- B[, "upper"]

# Sensible k for the time smooth: scale with span of years, but capped.
k_time <- {
  ny <- dplyr::n_distinct(lubridate::year(cp_df$SAMPLE_DATE))
  min(60, max(35, round(0.6 * ny)))
}

# ---- fit censored Poisson GAMM ----
# Model: season intercepts + season-varying smooth of decimal time + site RE.
# Family cpois() handles interval-censored integer counts via (lower, upper).
cpois_model <- bam(
  cbind(lower, upper) ~
    season_f +
    s(decimal_date, by = season_f, k = k_time) +
    s(SITE_ID.F, bs = "re"),
  family   = cpois(),         # interval-censored Poisson likelihood
  data     = cp_df,
  method   = "fREML",         # fast REML estimation
  discrete = TRUE,            # large-data speed-up
  select   = TRUE,            # shrink unnecessary wiggle
  gamma    = 1.2              # mild extra penalty (guard against overfit)
)

# ---- minimal reporting ----
# Print model summary (parametric + smooth terms, deviance explained, etc.).
print(summary(cpois_model))

# Run gam.check with plots diverted to a temp PDF (so only text prints here),
# then delete the temp file. If gam.check succeeded, echo its (text) output.
tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
gc_out <- try(gam.check(cpois_model, rep = 0), silent = TRUE)
dev.off(); unlink(tmp)
if (!inherits(gc_out, "try-error")) print(gc_out)

# Optional quick metric line:
cat(sprintf("Deviance explained: %.1f%% | AIC: %.1f | N: %d\n",
            100*summary(cpois_model)$dev.expl, AIC(cpois_model), nobs(cpois_model)))

# Plot season-specific time smooths (spring/autumn) for quick visual inspection.
gratia::draw(cpois_model, select = which(grepl("decimal_date", smooths(cpois_model))))


      
# =========================================================
# CPOIS Model diagnostics — robust & compact
# Requires: cpois_model (bam fit), cp_df with lower/upper
# Purpose:
#   - Check basis adequacy (k-index), concurvity, and dispersion.
#   - Use randomized quantile residuals (RQR) tailored for *censored Poisson*.
#   - Provide quick visual checks: QQ, residual-vs-fitted, histogram, ACF.
#   - Report simple overdispersion proxies and headline fit metrics.
# Notes:
#   - This code assumes you already fitted `cpois_model` and built `cp_df`
#     containing numeric bounds `lower`/`upper` (non-integers OK).
# =========================================================

# --- Randomised quantile residuals (robust) ---
# For each observation i with integer bounds [li, ui], compute:
#   U_i ~ Uniform(F(li-1; μ_i), F(ui; μ_i)), where F is Poisson CDF with mean μ_i.
# Then map U_i through the standard normal inverse CDF to get N(0,1) under the model.
# * This properly accounts for censoring and discreteness.
cpois_rqr <- function(fit, lower, upper, eps = 1e-12, seed = 123) {
  n  <- length(lower)
  mu <- as.numeric(pmax(fitted(fit, type = "response"), eps))  # fitted λ (guard tiny)

  # integer bounds; enforce ui >= li and li >= 0
  li <- pmax(0L, ceiling(lower))   # move up to nearest integer
  ui <- floor(upper)               # move down to nearest integer
  ui <- pmax(ui, li)               # ensure non-empty interval

  Flo <- ppois(li - 1L, mu)        # F(lo-1) : mass strictly below li
  Fup <- ppois(ui,      mu)        # F(up)   : mass up to ui

  w <- pmax(Fup - Flo, 0)          # interval probability (non-negative)
  set.seed(seed)
  u <- Flo + w * runif(n)          # randomized PIT within the interval
  u <- pmin(pmax(u, eps), 1 - eps) # clamp away from 0/1 for stable qnorm

  qnorm(u)                         # RQR ~ N(0,1) if model well-specified
}

# --- Truncated Poisson conditional mean E[Y | li <= Y <= ui] ---
# Deterministic imputation used only for a Pearson-like dispersion proxy.
# If the interval probability is extremely small, fall back to a clamped round(λ).
tpois_mean <- function(li, ui, lambda, eps = 1e-12) {
  li <- as.integer(pmax(0L, li))
  ui <- as.integer(pmax(li, ui))           # ensure ui >= li
  # support
  xs <- li:ui
  # handle single-point interval quickly
  if (length(xs) == 1L) return(xs)
  # pmf on support
  pr <- dpois(xs, lambda)
  s  <- sum(pr)
  if (!is.finite(s) || s < eps) {
    # fallback: nearest integer to lambda within [li, ui]
    return(max(li, min(ui, round(lambda))))
  }
  sum(xs * pr) / s
}

diagnose_cpois <- function(fit, data) {
  # 1) Basis/smooth adequacy (suppress plots)
  #    - gam.check text includes k-index (should be ~1) and EDF vs k checks.
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  gc_out <- try(gam.check(fit, rep = 0), silent = TRUE)
  dev.off(); unlink(tmp)
  if (!inherits(gc_out, "try-error")) print(gc_out)

  # 2) Concurvity (safe)
  #    - High concurvity suggests redundant smooth structure (instability).
  con <- try(concurvity(fit, full = TRUE), silent = TRUE)
  if (!inherits(con, "try-error")) {
    cat("\nConcurvity (rounded):\n")
    if (is.list(con) && !is.null(con$estimate)) print(round(con$estimate, 3)) else print(con)
  }

  # 3) Randomised quantile residuals
  #    - Correct for censoring; should be ~N(0,1) if the model is well specified.
  rqr <- cpois_rqr(fit, data$lower, data$upper)

  # QQ plot for RQR
  par(mfrow = c(1,1))
  qqnorm(rqr, main = "QQ plot — RQ residuals"); qqline(rqr)

  # Residuals vs fitted (no smoother) + histogram
  mu <- fitted(fit, type = "response")
  dd <- tibble(mu = mu, rqr = rqr)

  print(
    ggplot(dd, aes(mu, rqr)) +
      geom_point(alpha = 0.15, size = 0.6) +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(title = "RQR vs Fitted mean (λ)", x = "Fitted λ", y = "RQ residual") +
      theme_minimal(base_size = 12)
  )
  print(
    ggplot(dd, aes(rqr)) +
      geom_histogram(bins = 40, fill = "#3D5A80", colour = "white", linewidth = 0.15) +
      labs(title = "Histogram of RQ residuals", x = "RQ residual", y = "Frequency") +
      theme_minimal(base_size = 12)
  )
  # Autocorrelation: substantial spikes may indicate temporal dependence not captured.
  acf(rqr[is.finite(rqr)], na.action = na.pass, main = "ACF of RQ residuals")

  # 4) Overdispersion proxies (no sampling)
  #    - SD of RQR: should be ~1 under correct variance.
  #    - Pearson-like: use truncated-Poisson E[Y|interval] as "observed" y.
  sd_rqr <- sd(rqr, na.rm = TRUE)

  li <- pmax(0L, ceiling(data$lower))
  ui <- pmax(li, floor(data$upper))
  y_bar <- mapply(tpois_mean, li, ui, mu)   # deterministic imputation on integer support
  pearson <- sum((y_bar - mu)^2 / pmax(mu, 1e-8))   # guard by tiny μ
  disp    <- pearson / fit$df.residual               # dispersion ≈ 1 is ideal

  cat(sprintf("\nOverdispersion — SD(RQR): %.2f (≈1 ideal) | Pearson proxy: %.2f\n",
              sd_rqr, disp))
}

# ---- Run on your model/data ----
# Ensure `cpois_model` and `cp_df` exist in the workspace before calling.
diagnose_cpois(cpois_model, cp_df)




#4.4 CENSORED NORMAL GAMM MODEL
# ================================================
# CENSORED NORMAL (cnorm) — one family + partial effects
# Goal:
#   Fit an interval-censored Normal GAMM to sqrt-transformed
#   count intervals for a single family, with:
#     • season main effect (spring vs autumn)
#     • a single smooth temporal trend
#     • site-level random intercepts
#   Then print a text summary and the partial time smooth.
# Assumes:
#   - `models_data[[fam_name]]$cnorm` exists from Phase 3 and
#     already contains count-scale bounds `lower`/`upper`,
#     `season_f`, `decimal_date`, `SITE_ID.F`, `SAMPLE_DATE`.
# ================================================

# ---- choose family ----
# Change the string to any of your four focus families as needed.
fam_name <- "Aphelocheiridae"   # change as needed

# ---- get modelling data built earlier (same as cpois) ----
# Pull the cnorm-ready table for the chosen family.
cn_df <- models_data[[fam_name]]$cnorm

# ---- build sqrt-intervals (avoid sqrt of negatives) ----
# Convert the count-scale interval bounds to sqrt(count) scale,
# clamping at 0 first to prevent sqrt of negative values.
to_sqrt_bounds <- function(l, u) {
  l2 <- pmax(as.numeric(l), 0)
  u2 <- pmax(as.numeric(u), 0)
  cbind(lower_t = sqrt(l2), upper_t = sqrt(u2))
}
SB <- to_sqrt_bounds(cn_df$lower, cn_df$upper)
cn_df$lower_t <- SB[, "lower_t"]
cn_df$upper_t <- SB[, "upper_t"]

# ---- sensible k from span of years ----
# Choose basis dimension k as a function of the number of sampled years:
# not too small (>=10) and capped at 30 to avoid overfitting.
k_time <- {
  ny <- dplyr::n_distinct(lubridate::year(cn_df$SAMPLE_DATE))
  min(30, max(10, round(0.6 * ny)))
}

# ---- fit cnorm: season main effect + single time smooth ----
# Interval response: cbind(lower_t, upper_t) on sqrt scale.
#   • season_f as a parametric factor (spring/autumn intercept shift)
#   • s(decimal_date) for a single smooth trend (same shape across seasons)
#   • s(SITE_ID.F, "re") as site random intercepts
# Family: cnorm() — censored Normal with identity link on sqrt scale.
# Fitting: bam(..., method="fREML", discrete=TRUE, select=TRUE, gamma=1.2)
cnorm_model <- bam(
  cbind(lower_t, upper_t) ~
    season_f +
    s(decimal_date, k = k_time) +
    s(SITE_ID.F, bs = "re"),
  family   = cnorm(),
  data     = cn_df,
  method   = "fREML",
  discrete = TRUE,
  select   = TRUE,
  gamma    = 1.2
)

# ---- summary ----
# Prints parametric terms, smooth EDFs, deviance explained, etc.
summary(cnorm_model)

# ---- partial-effect plots ----
# (1) time smooth: shows the marginal temporal trend on sqrt scale.
draw(cnorm_model, select = 1) +
  ggplot2::ggtitle(paste("Time smooth —", fam_name, "(cnorm on sqrt count)"))

# (2) optional: show site RE distribution (not a smooth, but useful)
# If you want to plot site random effects, uncomment the next line
# and ensure the object name matches (`cnorm_model2` here).
# gratia::draw(cnorm_model, select = 2)  # uncomment if you want the RE plot


      
# =========================================================
# CNORM Model Diagnostics 
# Requires: cnorm_model (bam fit), cn_df with lower_t/upper_t
# Purpose:
#   - Run lightweight but meaningful checks for an interval-censored
#     Normal GAMM fitted on sqrt(count) scale (mgcv::cnorm()).
#   - Provide: basis adequacy, concurvity, randomized-quantile residuals,
#     simple residual panels (QQ / vs-fitted / hist / ACF), dispersion cue,
#     and headline fit metrics. Optionally shows time smooth partial effects.
# Notes:
#   - All plots are cheap and reproducible.
#   - Expects cn_df to already have lower_t/upper_t (sqrt bounds). If not,
#     they are created from count-scale lower/upper via ensure_sqrt_bounds().
# =========================================================

# --- Ensure sqrt-bounds exist (build from count-scale bounds if needed)
# Input:  df with either {lower_t, upper_t} OR {lower, upper} on count scale
# Output: df guaranteed to include lower_t/upper_t on sqrt(count) scale
ensure_sqrt_bounds <- function(df) {
  if (!all(c("lower_t","upper_t") %in% names(df))) {
    stopifnot(all(c("lower","upper") %in% names(df)))
    df <- df %>%
      mutate(
        lower_t = sqrt(pmax(as.numeric(lower), 0)),  # clamp at 0 before sqrt
        upper_t = sqrt(pmax(as.numeric(upper), 0))
      )
  }
  df
}

# --- Randomised quantile residuals for censored Normal ---
# Constructs PIT values from the Normal CDF between [lt, ut] and randomises
# uniformly within the censoring interval; then maps to N(0,1) via qnorm.
# Arguments:
#   fit : mgcv::bam fit with family=cnorm()
#   lt,ut : sqrt-scale lower/upper bounds
#   eps : numerical guard against exact 0/1 PIT values
#   seed : RNG seed for reproducibility of the randomisation step
cnorm_rqr <- function(fit, lt, ut, eps = 1e-10, seed = 123) {
  mu <- as.numeric(fitted(fit, type = "response"))    # on sqrt scale
  sd <- sqrt(summary(fit)$scale)                      # sigma on sqrt scale

  Flo <- pnorm(lt, mean = mu, sd = sd)                # F(lower)
  Fup <- pnorm(ut, mean = mu, sd = sd)                # F(upper)

  w <- pmax(Fup - Flo, 0)                             # interval mass
  set.seed(seed)
  u <- Flo + w * runif(length(mu))                    # randomised PIT in [Flo,Fup]
  u <- pmin(pmax(u, eps), 1 - eps)                    # clamp away from 0/1
  qnorm(u)                                            # -> approx N(0,1) if model is correct
}

diagnose_cnorm <- function(fit, data) {
  # Standardise inputs: guarantee sqrt-scale bounds
  data <- ensure_sqrt_bounds(data)

  # 1) Basis/smooth adequacy (suppress 2×2 plot)
  # Uses gam.check() with rep=0 (no expensive re-smoothing); plots are diverted.
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  gc_out <- try(gam.check(fit, rep = 0), silent = TRUE)
  dev.off(); unlink(tmp)
  if (!inherits(gc_out, "try-error")) print(gc_out)

  # 2) Concurvity (safe)
  # Detects non-linear dependence among smooth terms; prints rounded estimates.
  con <- try(concurvity(fit, full = TRUE), silent = TRUE)
  if (!inherits(con, "try-error")) {
    cat("\nConcurvity (rounded):\n")
    if (is.list(con) && !is.null(con$estimate)) print(round(con$estimate, 3)) else print(con)
  }

  # 3) Randomised quantile residuals (Normal PIT on sqrt scale)
  # RQRs should be ~N(0,1) with minimal structure if model is appropriate.
  rqr <- cnorm_rqr(fit, data$lower_t, data$upper_t)

  # QQ plot (RQR) — quick Gaussianity check
  par(mfrow = c(1,1))
  qqnorm(rqr, main = "QQ plot — RQ residuals (cnorm)"); qqline(rqr)

  # Residuals vs fitted & histogram — check for mean-variance issues/outliers
  mu <- fitted(fit, type = "response")
  dd <- tibble(mu = mu, rqr = rqr)

  print(
    ggplot(dd, aes(mu, rqr)) +
      geom_point(alpha = 0.15, size = 0.6) +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(title = "RQR vs Fitted mean (sqrt scale)",
           x = "Fitted mean (sqrt count)", y = "Randomised quantile residual") +
      theme_minimal(base_size = 12)
  )

  print(
    ggplot(dd, aes(rqr)) +
      geom_histogram(bins = 40, fill = "#3D5A80", colour = "white", linewidth = 0.15) +
      labs(title = "Histogram of RQ residuals (cnorm)",
           x = "RQ residual", y = "Frequency") +
      theme_minimal(base_size = 12)
  )

  # ACF of RQR — temporal autocorrelation check (should be small if independent)
  acf(rqr[is.finite(rqr)], na.action = na.pass, main = "ACF of RQ residuals (cnorm)")

  # 4) Dispersion cue: SD(RQR) should be ~1 under correct specification
  sd_rqr <- sd(rqr, na.rm = TRUE)
  cat(sprintf("\nResidual spread — SD(RQR): %.2f (≈1 ideal)\n", sd_rqr))

  # 5) Time smooth (partial effect) — printed by season structure
  # Draw all smooths containing 'decimal_date' (single or by-season).
  sm_ids <- which(grepl("decimal_date", smooths(fit)))
  if (length(sm_ids)) print(draw(fit, select = sm_ids))

  # 6) Headline metrics — quick readout of fit quality
  cat(sprintf("Deviance explained: %.1f%% | AIC: %.1f | N: %d | sigma(sqrt-scale): %.3f\n",
              100 * summary(fit)$dev.expl, AIC(fit), nobs(fit), sqrt(summary(fit)$scale)))
}

# ---- Run on your cnorm model/data ----
# Expects:
#   - cnorm_model : mgcv::bam fit with family = cnorm()
#   - cn_df       : data.frame used to fit (has lower_t/upper_t or lower/upper)
diagnose_cnorm(cnorm_model, cn_df)


      

#4.4.1 CNORM MODEL(INTERVAL WEIGHTING + POWER TUNING with REPORTING_AREA COVARIATE)
# =========================================================
# Make cnorm tighter: interval weighting + power tuning
# =========================================================

# --- helper: RQ residuals for cnorm (same as before)
#     Given a fitted censored-normal model on the transformed scale, compute
#     randomized quantile residuals by sampling uniformly within each case’s
#     CDF interval [F(lower), F(upper)]. Ideal dispersion: SD ≈ 1.
cnorm_rqr <- function(fit, lt, ut, eps = 1e-10, seed = 123) {
  mu <- as.numeric(fitted(fit, type = "response"))  # fitted mean on transformed (e.g., sqrt) scale
  sd <- sqrt(summary(fit)$scale)                    # model scale (σ) on the same scale as mu
  Flo <- pnorm(lt, mean = mu, sd = sd)              # CDF at lower bound
  Fup <- pnorm(ut, mean = mu, sd = sd)              # CDF at upper bound
  w   <- pmax(Fup - Flo, 0)                         # interval probability (non-negative)
  set.seed(seed)
  u <- Flo + w * runif(length(mu))                  # uniform draw within [Flo, Fup]
  u <- pmin(pmax(u, eps), 1 - eps)                  # clamp away from 0/1 to avoid ±Inf
  qnorm(u)                                          # convert PIT to Normal residuals
}

# --- power-transform the bounds *after* you’ve computed correct 1-s.f. count-scale bounds
#     Input: df with count-scale 'lower'/'upper' (non-integer, 1-s.f. rounding).
#     Output: matching bounds on transformed scale: sqrt if p=0.5, or x^p otherwise.
transform_bounds_p <- function(df, p = 0.5) {
  stopifnot(all(c("lower","upper") %in% names(df)))
  lt <- pmax(as.numeric(df$lower), 0)               # guard against negative lower
  ut <- pmax(as.numeric(df$upper), 0)               # guard against negative upper
  if (p == 0.5) {
    data.frame(lower_t = sqrt(lt), upper_t = sqrt(ut))
  } else {
    data.frame(lower_t = lt^p,     upper_t = ut^p)
  }
}

# --- one fit with given p and weighting exponent alpha_w
#     - Transforms bounds with power p
#     - Weights each row by interval width^{-alpha_w} so narrower intervals get more weight
#     - Optionally adds an EA reporting-area random effect if 'sites_clean' exists
fit_cnorm_power <- function(df, p = 0.5, alpha_w = 1, add_area_re = TRUE) {
  B  <- transform_bounds_p(df, p)                    # transform bounds to 'lower_t'/'upper_t'
  df$lower_t <- B$lower_t; df$upper_t <- B$upper_t

  # weights: narrower intervals get higher weight
  width_t <- pmax(df$upper_t - df$lower_t, 1e-6)    # avoid division by ~0
  wts     <- (1 / width_t)^alpha_w                   # alpha_w=1 gives strong emphasis

  # add EA region RE if you have sites_clean
  # (joins REPORTING_AREA by SITE_ID; only if not already present and object exists)
  if (add_area_re && !"REPORTING_AREA" %in% names(df) && exists("sites_clean")) {
    df <- dplyr::left_join(
      df,
      dplyr::select(sites_clean, SITE_ID, REPORTING_AREA),  # qualify select to avoid conflicts
      by = "SITE_ID"
    )
  }
  if ("REPORTING_AREA" %in% names(df)) df$REPORTING_AREA <- factor(df$REPORTING_AREA)

  # sensible k for the time smooth: scales with number of unique years (cap 35)
  k_time <- {
    ny <- dplyr::n_distinct(lubridate::year(df$SAMPLE_DATE))
    min(35, max(12, round(0.7 * ny)))
  }

  # Model formula:
  # - season_f fixed effect (spring/autumn intercept shift)
  # - single time smooth (no by=season to keep parsimony here)
  # - site random intercept; optional reporting-area RE if available
  fml <- if ("REPORTING_AREA" %in% names(df)) {
    cbind(lower_t, upper_t) ~ season_f + s(decimal_date, k = k_time) +
      s(SITE_ID.F, bs = "re") + s(REPORTING_AREA, bs = "re")
  } else {
    cbind(lower_t, upper_t) ~ season_f + s(decimal_date, k = k_time) +
      s(SITE_ID.F, bs = "re")
  }

  # Fit censored normal GAM (bam for speed):
  # - weights emphasize precise intervals
  # - select=TRUE + gamma=1.2 provide gentle extra regularization
  m <- bam(
    fml,
    family   = cnorm(),
    data     = df,
    weights  = wts,
    method   = "fREML",
    discrete = TRUE,
    select   = TRUE,
    gamma    = 1.2
  )

  # dispersion cue via RQ residuals (should be ~1 if well calibrated)
  rqr <- cnorm_rqr(m, df$lower_t, df$upper_t)
  sd_rqr <- sd(rqr, na.rm = TRUE)

  list(model = m, data = df, p = p, alpha_w = alpha_w, sd_rqr = sd_rqr)  # return fit + diagnostics
}

# ---- run a small grid over p and pick what tightens RQ spread the most
cn_df <- models_data[["Aphelocheiridae"]]$cnorm  # or your chosen family
# IMPORTANT: ensure cn_df$lower/upper are the *correct 1-s.f.* count-scale bounds

grid_p   <- c(0.35, 0.4, 0.45, 0.5, 0.6)   # √ is 0.5; lower p compresses high counts more
alpha_w  <- 1                               # try 0.5 and 1 to vary interval weighting strength

# Fit across the grid of p, holding alpha_w fixed; collect SD of RQ residuals
fits <- lapply(grid_p, function(pp) fit_cnorm_power(cn_df, p = pp, alpha_w = alpha_w))
sd_tab <- data.frame(
  p = sapply(fits, `[[`, "p"),
  SD_RQR = sapply(fits, `[[`, "sd_rqr")
)
print(sd_tab)  # inspect which p gives SD_RQR closest to 1 and/or smallest

# choose the best (smallest SD_RQR)
best_idx <- which.min(sd_tab$SD_RQR)
best_fit <- fits[[best_idx]]
cat(sprintf("\nChosen p = %.2f | SD(RQR) = %.2f (alpha_w = %.2f)\n",
            best_fit$p, best_fit$sd_rqr, best_fit$alpha_w))

# inspect k-index and smooth quickly (plots suppressed by sending to a temp PDF)
tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
print(try(gam.check(best_fit$model, rep = 0), silent = TRUE))
dev.off(); unlink(tmp)

# optional: draw the time smooth (partial effect) for the chosen p
draw(best_fit$model, select = 1) + ggtitle(sprintf("Time smooth — p=%.2f (weighted cnorm)", best_fit$p))


               

#4.4.2 CNORM MODEL( P=O.35 WITH NO REPORTING_AREA COVARIATE)
# =========================================================
# CENSORED-NORMAL with p = 0.35 (no REPORTING_AREA term)
# Keeps everything else the same: 1-s.f. intervals, weights,
# season main effect, single time smooth, site RE.
# ---------------------------------------------------------
# Purpose of this script:
#   • Fit a censored Normal GAM on power-transformed bounds (p = 0.35).
#   • Use interval-width weights (narrower intervals = higher weight).
#   • Include fixed season effect, one smooth trend over decimal_date,
#     and a random intercept for SITE_ID.F.
# Inputs expected (already prepared upstream):
#   • models_data[[fam_name]]$cnorm with columns:
#       lower/upper (1-s.f. count-scale bounds), SAMPLE_DATE, decimal_date,
#       season_f (or SAMPLE_DATE to derive it), SITE_ID or SITE_ID.F, TOTAL_NUMBER (optional).
# Notes:
#   • This version intentionally omits a REPORTING_AREA random effect.
#   • Only prints model summary and text-only gam.check output (no plots saved).
# =========================================================

# --- helper: RQ residuals for cnorm ---
# Computes randomized quantile residuals for the censored Normal fit:
#   1) Get fitted mean (mu) and residual scale (sd) on the transformed scale.
#   2) Compute Normal CDF at lower/upper bounds.
#   3) Draw a uniform random value within [F(lower), F(upper)] to respect censoring.
#   4) Map back to z-scale via qnorm for approximately N(0,1) residuals.
cnorm_rqr <- function(fit, lt, ut, eps = 1e-10, seed = 123) {
  mu <- as.numeric(fitted(fit, type = "response"))
  sd <- sqrt(summary(fit)$scale)
  Flo <- pnorm(lt, mean = mu, sd = sd)
  Fup <- pnorm(ut, mean = mu, sd = sd)
  w   <- pmax(Fup - Flo, 0)
  set.seed(seed)
  u <- Flo + w * runif(length(mu))
  u <- pmin(pmax(u, eps), 1 - eps)
  qnorm(u)
}

# --- power-transform bounds AFTER you’ve built correct 1-s.f. count-scale bounds
# Transforms raw count-scale bounds to the p-power scale used by the cnorm family.
# Here p = 0.35 (from your tuning). Assumes non-negative bounds.
transform_bounds_p <- function(df, p = 0.35) {
  stopifnot(all(c("lower","upper") %in% names(df)))
  lt <- pmax(as.numeric(df$lower), 0)
  ut <- pmax(as.numeric(df$upper), 0)
  data.frame(lower_t = lt^p, upper_t = ut^p)
}

# --- fit cnorm with p=0.35; NO REPORTING_AREA term ---
# Wraps the entire fitting pipeline:
#   • Builds transformed bounds lower_t/upper_t with p = 0.35.
#   • Computes precision weights = 1 / (interval width)^alpha_w.
#   • Chooses k for the time smooth based on year span.
#   • Fits bam() with family = cnorm().
# Returns a small list with model, data used, SD of RQ residuals, and k_time.
fit_cnorm_p035 <- function(df, alpha_w = 1) {
  B  <- transform_bounds_p(df, p = 0.35)
  df$lower_t <- B$lower_t
  df$upper_t <- B$upper_t

  # precision weights: narrower transformed intervals -> higher weight
  width_t <- pmax(df$upper_t - df$lower_t, 1e-6)
  wts     <- (1 / width_t)^alpha_w

  # sensible k for time smooth
  k_time <- {
    ny <- dplyr::n_distinct(lubridate::year(df$SAMPLE_DATE))
    min(35, max(12, round(0.7 * ny)))
  }

  m <- bam(
    cbind(lower_t, upper_t) ~
      season_f +
      s(decimal_date, k = k_time) +
      s(SITE_ID.F, bs = "re"),
    family   = cnorm(),
    data     = df,
    weights  = wts,
    method   = "fREML",
    discrete = TRUE,
    select   = TRUE,
    gamma    = 1.2
  )

  # quick dispersion cue
  rqr <- cnorm_rqr(m, df$lower_t, df$upper_t)
  sd_rqr <- sd(rqr, na.rm = TRUE)

  list(model = m, data = df, sd_rqr = sd_rqr, k_time = k_time)
}

# ----------------- RUN -----------------
# Choose the family to fit:
#   • Ensure models_data[[fam_name]]$cnorm exists and has correct lower/upper.
fam_name <- "Aphelocheiridae"  # change as needed
cn_df    <- models_data[[fam_name]]$cnorm   # assumes lower/upper are correct 1-s.f. bounds

# Fit the weighted cnorm with p = 0.35 (no area RE). Keep alpha_w = 1 unless you have reason to tweak.
res_p035 <- fit_cnorm_p035(cn_df, alpha_w = 1)

# summaries (no extra plots)
#   • summary(): parametric terms, smooth EDFs, deviance explained, scale, etc.
print(summary(res_p035$model))

#   • gam.check(): text output only (plots diverted to and discarded from a temp PDF)
tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
print(try(gam.check(res_p035$model, rep = 0), silent = TRUE))
dev.off(); unlink(tmp)

#   • Compact headline metrics for your log
cat(sprintf("SD(RQ residuals): %.2f | Dev.expl: %.1f%% | AIC: %.1f | N: %d\n",
            res_p035$sd_rqr, 100 * summary(res_p035$model)$dev.expl,
            AIC(res_p035$model), nobs(res_p035$model)))

# optional: partial effect of time
#   • Visualises the single time smooth on the transformed scale (no REs shown).
draw(res_p035$model, select = 1) +
  ggtitle(paste0("Time smooth — ", fam_name, " (cnorm, p=0.35, weighted)"))




# =========================================================
# Model Diagnostics for the p = 0.35 cnorm model (no REPORTING_AREA)
# Works with objects returned by `res_p035 <- fit_cnorm_p035(...)`
# ---------------------------------------------------------
# What this script does:
#   • Accepts a fitted censored-Normal GAM (on a p-power scale) and the data used.
#   • Rebuilds transformed bounds if missing, using the same power p (default 0.35).
#   • Runs text-only k-adequacy and concurvity checks.
#   • Computes randomized quantile residuals (RQR) for censored Normal and
#     visualizes them (QQ, residuals vs fitted, histogram, ACF).
#   • Prints headline fit metrics (DevExpl, AIC, N, sigma on p-scale).
# Assumptions:
#   • `data` has raw count-scale bounds `lower`/`upper`; or already has `lower_t`/`upper_t`.
#   • The model `fit` is a `bam()` fit with `family = cnorm()`.
#   • The same p used at fit time (0.35 here) is passed to diagnostics for consistency.
# =========================================================


# --- helper: build transformed bounds if needed (p-scale) ---
# Takes raw count-scale bounds (lower/upper) and applies the p-power transform,
# returning a data.frame with columns `lower_t` and `upper_t` (both on the p-scale).
transform_bounds_p <- function(df, p = 0.35) {
  stopifnot(all(c("lower","upper") %in% names(df)))
  lt <- pmax(as.numeric(df$lower), 0)
  ut <- pmax(as.numeric(df$upper), 0)
  data.frame(lower_t = lt^p, upper_t = ut^p)
}

# --- Randomised quantile residuals for censored Normal (generic p-scale) ---
# For each observation, draws a uniform value between F(lower_t) and F(upper_t)
# under N(mu, sd^2) on the p-scale, then converts to a standard normal score.
# This produces residuals that should be ~N(0,1) if the model/distribution is correct.
cnorm_rqr <- function(fit, lt, ut, eps = 1e-10, seed = 123) {
  mu <- as.numeric(fitted(fit, type = "response"))
  sd <- sqrt(summary(fit)$scale)
  Flo <- pnorm(lt, mean = mu, sd = sd)
  Fup <- pnorm(ut, mean = mu, sd = sd)
  w   <- pmax(Fup - Flo, 0)
  set.seed(seed)
  u <- Flo + w * runif(length(mu))
  u <- pmin(pmax(u, eps), 1 - eps)
  qnorm(u)
}

diagnose_cnorm_p <- function(fit, data, p = 0.35) {
  # Ensure transformed bounds exist
  if (!all(c("lower_t","upper_t") %in% names(data))) {
    B <- transform_bounds_p(data, p = p)
    data$lower_t <- B$lower_t; data$upper_t <- B$upper_t
  }

  # 1) k-adequacy (no 2x2 plot)
  # Diverts any plots into a temp device to avoid clutter; prints text output if available.
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  print(try(gam.check(fit, rep = 0), silent = TRUE))
  dev.off(); unlink(tmp)

  # 2) Concurvity
  # Safe wrapper: if concurvity computation fails, we just skip it.
  con <- try(concurvity(fit, full = TRUE), silent = TRUE)
  if (!inherits(con, "try-error")) {
    cat("\nConcurvity (rounded):\n")
    if (is.list(con) && !is.null(con$estimate)) print(round(con$estimate, 3)) else print(con)
  }

  # 3) Randomised quantile residuals
  # Build RQR on the p-scale bounds then visualize basic diagnostics.
  rqr <- cnorm_rqr(fit, data$lower_t, data$upper_t)

  par(mfrow = c(1,1))
  qqnorm(rqr, main = "QQ plot — RQ residuals (cnorm, p=0.35)"); qqline(rqr)

  mu <- fitted(fit, type = "response")
  dd <- tibble(mu = mu, rqr = rqr)

  print(
    ggplot(dd, aes(mu, rqr)) +
      geom_point(alpha = 0.15, size = 0.6) +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(title = "RQR vs Fitted mean (p-scale)",
           x = "Fitted mean on p-scale", y = "Randomised quantile residual") +
      theme_minimal(base_size = 12)
  )

  print(
    ggplot(dd, aes(rqr)) +
      geom_histogram(bins = 40, fill = "#3D5A80", colour = "white", linewidth = 0.15) +
      labs(title = "Histogram of RQ residuals (cnorm, p=0.35)",
           x = "RQ residual", y = "Frequency") +
      theme_minimal(base_size = 12)
  )

  acf(rqr[is.finite(rqr)], na.action = na.pass,
      main = "ACF of RQ residuals (cnorm, p=0.35)")

  # 4) Dispersion cue
  # SD of RQR should be close to 1. Values <<1 suggest underdispersion;
  # values >>1 suggest overdispersion or model misspecification.
  sd_rqr <- sd(rqr, na.rm = TRUE)
  cat(sprintf("\nResidual spread — SD(RQR): %.2f (≈1 ideal)\n", sd_rqr))

  # 5) Time smooth
  # Prints the partial effect(s) for the smooth(s) involving decimal_date (if present).
  sm_ids <- which(grepl("decimal_date", smooths(fit)))
  if (length(sm_ids)) print(draw(fit, select = sm_ids))

  # 6) Headline metrics
  # Quick, text-only summary to add to logs.
  cat(sprintf("Deviance explained: %.1f%% | AIC: %.1f | N: %d | sigma(p-scale): %.3f\n",
              100 * summary(fit)$dev.expl, AIC(fit), nobs(fit), sqrt(summary(fit)$scale)))
}

# ---- Run on your fitted model/data (from previous step) ----
# `res_p035` is expected to be created by fit_cnorm_p035(); we reuse the same p (=0.35).
diagnose_cnorm_p(res_p035$model, res_p035$data, p = 0.35)



               
# ================================
# Minor checks for p = 0.35 cnorm
# --------------------------------
# Purpose of this script:
#   • Recreate the per-row precision weights used in the fitted cnorm model.
#   • Check whether those weights are stable over time (summary + quick plot).
#   • Run a leave-one-site-out (for top 10 sites) influence check on the
#     population-level time curve (with site RE excluded in prediction).
#   • Compare the baseline single time smooth vs a season-varying time smooth.
# Notes:
#   • Assumes `res_p035` exists and was returned by `fit_cnorm_p035(...)`.
#   • Assumes `models_data` preprocessing created season_f, decimal_date, SITE_ID.F.
#   • All modelling settings are kept exactly as in your fit; no code is changed.
# ================================

# Extract fitted model and data from the p=0.35 cnorm result object
m  <- res_p035$model
df <- res_p035$data

# Recreate the analysis weights used in the fit
# - width_t: width of the transformed (p-scale) interval = upper_t - lower_t
# - wt: precision weight = 1 / width_t (same exponent as fit: alpha_w = 1)
# - year: helper for time summaries
df <- df %>%
  mutate(width_t = pmax(upper_t - lower_t, 1e-6),
         wt      = 1 / width_t,
         year    = year(SAMPLE_DATE))

# -------------------------------------------------
# 1) Weight stability through time (quick summary)
# -------------------------------------------------
# For each year: number of rows, median interval width, and 10th/median/90th
# percentiles of the precision weights. Large temporal drift could indicate
# a change in interval construction or data quality through time.
w_summary <- df %>%
  group_by(year) %>%
  summarise(n = n(),
            med_width = median(width_t),
            p10_wt    = quantile(wt, 0.10),
            med_wt    = median(wt),
            p90_wt    = quantile(wt, 0.90),
            .groups = "drop")
print(w_summary, n = 10)

# Optional tiny diagnostic plot (comment out if not needed)
# Median weight per year with an uncertainty ribbon (10–90% quantiles).
ggplot(w_summary, aes(year, med_wt)) +
  geom_ribbon(aes(ymin = p10_wt, ymax = p90_wt), alpha = 0.15) +
  geom_line() +
  labs(title = "Precision weight stability over time",
       x = "Year", y = "Median weight (ribbon = 10–90% range)") +
  theme_minimal(base_size = 12)

# Spearman correlation between weight and time (should be small)
# If |rho| is large, weights systematically shift over years, which you may want to note.
cor_w_time <- with(df, cor(wt, year, method = "spearman"))
cat(sprintf("\nSpearman(weight, year) = %.3f (|rho| small is good)\n", cor_w_time))

# -------------------------------------------------
# 2) Leave-one-site-out (top 10 sites) influence
#     on population time curve (site RE excluded)
# -------------------------------------------------
# Build a population-level prediction grid across time:
#  - Fixed season_f = "spring" for a single representative trend;
#  - SITE_ID.F is set to any valid level (random effect is excluded in predict()).
tgrid <- data.frame(
  decimal_date = seq(min(df$decimal_date), max(df$decimal_date), length.out = 200),
  season_f     = factor("spring", levels = levels(df$season_f)),
  SITE_ID.F    = df$SITE_ID.F[1]   # placeholder (RE excluded anyway)
)
# Baseline population curve from the full model (exclude site RE)
pred0 <- as.numeric(predict(m, newdata = tgrid, type = "response",
                            exclude = "s(SITE_ID.F)"))

# Identify the 10 sites with the largest number of observations
top_sites <- df %>% count(SITE_ID.F, sort = TRUE) %>% slice_head(n = 10) %>% pull(SITE_ID.F)

# For each top site:
#   • Refit the model dropping that site's data.
#   • Predict the population curve on the same tgrid (exclude site RE).
#   • Compute max/median absolute differences relative to the baseline curve.
loso_tab <- map_df(top_sites, function(sid) {
  d2 <- df %>% filter(SITE_ID.F != sid)
  ny <- n_distinct(year(d2$SAMPLE_DATE))
  k_time <- min(35, max(12, round(0.7 * ny)))
  m2 <- bam(
    cbind(lower_t, upper_t) ~ season_f + s(decimal_date, k = k_time) + s(SITE_ID.F, bs = "re"),
    family   = cnorm(),
    data     = d2,
    weights  = (1 / pmax(d2$upper_t - d2$lower_t, 1e-6)),
    method   = "fREML",
    discrete = TRUE,
    select   = TRUE,
    gamma    = 1.2
  )
  p2 <- as.numeric(predict(m2, newdata = tgrid, type = "response",
                           exclude = "s(SITE_ID.F)"))
  tibble(site = as.character(sid),
         max_abs_diff = max(abs(p2 - pred0)),
         med_abs_diff = median(abs(p2 - pred0)))
})
print(loso_tab)

# -------------------------------------------------
# 3) Season-varying trend vs single trend (model test)
# -------------------------------------------------
# Compare:
#   Base model  : single smooth over time (same as in res_p035 fit)
#   Alternate   : season-varying smooth s(decimal_date, by = season_f)
# Report AIC and deviance explained; ΔAIC indicates which is preferred.
ny <- n_distinct(df$year)
k_time <- min(35, max(12, round(0.7 * ny)))

m_alt <- bam(
  cbind(lower_t, upper_t) ~ season_f +
    s(decimal_date, by = season_f, k = k_time) +
    s(SITE_ID.F, bs = "re"),
  family   = cnorm(),
  data     = df,
  weights  = (1 / pmax(df$upper_t - df$lower_t, 1e-6)),
  method   = "fREML",
  discrete = TRUE,
  select   = TRUE,
  gamma    = 1.2
)

# Gather headline metrics for both models and print a concise comparison
base_AIC <- AIC(m); alt_AIC <- AIC(m_alt)
base_dev <- summary(m)$dev.expl; alt_dev <- summary(m_alt)$dev.expl

cat(sprintf("\nModel comparison (cnorm p=0.35)\n  Base (single time smooth):   AIC = %.1f, DevExpl = %.1f%%\n  Alt  (season-varying time): AIC = %.1f, DevExpl = %.1f%%\n  ΔAIC (alt - base) = %.1f  -> prefer %s\n",
            base_AIC, 100*base_dev, alt_AIC, 100*alt_dev, alt_AIC - base_AIC,
            ifelse(alt_AIC < base_AIC, "season-varying", "single-smooth")))



               
#4.4.3 CNORM MODEL( P=O.35 WITH SEASON VARYING TREND)

# ================================================
# Cnorm (p = 0.35) with season-varying trend
# ------------------------------------------------
# What this script does:
#   • Fits a censored-normal GAMM on a power-transformed scale (p = 0.35).
#   • Uses precision weights based on interval width (narrower = higher weight).
#   • Includes a season main effect and a season-varying smooth trend over time.
#   • Adds a site random intercept (s(SITE_ID.F, bs="re")).
#   • Runs diagnostics (via your existing diagnose_cnorm()).
#   • Produces population-level trend curves (random effects excluded) and a plot.
#
# Expected columns in `cn_df`:
#   lower, upper          -> raw 1-s.f. bounds on the count scale
#   lower_t, upper_t      -> transformed bounds (if already present)
#   season_f              -> factor with levels c("spring","autumn")
#   decimal_date          -> numeric time (e.g., 1995.5)
#   SITE_ID.F             -> site factor used for the random effect
#   SAMPLE_DATE           -> Date, used to size the smooth basis (k)
# ================================================

# cn_df: your family-specific censored-normal table with
# columns lower, upper, lower_t, upper_t, season_f, decimal_date, SITE_ID.F, SAMPLE_DATE

fit_cnorm_seasonvary <- function(df, p = 0.35, alpha_w = 1) {
  # Ensure transformed bounds exist on the p-scale:
  # If lower_t/upper_t were not precomputed, create them from raw bounds.
  if (!all(c("lower_t","upper_t") %in% names(df))) {
    df <- df %>% mutate(
      lower_t = pmax(as.numeric(lower), 0)^p,  # guard against negatives before transform
      upper_t = pmax(as.numeric(upper), 0)^p
    )
  }
  # Precision weights: intervals with smaller width on the p-scale carry more weight.
  # Raising by alpha_w lets you increase/decrease emphasis on narrow intervals.
  wts <- (1 / pmax(df$upper_t - df$lower_t, 1e-6))^alpha_w

  # Choose the size of the time smooth (k) from the span of years in the data.
  k_time <- {
    ny <- n_distinct(lubridate::year(df$SAMPLE_DATE))
    min(35, max(12, round(0.7 * ny)))
  }

  # Fit the cnorm BAM:
  #   • season_f main effect
  #   • season-varying smooth of decimal_date
  #   • site random intercept
  #   • interval-width weights
  m <- bam(
    cbind(lower_t, upper_t) ~ season_f +
      s(decimal_date, by = season_f, k = k_time) +
      s(SITE_ID.F, bs = "re"),
    family   = cnorm(),
    data     = df,
    weights  = wts,
    method   = "fREML",
    discrete = TRUE,
    select   = TRUE,
    gamma    = 1.2
  )

  # Return the fitted model, data (with transformed bounds), and meta info
  list(model = m, data = df, p = p, k_time = k_time)
}

# ---- Fit the final model for the current family
# `cn_df` must already be defined in your environment as described above.
final <- fit_cnorm_seasonvary(cn_df, p = 0.35)

# ---- Run your existing cnorm diagnostics on the final model
# Uses your previously defined diagnose_cnorm(fit, data) helper.
diagnose_cnorm(final$model, final$data)

# ---- Season-specific trend curves (population-level; RE excluded)
# Builds a time grid for each season, makes predictions on the p-scale with
# the site random effect excluded, then back-transforms to a relative index.
make_trend_curves <- function(fit, df, p = 0.35) {
  tgrid <- expand.grid(
    decimal_date = seq(min(df$decimal_date), max(df$decimal_date), length.out = 250),
    season_f     = levels(df$season_f)
  )
  # Dummy site factor value (any valid level is fine since we exclude the RE)
  tgrid$SITE_ID.F <- df$SITE_ID.F[1]

  # Predict on the p-scale; exclude s(SITE_ID.F) to remove random-effect influence.
  pred_p <- as.numeric(predict(fit, newdata = tgrid, type = "response",
                               exclude = "s(SITE_ID.F)"))
  # Back-transform to an index on the original count scale.
  # Note: (E[X^p])^(1/p) is a comparable index, not an unbiased mean count.
  tgrid$index_count <- pmax(pred_p, 0)^(1/p)
  tgrid
}

# Build the trend dataframe and plot the season-specific indices through time
trend_df <- make_trend_curves(final$model, final$data, p = final$p)

ggplot(trend_df, aes(decimal_date, index_count, colour = season_f)) +
  geom_line(linewidth = 0.9) +
  labs(title = "Season-specific trend (relative abundance index)",
       x = "Year", y = "Index (back-transformed from p-scale)",
       colour = "Season") +
  theme_minimal(base_size = 12)

# ---- Quick headline metrics (concise model summary)
sum_final <- summary(final$model)
cat(sprintf("\nFINAL (season-varying) — Dev.expl: %.1f%% | AIC: %.1f | k_time: %d\n",
            100*sum_final$dev.expl, AIC(final$model), final$k_time))



               
# ================================================================
# Batch cnorm(p = 0.35) for 4 families
#  - Prints full summary() and gam.check() text for each family
#  - No diagnostic plots are shown
#  - Returns/prints a compact summary table at the end
# Requirements:
#   models_data[["<family>"]]$cnorm  (with lower/upper bounds, season_f,
#   decimal_date, SITE_ID.F, SAMPLE_DATE)
# ------------------------------------------------
# What this script provides:
#   • A helper to power-transform interval bounds and build precision weights.
#   • A helper to compute randomised-quantile residuals for cnorm fits.
#   • A family fitter that:
#       – guards for missing season/site factors,
#       – sizes the time-smooth basis from the year span,
#       – fits a season-varying time smooth + site RE with weights,
#       – prints textual diagnostics (summary + gam.check text),
#       – returns a compact list for tabulation.
#   • A loop over your four families and a tidy summary table.
# Notes:
#   • No figures are produced here (gam.check plots are diverted to a temp PDF).
#   • The weights favour narrower transformed intervals (1/width on p-scale).
#   • The “index of fit tightness” we print is SD of RQ residuals (ideal ≈ 1).
# ================================================================

FAMILIES <- c("Aphelocheiridae","Brachycentridae","Odontoceridae","Cordulegastridae")

# --- helpers ----------------------------------------------------

# Ensure p-power transformed bounds & precision weights
# ----------------------------------------------------
# Expects raw count-scale interval columns `lower` and `upper`.
# Produces:
#   • lower_t, upper_t  -> bounds on the p-scale (p=0.35),
#   • width_t           -> interval width on p-scale (clamped at 1e-6),
#   • wt                -> precision weight = 1/width_t,
#   • year              -> convenience field for choosing k.
ensure_bounds_p <- function(df, p = 0.35) {
  stopifnot(all(c("lower","upper") %in% names(df)))
  df %>%
    mutate(
      lower_t = pmax(as.numeric(lower), 0)^p,
      upper_t = pmax(as.numeric(upper), 0)^p,
      width_t = pmax(upper_t - lower_t, 1e-6),
      wt      = 1 / width_t,                 # narrower intervals => larger weight
      year    = year(SAMPLE_DATE)
    )
}

# Randomised-quantile residuals for censored Normal (on p-scale)
# --------------------------------------------------------------
# Given a fitted cnorm model and transformed bounds (lt/ut),
# draw a uniform from the conditional CDF slice and map via qnorm.
# Use as a quick dispersion/shape diagnostic (SD ≈ 1 is ideal).
cnorm_rqr <- function(fit, lt, ut, eps = 1e-10, seed = 123) {
  mu <- as.numeric(fitted(fit, type = "response"))
  sd <- sqrt(summary(fit)$scale)
  Flo <- pnorm(lt, mean = mu, sd = sd)
  Fup <- pnorm(ut, mean = mu, sd = sd)
  set.seed(seed)
  u <- pmin(pmax(Flo + pmax(Fup - Flo, 0) * runif(length(mu)), eps), 1 - eps)
  qnorm(u)
}

# Fit one family and PRINT summary + gam.check (no plots)
# -------------------------------------------------------
# Inputs:
#   fam     -> family name key in models_data
#   md      -> list-like container with models_data[[fam]]$cnorm
#   p       -> power for transformation (fixed to 0.35 as per design)
# Side effects:
#   • Prints model summary and textual gam.check output.
# Returns:
#   • Compact list (model, data, diagnostics) for tabulation.
fit_cnorm_family_print <- function(fam, md = models_data, p = 0.35) {
  stopifnot(!is.null(md[[fam]]), "cnorm" %in% names(md[[fam]]))
  df <- md[[fam]]$cnorm %>% ensure_bounds_p(p)

  # guard rails: season factor & site factor
  # ------------------------------------------------
  # If season_f was not created upstream, derive spring/autumn from months.
  # Ensure site IDs are a proper factor for the random intercept term.
  if (!("season_f" %in% names(df))) {
    df <- df %>%
      mutate(m = month(SAMPLE_DATE),
             season_f = factor(if_else(m %in% 3:5, "spring", "autumn"),
                               levels = c("spring","autumn")))
  }
  if (!is.factor(df$SITE_ID.F)) df$SITE_ID.F <- factor(df$SITE_ID)

  # choose k from span of years
  # ------------------------------------------------
  # Basis dimension increases with the number of distinct sampling years,
  # but is clamped to [12, 35] to avoid under/over-fitting extremes.
  k_time <- {
    ny <- n_distinct(df$year)
    min(35, max(12, round(0.7 * ny)))
  }

  # season-varying time smooth + site RE; precision weights
  # -------------------------------------------------------
  # The cnorm family takes cbind(lower_t, upper_t) as the response.
  # We include:
  #   • season_f (main effect),
  #   • s(decimal_date, by = season_f, k = k_time) for different trends by season,
  #   • s(SITE_ID.F, bs="re") for a site random intercept,
  #   • df$wt to emphasize tighter intervals on the p-scale.
  m <- bam(
    cbind(lower_t, upper_t) ~ season_f +
      s(decimal_date, by = season_f, k = k_time) +
      s(SITE_ID.F, bs = "re"),
    family   = cnorm(),
    data     = df,
    weights  = df$wt,
    method   = "fREML",
    discrete = TRUE,
    select   = TRUE,
    gamma    = 1.2
  )

  # ---- PRINTED OUTPUTS (no plots) ----
  cat("\n====================================================\n")
  cat(sprintf("Family: %s  |  cnorm(p = %.2f)\n", fam, p))
  cat("====================================================\n")
  print(summary(m))

  # gam.check text without plotting the 2x2 panel
  # ------------------------------------------------
  # We divert any plots into a temporary PDF (and delete it),
  # keeping only the textual adequacy checks in the console.
  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
  gc_out <- try(gam.check(m, rep = 0), silent = TRUE)
  dev.off(); unlink(tmp)
  if (!inherits(gc_out, "try-error")) {
    cat("\n--- gam.check (text) ---\n")
    print(gc_out)
  }

  # simple numeric diagnostics
  # ------------------------------------------------
  # SD of RQ residuals (ideal ≈ 1), deviance explained, AIC, N, and k used.
  rqr_sd <- sd(cnorm_rqr(m, df$lower_t, df$upper_t), na.rm = TRUE)
  devex  <- summary(m)$dev.expl
  aic    <- AIC(m)

  cat(sprintf("\nSD(RQ residuals): %.2f | Dev.expl: %.1f%% | AIC: %.1f | N: %d | k_time: %d\n",
              rqr_sd, 100*devex, aic, nobs(m), k_time))

  # return compact object for summary table
  list(
    family   = fam,
    model    = m,
    data     = df,
    rqr_sd   = rqr_sd,
    dev_expl = devex,
    AIC      = aic,
    n        = nobs(m),
    k_time   = k_time
  )
}

# --- RUN for all families --------------------------------------
# Map the fitter across the four families and name the output list
# for convenient downstream access (cn_results[["Aphelocheiridae"]] etc.).
cn_results <- map(FAMILIES, fit_cnorm_family_print) %>% set_names(FAMILIES)

# --- SUMMARY TABLE ---------------------------------------------
# Bind the compact diagnostics from each family into a single table and
# print a clean, rounded view suitable for a log or appendix.
cn_summary <- bind_rows(lapply(cn_results, function(x) {
  data.frame(
    family   = x$family,
    n        = x$n,
    k_time   = x$k_time,
    rqr_sd   = round(x$rqr_sd, 2),
    dev_expl = round(100 * x$dev_expl, 1),
    AIC      = round(x$AIC, 1),
    stringsAsFactors = FALSE
  )
}))
cat("\n==================== Summary table (cnorm, p=0.35) ====================\n")
print(cn_summary, row.names = FALSE)



               
# ===============================================================
# Quick check: cnorm with p = 0.35 and LIGHTER interval weighting
# ---------------------------------------------------------------
# What this script does
#   • Keeps your chosen power p = 0.35 (transforming count-scale bounds).
#   • Uses precision weights raised to alpha_w = 0.5 (so wide intervals are
#     still down-weighted, but less aggressively than alpha_w = 1).
#   • Fits the same cnorm GAMM form as before (season-varying smooth + site RE).
#   • Reports SD of randomised-quantile residuals (≈1 ideal), AIC, deviance
#     explained and k used for the time smooth, for each family.
#
# Inputs expected (built earlier in your pipeline)
#   • models_data[["<family>"]]$cnorm with: lower/upper (count-scale bounds),
#     season_f, decimal_date, SITE_ID.F, SAMPLE_DATE (and optionally lower_t/upper_t).
# ===============================================================

# Families to evaluate (same set you’ve been using)
families <- c("Aphelocheiridae","Brachycentridae","Odontoceridae","Cordulegastridae")

# ---------------------------------------------------------------
# Helper: extract sigma from a fitted cnorm model
#   • mgcv prints the family as "cnorm(<sigma>)" in the family string.
#   • We parse that out so we can compute RQ residuals with the correct SD.
# ---------------------------------------------------------------
get_sigma <- function(fit) {
  fam_str <- fit$family$family
  as.numeric(sub(".*cnorm\\(([^)]+)\\).*", "\\1", fam_str))
}

# ---------------------------------------------------------------
# Helper: SD of randomised-quantile residuals on the p-scale
#   • Computes PIT for the censored Normal on transformed scale using
#     the model’s fitted mean (mu) and sigma from get_sigma().
#   • Returns the SD of qnorm(U); ideal value is ~1 if dispersion is OK.
# ---------------------------------------------------------------
rqr_sd_corrected <- function(fit, lt, ut) {
  mu <- as.numeric(fitted(fit, type = "response"))
  sd <- get_sigma(fit)
  u  <- pmin(pmax(pnorm(lt, mu, sd) + (pnorm(ut, mu, sd) - pnorm(lt, mu, sd)) * runif(length(mu)), 1e-10), 1 - 1e-10)
  sd(qnorm(u), na.rm = TRUE)
}

# ---------------------------------------------------------------
# Core fitter for one dataset with alpha_w = 0.5
#   Arguments:
#     df       : family-specific cnorm table (from models_data[[fam]]$cnorm)
#     p        : power for transforming raw bounds (default 0.35)
#     alpha_w  : exponent for precision weights (default 0.5 here)
#
#   Steps:
#     1) Ensure transformed bounds exist (lower_t/upper_t on p-scale).
#     2) Build precision weights = (1/width)^alpha_w on p-scale.
#     3) Choose k for time smooth from span of years.
#     4) Fit cnorm with season main effect, season-varying time smooth,
#        and site random intercept; use weights.
#     5) Return compact metrics for quick comparison.
# ---------------------------------------------------------------
fit_cnorm_alpha <- function(df, p = 0.35, alpha_w = 0.5) {
  if (!all(c("lower_t","upper_t") %in% names(df))) {
    df <- df %>% mutate(lower_t = pmax(as.numeric(lower),0)^p,
                        upper_t = pmax(as.numeric(upper),0)^p)
  }
  df <- df %>% mutate(wt = (1 / pmax(upper_t - lower_t, 1e-6))^alpha_w,
                      year = year(SAMPLE_DATE))
  ny <- n_distinct(df$year)
  k_time <- min(35, max(12, round(0.7 * ny)))
  m <- bam(
    cbind(lower_t, upper_t) ~ season_f + s(decimal_date, by = season_f, k = k_time) + s(SITE_ID.F, bs = "re"),
    family   = cnorm(), data = df, weights  = df$wt,
    method   = "fREML", discrete = TRUE, select = TRUE, gamma = 1.2
  )
  data.frame(
    SD_RQR = round(rqr_sd_corrected(m, df$lower_t, df$upper_t), 2),
    AIC    = round(AIC(m), 1),
    DevExpl = round(100*summary(m)$dev.expl, 1),
    k_time = k_time
  )
}

# ---------------------------------------------------------------
# Run for all four families and print a tidy comparison table
#   • Each row: family, SD(RQR), AIC, % deviance explained, k_time.
# ---------------------------------------------------------------
alpha05_results <- map_dfr(families, function(fam) {
  df <- models_data[[fam]]$cnorm
  out <- fit_cnorm_alpha(df, p = 0.35, alpha_w = 0.5)
  cbind(family = fam, out)
})
print(alpha05_results, row.names = FALSE)



               
# ================================================================
# Tune p (common grid) and weight exponent α for all families
# - Uses models_data[[fam]]$cnorm (with lower/upper, season_f, decimal_date, SITE_ID.F, SAMPLE_DATE)
# - Tries p ∈ {0.35, 0.40, 0.45} and α ∈ {0, 0.5}
# - Reports SD(RQR) (using fitted sigma), AIC, DevExpl, and suggests best
# ================================================================

FAMILIES <- c("Aphelocheiridae","Brachycentridae","Odontoceridae","Cordulegastridae")
P_GRID   <- c(0.35, 0.40, 0.45)
ALPHAS   <- c(0.0, 0.5)

# --- helpers ---
# Extract σ from a fitted cnorm() family string, e.g. "cnorm(0.89)" -> 0.89
get_sigma <- function(fit) {
  fam_str <- fit$family$family
  as.numeric(sub(".*cnorm\\(([^)]+)\\).*", "\\1", fam_str))
}

# Compute SD of randomised-quantile residuals for censored Normal on the transformed scale.
# Uses fitted μ and σ from the model; clamps PIT away from {0,1}; SD ≈ 1 indicates good dispersion.
rqr_sd_corrected <- function(fit, lt, ut) {
  mu <- as.numeric(fitted(fit, type = "response"))
  sd <- get_sigma(fit)
  u  <- pmin(pmax(pnorm(lt, mu, sd) + (pnorm(ut, mu, sd) - pnorm(lt, mu, sd)) * runif(length(mu)), 1e-10), 1 - 1e-10)
  sd(qnorm(u), na.rm = TRUE)
}

# Prepare bounds and weights for a given power p and weight exponent alpha_w:
# - (Re)build transformed bounds lower_t/upper_t = lower^p / upper^p (lower/upper on count-scale)
# - width_t = interval width on transformed scale
# - wt_raw = 1/width_t; wt = wt_raw^alpha_w, lightly capped at 95th percentile to avoid extremes
# - year column used to set k for time smooth
prep_bounds <- function(df, p, alpha_w) {
  if (!all(c("lower_t","upper_t") %in% names(df))) {
    df <- df %>% mutate(lower_t = pmax(as.numeric(lower),0)^p,
                        upper_t = pmax(as.numeric(upper),0)^p)
  } else {
    # If lower_t/upper_t exist but for a different p, recompute for current p
    df <- df %>% mutate(lower_t = pmax(as.numeric(lower),0)^p,
                        upper_t = pmax(as.numeric(upper),0)^p)
  }
  df %>%
    mutate(
      width_t = pmax(upper_t - lower_t, 1e-6),
      wt_raw  = 1 / width_t,
      # weight exponent and gentle capping to avoid a few huge weights dominating
      wt      = (wt_raw^alpha_w) %>% pmin(quantile(., 0.95, na.rm = TRUE)),
      year    = year(SAMPLE_DATE)
    )
}

# Fit one cnorm model at a specific (p, alpha_w) setting:
# - Ensures season_f (spring/autumn) and SITE_ID.F exist
# - Chooses k_time from span of years (slightly generous)
# - Fits: season main effect + season-varying time smooth + site RE, with precision weights
# - Returns a one-row data.frame of metrics (SD_RQR, AIC, DevExpl, k_time, n)
fit_once <- function(df, p, alpha_w) {
  d <- prep_bounds(df, p, alpha_w)
  if (!("season_f" %in% names(d))) {
    d <- d %>% mutate(m = month(SAMPLE_DATE),
                      season_f = factor(if_else(m %in% 3:5, "spring", "autumn"),
                                        levels = c("spring","autumn")))
  }
  if (!is.factor(d$SITE_ID.F)) d$SITE_ID.F <- factor(d$SITE_ID)

  ny <- n_distinct(d$year)
  k_time <- min(40, max(14, round(0.8 * ny)))   # slightly larger k to appease low k-index families

  m <- bam(
    cbind(lower_t, upper_t) ~ season_f +
      s(decimal_date, by = season_f, k = k_time) +
      s(SITE_ID.F, bs = "re"),
    family   = cnorm(),
    data     = d,
    weights  = d$wt,
    method   = "fREML",
    discrete = TRUE,
    select   = TRUE,
    gamma    = 1.2
  )

  data.frame(
    p        = p,
    alpha_w  = alpha_w,
    SD_RQR   = round(rqr_sd_corrected(m, d$lower_t, d$upper_t), 2),
    AIC      = round(AIC(m), 1),
    DevExpl  = round(100*summary(m)$dev.expl, 1),
    k_time   = k_time,
    n        = nobs(m),
    stringsAsFactors = FALSE
  )
}

# Grid-search tuner for one family:
# - Evaluates all p×α combinations
# - Picks the setting with SD_RQR closest to 1, breaking ties by lowest AIC
# - Returns the full table and the chosen 'best' row
tune_family <- function(fam) {
  base <- models_data[[fam]]$cnorm
  grid <- expand.grid(p = P_GRID, alpha_w = ALPHAS, KEEP.OUT.ATTRS = FALSE)
  res  <- purrr::pmap_dfr(grid, ~ fit_once(base, ..1, ..2))
  # pick SD_RQR closest to 1 (tie-break by lowest AIC)
  res <- res %>% mutate(dist = abs(SD_RQR - 1))
  best <- res %>% filter(dist == min(dist)) %>% slice_min(AIC, n = 1)
  list(table = res, best = best)
}

# ---- run
# Iterate families, store the tuning table + best choice per family; attach family name
tuned <- lapply(FAMILIES, function(f) { out <- tune_family(f); out$family <- f; out })
names(tuned) <- FAMILIES

# ---- print per-family tables and choices
# For each family: print the full results table, then a concise "Recommended" line
for (f in FAMILIES) {
  cat("\n================ ", f, " ================\n", sep = "")
  print(tuned[[f]]$table, row.names = FALSE)
  cat("→ Recommended: ",
      with(tuned[[f]]$best, sprintf("p=%.2f, α=%.1f (SD_RQR=%0.2f, AIC=%0.1f, Dev=%.1f%%, k=%d, n=%d)",
                                    p, alpha_w, SD_RQR, AIC, DevExpl, k_time, n)),
      "\n", sep = "")
}

# ---- compact summary of chosen settings
# Bind all 'best' rows into one data.frame for quick comparison across families
chosen <- do.call(rbind, lapply(tuned, function(x) cbind(family = x$family, x$best)))
cat("\n========= Chosen settings per family =========\n")
print(chosen, row.names = FALSE)



#4.4.4 FINAL CNORM MODEL FOR ALL FAMILIES WITH NO INTERVAL WEIGHTING
# ================================
# Final cnorm fits per family (no plots)
# p per family; alpha = 0 (no interval weighting)
# Prints summary() and textual gam.check()
# ================================

# Family-specific p (from your tuning)
final_p <- tibble::tribble(
  ~family,             ~p,
  "Aphelocheiridae",   0.35,
  "Brachycentridae",   0.35,
  "Odontoceridae",     0.40,
  "Cordulegastridae",  0.40
)

# Helper: build p-transformed bounds, no weights, guard season/site factors
prep_cnorm_df <- function(df, p) {
  stopifnot(all(c("lower","upper") %in% names(df)))
  df %>%
    mutate(
      lower_t = pmax(as.numeric(lower), 0)^p,        # transform count-scale lower bound to p-scale
      upper_t = pmax(as.numeric(upper), 0)^p         # transform count-scale upper bound to p-scale
    ) %>%
    { if (!("season_f" %in% names(.)))               # if season factor missing, build spring/autumn from month
        mutate(., m = month(SAMPLE_DATE),
                  season_f = factor(if_else(m %in% 3:5, "spring", "autumn"),
                                    levels = c("spring","autumn")))
      else . } %>%
    { if (!is.factor(.$SITE_ID.F)) mutate(., SITE_ID.F = factor(SITE_ID)) else . } %>%  # ensure site RE factor exists
    mutate(year = year(SAMPLE_DATE))                  # convenience: calendar year for k selection
}

# Correct-sigma RQ residual SD
get_sigma <- function(fit) {
  fam_str <- fit$family$family
  as.numeric(sub(".*cnorm\\(([^)]+)\\).*", "\\1", fam_str))  # parse sigma from cnorm(<sigma>)
}
sd_rqr_corrected <- function(fit, lt, ut) {
  mu <- as.numeric(fitted(fit, type = "response"))           # fitted mean on p-scale
  sd <- get_sigma(fit)                                       # model sigma on p-scale
  u  <- pmin(pmax(pnorm(lt, mu, sd) + (pnorm(ut, mu, sd) - pnorm(lt, mu, sd)) * runif(length(mu)), 1e-10), 1 - 1e-10)
  sd(qnorm(u), na.rm = TRUE)                                 # SD of PIT→Normal = dispersion cue (~1 ideal)
}

fit_cnorm_final <- function(fam, p, md = models_data) {
  stopifnot(!is.null(md[[fam]]), "cnorm" %in% names(md[[fam]]))
  df <- prep_cnorm_df(md[[fam]]$cnorm, p)                    # prepare df (bounds→p-scale, season/site guards)

  # k chosen from span of years; a touch generous helps low k-index cases
  ny <- dplyr::n_distinct(df$year)
  k_time <- min(40, max(14, round(0.8 * ny)))

  m <- bam(
    cbind(lower_t, upper_t) ~ season_f +
      s(decimal_date, by = season_f, k = k_time) +          # season-varying smooth over time
      s(SITE_ID.F, bs = "re"),                              # site random intercept
    family   = cnorm(),
    data     = df,
    method   = "fREML",
    discrete = TRUE,
    select   = TRUE,
    gamma    = 1.2
  )

  # ---- print full outputs (no plots) ----
  cat("\n====================================================\n")
  cat(sprintf("Family: %s  |  cnorm(p = %.2f)\n", fam, p))
  cat("====================================================\n")
  print(summary(m))                                          # detailed param/smooth summary

  tmp <- tempfile(fileext = ".pdf"); pdf(tmp)                # suppress 2×2 gam.check panel to text-only
  gc_out <- try(gam.check(m, rep = 0), silent = TRUE)
  dev.off(); unlink(tmp)
  if (!inherits(gc_out, "try-error")) {
    cat("\n--- gam.check (text) ---\n")
    print(gc_out)
  }

  # compact metrics
  rqr_sd <- sd_rqr_corrected(m, df$lower_t, df$upper_t)      # dispersion ≈1 good
  aic    <- AIC(m); devx <- summary(m)$dev.expl               # parsimony + fit
  cat(sprintf("\nSD(RQR): %.2f | Dev.expl: %.1f%% | AIC: %.1f | N: %d | k_time: %d\n",
              rqr_sd, 100*devx, aic, nobs(m), k_time))

  list(family = fam, p = p, model = m, data = df,
       rqr_sd = rqr_sd, AIC = aic, dev_expl = devx, k_time = k_time, n = nobs(m))
}

# Run all four and build a summary table
final_results <- pmap(final_p, ~ fit_cnorm_final(..1, ..2))  # iterate families with their chosen p
names(final_results) <- final_p$family

final_summary <- dplyr::bind_rows(lapply(final_results, function(x) {
  data.frame(
    family   = x$family,
    p        = x$p,
    n        = x$n,
    k_time   = x$k_time,
    SD_RQR   = round(x$rqr_sd, 2),
    DevExpl  = round(100*x$dev_expl, 1),
    AIC      = round(x$AIC, 1),
    stringsAsFactors = FALSE
  )
}))
cat("\n==================== Final cnorm summary ====================\n")
print(final_summary, row.names = FALSE)                      # compact per-family comparison



  
# ================================================================
# CNORM — per-family Model Diagnostics + partial-effect plots
# Assumes: final_results (named list) with elements containing:
#   $family, $p, $model (bam fit), $data (with lower_t, upper_t)
# ----------------------------------------------------------------
# What this script does:
# - Defines helpers to extract the cnorm sigma and to compute
#   randomised-quantile residuals (RQR) on the p-scale.
# - For each fitted family model in `final_results`, it:
#     * draws a QQ plot of RQRs,
#     * prints Residuals vs Fitted and Residual histogram panels,
#     * shows residual ACF,
#     * prints season-specific time smooth partial effects,
#     * and prints a compact console footer with SD(RQR), DevExpl, AIC, N.
# Notes:
# - No model refitting happens here; we only diagnose fitted objects.
# - Random effects are not excluded in partial-effect plots (these are
#   standard gratia::draw() outputs for the smooth terms).
# ================================================================

# ---- helpers ----------------------------------------------------
# Purpose: parse the fitted Normal sigma from mgcv's cnorm family string,
# e.g. "cnorm(0.612)". Falls back to 1.0 if parsing fails.
get_cnorm_sigma <- function(fit) {
  fam_str <- fit$family$family
  m <- regexpr("cnorm\\(([^)]+)\\)", fam_str)
  if (m > 0) as.numeric(sub("cnorm\\(([^)]+)\\)", "\\1", regmatches(fam_str, m))) else 1.0
}

# Purpose: compute RQRs using the *fitted* sigma on the p-scale.
# Inputs are transformed bounds (lower_t/upper_t) on the same p-scale as the model.
cnorm_rqr <- function(fit, lower_t, upper_t, eps = 1e-10, seed = 123) {
  mu <- as.numeric(fitted(fit, type = "response"))
  sd <- get_cnorm_sigma(fit)
  Flo <- pnorm(lower_t, mean = mu, sd = sd)
  Fup <- pnorm(upper_t, mean = mu, sd = sd)
  set.seed(seed)
  u <- pmin(pmax(Flo + pmax(Fup - Flo, 0) * runif(length(mu)), eps), 1 - eps)
  qnorm(u)
}

# Make a single family's full diagnostic + partial-effect plots
# Expects `res` to be one element of final_results with $model, $data, $family, $p.
plot_cnorm_family <- function(res) {
  fit <- res$model
  df  <- res$data
  fam <- res$family
  pwr <- res$p

  # --- RQ residuals (correct sigma) ---
  rqr <- cnorm_rqr(fit, df$lower_t, df$upper_t)
  mu  <- fitted(fit, type = "response")

  # 1) QQ plot (base graphics) — visual check for ~Normal RQR
  par(mfrow = c(1,1))
  qqnorm(rqr, main = paste0("QQ plot — RQ residuals (", fam, ", p=", pwr, ")"))
  qqline(rqr)

  # Prep a small tibble for ggplot panels
  dd <- data.frame(mu = mu, rqr = rqr)

  # 2) Residuals vs Fitted (no smooth line) — check mean/variance structure
  print(
    ggplot(dd, aes(mu, rqr)) +
      geom_point(alpha = 0.15, size = 0.6) +
      geom_hline(yintercept = 0, linetype = 2) +
      labs(
        title = paste0("RQR vs Fitted (sqrt/p-scale) — ", fam, " (p=", pwr, ")"),
        x = "Fitted mean (sqrt/p-scale)", y = "Randomised quantile residual"
      ) +
      theme_minimal(base_size = 12)
  )

  # 3) Histogram of residuals — check approximate symmetry/scale
  print(
    ggplot(dd, aes(rqr)) +
      geom_histogram(bins = 40, fill = "#3D5A80", colour = "white", linewidth = 0.15) +
      labs(
        title = paste0("Histogram of RQ residuals — ", fam, " (p=", pwr, ")"),
        x = "RQ residual", y = "Frequency"
      ) +
      theme_minimal(base_size = 12)
  )

  # 4) ACF of residuals — check serial dependence
  acf(rqr[is.finite(rqr)], na.action = na.pass,
      main = paste0("ACF of RQ residuals — ", fam, " (p=", pwr, ")"))

  # 5) Partial effects — season-specific time smooths (printed separately)
  sm_ids   <- which(grepl("decimal_date", smooths(fit)))
  sm_names <- smooths(fit)[sm_ids]
  for (j in seq_along(sm_ids)) {
    season_lab <- if (grepl("spring", sm_names[j], ignore.case = TRUE)) "spring"
      else if (grepl("autumn", sm_names[j], ignore.case = TRUE)) "autumn"
      else sm_names[j]
    p_sm <- draw(fit, select = sm_ids[j]) +
      ggtitle(paste0("Time smooth (", season_lab, ") — ", fam, " (p=", pwr, ")"))
    print(p_sm)
  }

  # Console footer: quick calibration/fit cues for the family
  cat(sprintf("\n[%s] SD(RQR)=%.2f | Dev.expl=%.1f%% | AIC=%.1f | N=%d\n",
              fam, sd(rqr, na.rm = TRUE), 100*summary(fit)$dev.expl,
              AIC(fit), nobs(fit)))
  invisible(NULL)
}

# ---- RUN for all families --------------------------------------
# Guard: ensure final_results exists and is a list; then loop and plot for each family.
stopifnot(exists("final_results"), is.list(final_results))
invisible(lapply(final_results, plot_cnorm_family))



  
# ================================================================
# National marginal trend prediction from final cnorm models
# ---------------------------------------------------------------
# What this script does
# 1) Takes each fitted cnorm model in `final_results` (one per family).
# 2) Builds a date grid and generates *marginal* predictions:
#      - excludes the site random effect via `exclude = "s(SITE_ID.F)"`
#      - predicts separately for spring and autumn, plus a 50/50 seasonal average
# 3) Simulates parameter uncertainty from the model’s covariance to produce
#    95% intervals on the p-scale, and back-transforms to a count-scale index.
# 4) Binds all families together, writes a CSV, and (optionally) plots.
#
# Notes
# - The “count-scale” here is a back-transformed index: (E[Y^p])^(1/p).
#   It is intended for relative comparisons over time; it is not an unbiased
#   estimator of E[Y] on the raw-count scale.
# - We set the site factor to a *real observed* level then exclude the RE,
#   which avoids “new level” warnings while still giving marginal curves.
# - Plots are optional; set MAKE_PLOTS <- FALSE to skip.
# ================================================================


MAKE_PLOTS <- TRUE

# --- helpers ----------------------------------------------------
# Extract the fitted Normal sigma used by mgcv's cnorm() family from the model object.
# (Kept here for completeness; this script works off the linear predictor draws.)
get_sigma_cnorm <- function(fit) {
  fam_str <- fit$family$family
  m <- regexpr("cnorm\\(([^)]+)\\)", fam_str)
  if (m > 0) as.numeric(sub("cnorm\\(([^)]+)\\)", "\\1", regmatches(fam_str, m))) else 1.0
}

# Draw from the posterior of the linear predictor (η) at `newdata`
# by simulating coefficients ~ N(b̂, V̂) and multiplying by the model matrix.
# - `exclude` can drop specific smooths (e.g., the site RE) from the lpmatrix.
# - Returns a matrix [nsim × n_new] of linear predictor draws on the model's p-scale.
posterior_mu <- function(fit, newdata, nsim = 500, exclude = NULL, seed = 123, unconditional = FALSE) {
  X <- predict(fit, newdata = newdata, type = "lpmatrix", exclude = exclude)
  b <- coef(fit)
  V <- vcov(fit, unconditional = unconditional)
  set.seed(seed)
  B <- MASS::mvrnorm(n = nsim, mu = b, Sigma = V)
  as.matrix(B %*% t(X))
}

# Create a regular date grid spanning the observed sample dates in a family dataset.
# Adds decimal_date for the time smooth.
make_grid <- function(df, n = 200) {
  rng <- range(df$SAMPLE_DATE, na.rm = TRUE)
  dates <- seq(rng[1], rng[2], length.out = n)
  tibble(
    SAMPLE_DATE = dates,
    decimal_date = lubridate::decimal_date(dates)
  )
}

# Summarise a matrix of draws into mean and 95% CI on the p-scale,
# then back-transform to a count-scale *index* (^(1/p)).
# Includes the family/season/date fields for plotting/aggregation.
summarise_draws <- function(mu_draws, family, season, dates, p) {
  qs <- t(apply(mu_draws, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
  mu <- colMeans(mu_draws, na.rm = TRUE)
  tibble::tibble(
    family = family,
    season = season,
    SAMPLE_DATE = dates,
    decimal_date = lubridate::decimal_date(dates),
    mu_p  = mu,            # mean on p-scale
    lo_p  = qs[,1],        # 2.5% p-scale
    med_p = qs[,2],        # 50%  p-scale
    hi_p  = qs[,3],        # 97.5% p-scale
    mu_count = pmax(mu, 0)^(1/p),       # back-transform index
    lo_count = pmax(qs[,1], 0)^(1/p),   # back-transform index (lower)
    hi_count = pmax(qs[,3], 0)^(1/p)    # back-transform index (upper)
  )
}

# For one family's fitted model:
# - Build a time grid, predict spring & autumn marginal curves (exclude site RE),
# - Optionally compute a season-average curve by averaging spring and autumn draws,
# - Return a tidy tibble of summaries on both p-scale and back-transformed index.
predict_marginal <- function(res, nsim = 500, seed = 123, make_avg = TRUE) {
  fit <- res$model
  df  <- res$data
  fam <- res$family
  pwr <- res$p

  grid <- make_grid(df, n = 200)

  # Choose a *real* site level to define SITE_ID.F in newdata, then exclude the RE.
  # This avoids factor-level warnings while producing marginal predictions.
  first_site <- as.character(df$SITE_ID.F[1])
  site_levels <- levels(df$SITE_ID.F)
  if (!first_site %in% site_levels) first_site <- site_levels[1]

  seasons <- c("spring","autumn")
  out_list <- vector("list", length(seasons))

  for (i in seq_along(seasons)) {
    newd <- grid %>%
      dplyr::mutate(
        season_f = factor(seasons[i], levels = c("spring","autumn")),
        SITE_ID.F = factor(first_site, levels = site_levels)
      )
    # Exclude the site random effect to get population-level (marginal) curves.
    draws <- posterior_mu(fit, newd, nsim = nsim, exclude = "s(SITE_ID.F)", seed = seed + i)
    out_list[[i]] <- summarise_draws(draws, fam, seasons[i], newd$SAMPLE_DATE, pwr)
  }

  res_seasons <- dplyr::bind_rows(out_list)

  if (!make_avg) return(res_seasons)

  # Make a season-average by averaging draws (preserves uncertainty properly).
  newd_spr <- grid %>%
    dplyr::mutate(season_f = factor("spring", levels = c("spring","autumn")),
                  SITE_ID.F = factor(first_site, levels = site_levels))
  newd_aut <- grid %>%
    dplyr::mutate(season_f = factor("autumn", levels = c("spring","autumn")),
                  SITE_ID.F = factor(first_site, levels = site_levels))

  Dspr <- posterior_mu(fit, newd_spr, nsim = nsim, exclude = "s(SITE_ID.F)", seed = seed + 101)
  Daut <- posterior_mu(fit, newd_aut, nsim = nsim, exclude = "s(SITE_ID.F)", seed = seed + 102)
  Davg <- (Dspr + Daut) / 2

  res_avg <- summarise_draws(Davg, fam, "avg", grid$SAMPLE_DATE, pwr)

  dplyr::bind_rows(res_seasons, res_avg)
}

# ---- run for all families --------------------------------------
# Safety: require `final_results` to exist and be a list of fitted family objects.
stopifnot(exists("final_results"), is.list(final_results))

preds_marg_list <- lapply(final_results, function(res) {
  # Wrap in try() so one failing family does not stop the rest; keeps family name in error.
  try(predict_marginal(res, nsim = 600, seed = 2025), silent = TRUE)
})

# Drop any families that failed (e.g., due to missing fields)
preds_marg_list <- purrr::compact(preds_marg_list)

# Bind all families together and retain the key columns for saving/plotting
preds_marg_tbl <- dplyr::bind_rows(preds_marg_list) %>%
  dplyr::select(family, season, SAMPLE_DATE, mu_p, lo_p, hi_p, mu_count, lo_count, hi_count)

# Persist results for the report
dir.create("outputs", showWarnings = FALSE, recursive = TRUE)
readr::write_csv(preds_marg_tbl, file.path("outputs", "cnorm_predictions_marginal.csv"))

# ---- minimal plots (optional) ----------------------------------
# For each family:
#   1) Show spring & autumn marginal curves with 95% ribbons on the p-scale.
#   2) If the seasonal average was computed, show that as a single curve.
if (MAKE_PLOTS) {
  fams <- unique(preds_marg_tbl$family)
  for (fam in fams) {
    dfp <- preds_marg_tbl %>% dplyr::filter(family == fam)

    # Spring & autumn curves (p-scale)
    p1 <- ggplot(dfp %>% dplyr::filter(season %in% c("spring","autumn")),
                 aes(SAMPLE_DATE, mu_p, colour = season, fill = season)) +
      geom_ribbon(aes(ymin = lo_p, ymax = hi_p), alpha = 0.15, colour = NA) +
      geom_line(linewidth = 0.6) +
      labs(title = paste0("National trend (p-scale) — ", fam),
           x = NULL, y = "Mean abundance (p-scale)") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
    print(p1)

    # Season-average (if present)
    if ("avg" %in% unique(dfp$season)) {
      p2 <- ggplot(dfp %>% dplyr::filter(season == "avg"),
                   aes(SAMPLE_DATE, mu_p)) +
        geom_ribbon(aes(ymin = lo_p, ymax = hi_p), alpha = 0.18) +
        geom_line(linewidth = 0.7) +
        labs(title = paste0("National trend (season-average, p-scale) — ", fam),
             x = NULL, y = "Mean abundance (p-scale)") +
        theme_minimal(base_size = 12)
      print(p2)
    }
  }
}

