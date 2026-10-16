# =============================================================================
# PREDICT2 — Figures
#
# Produces the five figures used in the paper, as PDFs:
#   1. Fig1_postprand_combined.pdf  — postprandial response spaghetti plots
#   2. Fig2_rsq_combined.pdf        — bootstrapped adjusted-R² by feature domain
#   3. Fig3_correlation_heatmap.pdf — Pearson correlations between markers
#   4. Fig4_ML_combined.pdf         — ML predicted-vs-measured scatter plots
#   5. Fig5_ROC_combined.pdf        — ASCVD ROC curves (fasting vs post-prandial)
#
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(gridExtra)
  library(ranger)
  library(xgboost)
  library(pROC)
  library(ggcorrplot)
})

# ── constants ────────────────────────────────────────────────────────────────
CPEP_CONV <- 331.094   # Quest CPEP (ng/mL) → clinic c_peptide (pmol/L)
N_BOOT    <- 1000

# ── paths (resolved relative to this script, so any working dir works) ───────
get_script_dir <- function() {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg))
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable())
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  getwd()
}
SCRIPT_DIR <- get_script_dir()                 # folder this script lives in
DATA <- file.path(SCRIPT_DIR, "input")         # bundled inputs: input/figureN/...
OUT  <- SCRIPT_DIR                             # outputs land next to this script
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── shared helpers ───────────────────────────────────────────────────────────

# Parse microbiome sample id → username (mirrors the Python notebooks exactly:
#   idx.split('__')[0].replace('s','predict').replace('15','').replace('profile','')[:-1])
parse_mb_username <- function(idx) {
  x <- sub("__.*$", "", idx)         # split('__')[0]
  x <- gsub("s", "predict", x)       # .replace('s','predict')
  x <- gsub("15", "", x)             # .replace('15','')
  x <- gsub("profile", "", x)        # .replace('profile','')
  x <- substr(x, 1L, nchar(x) - 1L)  # [:-1]
  paste0("predict2_000", x)
}

# FFQ column mapping (shared by figures 2 and 4)
FFQ_COLS <- c("total_sugar_g","protein_g","fat_g","carbohydrate_g","starch_g",
              "glucose_g","galactose_g","fructose_g","sucrose_g","maltose_g",
              "lactose_g","AOAC_g","sat_fa_g","trans_g","cholesterol_mg",
              "mono_fa_g","poly_fa_g","energy_kcal")
FFQ_RENAME <- c(
  "username"       = "predict2_id",
  "total_sugar_g"  = "X.Total.sugars..g.",
  "protein_g"      = "X.Total.protein..g.",
  "carbohydrate_g" = "Carbohydrate..g.",
  "fat_g"          = "Total.fat..g.",
  "starch_g"       = "X.Starch..g.",
  "glucose_g"      = "X.Glucose..g.",
  "galactose_g"    = "X.Galactose..g.",
  "fructose_g"     = "X.Fructose..g.",
  "sucrose_g"      = "X.Sucrose..g.",
  "maltose_g"      = "X.Maltose..g.",
  "lactose_g"      = "X.Lactose..g.",
  "AOAC_g"         = "X.Total.dietary.fiber..g.",
  "sat_fa_g"       = "Total.saturated.fatty.acids..g.",
  "trans_g"        = "X.Total.trans.fatty.acitds..g.",
  "cholesterol_mg" = "Cholesterol..mg.",
  "mono_fa_g"      = "Total.monounsaturated.fatty.acids..g.",
  "poly_fa_g"      = "Total.polyunsaturated.fatty.acids..g.",
  "energy_kcal"    = "Energy..kcal."
)
load_ffq <- function(path) {
  readxl::read_excel(path) |>
    rename(any_of(FFQ_RENAME)) |>
    select(username, any_of(FFQ_COLS))
}

# =============================================================================
# FIGURE 1 — postprandial response spaghetti plots
# =============================================================================
build_figure1 <- function() {
  D <- file.path(DATA, "figure1")

  person_df    <- read_csv(file.path(D, "person.csv"),       show_col_types = FALSE)
  bio_df       <- read_csv(file.path(D, "bio.csv"),          show_col_types = FALSE)
  base_dataset <- read_csv(file.path(D, "base_dataset.csv"), show_col_types = FALSE)

  bio_wide <- bio_df |>
    filter(biomarker == "CPEP") |>
    pivot_wider(id_cols = username, names_from = biomarker,
                values_from = result_value, values_fn = first)

  df <- base_dataset |>
    left_join(person_df, by = "username") |>
    left_join(bio_wide,  by = "username") |>
    mutate(c_peptide_0 = coalesce(CPEP * CPEP_CONV, c_peptide_0)) |>
    filter(!is.na(person_md_sex), !is.na(person_clinic_bmi))

  cat("Fig1 — N participants:", n_distinct(df$username), "\n")

  make_spaghetti <- function(df_long, x_var, y_var, x_lab, y_lab, title_n) {
    mean_med <- df_long |>
      group_by({{ x_var }}) |>
      summarise(mean_val   = mean({{ y_var }}, na.rm = TRUE),
                median_val = median({{ y_var }}, na.rm = TRUE),
                .groups = "drop")

    ggplot() +
      geom_line(data = df_long,
                aes(x = {{ x_var }}, y = {{ y_var }}, group = username),
                colour = "grey60", alpha = 0.35, linewidth = 0.5) +
      geom_point(data = df_long,
                 aes(x = {{ x_var }}, y = {{ y_var }}, group = username),
                 colour = "black", alpha = 0.5, size = 0.8) +
      geom_line(data = mean_med,
                aes(x = {{ x_var }}, y = mean_val, colour = "Mean"), linewidth = 1) +
      geom_point(data = mean_med,
                 aes(x = {{ x_var }}, y = mean_val, colour = "Mean"), size = 2.5) +
      geom_line(data = mean_med,
                aes(x = {{ x_var }}, y = median_val, colour = "Median"), linewidth = 1) +
      geom_point(data = mean_med,
                 aes(x = {{ x_var }}, y = median_val, colour = "Median"), size = 2.5) +
      scale_colour_manual(values = c(Mean = "royalblue", Median = "red2"),
                          name = paste0("n = ", title_n, "   ")) +
      labs(x = x_lab, y = y_lab) +
      theme_bw(base_size = 16) +
      theme(legend.position  = "top",
            legend.direction = "horizontal",
            legend.title     = element_text(size = 18),
            legend.text      = element_text(size = 18),
            legend.key.size  = unit(1.4, "lines"),
            axis.title       = element_text(size = 17),
            axis.text        = element_text(size = 14),
            panel.grid       = element_blank())
  }

  make_panel <- function(cols, times, y_lab) {
    long <- df |>
      select(username, all_of(cols)) |>
      drop_na() |>
      pivot_longer(-username, names_to = "timepoint", values_to = "value") |>
      mutate(time_min = times[match(timepoint, cols)])
    make_spaghetti(long, time_min, value,
                   "Time (minutes)", y_lab, n_distinct(long$username)) +
      scale_x_continuous(breaks = times)
  }

  p1a <- make_panel(c("c_peptide_0", "c_peptide_120"), c(0, 120), "C-Peptide (pmol/L)")
  p1b <- make_panel(c("trig_0", "trig_360"),           c(0, 360), "TG (mmol/L)")
  p1c <- make_panel(c("glucose_0","glucose_15","glucose_30","glucose_60",
                      "glucose_120","glucose_180","glucose_240","glucose_360"),
                    c(0, 15, 30, 60, 120, 180, 240, 360), "Glucose (mmol/L)")

  combined <- arrangeGrob(p1c, p1b, p1a, ncol = 3)  # Glucose, TG, C-Pep
  ggsave(file.path(OUT, "Fig1_postprand_combined.pdf"), plot = combined, width = 18, height = 7)
  cat("Figure 1 saved.\n")
}

# =============================================================================
# FIGURE 2 — bootstrapped adjusted-R² by feature domain
# =============================================================================
# Build the analysis dataset and bootstrapped adjusted-R² results.
compute_fig2_rsq <- function() {
  set.seed(42)
  D  <- file.path(DATA, "figure2")
  D4 <- file.path(DATA, "figure4")   # FFQ xlsx lives here

  raw_meals    <- read_csv(file.path(D, "raw_meals.csv"),     show_col_types = FALSE)
  context_data <- read_csv(file.path(D, "outputs_meals.csv"), show_col_types = FALSE)
  base_dataset <- read_csv(file.path(D, "base_dataset.csv"),  show_col_types = FALSE)
  ffq_df       <- load_ffq(file.path(D4, "DIET_DHQ3_851.xlsx"))

  # ── microbiome raw abundances → per-taxon features (PCA is applied later,
  #    inside boot_rsq, matching Python figure2.ipynb) ──────────────────────
  mb_raw <- read_tsv(file.path(D, "microbiome PCs bray.tsv"), show_col_types = FALSE)
  mb_mat <- mb_raw |>
    column_to_rownames("clade_name") |>
    select(-NCBI_tax_id) |>
    t() |>
    as.data.frame()

  mb_usernames_all <- parse_mb_username(rownames(mb_mat))

  # mean-impute, then drop zero-variance taxa (matches SimpleImputer + StandardScaler)
  mb_numeric <- as.matrix(mb_mat)
  storage.mode(mb_numeric) <- "numeric"
  col_means_mb <- colMeans(mb_numeric, na.rm = TRUE)
  for (j in seq_len(ncol(mb_numeric))) {
    nas <- is.na(mb_numeric[, j])
    if (any(nas)) {
      fill <- if (is.na(col_means_mb[j]) || is.nan(col_means_mb[j])) 0 else col_means_mb[j]
      mb_numeric[nas, j] <- fill
    }
  }
  col_sds_mb <- apply(mb_numeric, 2, sd)
  mb_numeric <- mb_numeric[, col_sds_mb > 0, drop = FALSE]

  mb_pca_df <- as.data.frame(mb_numeric)
  colnames(mb_pca_df) <- paste0("mb_taxa_", seq_len(ncol(mb_pca_df)))
  mb_pca_df$username  <- mb_usernames_all

  # ── build analysis dataset ─────────────────────────────────────────────────
  df <- base_dataset |>
    left_join(raw_meals, by = "username") |>
    mutate(
      c_peptide_0          = coalesce(CPEP * CPEP_CONV, c_peptide_0),
      person_md_sex_cat    = as.integer(as.factor(person_md_sex)) - 1L,
      person_cpep_120_rise = c_peptide_120 - c_peptide_0,
      person_trig_360_rise = trig_360 - trig_0,
      meal_fat2carbs       = meal_fat / meal_carbohydrate
    )

  ctx_cols <- c("meal_baseline_difference","previous_meal_carbohydrate",
                "previous_meal_time_since","meal_carb_3hr_bef","meal_carb_6hr_bef",
                "meal_carb_9hr_bef","meal_carb_12hr_bef")
  ctx_need <- intersect(ctx_cols, names(context_data))
  if (length(ctx_need) > 0) {
    ctx_sub <- context_data |>
      select(meal_username, meal_set_meal, meal_start_at, any_of(ctx_need)) |>
      rename(username = meal_username)
    df <- df |>
      left_join(ctx_sub, by = c("username","meal_set_meal","meal_start_at"),
                suffix = c("", ".ctx"))
  }

  df <- df |>
    left_join(mb_pca_df, by = "username") |>
    left_join(ffq_df,    by = "username")

  # ── feature groups ─────────────────────────────────────────────────────────
  feature_groups <- list(
    Age               = "person_md_age",
    Sex               = "person_md_sex_cat",
    Anthropometry     = c("person_clinic_weight","person_clinic_height","person_clinic_bmi",
                          "person_clinic_waist","person_clinic_hip","person_clinic_waist2hip"),
    Lipids            = intersect(c("CHOL","HDL","LDL","Non-HDL","THR"), names(df)),
    Cellular          = intersect(c("Lymphocytes","MCH","MCV","Monocytes","Neutrophils",
                                    "PCV","PLT","RBC","RDW","WBC","Basophils","Eosinophils",
                                    "ALT","HB"), names(df)),
    Serum             = intersect(c("PGLU","CPEP","HbA1c"), names(df)),
    Microbiome        = grep("^mb_taxa_", names(df), value = TRUE),
    `Meal Content`    = intersect(c("meal_fat2carbs","meal_calories","meal_carbohydrate",
                                    "meal_fat","meal_protein","meal_fibre","meal_sugar"), names(df)),
    `Meal Context`    = intersect(ctx_cols, names(df)),
    FFQ               = intersect(FFQ_COLS, names(df)),
    `Fasting TG`      = "trig_0",
    `Fasting Glucose` = "glucose_0",
    `Fasting C-pep`   = "c_peptide_0"
  )

  outcomes <- c(
    "meal_iauc"            = "GLU 2h iAUC",
    "person_trig_360_rise" = "TG 6h Rise",
    "person_cpep_120_rise" = "C-PEP 2h Rise"
  )

  # PCA n_components per group (matching Python sub_pipeline); NULL = scale only
  pca_k <- list(
    Age = NULL, Sex = NULL, Anthropometry = 2L, Lipids = 2L, Cellular = 2L,
    Serum = 2L, Microbiome = 40L, `Meal Content` = 5L, `Meal Context` = 5L,
    FFQ = 5L, `Fasting TG` = NULL, `Fasting Glucose` = NULL, `Fasting C-pep` = NULL
  )

  # ── bootstrapped adjusted-R² (matches Python boostrapped_rsquared) ──────────
  #   impute → scale → optional fixed-k PCA on full data, then bootstrap OLS.
  #   point estimate = mean(boot adj R²); CI = 2.5/97.5 quantiles;
  #   p-value = one-sample t-test of boot values vs 0.
  boot_rsq <- function(y, X, n_boot = N_BOOT, n_pca = NULL) {
    valid_y <- !is.na(y)
    y <- y[valid_y]; X <- X[valid_y, , drop = FALSE]
    n <- nrow(X)
    if (n < 5 || ncol(X) == 0)
      return(c(rsq = NA, lo = NA, hi = NA, n = n, pval = NA))

    for (j in seq_len(ncol(X))) {
      nas <- is.na(X[, j])
      if (any(nas) && !all(nas)) X[nas, j] <- mean(X[!nas, j])
    }
    all_na <- apply(X, 2, function(col) all(is.na(col)))
    if (any(all_na)) X <- X[, !all_na, drop = FALSE]
    if (ncol(X) == 0) return(c(rsq = NA, lo = NA, hi = NA, n = n, pval = NA))

    X <- scale(X)
    X[!is.finite(X)] <- 0
    nonzero_cols <- apply(X, 2, function(x) any(x != 0))
    X <- X[, nonzero_cols, drop = FALSE]
    if (ncol(X) == 0) return(c(rsq = NA, lo = NA, hi = NA, n = n, pval = NA))

    if (!is.null(n_pca) && ncol(X) > 1) {
      k <- min(n_pca, ncol(X), n - 1L)
      X <- prcomp(X, center = FALSE, scale. = FALSE)$x[, 1:k, drop = FALSE]
    }

    boot_vals <- replicate(n_boot, {
      idx <- sample(n, replace = TRUE)
      summary(lm(y[idx] ~ X[idx, , drop = FALSE]))$adj.r.squared * 100
    })

    rsq_mean <- mean(boot_vals)
    pval <- if (sd(boot_vals) > 0)
              tryCatch(t.test(boot_vals, mu = 0)$p.value, error = function(e) NA)
            else if (rsq_mean > 0) 0 else 1

    c(rsq  = rsq_mean,
      lo   = as.numeric(quantile(boot_vals, 0.025)),
      hi   = as.numeric(quantile(boot_vals, 0.975)),
      n    = n,
      pval = pval)
  }

  # ── Python-matched data selection for biological groups ────────────────────
  df_bio_glu  <- df |> filter(meal_set_meal == "Meal JJ")
  df_bio_tg   <- df |> filter(meal_set_meal == "Meal JJ") |>
                   distinct(person_trig_360_rise, username, .keep_all = TRUE)
  df_bio_cpep <- df |> filter(meal_set_meal == "Meal JJ") |>
                   distinct(person_cpep_120_rise, username, .keep_all = TRUE)

  results <- map_dfr(names(outcomes), function(out_col) {
    out_name <- outcomes[[out_col]]
    map_dfr(names(feature_groups), function(grp) {
      feats <- intersect(feature_groups[[grp]], names(df))
      if (length(feats) == 0) return(NULL)

      if (grp %in% c("Meal Content", "Meal Context")) {
        if (out_col == "meal_iauc") {
          data_sub <- df
        } else if (grp == "Meal Content") {
          data_sub <- df |>
            filter(meal_set_meal %in% c("Meal JJ", "Meal HH")) |>
            distinct(.data[[out_col]], .keep_all = TRUE)
        } else {
          data_sub <- df |>
            filter(meal_set_meal == "Meal JJ") |>
            distinct(.data[[out_col]], .keep_all = TRUE)
        }
      } else {
        data_sub <- switch(out_col,
          meal_iauc            = df_bio_glu,
          person_trig_360_rise = df_bio_tg,
          person_cpep_120_rise = df_bio_cpep
        )
      }

      data_sub <- data_sub |> filter(!is.na(.data[[out_col]]))

      # microbiome: drop mean-imputed rows so full-data adj R² stays inside the CI
      if (grp == "Microbiome" && length(feats) > 0)
        data_sub <- data_sub |> filter(!is.na(.data[[feats[1]]]))

      y <- data_sub[[out_col]]
      X <- as.matrix(data_sub[, feats, drop = FALSE])

      res <- boot_rsq(y, X, n_pca = pca_k[[grp]])
      # [[ ]] drops the element names so the columns are plain (unnamed) numerics
      tibble(outcome = out_name, group = grp,
             rsq = res[["rsq"]], lo = res[["lo"]], hi = res[["hi"]],
             n = res[["n"]], pval = res[["pval"]])
    })
  })

  # clamp negative adj R² to 0 (penalty artifact — no negative variance explained)
  results <- results |> mutate(rsq = pmax(rsq, 0), lo = pmax(lo, 0))

  cat("\nR² results:\n"); print(results, n = 40)
  list(results = results, df = df)
}

# Figure 2 plot + supporting statistics (ICC, CV, Levene's test)
build_figure2 <- function() {
  fit     <- compute_fig2_rsq()
  results <- fit$results
  df      <- fit$df

  # colors matched to PREDICT 1 (Berry et al. 2020 Nat Med)
  domain_colors <- c(
    "Age" = "darkolivegreen", "Sex" = "blue", "Anthropometry" = "deeppink",
    "Lipids" = "darkorange", "Cellular" = "steelblue", "Serum" = "deepskyblue",
    "Microbiome" = "lightpink", "Meal Content" = "gray30", "Meal Context" = "purple3",
    "FFQ" = "yellow2", "Fasting TG" = "#499894", "Fasting Glucose" = "#499894",
    "Fasting C-pep" = "#499894"
  )

  # ── plot ───────────────────────────────────────────────────────────────────
  plot_rsq <- function(df_res, outcome_label, exclude_groups = NULL,
                       title_label = outcome_label) {
    d <- df_res |>
      filter(outcome == outcome_label, !is.na(rsq)) |>
      (\(x) if (!is.null(exclude_groups)) filter(x, !group %in% exclude_groups) else x)() |>
      arrange(rsq) |>
      mutate(group = fct_inorder(group))

    ggplot(d, aes(x = rsq, y = group, fill = group)) +
      geom_col(show.legend = FALSE) +
      geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.3, orientation = "y") +
      scale_fill_manual(values = domain_colors) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.06))) +
      labs(title = title_label, x = "Adjusted R² (%)", y = NULL) +
      theme_bw(base_size = 13) +
      theme(plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 16),
            axis.text.x = element_text(size = 13),
            axis.title.x = element_text(size = 15),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank())
  }

  p_iauc <- plot_rsq(results, "GLU 2h iAUC",
                     exclude_groups = c("Fasting TG", "Fasting C-pep"),
                     title_label = "Glucose 2h iAUC")
  p_trig <- plot_rsq(results, "TG 6h Rise",
                     exclude_groups = c("Fasting Glucose", "Fasting C-pep"))
  p_cpep <- plot_rsq(results, "C-PEP 2h Rise",
                     exclude_groups = c("Fasting TG", "Fasting Glucose"),
                     title_label = "C-Peptide 2h Rise")

  combined <- arrangeGrob(p_iauc, p_trig, p_cpep, ncol = 3)
  ggsave(file.path(OUT, "Fig2_rsq_combined.pdf"), plot = combined, width = 18, height = 6)

  # ── ICC for repeat Meal JJ measurements ────────────────────────────────────
  mcm_df <- df |>
    filter(meal_set_meal == "Meal JJ") |>
    arrange(username, meal_start_at) |>
    group_by(username) |>
    filter(n() >= 2) |>
    slice(1:2) |>
    mutate(rep = c("first","second")) |>
    ungroup()

  make_icc_data <- function(var) {
    mcm_df |>
      select(username, rep, val = all_of(var)) |>
      pivot_wider(names_from = rep, values_from = val) |>
      drop_na()
  }
  icc_result <- function(var, label) {
    d <- make_icc_data(var)
    if (nrow(d) < 5) { cat(label, ": too few repeats\n"); return(invisible(NULL)) }
    x_bar   <- mean(c(d$first, d$second))
    s2      <- (sum((d$first - x_bar)^2) + sum((d$second - x_bar)^2)) / (2 * nrow(d))
    icc_val <- mean((d$first - x_bar) * (d$second - x_bar)) / (nrow(d) * s2 / nrow(d)) * 100
    boot_icc <- replicate(1000, {
      idx <- sample(nrow(d), replace = TRUE); bd <- d[idx, ]
      xb  <- mean(c(bd$first, bd$second))
      s2b <- (sum((bd$first - xb)^2) + sum((bd$second - xb)^2)) / (2 * nrow(bd))
      mean((bd$first - xb) * (bd$second - xb)) / (nrow(bd) * s2b / nrow(bd)) * 100
    })
    cat(sprintf("ICC %s: %.1f%% [%.1f%%, %.1f%%], n=%d\n", label, icc_val,
                quantile(boot_icc, 0.025), quantile(boot_icc, 0.975), nrow(d)))
  }
  icc_result("meal_iauc",            "GLU 2h iAUC")
  icc_result("person_trig_360_rise", "TG 6h Rise")
  icc_result("person_cpep_120_rise", "C-PEP 2h Rise")

  # ── inter-individual variation: CV and Levene's test ───────────────────────
  person_df  <- read_csv(file.path(DATA, "figure1", "person.csv"),    show_col_types = FALSE)
  biomarkers <- read_csv(file.path(DATA, "figure2", "biomarkers.csv"), show_col_types = FALSE)
  cpep_wide  <- biomarkers |>
    filter(biomarker == "CPEP") |>
    pivot_wider(id_cols = username, names_from = biomarker,
                values_from = result_value, values_fn = first)
  person_df <- person_df |>
    left_join(cpep_wide, by = "username") |>
    mutate(c_peptide_0 = coalesce(CPEP * CPEP_CONV, c_peptide_0),
           trig_rise   = trig_360 - trig_0,
           cpep_rise   = c_peptide_120 - c_peptide_0)

  cv_pct <- function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
  levene_test <- function(groups) {
    all_dev  <- unlist(lapply(groups, function(x) abs(x - median(x, na.rm = TRUE))))
    group_id <- unlist(lapply(seq_along(groups), function(i) rep(i, length(groups[[i]]))))
    res      <- summary(aov(all_dev ~ factor(group_id)))[[1]]
    list(F = res$`F value`[1], p = res$`Pr(>F)`[1])
  }

  cat("\n=== Inter-Individual Variation ===\n")
  cat(sprintf("CV fasting glucose:        %.1f%%  (n=%d)\n",
              cv_pct(person_df$glucose_0),   sum(!is.na(person_df$glucose_0))))
  cat(sprintf("CV fasting TG:             %.1f%%  (n=%d)\n",
              cv_pct(person_df$trig_0),      sum(!is.na(person_df$trig_0))))
  cat(sprintf("CV fasting C-peptide:      %.1f%%  (n=%d)\n",
              cv_pct(person_df$c_peptide_0), sum(!is.na(person_df$c_peptide_0))))
  cat(sprintf("CV glucose iAUC 0-2h:      %.1f%%  (n=%d)\n",
              cv_pct(df$meal_iauc),          sum(!is.na(df$meal_iauc))))
  cat(sprintf("CV TG 6h-rise:             %.1f%%  (n=%d)\n",
              cv_pct(person_df$trig_rise),   sum(!is.na(person_df$trig_rise))))
  cat(sprintf("CV C-peptide 2h-rise:      %.1f%%  (n=%d)\n",
              cv_pct(person_df$cpep_rise),   sum(!is.na(person_df$cpep_rise))))

  glu_cols <- c("glucose_0","glucose_15","glucose_30","glucose_60",
                "glucose_120","glucose_180","glucose_240","glucose_360")
  lev_tg   <- levene_test(list(na.omit(person_df$trig_0), na.omit(person_df$trig_360)))
  lev_glu  <- levene_test(lapply(glu_cols, function(col) na.omit(person_df[[col]])))
  lev_cpep <- levene_test(list(na.omit(person_df$c_peptide_0), na.omit(person_df$c_peptide_120)))
  cat(sprintf("Levene TG (across timepoints):       F=%.2f, p=%.3e\n", lev_tg$F,   lev_tg$p))
  cat(sprintf("Levene glucose (across timepoints):  F=%.2f, p=%.3e\n", lev_glu$F,  lev_glu$p))
  cat(sprintf("Levene C-peptide (across timepoints):F=%.2f, p=%.3e\n", lev_cpep$F, lev_cpep$p))

  cat("Figure 2 saved.\n")
}

# =============================================================================
# FIGURE 3 — Pearson correlation heatmap
# =============================================================================
build_figure3 <- function() {
  D <- file.path(DATA, "figure3")

  biometrics   <- read_csv(file.path(D, "biometrics.csv"),   show_col_types = FALSE)
  biomark      <- read_csv(file.path(D, "biomark.csv"),      show_col_types = FALSE)
  base_dataset <- read_csv(file.path(D, "base_dataset.csv"), show_col_types = FALSE)

  bio_wide <- biomark |>
    filter(!is.na(result_value)) |>
    pivot_wider(id_cols = username, names_from = biomarker,
                values_from = result_value, values_fn = first)

  df <- base_dataset |>
    left_join(biometrics, by = "username") |>
    left_join(bio_wide,   by = "username") |>
    mutate(c_peptide_0 = coalesce(CPEP * CPEP_CONV, c_peptide_0),
           glucose_0   = as.numeric(glucose_0))

  # participants with exactly 4 Meal HH/JJ records (2 HH + 2 JJ)
  df_meals    <- df |> filter(meal_set_meal %in% c("Meal HH", "Meal JJ"))
  valid_users <- df_meals |>
    group_by(username) |>
    summarise(n = n(), .groups = "drop") |>
    filter(n == 4) |>
    pull(username)
  df_filtered <- df_meals |> filter(username %in% valid_users)
  cat("Fig3 — N participants (4 Meal HH/JJ records):", n_distinct(df_filtered$username), "\n")

  corr_df <- df_filtered |>
    select(glucose_0, meal_iauc, trig_0, trig_360, c_peptide_0, c_peptide_120) |>
    drop_na() |>
    distinct() |>
    mutate(trig_rise = trig_360 - trig_0,
           cpep_rise = c_peptide_120 - c_peptide_0) |>
    select(glucose_0, meal_iauc, trig_0, trig_rise, c_peptide_0, cpep_rise)

  corr_labels <- c(
    "glucose_0"   = "Fasting\nGlucose",
    "meal_iauc"   = "Glucose\niAUC 0-2h",
    "trig_0"      = "Fasting\nTG",
    "trig_rise"   = "TG\n6h-rise",
    "c_peptide_0" = "Fasting\nC-Peptide",
    "cpep_rise"   = "C-Peptide\n2h-rise"
  )
  cat("Fig3 — N observations for correlation:", nrow(corr_df), "\n")

  corr_df_labeled <- corr_df
  colnames(corr_df_labeled) <- corr_labels[colnames(corr_df)]

  cor_mat  <- cor(corr_df_labeled, use = "pairwise.complete.obs", method = "pearson")
  cor_pmat <- ggcorrplot::cor_pmat(corr_df_labeled)

  p3 <- ggcorrplot(
    cor_mat, hc.order = FALSE, type = "lower", lab = TRUE, lab_size = 4,
    p.mat = cor_pmat, sig.level = 0.05, insig = "blank",
    colors = c("#2166ac", "white", "#d6604d"), outline.col = "white",
    ggtheme = theme_bw(base_size = 14)
  ) +
    labs(title = "Pearson correlations between postprandial markers",
         caption = paste0("n = ", nrow(corr_df), " observations")) +
    theme(plot.title  = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(OUT, "Fig3_correlation_heatmap.pdf"), plot = p3, width = 8, height = 7)
  cat("Figure 3 saved.\n")
}

# =============================================================================
# FIGURE 4 — ML predicted-vs-measured scatter plots
# =============================================================================
build_figure4 <- function() {
  set.seed(42)
  D <- file.path(DATA, "figure4")

  biometrics   <- read_csv(file.path(D, "biometrics.csv"),    show_col_types = FALSE)
  context_data <- read_csv(file.path(D, "outputs_meals.csv"), show_col_types = FALSE)
  base_dataset <- read_csv(file.path(D, "base_dataset.csv"),  show_col_types = FALSE)
  ffq_df       <- load_ffq(file.path(D, "DIET_DHQ3_851.xlsx"))

  # ── microbiome PCA (10 components) ─────────────────────────────────────────
  mb_raw <- read_tsv(file.path(D, "microbiome PCs bray.tsv"), show_col_types = FALSE)
  mb_mat <- mb_raw |>
    column_to_rownames("clade_name") |>
    select(-NCBI_tax_id) |>
    t() |>
    as.data.frame()
  mb_mat$username <- parse_mb_username(rownames(mb_mat))

  mb_numeric <- mb_mat |> select(username, where(is.numeric))
  mb_feat    <- mb_numeric |> select(-username)
  mb_feat    <- mb_feat[, apply(mb_feat, 2, function(x) !all(is.na(x))), drop = FALSE]
  complete_rows <- rowSums(is.na(mb_feat)) == 0
  mb_pca     <- prcomp(mb_feat[complete_rows, ], center = TRUE, scale. = TRUE)
  mb_pca_df  <- as.data.frame(mb_pca$x[, 1:min(10, ncol(mb_pca$x))])
  colnames(mb_pca_df) <- paste0("mb_pca_", seq_len(ncol(mb_pca_df)))
  mb_pca_df$username  <- mb_numeric$username[complete_rows]

  # ── meal context ───────────────────────────────────────────────────────────
  MEAL_CONTEXT <- c("meal_baseline_difference","previous_meal_carbohydrate",
                    "previous_meal_time_since","meal_carb_3hr_bef","meal_carb_6hr_bef",
                    "meal_carb_9hr_bef","meal_carb_12hr_bef")
  ctx_need <- intersect(MEAL_CONTEXT, names(context_data))
  ctx_sub  <- context_data |>
    select(meal_username, meal_set_meal, meal_start_at, any_of(ctx_need)) |>
    rename(username = meal_username)

  # ── build dataset ──────────────────────────────────────────────────────────
  df <- base_dataset |>
    left_join(biometrics, by = "username") |>
    left_join(mb_pca_df,  by = "username") |>
    left_join(ffq_df,     by = "username") |>
    left_join(ctx_sub, by = c("username","meal_set_meal","meal_start_at"),
              suffix = c("", ".ctx")) |>
    mutate(
      c_peptide_0_pmol = coalesce(as.numeric(person_cpep) * CPEP_CONV, as.numeric(c_peptide_0)),
      trig_rise_6h     = as.numeric(trig_360) - as.numeric(trig_0),
      cpep_rise_2h     = as.numeric(c_peptide_120) - c_peptide_0_pmol,
      meal_fat2carbs   = as.numeric(meal_fat) / as.numeric(meal_carbohydrate),
      person_md_sex    = as.integer(as.factor(person_md_sex)) - 1L,
      person_md_age    = as.numeric(person_md_age),
      meal_iauc        = as.numeric(meal_iauc)
    )

  ANTHROPOMETRY <- c("person_clinic_weight","person_clinic_height","person_clinic_bmi",
                     "person_clinic_waist","person_clinic_hip","person_clinic_waist2hip")
  CLINICAL_BIO  <- c("person_clinic_pulse_1","person_clinic_pulse_2","person_clinic_pulse_3",
                     "person_clinic_pulse_1_diastolic","person_clinic_pulse_2_diastolic",
                     "person_clinic_pulse_3_diastolic","person_clinic_pulse_1_systolic",
                     "person_clinic_pulse_2_systolic","person_clinic_pulse_3_systolic",
                     "person_affinity_il6","person_affinity_alt","person_affinity_chol",
                     "person_affinity_hdl","person_affinity_ldl","person_affinity_non_hdl",
                     "person_affinity_thr","person_affinity_tsh","person_pglu","person_trig",
                     "person_ins","person_cpep","person_hba1c","person_basophils",
                     "person_eosinophils","person_hb","person_lymphocytes","person_mch",
                     "person_mcv","person_monocytes","person_neutrophils","person_nrbc",
                     "person_pcv","person_plt","person_rbc","person_rdw","person_wbc")
  MB_FEATURES  <- grep("^mb_pca_", names(df), value = TRUE)
  MEAL_CONTENT <- c("meal_fat2carbs","meal_calories","meal_carbohydrate",
                    "meal_fat","meal_protein","meal_fibre","meal_sugar")

  features_all <- c(ANTHROPOMETRY, CLINICAL_BIO, "person_md_age",
                    MB_FEATURES, MEAL_CONTEXT, MEAL_CONTENT, FFQ_COLS, "person_md_sex")
  features_all <- intersect(features_all, names(df))
  cat("Fig4 — features available:", length(features_all), "\n")

  # ── helpers ────────────────────────────────────────────────────────────────
  pearson_r   <- function(actual, predicted) cor(actual, predicted, use = "complete.obs")
  impute_mean <- function(mat) {
    col_means <- colMeans(mat, na.rm = TRUE)
    for (j in seq_len(ncol(mat))) mat[is.na(mat[, j]), j] <- col_means[j]
    mat
  }
  prepare_X <- function(data, feats) {
    X <- as.matrix(data[, intersect(feats, names(data)), drop = FALSE])
    mode(X) <- "numeric"
    X[!is.finite(X)] <- NA
    impute_mean(X)
  }

  kfold_cv <- function(df_train, feats, target, n_folds = 5, model_type = "rf") {
    df_sub <- df_train |> select(all_of(c(feats, target))) |> drop_na(all_of(target))
    df_sub <- df_sub[is.finite(df_sub[[target]]), ]
    n <- nrow(df_sub)
    fold_ids <- sample(rep(1:n_folds, length.out = n))

    all_y <- numeric(0); all_yp <- numeric(0); rs <- numeric(n_folds)
    for (k in seq_len(n_folds)) {
      tr <- df_sub[fold_ids != k, ]; te <- df_sub[fold_ids == k, ]
      X_tr <- prepare_X(tr, feats); y_tr <- tr[[target]]
      X_te <- prepare_X(te, feats); y_te <- te[[target]]

      if (model_type == "rf") {
        m  <- ranger(y ~ ., data = as.data.frame(cbind(y = y_tr, X_tr)),
                     num.trees = 500, num.threads = 1, seed = 42)
        yp <- predict(m, data = as.data.frame(X_te))$predictions
      } else {
        m  <- xgb.train(params = list(objective = "reg:squarederror",
                                      max_depth = 2, lambda = 10, eta = 0.1, nthread = 1),
                        data = xgb.DMatrix(X_tr, label = y_tr), nrounds = 100, verbose = 0)
        yp <- predict(m, xgb.DMatrix(X_te))
      }
      rs[k]  <- pearson_r(y_te, yp)
      all_y  <- c(all_y, y_te); all_yp <- c(all_yp, yp)
    }
    list(scores = rs, y = all_y, yp = all_yp)
  }

  plot_scatter <- function(y_meas, y_pred, n, title_str, x_lab, y_lab, lims = NULL) {
    df_p <- tibble(measured = y_meas, predicted = y_pred)
    r    <- pearson_r(y_meas, y_pred)
    pv   <- cor.test(y_meas, y_pred)$p.value
    if (is.null(lims)) lims <- range(c(y_meas, y_pred), na.rm = TRUE)
    ggplot(df_p, aes(x = measured, y = predicted)) +
      geom_point(alpha = 0.5, size = 1.5) +
      geom_smooth(method = "lm", se = FALSE, colour = "#2166ac") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey50") +
      coord_cartesian(xlim = lims, ylim = lims) +
      annotate("text", x = lims[1] + 0.05 * diff(lims), y = lims[2] - 0.05 * diff(lims),
               label = sprintf("R = %.2f\np = %.2e\nn = %d", r, pv, n),
               hjust = 0, vjust = 1, size = 5) +
      labs(title = title_str, x = x_lab, y = y_lab) +
      theme_bw(base_size = 16) +
      theme(aspect.ratio = 1)
  }

  # ── Model 1: TG 6h Rise (Random Forest) ────────────────────────────────────
  target1 <- "trig_rise_6h"
  df1 <- df |>
    filter(meal_set_meal == "Meal JJ") |>
    distinct(username, .keep_all = TRUE) |>
    select(username, all_of(features_all), all_of(target1)) |>
    filter(!is.na(.data[[target1]]), is.finite(.data[[target1]]))

  val_ids1 <- df1 |> distinct(username) |> slice_sample(prop = 0.2) |> pull(username)
  df1_tr   <- df1 |> filter(!username %in% val_ids1)
  df1_val  <- df1 |> filter(username %in% val_ids1)
  cat(sprintf("TG: train=%d, val=%d\n", nrow(df1_tr), nrow(df1_val)))

  cv1 <- kfold_cv(df1_tr, features_all, target1, model_type = "rf")
  cat(sprintf("TG CV R (mean ± sd): %.3f ± %.3f\n", mean(cv1$scores), sd(cv1$scores)))

  X1_tr <- prepare_X(df1_tr, features_all);  y1_tr  <- df1_tr[[target1]]
  X1_val<- prepare_X(df1_val, features_all); y1_val <- df1_val[[target1]]
  m1 <- ranger(y ~ ., data = as.data.frame(cbind(y = y1_tr, X1_tr)),
               num.trees = 1000, num.threads = 1, seed = 42)
  y1_pred_val <- predict(m1, data = as.data.frame(X1_val))$predictions
  cat(sprintf("TG Val R: %.3f\n", pearson_r(y1_val, y1_pred_val)))

  lims1 <- range(c(cv1$y, cv1$yp, y1_val, y1_pred_val), na.rm = TRUE)
  p4a_train <- plot_scatter(cv1$y, cv1$yp, nrow(df1_tr), "TG 6h Rise – Train",
                            "Measured TG Rise", "Predicted TG Rise", lims = lims1)
  p4a_val   <- plot_scatter(y1_val, y1_pred_val, nrow(df1_val), "TG 6h Rise – Test",
                            "Measured TG Rise", "Predicted TG Rise", lims = lims1)

  # ── Model 2: C-PEP 2h Rise (XGBoost) ───────────────────────────────────────
  target2 <- "cpep_rise_2h"
  df2 <- df |>
    filter(meal_set_meal == "Meal JJ") |>
    distinct(username, .keep_all = TRUE) |>
    select(username, all_of(features_all), all_of(target2)) |>
    filter(!is.na(.data[[target2]]), is.finite(.data[[target2]]))

  val_ids2 <- df2 |> distinct(username) |> slice_sample(prop = 0.2) |> pull(username)
  df2_tr   <- df2 |> filter(!username %in% val_ids2)
  df2_val  <- df2 |> filter(username %in% val_ids2)
  cat(sprintf("CPEP: train=%d, val=%d\n", nrow(df2_tr), nrow(df2_val)))

  cv2 <- kfold_cv(df2_tr, features_all, target2, model_type = "xgb")
  cat(sprintf("CPEP CV R (mean ± sd): %.3f ± %.3f\n", mean(cv2$scores), sd(cv2$scores)))

  X2_tr  <- prepare_X(df2_tr, features_all);  y2_tr  <- df2_tr[[target2]]
  X2_val <- prepare_X(df2_val, features_all); y2_val <- df2_val[[target2]]
  m2 <- xgb.train(params = list(objective = "reg:squarederror",
                                max_depth = 2, lambda = 10, eta = 0.1, nthread = 1),
                  data = xgb.DMatrix(X2_tr, label = y2_tr), nrounds = 100, verbose = 0)
  y2_pred_val <- predict(m2, xgb.DMatrix(X2_val))
  cat(sprintf("CPEP Val R: %.3f\n", pearson_r(y2_val, y2_pred_val)))

  lims2 <- range(c(cv2$y, cv2$yp, y2_val, y2_pred_val), na.rm = TRUE)
  p4b_train <- plot_scatter(cv2$y, cv2$yp, nrow(df2_tr), "C-Peptide 2h Rise – Train",
                            "Measured C-Peptide Rise", "Predicted C-Peptide Rise", lims = lims2)
  p4b_val   <- plot_scatter(y2_val, y2_pred_val, nrow(df2_val), "C-Peptide 2h Rise – Test",
                            "Measured C-Peptide Rise", "Predicted C-Peptide Rise", lims = lims2)

  # ── Model 3: Glucose iAUC (Random Forest) ──────────────────────────────────
  target3 <- "meal_iauc"
  df_cal <- df |>
    group_by(username) |>
    summarise(cal_iauc = mean(meal_iauc[meal_set_meal == "Meal CC"], na.rm = TRUE),
              .groups = "drop")

  df3 <- df |>
    left_join(df_cal, by = "username") |>
    select(username, all_of(features_all), cal_iauc, all_of(target3)) |>
    filter(!is.na(.data[[target3]]), is.finite(.data[[target3]]))

  feats3 <- intersect(c(features_all, "cal_iauc"), names(df3))

  val_ids3 <- sample(seq_len(nrow(df3)), round(nrow(df3) * 0.2))
  val_mask <- seq_len(nrow(df3)) %in% val_ids3
  df3_tr   <- df3[!val_mask, ]; df3_val <- df3[val_mask, ]
  cat(sprintf("iAUC: train=%d, val=%d\n", nrow(df3_tr), nrow(df3_val)))

  cv3 <- kfold_cv(df3_tr, feats3, target3, model_type = "rf")
  cat(sprintf("iAUC CV R (mean ± sd): %.3f ± %.3f\n", mean(cv3$scores), sd(cv3$scores)))

  X3_tr  <- prepare_X(df3_tr, feats3);  y3_tr  <- df3_tr[[target3]]
  X3_val <- prepare_X(df3_val, feats3); y3_val <- df3_val[[target3]]
  m3 <- ranger(y ~ ., data = as.data.frame(cbind(y = y3_tr, X3_tr)),
               num.trees = 1000, num.threads = 1, seed = 42)
  y3_pred_val <- predict(m3, data = as.data.frame(X3_val))$predictions
  cat(sprintf("iAUC Val R: %.3f\n", pearson_r(y3_val, y3_pred_val)))

  lims3 <- range(c(cv3$y, cv3$yp, y3_val, y3_pred_val), na.rm = TRUE)
  p4c_train <- plot_scatter(cv3$y, cv3$yp, nrow(df3_tr), "Glucose iAUC – Train",
                            "Measured iAUC", "Predicted iAUC", lims = lims3)
  p4c_val   <- plot_scatter(y3_val, y3_pred_val, nrow(df3_val), "Glucose iAUC – Test",
                            "Measured iAUC", "Predicted iAUC", lims = lims3)

  all_fig4 <- arrangeGrob(p4c_train, p4c_val,
                          p4a_train, p4a_val,
                          p4b_train, p4b_val, nrow = 3, ncol = 2)
  ggsave(file.path(OUT, "Fig4_ML_combined.pdf"), plot = all_fig4, width = 10, height = 15)
  cat("Figure 4 saved.\n")
}

# =============================================================================
# FIGURE 5 — ASCVD ROC curves
# =============================================================================
build_figure5 <- function() {
  D  <- file.path(DATA, "figure5")
  D2 <- file.path(DATA, "figure2")   # shared base_dataset

  meals_df     <- read_csv(file.path(D,  "meals.csv"),        show_col_types = FALSE)
  biomark      <- read_csv(file.path(D,  "biomark.csv"),      show_col_types = FALSE)
  base_dataset <- read_csv(file.path(D2, "base_dataset.csv"), show_col_types = FALSE)

  bio_wide <- biomark |>
    filter(!is.na(result_value)) |>
    pivot_wider(id_cols = username, names_from = biomarker,
                values_from = result_value, values_fn = first)

  df <- base_dataset |>
    left_join(meals_df, by = "username") |>
    left_join(bio_wide, by = "username") |>
    mutate(c_peptide_0 = coalesce(CPEP * CPEP_CONV, c_peptide_0),
           trig_rise   = trig_360 - trig_0,
           cpep_rise   = c_peptide_120 - c_peptide_0)
  cat("Fig5 — N participants:", n_distinct(df$username), "\n")

  # ASCVD binary outcome: top 30th percentile = high risk
  person_ascvd <- df |>
    select(username, person_ascvd_risk) |>
    distinct() |>
    filter(!is.na(person_ascvd_risk))
  q70 <- quantile(person_ascvd$person_ascvd_risk, 0.70, na.rm = TRUE)
  person_ascvd <- person_ascvd |> mutate(ascvd_binary = as.integer(person_ascvd_risk >= q70))
  df <- df |> left_join(person_ascvd |> select(username, ascvd_binary), by = "username")
  cat("Fig5 — ASCVD binary (top 30%): n high =", sum(df$ascvd_binary == 1, na.rm = TRUE),
      ", n low =", sum(df$ascvd_binary == 0, na.rm = TRUE), "\n")

  make_roc_plot <- function(df_sub, fasting_var, pp_var, title, combine_fn = NULL) {
    d_fast <- df_sub |> select(username, all_of(c(fasting_var, "ascvd_binary"))) |>
      drop_na() |> distinct()
    d_pp   <- df_sub |> select(username, all_of(c(pp_var, "ascvd_binary"))) |>
      drop_na() |> distinct()
    d_both <- df_sub |> select(username, all_of(c(fasting_var, pp_var, "ascvd_binary"))) |>
      drop_na() |> distinct()

    roc_fast <- roc(d_fast$ascvd_binary, d_fast[[fasting_var]], quiet = TRUE)
    roc_pp   <- roc(d_pp$ascvd_binary,   d_pp[[pp_var]],        quiet = TRUE)
    combined_score <- if (is.null(combine_fn))
      d_both[[fasting_var]] + d_both[[pp_var]]           # raw sum (TG, C-peptide)
    else
      combine_fn(d_both[[fasting_var]], d_both[[pp_var]])
    roc_both <- roc(d_both$ascvd_binary, combined_score, quiet = TRUE)

    df_roc <- bind_rows(
      tibble(fpr = 1 - roc_fast$specificities, tpr = roc_fast$sensitivities,
             label = sprintf("Fasting (AUC=%.2f)", auc(roc_fast))),
      tibble(fpr = 1 - roc_pp$specificities, tpr = roc_pp$sensitivities,
             label = sprintf("Post-prandial (AUC=%.2f)", auc(roc_pp))),
      tibble(fpr = 1 - roc_both$specificities, tpr = roc_both$sensitivities,
             label = sprintf("Combined (AUC=%.2f)", auc(roc_both)))
    )

    ggplot(df_roc, aes(x = fpr, y = tpr, colour = label)) +
      geom_line(linewidth = 0.8) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
      scale_colour_brewer(palette = "Set1", name = NULL) +
      labs(title = title, x = "1 − Specificity", y = "Sensitivity") +
      theme_bw(base_size = 16) +
      theme(legend.position   = c(0.62, 0.18),
            legend.text       = element_text(size = 16),
            legend.key.size   = unit(1.5, "lines"),
            legend.margin     = margin(8, 10, 8, 10),
            legend.background = element_rect(fill = alpha("white", 0.9), colour = "grey80"))
  }

  # Glucose combined: z(fasting)/2 + log(z(postprandial))/2, inf/nan → 0
  glu_combine <- function(fasting, pp) {
    z_fast   <- (fasting - mean(fasting, na.rm = TRUE)) / sd(fasting, na.rm = TRUE)
    z_pp     <- (pp - mean(pp, na.rm = TRUE)) / sd(pp, na.rm = TRUE)
    log_z_pp <- log(z_pp)
    log_z_pp[!is.finite(log_z_pp)] <- 0
    z_fast / 2 + log_z_pp / 2
  }

  p5a <- make_roc_plot(df, "glucose_0",   "meal_iauc", "ASCVD – Glucose",       glu_combine)
  p5b <- make_roc_plot(df, "trig_0",      "trig_rise", "ASCVD – Triglycerides")
  p5c <- make_roc_plot(df, "c_peptide_0", "cpep_rise", "ASCVD – C-Peptide")

  combined <- arrangeGrob(p5a, p5b, p5c, nrow = 3)
  ggsave(file.path(OUT, "Fig5_ROC_combined.pdf"), plot = combined, width = 6, height = 18)
  cat("Figure 5 saved.\n")
}

# =============================================================================
# run all — only when executed directly (skipped when sourced, e.g. by the
# supplementary-tables script, which reuses compute_fig2_rsq() and the helpers)
# =============================================================================
if (sys.nframe() == 0L) {
  # Route incidental grid drawing to a null device so no stray Rplots.pdf is
  # created under Rscript (ggsave still opens/closes its own device on top).
  if (!interactive()) grDevices::pdf(NULL)

  build_figure1()
  build_figure2()
  build_figure3()
  build_figure4()
  build_figure5()

  if (!interactive() && !is.null(grDevices::dev.list())) invisible(grDevices::dev.off())
  cat("\nAll figures written to:", OUT, "\n")
}
