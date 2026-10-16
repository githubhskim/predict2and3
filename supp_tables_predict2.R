# =============================================================================
# PREDICT2 — Supplementary tables
#
# Writes supp_tables_predict2.xlsx (next to this script) with three sheets:
#   - S Table Baseline    — baseline participant characteristics
#   - S Table Variables   — feature list per domain (one column per group)
#   - S Table Regression  — bootstrapped adjusted-R² per domain and outcome
#
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(writexl)
})

# ── source the figures script for shared paths/helpers + compute_fig2_rsq() ──
this_dir <- {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg)))
  else if (interactive() && requireNamespace("rstudioapi", quietly = TRUE) &&
           rstudioapi::isAvailable())
    dirname(rstudioapi::getActiveDocumentContext()$path)
  else getwd()
}
source(file.path(this_dir, "predict2_figures.R"))   # defines DATA, OUT, CPEP_CONV, compute_fig2_rsq()

# =============================================================================
# S Table Baseline — participant characteristics
# =============================================================================
D  <- file.path(DATA, "figure2")
D1 <- file.path(DATA, "figure1")
D4 <- file.path(DATA, "figure4")

raw_meals    <- read_csv(file.path(D,  "raw_meals.csv"),    show_col_types = FALSE)
biomarkers   <- read_csv(file.path(D,  "biomarkers.csv"),   show_col_types = FALSE)
base_dataset <- read_csv(file.path(D,  "base_dataset.csv"), show_col_types = FALSE)
person_fig1  <- read_csv(file.path(D1, "person.csv"),       show_col_types = FALSE)
bio4         <- read_csv(file.path(D4, "biometrics.csv"),   show_col_types = FALSE)

# Sex from figure1/person.csv (complete for all base_dataset users)
sex_df <- person_fig1 |>
  group_by(username) |>
  summarise(person_md_sex = first(person_md_sex), .groups = "drop")

# Demographics from figure4/biometrics.csv; pulse/BP = mean of measurements 2 & 3
demo_df <- bio4 |>
  group_by(username) |>
  summarise(
    person_md_age        = first(na.omit(person_md_age)),
    person_clinic_weight = first(na.omit(person_clinic_weight)),
    person_clinic_height = first(na.omit(person_clinic_height)),
    person_clinic_bmi    = first(na.omit(person_clinic_bmi)),
    person_clinic_waist  = first(na.omit(person_clinic_waist)),
    person_clinic_hip    = first(na.omit(person_clinic_hip)),
    waist_hip_ratio      = first(na.omit(waist_hip_ratio)),
    pulse                = mean(c(first(na.omit(person_clinic_pulse_2)),
                                  first(na.omit(person_clinic_pulse_3))), na.rm = TRUE),
    diastolic            = mean(c(first(na.omit(person_clinic_pulse_2_diastolic)),
                                  first(na.omit(person_clinic_pulse_3_diastolic))), na.rm = TRUE),
    systolic             = mean(c(first(na.omit(person_clinic_pulse_2_systolic)),
                                  first(na.omit(person_clinic_pulse_3_systolic))), na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(across(c(pulse, diastolic, systolic), ~ ifelse(is.nan(.), NA, .)))

# Biomarkers wide (Quest lab values)
bio_wide <- biomarkers |>
  filter(!is.na(result_value)) |>
  pivot_wider(id_cols = username, names_from = biomarker,
              values_from = result_value, values_fn = first)

# Postprandial outcomes (sources mirror predict2_figures.R for consistency)
iauc_outcomes <- raw_meals |>
  group_by(username) |>
  summarise(meal_iauc = mean(meal_iauc, na.rm = TRUE), .groups = "drop") |>
  mutate(meal_iauc = ifelse(is.nan(meal_iauc), NA, meal_iauc))

cpep_wide <- biomarkers |>
  filter(!is.na(result_value), biomarker == "CPEP") |>
  pivot_wider(id_cols = username, names_from = biomarker,
              values_from = result_value, values_fn = first)

bio_outcomes <- person_fig1 |>
  left_join(cpep_wide, by = "username") |>
  mutate(c_peptide_0 = coalesce(CPEP * CPEP_CONV, c_peptide_0),
         trig_rise   = trig_360 - trig_0,
         cpep_rise   = c_peptide_120 - c_peptide_0) |>
  select(username, trig_rise, cpep_rise)

df_base <- base_dataset |>
  left_join(sex_df,        by = "username") |>
  left_join(demo_df,       by = "username") |>
  left_join(bio_wide,      by = "username") |>
  left_join(iauc_outcomes, by = "username") |>
  left_join(bio_outcomes,  by = "username")

# Summary-row helpers
smry <- function(x, label, varname) {
  tibble(
    Variable                    = label,
    Note                        = "",
    `Variable name in datafile` = varname,
    n        = sum(!is.na(x)),
    Mean     = round(mean(x,   na.rm = TRUE), 2),
    SD       = round(sd(x,     na.rm = TRUE), 2),
    Median   = round(median(x, na.rm = TRUE), 2),
    `25th %` = round(quantile(x, 0.25, na.rm = TRUE), 2),
    `75th %` = round(quantile(x, 0.75, na.rm = TRUE), 2)
  )
}
missing_row <- function(label, varname = "") {
  tibble(
    Variable = label, Note = "missing", `Variable name in datafile` = varname,
    n = NA_integer_, Mean = NA_real_, SD = NA_real_, Median = NA_real_,
    `25th %` = NA_real_, `75th %` = NA_real_
  )
}

sex_n  <- sum(!is.na(df_base$person_md_sex))
sex_fm <- paste0(sum(df_base$person_md_sex == "F", na.rm = TRUE), "/",
                 sum(df_base$person_md_sex == "M", na.rm = TRUE))

t_baseline <- bind_rows(
  tibble(
    Variable = "Sex (F/M)", Note = "", `Variable name in datafile` = "person_md_sex",
    n = sex_n, Mean = NA_real_, SD = NA_real_, Median = NA_real_,
    `25th %` = NA_real_, `75th %` = NA_real_, `F/M` = sex_fm
  ),
  smry(df_base$person_md_age,        "Age (years)",              "person_md_age"),
  smry(df_base$person_clinic_weight, "Weight (kg)",              "person_clinic_weight"),
  smry(df_base$person_clinic_height, "Height (cm)",              "person_clinic_height"),
  smry(df_base$person_clinic_bmi,    "BMI (kg/m²)",              "person_clinic_bmi"),
  smry(df_base$person_clinic_waist,  "Waist circumference (cm)", "person_clinic_waist"),
  smry(df_base$person_clinic_hip,    "Hip circumference (cm)",   "person_clinic_hip"),
  smry(df_base$waist_hip_ratio,      "Waist/hip ratio",          "waist_hip_ratio"),
  missing_row("Pulse (bpm)"),
  missing_row("Diastolic BP (mmHg)"),
  missing_row("Systolic BP (mmHg)"),
  missing_row("IL-6 (ng/L)"),
  smry(df_base$ALT,          "Alanine aminotransferase ALT (U/L)",   "ALT_IU/L_Wholeblood_Hospital"),
  smry(df_base$CHOL,         "Total cholesterol (mmol/L)",           "CHOL_mmol/L_Wholeblood_Hospital"),
  smry(df_base$HDL,          "HDL (mmol/L)",                         "HDL_mmol/L_Wholeblood_Hospital"),
  smry(df_base$LDL,          "LDL (mmol/L)",                         "LDL_mmol/L_Wholeblood_Hospital"),
  smry(df_base[["Non-HDL"]], "non-HDL (mmol/L)",                     "Non-HDL_mmol/L_Wholeblood_Hospital"),
  smry(df_base$THR,          "Total cholesterol to HDL ratio (THR)", "THR__Wholeblood_Hospital"),
  missing_row("Thyroid stimulating hormone TSH (mU/L)"),
  smry(df_base$PGLU,         "Plasma glucose (mmol/L)",              "PGLU_mmol/L_Wholeblood_Hospital"),
  smry(df_base$TRIG,         "Triglycerides (mmol/L)",               "TRIG"),
  missing_row("Insulin (mU/L)"),
  smry(df_base$CPEP,         "C-peptide (ng/mL)",                    "CPEP_ug/L_Wholeblood_Hospital"),
  smry(df_base$HbA1c,        "HbA1c (%)",                            "HbA1c_%_Wholeblood_Hospital"),
  smry(df_base$Basophils,    "Basophils (×10⁹/L)",                   "Basophils_x10^9_Wholeblood_Hospital"),
  smry(df_base$Eosinophils,  "Eosinophils (×10⁹/L)",                 "Eosinophils_x10^9_Wholeblood_Hospital"),
  smry(df_base$HB,           "Hemoglobin (g/L)",                     "HB_g/L_Wholeblood_Hospital"),
  smry(df_base$Lymphocytes,  "Lymphocytes (×10⁹/L)",                 "Lymphocytes_x10^9_Wholeblood_Hospital"),
  smry(df_base$MCH,          "Mean Corpuscular Haemoglobin MCH (pg)","MCH_pg_Wholeblood_Hospital"),
  missing_row("Mean Corpuscular Volume MCV (fL)"),
  smry(df_base$Monocytes,    "Monocytes (×10⁹/L)",                   "Monocytes_x10^9_Wholeblood_Hospital"),
  smry(df_base$Neutrophils,  "Neutrophils (×10⁹/L)",                 "Neutrophils_x10^9_Wholeblood_Hospital"),
  smry(df_base$RBC,          "Nucleated Red Blood Cells NRBC (×10⁹/L)", "RBC_x10^12_Wholeblood_Hospital"),
  smry(df_base$PCV,          "Packed Cell Volume PCV (L/L)",         "PCV_L/L_Wholeblood_Hospital"),
  smry(df_base$PLT,          "Platelets PLT (×10⁹/L)",               "PLT_x10^9_Wholeblood_Hospital"),
  smry(df_base$RBC,          "Red Blood Cells RBC (×10¹²/L)",        "RBC_x10^12_Wholeblood_Hospital"),
  smry(df_base$RDW,          "Red Cell Distribution Width RDW (%)",  "RDW__Wholeblood_Hospital"),
  smry(df_base$WBC,          "White blood cells WBC (×10⁹/L)",       "WBC_x10^9_Wholeblood_Hospital"),
  smry(df_base$meal_iauc,    "Glucose iAUC 0-2h",                    "MEAL_iAUC"),
  smry(df_base$trig_rise,    "Triglycerides 6h-rise (mmol/L)",       "TRIG_RISE"),
  smry(df_base$cpep_rise,    "C-peptide 2h-rise (pmol/L)",           "CPEP_RISE")
)

# =============================================================================
# S Table Variables — feature list, one column per domain
# =============================================================================
var_lists <- list(
  "Demographic" = c("Age (years)", "Sex (M/F)"),
  "Anthropometric" = c("BMI (kg/m²)", "Height (cm)", "Weight (kg)",
                       "Waist circumference (cm)", "Hip circumference (cm)", "Waist/hip ratio"),
  "Heart Rate & Blood Pressure" = c("(not collected in PREDICT2 US)"),
  "Serum lipid markers" = c("Total cholesterol (mmol/L)", "HDL cholesterol (mmol/L)",
                            "LDL cholesterol (mmol/L)", "non-HDL cholesterol (mmol/L)",
                            "Total cholesterol/HDL ratio (THR)"),
  "Other serum biochemistry" = c("ALT (U/L)", "Basophils (×10⁹/L)", "Eosinophils (×10⁹/L)",
                                 "Haemoglobin (g/L)", "Lymphocytes (×10⁹/L)", "MCH (pg)",
                                 "MCV (fL)", "Monocytes (×10⁹/L)", "Neutrophils (×10⁹/L)",
                                 "PCV (L/L)", "Platelets (×10⁹/L)", "RBC (×10¹²/L)",
                                 "RDW (%)", "WBC (×10⁹/L)"),
  "Serum glycemic markers" = c("Plasma glucose (mmol/L)", "C-peptide (ng/mL)", "HbA1c (%)"),
  "Meal Composition" = c("Fat/carbohydrate ratio", "Total calories (kcal)", "Carbohydrate (g)",
                         "Fat (g)", "Protein (g)", "Fibre (g)", "Sugar (g)"),
  "Context (activity)" = c("(not available in PREDICT2 US)"),
  "Context (sleep)" = c("(not available in PREDICT2 US)"),
  "Context (timing)" = c("Baseline glycaemic difference", "Previous meal carbohydrate (g)",
                         "Time since previous meal (min)", "Carbohydrate 3h before meal (g)",
                         "Carbohydrate 6h before meal (g)", "Carbohydrate 9h before meal (g)",
                         "Carbohydrate 12h before meal (g)"),
  "FFQ" = c("Energy (kcal)", "Protein (g)", "Fat (g)", "Carbohydrate (g)", "Sugar (g)",
            "Starch (g)", "Dietary fiber (g)", "Glucose (g)", "Fructose (g)", "Galactose (g)",
            "Sucrose (g)", "Maltose (g)", "Lactose (g)", "Saturated FA (g)", "Trans FA (g)",
            "Monounsaturated FA (g)", "Polyunsaturated FA (g)", "Cholesterol (mg)")
)

max_rows  <- max(lengths(var_lists))
t_vars <- as.data.frame(
  lapply(var_lists, function(v) c(v, rep(NA_character_, max_rows - length(v)))),
  stringsAsFactors = FALSE, check.names = FALSE
)

# =============================================================================
# S Table Regression — bootstrapped adjusted-R² per domain and outcome
# =============================================================================
results <- compute_fig2_rsq()$results

# Map internal R group names → display names
group_display <- c(
  "Age" = "Age", "Sex" = "Sex", "Anthropometry" = "Anthropometry",
  "Lipids" = "Serum lipids markers", "Cellular" = "Other serum markers",
  "Serum" = "Serum glycemic markers", "Microbiome" = "Microbiome",
  "Meal Content" = "Meal Composition", "Meal Context" = "Meal Context", "FFQ" = "FFQ",
  "Fasting TG" = "Fasting triglycerides", "Fasting Glucose" = "Fasting Glucose",
  "Fasting C-pep" = "Fasting Cpeptide"
)

# Row builders (numeric columns only — significance-flag columns removed)
data_row <- function(meals_lbl, grp, rsq, hi, lo, pval, dash = FALSE) {
  if (dash) {
    return(tibble(`Meals Analysed` = "", Group = grp, `Adjusted R²` = "-",
                  `95% CI lower` = "-", `95% CI upper` = "-", `p-value` = "-"))
  }
  tibble(
    `Meals Analysed` = meals_lbl,
    Group            = grp,
    `Adjusted R²`    = ifelse(is.na(rsq),  "-", as.character(round(rsq, 2))),
    `95% CI lower`   = ifelse(is.na(lo),   "-", as.character(round(lo,  2))),
    `95% CI upper`   = ifelse(is.na(hi),   "-", as.character(round(hi,  2))),
    `p-value`        = ifelse(is.na(pval), "-", as.character(signif(pval, 3)))
  )
}
blank_row <- tibble(`Meals Analysed` = "", Group = "", `Adjusted R²` = "",
                    `95% CI lower` = "", `95% CI upper` = "", `p-value` = "")
title_row <- function(txt) mutate(blank_row, `Meals Analysed` = txt)

make_reg_block <- function(out_name) {
  res <- results |>
    filter(outcome == out_name) |>
    mutate(Group = recode(group, !!!group_display))

  rv <- function(disp_grp) {
    r <- res |> filter(Group == disp_grp)
    if (nrow(r) == 0) return(c(rsq = NA, hi = NA, lo = NA, pval = NA))
    c(rsq = r$rsq[1], hi = r$hi[1], lo = r$lo[1], pval = r$pval[1])
  }

  n_meal1  <- res |> filter(group == "Age")          |> pull(n) |> first()
  n_meal12 <- res |> filter(group == "Meal Content") |> pull(n) |> first()
  label_meal1  <- paste0("Meal 1 (m = ", n_meal1, " meals)")
  label_meal12 <- if (out_name == "GLU 2h iAUC")
    paste0("all meals (m = ", n_meal12, " meals)")
  else
    paste0("Meal 1, 2 (m = ", n_meal12, " meals)")

  fspec <- list(
    "TG 6h Rise"    = list(grp = "Fasting triglycerides", insert_after = "Serum lipids markers"),
    "GLU 2h iAUC"   = list(grp = "Fasting Glucose",       insert_after = "Serum glycemic markers"),
    "C-PEP 2h Rise" = list(grp = "Fasting Cpeptide",      insert_after = "Serum glycemic markers")
  )[[out_name]]

  meal1_seq <- c("Sex", "Age", "Anthropometry", "Blood pressure",
                 "Serum lipids markers", "Serum glycemic markers",
                 "Other serum markers", "Microbiome", "FFQ")

  rows1  <- list()
  ml_lbl <- label_meal1
  for (g in meal1_seq) {
    if (g == "Blood pressure") {
      rows1[[length(rows1) + 1L]] <- data_row("", g, NA, NA, NA, NA, dash = TRUE)
    } else {
      r <- rv(g)
      rows1[[length(rows1) + 1L]] <- data_row(ml_lbl, g, r["rsq"], r["hi"], r["lo"], r["pval"])
      ml_lbl <- ""
    }
    if (g == fspec$insert_after) {
      fr <- rv(fspec$grp)
      rows1[[length(rows1) + 1L]] <- data_row("", fspec$grp, fr["rsq"], fr["hi"], fr["lo"], fr["pval"])
    }
  }

  rows2   <- list()
  ml_lbl2 <- label_meal12
  for (g in c("Meal Composition", "Meal Context")) {
    r <- rv(g)
    rows2[[length(rows2) + 1L]] <- data_row(ml_lbl2, g, r["rsq"], r["hi"], r["lo"], r["pval"])
    ml_lbl2 <- ""
  }

  bind_rows(c(rows1, rows2))
}

outcome_display <- c(
  "TG 6h Rise"    = "Triglycerides 6h-rise",
  "GLU 2h iAUC"   = "Glucose 2h iAUC",
  "C-PEP 2h Rise" = "C-peptide 2h-rise"
)

t_reg <- bind_rows(lapply(names(outcome_display), function(out) {
  bind_rows(title_row(outcome_display[[out]]), make_reg_block(out), blank_row)
}))

# =============================================================================
# save
# =============================================================================
write_xlsx(
  list(
    `S Table Baseline`   = t_baseline,
    `S Table Variables`  = t_vars,
    `S Table Regression` = t_reg
  ),
  file.path(OUT, "supp_tables_predict2.xlsx")
)

cat("\nSupplementary tables saved to:", file.path(OUT, "supp_tables_predict2.xlsx"), "\n")
