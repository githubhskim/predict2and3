# ============================================================
# PREDICT3 CGM Exposome — Common Setup
# ------------------------------------------------------------
# Shared data prep used by every figure script:
#   1. Load predictor + raw CGM data (Cycle 1, PREDICT 3.0 UK/US)
#   2. Compute per-person CGM curve-shape metrics
#   3. Build the analysis dataset (one row per Meal 55 participant)
#   4. Recode ordinal / categorical variables
#   5. Add grouping variables used by the deep-dive figures
#   6. Define predictor groups, publication themes, and palettes
#
# Each fig*.R script sources this once (guarded), then draws its figure.
# Outcome: Meal 55 incremental AUC (iAUC). Cross-sectional, one row per person.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

# Install + load figure/stat helpers if needed
for (pkg in c("ggpubr", "rstatix", "ggtext", "ggh4x", "ggnewscale")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cran.r-project.org")
}
suppressPackageStartupMessages({
  library(ggpubr)
  library(rstatix)
  library(ggtext)
})

set.seed(42)

# ----------------------------------------------------------------
# 1.  LOAD & FILTER DATA
# ----------------------------------------------------------------

# Source CSVs live in this folder's input/ subdirectory (self-contained copies).
# find_data() prefers input/, then falls back to the working directory and its
# parent, so the scripts also run if launched from the project root.
find_data <- function(f) {
  candidates <- c(file.path("input", f), f, file.path("..", f))
  for (p in candidates) if (file.exists(p)) return(p)
  stop("Cannot find data file: ", f,
       "  (looked in input/, the working directory, and its parent)")
}

cat("Loading predictor data...\n")
pred_raw <- read.csv(find_data("1_zoe_harvard_p3.csv"), stringsAsFactors = FALSE)

cat("Loading raw CGM data...\n")
cgm_raw  <- read.csv(find_data("P3_CGM_glucose_meal55.csv"), stringsAsFactors = FALSE)

# Keep only testing cycle 1 and PREDICT 3.0 UK/US studies.
pred <- pred_raw %>%
  filter(testing_cycle == 1,
         study %in% c("PRODUCT UK 1.0 PREDICT 3.0",
                      "PRODUCT US 1.0 PREDICT 3.0"))

cat(sprintf("Cycle-1 PREDICT3 rows: %d  |  Unique users: %d\n",
            nrow(pred), n_distinct(pred$username)))

# ----------------------------------------------------------------
# 2.  COMPUTE PER-PERSON CGM CURVE-SHAPE METRICS
# ----------------------------------------------------------------

# Post-peak dip: minimum (glucose - baseline) AFTER the peak index.
# Negative = glucose dipped BELOW pre-meal baseline (reactive hypoglycemia proxy).
# iAUC cannot capture this because inc_glucose_level is floored at 0.
cgm_post_peak_dip_fn <- function(inc_gl, gl, baseline) {
  pk_idx <- which.max(inc_gl)
  if (length(pk_idx) == 0 || all(is.na(inc_gl))) return(NA_real_)
  post_vals <- (gl - baseline)[pk_idx:length(gl)]
  if (length(post_vals) == 0) return(NA_real_)
  min(post_vals, na.rm = TRUE)
}

# 2-hour incremental glucose: (glucose - baseline) at the timepoint
# closest to 120 min (7200 s from meal start). Negative = fell below baseline.
cgm_2h_inc_fn <- function(gl, baseline, time) {
  idx <- which.min(abs(time - 7200))
  if (length(idx) == 0) return(NA_real_)
  (gl - baseline)[idx[1]]
}

cat("Computing per-person CGM metrics...\n")
cgm_metrics <- cgm_raw %>%
  arrange(username, inc_measured_at) %>%
  group_by(username) %>%
  summarise(
    cgm_fasting_glucose  = first(glucose_baseline),
    cgm_peak_glucose     = max(glucose_level,     na.rm = TRUE),
    cgm_peak_inc         = max(inc_glucose_level, na.rm = TRUE),
    cgm_time_to_peak_min = inc_measured_at[which.max(inc_glucose_level)] / 60,
    cgm_post_peak_dip    = cgm_post_peak_dip_fn(inc_glucose_level,
                                                glucose_level, glucose_baseline),
    cgm_2h_inc           = cgm_2h_inc_fn(glucose_level, glucose_baseline,
                                         inc_measured_at),
    .groups = "drop"
  )

# ----------------------------------------------------------------
# 3.  BUILD ANALYSIS DATASET
# ----------------------------------------------------------------

df <- pred %>%
  left_join(cgm_metrics, by = "username") %>%
  filter(!is.na(meal_55_meal_iauc))   # restrict to Meal 55 participants

cat(sprintf("Analysis N: %d\n", nrow(df)))

# ----------------------------------------------------------------
# 4.  RECODE ORDINAL / CATEGORICAL VARIABLES
# ----------------------------------------------------------------

likert5 <- c("Strongly disagree" = 1, "Disagree" = 2, "Neutral" = 3,
             "Agree" = 4, "Strongly agree" = 5)

severity7 <- c("No discomfort."                = 0, "Minor discomfort."            = 1,
               "Mild discomfort."              = 2, "Moderate discomfort."         = 3,
               "Moderately severe discomfort." = 4, "Severe discomfort."           = 5,
               "Very severe discomfort."       = 6)

interference7 <- c("No interference."                = 0, "Minor interference."          = 1,
                   "Mild interference."              = 2, "Moderate interference."       = 3,
                   "Moderately severe interference." = 4, "Severe interference."         = 5,
                   "Very severe interference."       = 6)

inconvenience7 <- c("No symptoms or inconvenience at all."           = 0,
                    "Minor GI symptoms or inconvenience."             = 1,
                    "Mild GI symptoms or inconvenience."              = 2,
                    "Moderate GI symptoms or inconvenience."          = 3,
                    "Moderately severe GI symptoms or inconvenience." = 4,
                    "Severe GI symptoms or inconvenience."            = 5,
                    "Very severe GI symptoms or inconvenience."       = 6)

sleep_qual_map <- c("Very bad" = 1, "Fairly bad" = 2, "Fairly good" = 3, "Very good" = 4)

sleep_dur_map  <- c("Less than 5" = 4.5, "5-6" = 5.5, "6-7" = 6.5, "7-8" = 7.5,
                    "8-9" = 8.5, "9-10" = 9.5, "10-11" = 10.5, "More than 12" = 12.5)

gut_freq_ord   <- c("Once a week or less"        = 1, "Twice a week"             = 2,
                    "Three or four times a week" = 3, "Five or six times a week" = 4,
                    "Once a day"                 = 5, "Twice a day"              = 6,
                    "Three times a day"          = 7, "Four times a day"         = 8,
                    "Five or more times a day"   = 9)

df <- df %>%
  mutate(
    # Demographics
    sex_male = as.integer(sex == "male"),

    # Metabolic — cap physiologically implausible outliers
    hba1c           = ifelse(self_reported_hba1c > 15 | self_reported_hba1c < 3,
                             NA, self_reported_hba1c),
    fasting_glucose = ifelse(self_reported_pglu > 20 | self_reported_pglu < 2,
                             NA, self_reported_pglu),
    triglycerides   = ifelse(self_reported_trig > 15, NA, self_reported_trig),

    # Blood pressure — cap at physiologically plausible range
    systolic_clean  = ifelse(systolic  < 70 | systolic  > 220, NA, systolic),
    diastolic_clean = ifelse(diastolic < 40 | diastolic > 140, NA, diastolic),

    # Sleep
    sleep_duration = sleep_dur_map[sleep_h_per_day],

    # Hunger / satiety (Agree-Disagree scale)
    hunger_freq    = likert5[hunger_frequency],
    hunger_immed   = likert5[hunger_immediacy],
    hunger_rumble  = likert5[hunger_stomach_rumbling],
    hunger_irritab = likert5[hunger_irritability],
    hunger_lighthd = likert5[hunger_lightheaded],

    # Gut health
    gut_bloat        = severity7[gut_symp_bloated_distended],
    gut_diarrhea     = severity7[gut_symp_diarrhea],
    gut_acidreflux   = severity7[gut_symp_acid_reflux],
    gut_heartburn    = severity7[gut_symp_heartburn],
    gut_loosestools  = severity7[gut_symp_loose_stools],
    gut_daily_interf = interference7[gut_daily_intereference],
    gut_inconvenienc = inconvenience7[gut_inconvenience],
    gut_bristol      = as.integer(sub("Bristol ", "", gut_bristol_stool)),
    gut_freq_num     = gut_freq_ord[gut_frequency],

    # Physical activity — cap 99th pctile for MVPA; cap sitting & resistance hrs
    pa_mvpa_cap    = pmin(pa_mvpa_min,       quantile(pa_mvpa_min,       0.99, na.rm = TRUE)),
    pa_sitting_cap = pmin(pa_sitting_hrs,    quantile(pa_sitting_hrs,    0.99, na.rm = TRUE)),
    pa_resist_cap  = pmin(pa_resistance_hrs, quantile(pa_resistance_hrs, 0.99, na.rm = TRUE)),

    # Lifestyle
    smoking = case_when(
      smoking_status == "Never"                       ~ 0,
      smoking_status == "Not currently, but used to"  ~ 1,
      smoking_status == "Yes"                         ~ 2,
      TRUE ~ NA_real_
    ),

    # Ethnicity — free-text encoded as White (1) vs non-White (0)
    ethnicity_white = {
      eth_lo <- tolower(trimws(ifelse(is.na(ethnicity), "", ethnicity)))
      white_pat <- paste(c(
        "white","caucasian","british","irish","english","scottish","welsh",
        "european","german","french","italian","spanish","portuguese","greek",
        "nordic","scandinavian","dutch","polish","czech","slovak","hungarian",
        "romanian","austrian","swiss","belgian","swedish","norwegian","danish",
        "finnish","australian","new zealand","canadian","american"
      ), collapse = "|")
      dplyr::case_when(
        is.na(ethnicity) | ethnicity == ""     ~ NA_real_,
        stringr::str_detect(eth_lo, white_pat) ~ 1,
        TRUE                                   ~ 0
      )
    },

    # Sleep quality from FFQ (well-powered: N~26K vs N~600 for the survey version)
    sleep_quality_ffq = sleep_qual_map[ffq_answer_sleep_quality],

    # Diet preference — questionnaire, ordered most-to-least plant-based
    diet_pref = case_when(
      diet_preference_questionnaire == "Vegan"       ~ "Vegan",
      diet_preference_questionnaire == "Vegetarian"  ~ "Vegetarian",
      diet_preference_questionnaire == "Pescatarian" ~ "Pescatarian",
      diet_preference_questionnaire == "Flexitarian" ~ "Flexitarian",
      diet_preference_questionnaire == "Mixed diet"  ~ "Mixed diet",
      diet_preference_questionnaire %in%
        c("Other diet pattern","None of the above")  ~ "Other",
      TRUE ~ NA_character_
    ),
    diet_pref = factor(diet_pref,
                       levels = c("Vegan","Vegetarian","Pescatarian",
                                  "Flexitarian","Mixed diet","Other")),
    # App-based diet preference (parallel measure; 4 categories, no Flexitarian)
    diet_app = factor(case_when(
      diet_preference_app == "Vegan"       ~ "Vegan",
      diet_preference_app == "Vegetarian"  ~ "Vegetarian",
      diet_preference_app == "Pescatarian" ~ "Pescatarian",
      diet_preference_app == "Mixed diet"  ~ "Mixed diet",
      TRUE ~ NA_character_),
      levels = c("Vegan","Vegetarian","Pescatarian","Mixed diet")),

    # Menopause status
    menopause = case_when(
      menopause_status == "pre-menopausal"  ~ "Pre",
      menopause_status == "peri-menopausal" ~ "Peri",
      menopause_status == "post-menopausal" ~ "Post",
      TRUE ~ NA_character_
    ),
    menopause = factor(menopause, levels = c("Pre","Peri","Post")),

    # MHT / HRT — ordered treatment status (current=1, never=2, past=3)
    mht_status = case_when(
      menopause_mht_taking == "I am currently taking MHT"    ~ "Currently taking",
      menopause_mht_taking == "I have never taken MHT"       ~ "Never taken",
      menopause_mht_taking == "I have taken MHT in the past" ~ "Taken in past",
      TRUE ~ NA_character_
    ),
    mht_status = factor(mht_status,
                        levels = c("Currently taking","Never taken","Taken in past")),
    mht_num = as.numeric(mht_status),

    # Urbanicity — linked LSOA data (RUC11), collapsed to 4 ordered levels
    urbanicity = dplyr::case_when(
      grepl("Urban major|Urban minor", ruc11nm_lsoa_england_wales) ~ "Urban\n(conurbation)",
      grepl("Urban city",              ruc11nm_lsoa_england_wales) ~ "Urban\n(town)",
      grepl("Rural town",              ruc11nm_lsoa_england_wales) ~ "Rural\ntown",
      grepl("Rural village",           ruc11nm_lsoa_england_wales) ~ "Rural\nvillage",
      TRUE ~ NA_character_
    ),
    urbanicity = factor(urbanicity,
                        levels = c("Urban\n(conurbation)","Urban\n(town)",
                                   "Rural\ntown","Rural\nvillage")),
    urbanicity_num = as.numeric(urbanicity),

    # Shift work (binary)
    shift_work = dplyr::case_when(
      ffq_answer_work_shifts == TRUE  ~ 1,
      ffq_answer_work_shifts == FALSE ~ 0,
      TRUE ~ NA_real_
    ),

    # Convert iAUC from mmol/L·s to mmol/L·min for display (R² is unit-invariant).
    meal_55_meal_iauc_min = meal_55_meal_iauc / 60
  )

# ----------------------------------------------------------------
# 5.  GROUPING VARIABLES FOR THE DEEP-DIVE FIGURES
# ----------------------------------------------------------------

df <- df %>%
  mutate(
    smoking_cat = factor(case_when(
      smoking_status == "Never"                      ~ "Never",
      smoking_status == "Not currently, but used to" ~ "Former",
      smoking_status == "Yes"                        ~ "Current",
      TRUE ~ NA_character_),
      levels = c("Never", "Former", "Current")),

    mvpa_q = {
      q <- quantile(pa_mvpa_cap, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
      cut(pa_mvpa_cap, breaks = c(-Inf, q[1], q[2], q[3], Inf),
          labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"), include.lowest = TRUE)
    },
    sitting_q = {
      q <- quantile(pa_sitting_cap, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
      cut(pa_sitting_cap, breaks = c(-Inf, q[1], q[2], q[3], Inf),
          labels = c("Q1 (Low)", "Q2", "Q3", "Q4 (High)"), include.lowest = TRUE)
    },
    # Resistance training is zero-inflated → group by weekly hours bands
    resist_group = factor(case_when(
      pa_resist_cap == 0 ~ "0 h",
      pa_resist_cap <= 1 ~ "0-1 h",
      pa_resist_cap <= 3 ~ "1-3 h",
      pa_resist_cap >  3 ~ ">3 h",
      TRUE ~ NA_character_),
      levels = c("0 h", "0-1 h", "1-3 h", ">3 h")),

    sleep_cat = factor(ffq_answer_sleep_quality,
                       levels = c("Very bad","Fairly bad","Fairly good","Very good")),

    hunger_group = factor(case_when(
      hunger_freq <= 2 ~ "Low (1-2)",
      hunger_freq == 3 ~ "Mid (3)",
      hunger_freq >= 4 ~ "High (4-5)",
      TRUE ~ NA_character_),
      levels = c("Low (1-2)", "Mid (3)", "High (4-5)")),

    gut_group = factor(case_when(
      gut_inconvenienc <= 1 ~ "None/Minor",
      gut_inconvenienc <= 3 ~ "Mild/Mod",
      gut_inconvenienc >= 4 ~ "Mod/Severe",
      TRUE ~ NA_character_),
      levels = c("None/Minor", "Mild/Mod", "Mod/Severe"))
  )

# ----------------------------------------------------------------
# 6.  PREDICTOR GROUPS
#     Each element: list(vars = c(...), labels = c(...))
# ----------------------------------------------------------------

predictor_groups <- list(

  "Demographics &\nAnthropometric" = list(
    vars   = c("age_at_test","sex_male","bmi","waist_cm","ethnicity_white"),
    labels = c("Age","Sex (male)","BMI","Waist circumference","White ethnicity (vs other)")
  ),

  "Metabolic\nMarkers" = list(
    vars   = c("hba1c","fasting_glucose","triglycerides","systolic_clean","diastolic_clean"),
    labels = c("HbA1c","Fasting glucose","Triglycerides","Systolic BP","Diastolic BP")
  ),

  # Diet-quality predictors come from two sources for the same underlying
  # nutrients: the habitual FFQ and the real-time food diary. Both lists
  # below now share one nutrient stem ("Total <nutrient>") and differ only
  # by an explicit " - FFQ" / " - food diary" suffix, so paired variables
  # (e.g. fiber measured two ways) read as an obvious pair rather than as
  # similarly-spelled-but-distinct variables ("Dietary fiber" vs "Daily
  # fiber", as before).
  "Diet Quality\n(Habitual FFQ)" = list(
    vars   = c("ffq_opdi_score_zoe","ffq_hpdi_score_zoe","ffq_updi_score_zoe",
               "ffq_total_dietary_fiber_g","ffq_total_carbs_g","ffq_protein_g",
               "ffq_total_fat_g","ffq_total_sugars_g","ffq_energy_kcal","ffq_alcohol_g"),
    labels = c("Original Plant-based Diet Index (OPDI) - FFQ",
               "Healthy Plant-based Diet Index (HPDI) - FFQ",
               "Unhealthy Plant-based Diet Index (UPDI) - FFQ",
               "Total fiber - FFQ","Total carbohydrates - FFQ","Total protein - FFQ",
               "Total fat - FFQ","Total sugars - FFQ","Total energy (kcal) - FFQ",
               "Alcohol intake - FFQ")
  ),

  "Real-time\nFood Diary" = list(
    # _zoe versions throughout (N~76K vs N~15K for _trento; r > 0.94).
    vars   = c("food_log_opdi_score_zoe","food_log_hpdi_score_zoe","food_log_updi_score_zoe",
               "food_log_calories","food_log_carbohydrate","food_log_fibre",
               "food_log_fat","food_log_protein","food_log_sugar",
               "food_log_ultra_processed_calories"),
    labels = c("Original Plant-based Diet Index (OPDI) - food diary",
               "Healthy Plant-based Diet Index (HPDI) - food diary",
               "Unhealthy Plant-based Diet Index (UPDI) - food diary",
               "Total energy (kcal) - food diary","Total carbohydrates - food diary",
               "Total fiber - food diary","Total fat - food diary",
               "Total protein - food diary","Total sugars - food diary",
               "Ultra-processed food calories - food diary")
  ),

  "Sleep" = list(
    vars   = c("sleep_quality_ffq","sleep_duration"),
    labels = c("Sleep quality","Sleep duration (h)")
  ),

  "Hunger &\nSatiety" = list(
    vars   = c("hunger_freq","hunger_immed","hunger_rumble","hunger_irritab","hunger_lighthd"),
    labels = c("Hunger frequency","Hunger immediacy","Stomach rumbling when hungry",
               "Irritability when hungry","Lightheadedness when hungry")
  ),

  "Physical\nActivity" = list(
    vars   = c("pa_mvpa_cap","pa_sitting_cap","pa_resist_cap"),
    labels = c("MVPA (min/week)","Daily sitting time (h)","Resistance training (h/wk)")
  ),

  "Gut Health" = list(
    vars   = c("gut_bristol","gut_freq_num","gut_inconvenienc","gut_daily_interf",
               "gut_bloat","gut_diarrhea","gut_acidreflux","gut_heartburn","gut_loosestools"),
    labels = c("Bristol stool type","Bowel frequency","GI inconvenience (overall)",
               "GI daily interference","Bloating","Diarrhea","Acid reflux",
               "Heartburn","Loose stools")
  ),

  "Lifestyle &\nDiet Preference" = list(
    vars   = c("diet_pref","menopause","mht_num","smoking"),
    labels = c("Diet preference (questionnaire)","Menopause status (females)",
               "MHT/HRT taking status (females)","Smoking status")
  ),

  "Social &\nGeography" = list(
    vars   = c("urbanicity_num","shift_work","imd_decile"),
    labels = c("Urbanicity (urban to rural)","Shift work (yes vs no)",
               "Deprivation decile (IMD)")
  )
)

# ----------------------------------------------------------------
# 7.  PUBLICATION THEMES
# ----------------------------------------------------------------

# Journal theme: theme_bw skeleton with full panel border, faint grid, black
# text, no bold — matches the PREDICT2 house style.
theme_paper <- theme_bw(base_size = 11, base_family = "Helvetica") +
  theme(
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    panel.grid.major  = element_line(colour = "grey88", linewidth = 0.25),
    panel.grid.minor  = element_line(colour = "grey94", linewidth = 0.20),
    axis.text         = element_text(colour = "black"),
    axis.title        = element_text(colour = "black"),
    axis.ticks        = element_line(colour = "black", linewidth = 0.30),
    strip.background  = element_blank(),
    strip.text        = element_text(colour = "black"),
    plot.title        = element_text(hjust = 0, face = "plain"),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
    legend.title      = element_text(face = "plain")
  )

# Deep-dive panel theme (Figures 2 & supplementary)
theme_dive <- theme_paper +
  theme(
    axis.text       = element_text(size = 13),
    axis.title      = element_text(size = 14),
    plot.title      = element_text(size = 15, face = "plain", colour = "black"),
    panel.spacing.y = unit(0.1, "lines")
  )

# ----------------------------------------------------------------
# 8.  COLOUR PALETTE — ColorBrewer Set1 / Dark2
# ----------------------------------------------------------------
pal_blue   <- "#377EB8"   # Set1 blue
pal_brick  <- "#E41A1C"   # Set1 red
pal_green  <- "#4DAF4A"   # Set1 green
pal_amber  <- "#FF7F00"   # Set1 orange
pal_lblue  <- "#A6CEE3"   # Paired light blue
pal_teal   <- "#1B9E77"   # Dark2 teal
pal_purple <- "#984EA3"   # Set1 purple
pal_brown  <- "#A65628"   # Set1 brown
pal_plum   <- "#666666"   # Set1 grey
pal_mauve  <- "#7570B3"   # Dark2 slate

# Domain colours (10 domains), keyed by the two-line predictor-group names.
domain_colours <- c(
  "Demographics &\nAnthropometric" = "#E41A1C",  # red
  "Metabolic\nMarkers"             = "#377EB8",  # blue
  "Diet Quality\n(Habitual FFQ)"   = "#4DAF4A",  # green
  "Real-time\nFood Diary"          = "#984EA3",  # purple
  "Sleep"                          = "#FF7F00",  # orange
  "Hunger &\nSatiety"              = "#A65628",  # brown
  "Physical\nActivity"             = "#F781BF",  # pink
  "Gut Health"                     = "#1B9E77",  # teal
  "Lifestyle &\nDiet Preference"   = "#666666",  # dark grey
  "Social &\nGeography"            = "#7570B3"   # slate
)

# Single-line domain labels for legends (no \n breaks).
domain_labels_clean <- c(
  "Demographics &\nAnthropometric" = "Demographics & Anthropometric",
  "Metabolic\nMarkers"             = "Metabolic Markers",
  "Diet Quality\n(Habitual FFQ)"   = "Diet Quality (Habitual FFQ)",
  "Real-time\nFood Diary"          = "Real-time Food Diary",
  "Sleep"                          = "Sleep",
  "Hunger &\nSatiety"              = "Hunger & Satiety",
  "Physical\nActivity"             = "Physical Activity",
  "Gut Health"                     = "Gut Health",
  "Lifestyle &\nDiet Preference"   = "Lifestyle & Diet Preference",
  "Social &\nGeography"            = "Social & Geography"
)

# Minimum N to fit a bivariate model.
MIN_N <- 500

cat("Common setup complete.\n")
