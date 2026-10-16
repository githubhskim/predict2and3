# ============================================================================
# make_supp_tables_clean.R
# Final-review build of the supplementary workbook.
# Keeps only three sheets:
#   1. Predictor_R2  — Figure 1 predictor R² results
#   2. iAUC_by_Group — iAUC summary statistics by group
#   3. MultiOutcome  — Figure 4 multi-outcome heatmap results
# ============================================================================

suppressMessages({
  library(tidyverse)
  library(openxlsx)
})

# Source files live in this folder's input/ subdirectory (self-contained copies):
# the raw predictor CSV plus the two R² tables from the exposome analysis.
# find_data() prefers input/, then falls back to the working directory and its
# parent, so the script also runs if launched from the project root.
# Matches the convention in 00_common_setup.R.
find_data <- function(f) {
  candidates <- c(file.path("input", f), f, file.path("..", f))
  for (p in candidates) if (file.exists(p)) return(p)
  stop("Cannot find data file: ", f,
       "  (looked in input/, the working directory, and its parent)")
}

# ── Load data (same pipeline as main script) ───────────────────────────────
pred_raw <- read.csv(find_data("1_zoe_harvard_p3.csv"), stringsAsFactors = FALSE)
cat("Data loaded\n")

likert5 <- c("Strongly disagree"=1,"Disagree"=2,"Neutral"=3,"Agree"=4,"Strongly agree"=5)
inconvenience7 <- c("No symptoms or inconvenience at all."=0,
                    "Minor GI symptoms or inconvenience."=1,
                    "Mild GI symptoms or inconvenience."=2,
                    "Moderate GI symptoms or inconvenience."=3,
                    "Moderately severe GI symptoms or inconvenience."=4,
                    "Severe GI symptoms or inconvenience."=5,
                    "Very severe GI symptoms or inconvenience."=6)

df <- pred_raw %>%
  filter(testing_cycle == 1,
         study %in% c("PRODUCT UK 1.0 PREDICT 3.0", "PRODUCT US 1.0 PREDICT 3.0")) %>%
  filter(!is.na(meal_55_meal_iauc)) %>%
  mutate(
    hunger_freq      = likert5[hunger_frequency],
    gut_inconvenienc = inconvenience7[gut_inconvenience],
    iauc_min   = meal_55_meal_iauc / 60,
    fem        = sex == "female",
    mht_status = factor(case_when(
      menopause_mht_taking == "I am currently taking MHT"    ~ "Currently taking",
      menopause_mht_taking == "I have taken MHT in the past" ~ "Taken in past",
      menopause_mht_taking == "I have never taken MHT"       ~ "Never taken",
      TRUE ~ NA_character_),
      levels = c("Currently taking", "Never taken", "Taken in past")),
    menopause  = factor(case_when(
      menopause_status == "pre-menopausal"  ~ "Pre",
      menopause_status == "peri-menopausal" ~ "Peri",
      menopause_status == "post-menopausal" ~ "Post",
      TRUE ~ NA_character_), levels = c("Pre","Peri","Post")),
    smoking_cat = factor(case_when(
      smoking_status == "Never"                       ~ "Never",
      smoking_status == "Not currently, but used to" ~ "Former",
      smoking_status == "Yes"                         ~ "Current",
      TRUE ~ NA_character_), levels = c("Never","Former","Current")),
    sleep_cat = factor(ffq_answer_sleep_quality,
                       levels = c("Very bad","Fairly bad","Fairly good","Very good")),
    pa_mvpa_cap    = pmin(pa_mvpa_min,      quantile(pa_mvpa_min,      0.99, na.rm = TRUE)),
    pa_sitting_cap = pmin(pa_sitting_hrs,   quantile(pa_sitting_hrs,   0.99, na.rm = TRUE)),
    pa_resist_cap  = pmin(pa_resistance_hrs, quantile(pa_resistance_hrs, 0.99, na.rm = TRUE)),
    mvpa_q = {
      q <- quantile(pa_mvpa_cap, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
      cut(pa_mvpa_cap, breaks = c(-Inf, q[1], q[2], q[3], Inf),
          labels = c("Q1 (Low)","Q2","Q3","Q4 (High)"), include.lowest = TRUE)
    },
    sitting_q = {
      q <- quantile(pa_sitting_cap, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
      cut(pa_sitting_cap, breaks = c(-Inf, q[1], q[2], q[3], Inf),
          labels = c("Q1 (Low)","Q2","Q3","Q4 (High)"), include.lowest = TRUE)
    },
    resist_group = factor(case_when(
      pa_resist_cap == 0 ~ "0 h",
      pa_resist_cap <= 1 ~ "0-1 h",
      pa_resist_cap <= 3 ~ "1-3 h",
      pa_resist_cap >  3 ~ ">3 h",
      TRUE ~ NA_character_), levels = c("0 h","0-1 h","1-3 h",">3 h")),
    diet_pref = factor(case_when(
      diet_preference_questionnaire == "Vegan"       ~ "Vegan",
      diet_preference_questionnaire == "Vegetarian"  ~ "Vegetarian",
      diet_preference_questionnaire == "Pescatarian" ~ "Pescatarian",
      diet_preference_questionnaire == "Flexitarian" ~ "Flexitarian",
      diet_preference_questionnaire == "Mixed diet"  ~ "Mixed diet",
      TRUE ~ NA_character_),
      levels = c("Vegan","Vegetarian","Pescatarian","Flexitarian","Mixed diet"))
  )
cat("N =", nrow(df), "\n")

# ── Helpers ────────────────────────────────────────────────────────────────
grp_summary <- function(data, group_var, value_var = "iauc_min", panel_label = NULL) {
  data %>%
    filter(!is.na(.data[[group_var]]), !is.na(.data[[value_var]])) %>%
    group_by(Group = .data[[group_var]]) %>%
    summarise(
      N       = n(),
      Mean    = round(mean(.data[[value_var]]),  2),
      SD      = round(sd(.data[[value_var]]),    2),
      Median  = round(median(.data[[value_var]]),2),
      Q1      = round(quantile(.data[[value_var]], 0.25), 2),
      Q3      = round(quantile(.data[[value_var]], 0.75), 2),
      .groups = "drop") %>%
    mutate(Panel = panel_label, .before = 1)
}

# Ordinal line panels summarised on their raw score.
# code_labels: named vector mapping the numeric survey code (names) to the
# readable response text (values). Group is shown as the response text; the
# factor levels preserve ordinal scale order for the rows.
score_summary <- function(data, var, label, code_labels) {
  data %>%
    filter(!is.na(.data[[var]]), !is.na(iauc_min)) %>%
    mutate(Group = factor(as.character(.data[[var]]),
                          levels = names(code_labels),
                          labels = unname(code_labels))) %>%
    group_by(Group) %>%
    summarise(N = n(), Mean = round(mean(iauc_min), 2), SD = round(sd(iauc_min), 2),
              Median = round(median(iauc_min), 2),
              Q1 = round(quantile(iauc_min, 0.25), 2),
              Q3 = round(quantile(iauc_min, 0.75), 2), .groups = "drop") %>%
    arrange(Group) %>%
    mutate(Group = as.character(Group), Panel = label, .before = 1)
}

# Canonical domain order (matches predictor_groups in the main script).
# The analysis script flattens the "\n" in domain labels before writing its
# CSVs, so these levels must match the flattened text.
domain_levels <- c("Demographics & Anthropometric", "Metabolic Markers",
                   "Diet Quality (Habitual FFQ)", "Real-time Food Diary",
                   "Sleep", "Hunger & Satiety", "Physical Activity",
                   "Gut Health", "Lifestyle & Diet Preference", "Social & Geography")
clean_domain <- function(x) gsub("\n", " ", x)

# ── SHEET 1: Predictor R² results (Figure 1) ───────────────────────────────
cat("Building Sheet 1: Predictor R² results\n")
predictor_r2 <- read.csv(find_data("cgm_exposome_r2_table.csv"), stringsAsFactors = FALSE) %>%
  mutate(Domain = factor(Domain, levels = domain_levels)) %>%
  arrange(Domain, desc(Adj.R...)) %>%
  mutate(Domain = clean_domain(as.character(Domain))) %>%
  select(
    Domain       = Domain,
    Predictor    = Predictor,
    `Adj R² (%)` = Adj.R...,
    `P-value`    = P.value,
    `BH-FDR q`   = BH.FDR.q,
    N            = N
  )

# ── SHEET 2: iAUC summary stats by group ───────────────────────────────────
# Panels from the main and supplementary figures are pooled into one sheet:
# variables move between those figures, so splitting the tables by figure goes
# stale. Order below follows main-figure panels first, then supplementary.
cat("Building Sheet 2: iAUC summaries\n")
sum_meno   <- grp_summary(df %>% filter(fem), "menopause",    panel_label = "Menopausal stage (females only)")
sum_mht    <- grp_summary(df %>% filter(fem), "mht_status",   panel_label = "MHT/HRT status (females only)")
sum_diet   <- grp_summary(df,                 "diet_pref",    panel_label = "Diet preference (questionnaire)")
sum_smoke  <- grp_summary(df,                 "smoking_cat",  panel_label = "Smoking status")
sum_resist <- grp_summary(df,                 "resist_group", panel_label = "Resistance training (h/week)")
sum_sleep  <- grp_summary(df,                 "sleep_cat",    panel_label = "Sleep quality")
sum_mvpa   <- grp_summary(df,                 "mvpa_q",       panel_label = "MVPA quartile")
sum_sit    <- grp_summary(df,                 "sitting_q",    panel_label = "Sitting time quartile")

sum_hunger <- score_summary(df, "hunger_freq", "Hunger frequency (\"I feel hungry often\")",
                            code_labels = c("1" = "Strongly disagree", "2" = "Disagree",
                                            "3" = "Neutral", "4" = "Agree",
                                            "5" = "Strongly agree"))
sum_gi     <- score_summary(df, "gut_inconvenienc", "GI inconvenience (symptom severity)",
                            code_labels = c("0" = "No symptoms", "1" = "Minor", "2" = "Mild",
                                            "3" = "Moderate", "4" = "Moderately severe",
                                            "5" = "Severe", "6" = "Very severe"))

iauc_by_group <- bind_rows(
  # main figure
  sum_meno   %>% mutate(Group = as.character(Group)),
  sum_mht    %>% mutate(Group = as.character(Group)),
  sum_diet   %>% mutate(Group = as.character(Group)),
  sum_smoke  %>% mutate(Group = as.character(Group)),
  sum_resist %>% mutate(Group = as.character(Group)),
  # supplementary figure (physical activity first)
  sum_mvpa   %>% mutate(Group = as.character(Group)),
  sum_sit    %>% mutate(Group = as.character(Group)),
  sum_sleep  %>% mutate(Group = as.character(Group)),
  sum_hunger,
  sum_gi
)

# ── SHEET 3: multi-outcome heatmap (Figure 4) ──────────────────────────────
cat("Building Sheet 3: multi-outcome results\n")
multi_outcome <- read.csv(find_data("cgm_multi_outcome_r2_table.csv"), stringsAsFactors = FALSE) %>%
  mutate(Domain = factor(Domain, levels = domain_levels)) %>%
  arrange(CGM_Outcome, Domain, desc(Adj.R2.pct)) %>%
  mutate(Domain = clean_domain(as.character(Domain))) %>%
  select(
    `CGM Outcome`         = CGM_Label,
    Domain                = Domain,
    Predictor             = Predictor,
    # Signed: magnitude is the adjusted R², sign gives the direction of the
    # association (negative = predictor is associated with a lower outcome).
    # Rows stay sorted by |R²| so the strongest effects lead regardless of sign.
    `Signed adj R² (%)`   = Signed.R2.pct,
    `P-value`             = P.value,
    `BH-FDR q`            = P.FDR,
    N                     = N
  )

# ── Write Excel workbook ───────────────────────────────────────────────────
cat("Writing Excel file...\n")
wb <- createWorkbook()

hs <- createStyle(textDecoration = "bold", fgFill = "#D9E1F2",
                  border = "Bottom", borderColour = "#4472C4", wrapText = FALSE)

add_sheet <- function(wb, sheet_name, data, col_widths) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data, headerStyle = hs)
  setColWidths(wb, sheet_name, cols = seq_along(col_widths), widths = col_widths)
  nr <- nrow(data)
  if (nr > 1) {
    even_rows <- seq(3, nr + 1, by = 2)
    addStyle(wb, sheet_name,
             createStyle(fgFill = "#F7F9FD"),
             rows = even_rows, cols = seq_along(col_widths),
             gridExpand = TRUE, stack = TRUE)
  }
}

add_sheet(wb, "Predictor_R2",  predictor_r2,  col_widths = c(28, 45, 12, 14, 12, 10))
add_sheet(wb, "iAUC_by_Group", iauc_by_group, col_widths = c(56, 24, 8, 8, 8, 9, 8, 8))
add_sheet(wb, "MultiOutcome",  multi_outcome, col_widths = c(28, 26, 42, 18, 14, 12, 10))

saveWorkbook(wb, "supplementary_tables_final.xlsx", overwrite = TRUE)
cat("Saved: supplementary_tables_final.xlsx\n")
