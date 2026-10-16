# ============================================================
# Figure 2 — Deep dive: novel PREDICT3 lifestyle variables
# ------------------------------------------------------------
# Five lifestyle predictors that are BOTH significant (q < 0.05) and
# well-powered. Each row: iAUC distribution by group (violin + boxplot with
# pairwise Bonferroni Wilcoxon brackets) | mean glucose trajectory by group.
#
#   A  Menopause status   — hormonal transition changes insulin sensitivity
#   B  MHT / HRT status   — currently on MHT = lowest iAUC (estrogen protective)
#   C  Diet preference    — plant → animal iAUC gradient
#   D  Smoking status     — never / former / current
#   E  Resistance training — weekly-hours dose response
#
# Weak / underpowered / null exposures live in fig3_other_vars_deepdive.R.
#
# Run:  Rscript fig2_novel_vars_deepdive.R
# Outputs: fig2_novel_vars_deepdive_final.pdf
# ============================================================

if (!exists("df"))              source("00_common_setup.R")
if (!exists("cgm_traj_base"))   source("deepdive_helpers.R")

# ---- Row A: Menopause status ----
men_df  <- df %>% filter(!is.na(menopause), !is.na(meal_55_meal_iauc_min)) %>%
  mutate(menopause = droplevels(menopause))
men_pal <- c("Pre" = pal_green, "Peri" = pal_amber, "Post" = pal_brick)

p_men <- make_boxplot(
  men_df, "menopause",
  pal = men_pal, pw = pairwise_brackets(men_df, "menopause"),
  n_df = count_labels(men_df, "menopause"),
  title_str = "Menopausal stage", x_lab = "Menopause stage")
p_men_traj <- make_trajectory("menopause", c("Pre","Peri","Post"), men_pal,
  title_str = "Glucose: menopause", legend_title = "Menopause")

# ---- Row B: MHT / HRT status ----
mht_df  <- df %>% filter(!is.na(mht_status), !is.na(meal_55_meal_iauc_min)) %>%
  mutate(mht_status = droplevels(mht_status))
mht_pal <- c("Currently taking" = pal_blue, "Never taken" = pal_lblue,
             "Taken in past" = pal_brick)

p_mht <- make_boxplot(
  mht_df, "mht_status",
  pal = mht_pal, pw = pairwise_brackets(mht_df, "mht_status"),
  n_df = count_labels(mht_df, "mht_status"),
  title_str = "MHT/HRT status", x_lab = "MHT/HRT status", x_angle = 12)
p_mht_traj <- make_trajectory("mht_status",
  c("Currently taking","Never taken","Taken in past"), mht_pal,
  title_str = "Glucose: MHT", legend_title = "MHT status")

# ---- Row C: Diet preference (questionnaire), plant → animal gradient ----
diet_df  <- df %>%
  filter(!is.na(diet_pref), diet_pref != "Other", !is.na(meal_55_meal_iauc_min)) %>%
  mutate(diet_pref = droplevels(diet_pref))
diet_pal <- c("Vegan" = pal_green, "Vegetarian" = pal_teal, "Pescatarian" = pal_blue,
              "Flexitarian" = pal_amber, "Mixed diet" = pal_brick)

p_diet <- make_boxplot(
  diet_df, "diet_pref",
  pal = diet_pal, pw = pairwise_brackets(diet_df, "diet_pref", sig_only = TRUE),
  n_df = count_labels(diet_df, "diet_pref"),
  title_str = "Diet preference (questionnaire)",
  x_lab = "Self-reported diet pattern", x_angle = 18)
p_diet_traj <- make_trajectory("diet_pref",
  c("Vegan","Vegetarian","Pescatarian","Flexitarian","Mixed diet"), diet_pal,
  title_str = "Glucose: diet preference", legend_title = "Diet")

# ---- Row D: Smoking status ----
smoke_df  <- df %>% filter(!is.na(smoking_cat), !is.na(meal_55_meal_iauc_min))
smoke_pal <- c("Never" = pal_green, "Former" = pal_amber, "Current" = pal_brick)

p_smoke <- make_boxplot(
  smoke_df, "smoking_cat",
  pal = smoke_pal, pw = pairwise_brackets(smoke_df, "smoking_cat"),
  n_df = count_labels(smoke_df, "smoking_cat"),
  title_str = "Smoking status", x_lab = "Smoking status")
p_smoke_traj <- make_trajectory("smoking_cat", c("Never","Former","Current"),
  smoke_pal, title_str = "Glucose: smoking", legend_title = "Smoking")

# ---- Row E: Resistance training (weekly hours bands) ----
resist_df  <- df %>% filter(!is.na(resist_group), !is.na(meal_55_meal_iauc_min))
resist_pal <- c("0 h" = pal_brick, "0-1 h" = pal_amber,
                "1-3 h" = pal_blue, ">3 h" = pal_green)

p_resist <- make_boxplot(
  resist_df, "resist_group",
  pal = resist_pal, pw = pairwise_brackets(resist_df, "resist_group"),
  n_df = count_labels(resist_df, "resist_group"),
  title_str = "Physical activity (resistance training)",
  x_lab = "Resistance training (h/week)")
p_resist_traj <- make_trajectory("resist_group", c("0 h","0-1 h","1-3 h",">3 h"),
  resist_pal, title_str = "Glucose: resistance training",
  legend_title = "Resistance (h/wk)")

# ---- Assemble & save ----
fig2 <- (
  make_row(p_men,    p_men_traj)    /
  make_row(p_mht,    p_mht_traj)    /
  make_row(p_diet,   p_diet_traj)   /
  make_row(p_smoke,  p_smoke_traj)  /
  make_row(p_resist, p_resist_traj)
)

ggsave("fig2_novel_vars_deepdive_final.pdf", fig2, width = 8.5, height = 15, limitsize = FALSE)
cat("Figure 2 saved.\n")
