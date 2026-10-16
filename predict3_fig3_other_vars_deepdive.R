# ============================================================
# Figure 3 — Deep dive: other lifestyle variables
# ------------------------------------------------------------
# The counterpoint to Figure 2: exposures that are weak (sleep quality),
# underpowered (hunger, n≈986), or null (MVPA, sitting, GI). Physical-activity
# panels lead so they sit alongside the resistance-training panel promoted to
# Figure 2, and the specificity contrast reads directly.
#
#   A  MVPA (quartiled)          — boxplot
#   B  Sitting time (quartiled)  — boxplot
#   C  Sleep quality (FFQ)       — boxplot
#   D  Hunger frequency          — line (dose-response across Likert scale)
#   E  GI inconvenience          — line (dose-response across 7-level severity)
#
# Boxplots pair with a mean glucose trajectory; the two ordinal line panels
# show the iAUC dose-response trend with a P-trend from a linear slope test.
#
# Run:  Rscript fig3_other_vars_deepdive.R
# Outputs: fig3_other_vars_deepdive_final.pdf
# ============================================================

if (!exists("df"))            source("00_common_setup.R")
if (!exists("cgm_traj_base")) source("deepdive_helpers.R")

# Format a P-trend significance star for the line panels.
trend_star <- function(p) dplyr::case_when(
  p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ " (ns)")

# ---- Panel A: MVPA (quartiled) ----
mvpa_df  <- df %>% filter(!is.na(mvpa_q), !is.na(meal_55_meal_iauc_min))
mvpa_pal <- c("Q1 (Low)" = pal_lblue, "Q2" = pal_blue,
              "Q3" = pal_amber, "Q4 (High)" = pal_brick)

p_mvpa <- make_boxplot(
  mvpa_df, "mvpa_q",
  pal = mvpa_pal, pw = pairwise_brackets(mvpa_df, "mvpa_q"),
  n_df = count_labels(mvpa_df, "mvpa_q"),
  title_str = "Physical activity (MVPA)", x_lab = "MVPA quartile (min/week)")
p_mvpa_traj <- make_trajectory("mvpa_q", c("Q1 (Low)","Q2","Q3","Q4 (High)"),
  mvpa_pal, title_str = "Glucose: MVPA", legend_title = "MVPA")

# ---- Panel B: Sitting time (quartiled) ----
sitting_df  <- df %>% filter(!is.na(sitting_q), !is.na(meal_55_meal_iauc_min))
sitting_pal <- c("Q1 (Low)" = pal_green, "Q2" = pal_amber,
                 "Q3" = pal_blue, "Q4 (High)" = pal_brick)

p_sitting <- make_boxplot(
  sitting_df, "sitting_q",
  pal = sitting_pal, pw = pairwise_brackets(sitting_df, "sitting_q"),
  n_df = count_labels(sitting_df, "sitting_q"),
  title_str = "Physical activity (daily sitting time)",
  x_lab = "Sitting time quartile (h/day)")
p_sitting_traj <- make_trajectory("sitting_q", c("Q1 (Low)","Q2","Q3","Q4 (High)"),
  sitting_pal, title_str = "Glucose: sitting time", legend_title = "Sitting")

# ---- Panel C: Sleep quality (FFQ) ----
sleep_df  <- df %>% filter(!is.na(sleep_cat), !is.na(meal_55_meal_iauc_min))
sleep_pal <- c("Very bad" = pal_brick, "Fairly bad" = pal_amber,
               "Fairly good" = pal_green, "Very good" = pal_blue)

p_sleep <- make_boxplot(
  sleep_df, "sleep_cat",
  pal = sleep_pal, pw = pairwise_brackets(sleep_df, "sleep_cat"),
  n_df = count_labels(sleep_df, "sleep_cat"),
  title_str = "Sleep quality", x_lab = "Self-reported sleep quality", x_angle = 18)
p_sleep_traj <- make_trajectory("sleep_cat",
  c("Very bad","Fairly bad","Fairly good","Very good"), sleep_pal,
  title_str = "Glucose: sleep quality", legend_title = "Sleep quality")

# ---- Panel D: Hunger frequency (line — ordinal Likert scale) ----
hunger_summ <- df %>%
  filter(!is.na(hunger_freq), !is.na(meal_55_meal_iauc_min)) %>%
  group_by(hunger_score = hunger_freq) %>%
  summarise(mean_iauc = mean(meal_55_meal_iauc_min),
            se = sd(meal_55_meal_iauc_min) / sqrt(n()), n = n(), .groups = "drop")

hunger_p <- summary(lm(meal_55_meal_iauc_min ~ hunger_freq,
                       data = df %>% filter(!is.na(hunger_freq),
                                            !is.na(meal_55_meal_iauc_min)))
                    )$coefficients["hunger_freq", "Pr(>|t|)"]

p_hunger <- ggplot(hunger_summ, aes(x = hunger_score, y = mean_iauc)) +
  geom_ribbon(aes(ymin = mean_iauc - 1.96 * se, ymax = mean_iauc + 1.96 * se),
              fill = pal_amber, alpha = 0.18) +
  geom_line(colour = pal_amber, linewidth = 0.9) +
  geom_point(aes(size = n), colour = pal_amber, fill = "white",
             shape = 21, stroke = 1.1) +
  geom_text(aes(y = mean_iauc + 1.96 * se + 5, label = paste0("n=", scales::comma(n))),
            vjust = 0, size = 4.3, colour = "black") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("P-trend=", formatC(hunger_p, format = "e", digits = 1),
                          trend_star(hunger_p)),
           size = 4.3, hjust = 1.05, vjust = 1.6, colour = "black", fontface = "italic") +
  scale_x_continuous(breaks = 1:5,
    labels = c("Strongly\ndisagree","Disagree","Neutral","Agree","Strongly\nagree"),
    expand = expansion(mult = c(0.07, 0.07))) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.05, 0.28))) +
  scale_size_continuous(range = c(2.2, 5.5), guide = "none") +
  labs(title = "Hunger frequency", x = '"I feel hungry often"',
       y = "Mean iAUC (mmol/L·min)") +
  theme_dive

hunger_traj_pal <- c("Low (1-2)" = pal_green, "Mid (3)" = pal_amber,
                     "High (4-5)" = pal_brick)
p_hunger_traj <- make_trajectory("hunger_group", c("Low (1-2)","Mid (3)","High (4-5)"),
  hunger_traj_pal, title_str = "Glucose: hunger", legend_title = "Hunger level")

# ---- Panel E: GI inconvenience (line — 7-level severity scale) ----
gut_summ <- df %>%
  filter(!is.na(gut_inconvenienc), !is.na(meal_55_meal_iauc_min)) %>%
  group_by(gut_score = gut_inconvenienc) %>%
  summarise(mean_iauc = mean(meal_55_meal_iauc_min),
            se = sd(meal_55_meal_iauc_min) / sqrt(n()), n = n(), .groups = "drop")

gut_p <- summary(lm(meal_55_meal_iauc_min ~ gut_inconvenienc,
                    data = df %>% filter(!is.na(gut_inconvenienc),
                                         !is.na(meal_55_meal_iauc_min)))
                 )$coefficients["gut_inconvenienc", "Pr(>|t|)"]

p_gut <- ggplot(gut_summ, aes(x = gut_score, y = mean_iauc)) +
  geom_ribbon(aes(ymin = mean_iauc - 1.96 * se, ymax = mean_iauc + 1.96 * se),
              fill = pal_brown, alpha = 0.18) +
  geom_line(colour = pal_brown, linewidth = 0.9) +
  geom_point(aes(size = n), colour = pal_brown, fill = "white",
             shape = 21, stroke = 1.1) +
  geom_text(aes(y = mean_iauc + 1.96 * se + 4, label = paste0("n=", scales::comma(n))),
            vjust = 0, size = 3.2, colour = "black") +
  annotate("text", x = -Inf, y = Inf,
           label = paste0("P-trend=", formatC(gut_p, format = "e", digits = 1),
                          trend_star(gut_p)),
           size = 4.3, hjust = -0.05, vjust = 1.6, colour = "black", fontface = "italic") +
  # Extra left padding so the wide centred "None" (x = 0) label is not clipped.
  scale_x_continuous(breaks = 0:6,
    labels = c("None","Minor","Mild","Moderate","Mod-Sev","Severe","V.Severe"),
    expand = expansion(mult = c(0.10, 0.05))) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0.05, 0.18))) +
  scale_size_continuous(range = c(2.2, 5.5), guide = "none") +
  labs(title = "GI inconvenience", x = "Severity", y = "Mean iAUC (mmol/L·min)") +
  theme_dive +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

gut_traj_pal <- c("None/Minor" = pal_green, "Mild/Mod" = pal_amber,
                  "Mod/Severe" = pal_brick)
p_gut_traj <- make_trajectory("gut_group", c("None/Minor","Mild/Mod","Mod/Severe"),
  gut_traj_pal, title_str = "Glucose: GI inconvenience", legend_title = "GI severity")

# ---- Assemble & save ----
fig3 <- (
  make_row(p_mvpa,    p_mvpa_traj)    /
  make_row(p_sitting, p_sitting_traj) /
  make_row(p_sleep,   p_sleep_traj)   /
  make_row(p_hunger,  p_hunger_traj)  /
  make_row(p_gut,     p_gut_traj)
)

ggsave("fig3_other_vars_deepdive_final.pdf", fig3, width = 8.5, height = 15, limitsize = FALSE)
cat("Figure 3 (MVPA + sitting + sleep + hunger + GI) saved.\n")
