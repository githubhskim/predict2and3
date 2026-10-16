# ============================================================
# Deep-dive helpers — shared by Figure 2 and the supplementary figure
# ------------------------------------------------------------
# Provides:
#   cgm_traj_base   raw CGM trace (0-120 min) joined with grouping variables
#   make_trajectory()  mean glucose trajectory by group (line + 95% CI ribbon)
#   make_boxplot()     violin + boxplot of iAUC by group, with pairwise brackets
#   make_row()         patchwork row: iAUC panel (wider) | trajectory (narrower)
#
# Sourced (guarded) by fig2_novel_vars_deepdive.R and fig_supp_mvpa_gi.R.
# ============================================================

if (!exists("df")) source("00_common_setup.R")

# --- CGM trajectory base ------------------------------------------------------
# Snap CGM timestamps to a 15-min grid; keep the 0-120 min post-meal window.
cgm_traj_base <- cgm_raw %>%
  select(username, inc_measured_at, glucose_level, inc_glucose_level) %>%
  inner_join(
    df %>% select(username, menopause, mht_status, sleep_cat, smoking_cat,
                  mvpa_q, hunger_group, gut_group, sitting_q, resist_group,
                  diet_pref, diet_app),
    by = "username"
  ) %>%
  mutate(time_min = round(inc_measured_at / 60 / 15) * 15) %>%
  filter(time_min >= 0, time_min <= 120)

# --- Mean glucose trajectory (0-120 min) by group ----------------------------
make_trajectory <- function(group_var, group_levels, colours,
                            title_str = NULL, legend_title = NULL) {
  summ <- cgm_traj_base %>%
    filter(!is.na(.data[[group_var]])) %>%
    mutate(grp = factor(.data[[group_var]], levels = group_levels)) %>%
    filter(!is.na(grp)) %>%   # drop categories outside group_levels
    group_by(time_min, grp) %>%
    summarise(mean_gl = mean(glucose_level, na.rm = TRUE),
              se      = sd(glucose_level,   na.rm = TRUE) / sqrt(n()),
              .groups = "drop")

  ggplot(summ, aes(x = time_min, y = mean_gl, colour = grp, fill = grp)) +
    geom_ribbon(aes(ymin = mean_gl - 1.96 * se, ymax = mean_gl + 1.96 * se),
                alpha = 0.30, colour = NA) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = colours, name = legend_title) +
    scale_fill_manual(values = colours, guide = "none") +
    scale_x_continuous(breaks = seq(0, 120, 30),
                       expand = expansion(mult = c(0.02, 0.05))) +
    # Top headroom so the in-panel legend clears the descending limb of the curve
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.28))) +
    labs(title = title_str, x = "Time post-meal (min)", y = "Glucose (mmol/L)") +
    theme_dive +
    theme(legend.position      = c(0.98, 0.98),
          legend.justification = c(1, 1),
          legend.background    = element_rect(fill = "white", colour = "black",
                                              linewidth = 0.3),
          legend.margin        = margin(3, 5, 3, 5),
          legend.text          = element_text(size = 12),
          legend.key.width     = unit(0.75, "cm"),
          legend.key.height    = unit(0.45, "cm"),
          legend.title         = element_blank())
}

# --- Violin + boxplot of iAUC by group, with pairwise significance brackets ---
make_boxplot <- function(data, x_var, y_var = "meal_55_meal_iauc_min",
                         pal, pw, n_df, title_str,
                         x_lab = NULL, y_lab = "iAUC (mmol/L·min)",
                         x_angle = 0) {
  p <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]],
                        fill = .data[[x_var]])) +
    geom_violin(alpha = 0.35, colour = NA, trim = TRUE,
                scale = "width", width = 0.78) +
    geom_boxplot(width = 0.13, outlier.shape = NA, fill = "white",
                 colour = "black", linewidth = 0.45) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 2.6,
                 colour = "black") +
    ggpubr::stat_pvalue_manual(pw, label = "p.adj.signif", tip.length = 0.012,
                               size = 5.5, hide.ns = TRUE, step.increase = 0.07) +
    geom_text(data = n_df, aes(y = y_pos, label = lbl, fill = NULL),
              size = 3.4, colour = "black") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0.10, 0.30))) +
    labs(title = title_str, x = x_lab, y = y_lab) +
    theme_dive
  if (x_angle != 0)
    p <- p + theme(axis.text.x = element_text(angle = x_angle, hjust = 1))
  p
}

# --- Pairwise Bonferroni Wilcoxon test, formatted for stat_pvalue_manual ------
# Convenience wrapper: keeps the four boxplot panels from repeating the same
# rstatix pipeline. `sig_only` drops non-significant brackets (declutters
# panels with many categories, e.g. diet preference).
pairwise_brackets <- function(data, x_var, y_var = "meal_55_meal_iauc_min",
                              sig_only = FALSE) {
  pw <- data %>%
    rstatix::pairwise_wilcox_test(as.formula(paste(y_var, "~", x_var)),
                                  p.adjust.method = "bonferroni") %>%
    rstatix::add_xy_position(x = x_var) %>%
    mutate(p.adj.signif = dplyr::recode(p.adj.signif, "****" = "***"))
  if (sig_only) pw <- pw %>% filter(p.adj.signif != "ns")
  pw
}

# --- Per-group n labels placed below the boxes --------------------------------
count_labels <- function(data, x_var, y_pos = -25) {
  data %>% count(across(all_of(x_var))) %>%   # across() keeps the original column name
    mutate(y_pos = y_pos, lbl = paste0("n=", scales::comma(n)))
}

# --- patchwork row: iAUC panel (wider) | trajectory (narrower) ----------------
make_row <- function(l, r) (l | r) + plot_layout(widths = c(1.0, 0.8))
