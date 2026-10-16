# ============================================================
# Figure 4 — Multi-outcome CGM exposome heatmap
# ------------------------------------------------------------
# CGM response is multi-dimensional: beyond total iAUC, the SHAPE of the
# postprandial curve carries distinct biological signal. Each predictor is
# fit against five CGM outcomes; the heatmap shows SIGNED adjusted R²
# (sign = direction of the predictor's effect).
#
#   iAUC          total glucose exposure above baseline (primary)
#   Peak rise     how high glucose spikes (max Δglucose)
#   Time to peak  how fast it peaks (gastric emptying / IR kinetics)
#   2h Δglucose   where glucose sits at 2 h — sustained or resolved
#   Post-peak dip how far glucose falls BELOW baseline (iAUC is blind to this)
#
# Layout: [predictor labels] | [domain colour strip] | [heatmap cells].
# Rows are the same non-zero-iAUC-R² predictors shown in Figure 1, in the same
# order (domain, then descending iAUC R²).
#
# Run:  Rscript fig4_cgm_multi_outcome_heatmap.R
# Outputs: fig4_cgm_multi_outcome_heatmap_final.pdf, input/cgm_multi_outcome_r2_table.csv
# ============================================================

if (!exists("df")) source("00_common_setup.R")

# ----------------------------------------------------------------
# 1.  FIT EACH PREDICTOR AGAINST EACH CGM OUTCOME
# ----------------------------------------------------------------

cgm_outcomes <- c(
  "meal_55_meal_iauc"    = "iAUC",
  "cgm_peak_inc"         = "Peak Rise",
  "cgm_time_to_peak_min" = "Time to Peak",
  "cgm_2h_inc"           = "2-hr Glucose",
  "cgm_post_peak_dip"    = "Post-peak Dip"
)

cat("Fitting multi-outcome CGM models...\n")
results_multi_list <- list()

for (outcome_nm in names(cgm_outcomes)) {
  if (!outcome_nm %in% names(df)) {
    cat(sprintf("  [OUTCOME MISSING] %s — skipping\n", outcome_nm)); next
  }

  for (grp in names(predictor_groups)) {
    grp_info <- predictor_groups[[grp]]

    for (i in seq_along(grp_info$vars)) {
      var   <- grp_info$vars[i]
      label <- grp_info$labels[i]
      if (!var %in% names(df)) next

      mod_data <- df %>%
        select(outcome = all_of(outcome_nm), predictor = all_of(var)) %>%
        mutate(predictor = suppressWarnings(as.numeric(predictor))) %>%
        filter(!is.na(outcome), !is.na(predictor))

      if (nrow(mod_data) < MIN_N) next

      tryCatch({
        m         <- lm(outcome ~ predictor, data = mod_data)
        sm        <- summary(m)
        r2        <- pmax(sm$adj.r.squared, 0)
        coef_sign <- sign(coef(m)["predictor"])
        fst       <- sm$fstatistic
        pval      <- pf(fst[1], fst[2], fst[3], lower.tail = FALSE)

        results_multi_list[[length(results_multi_list) + 1]] <- tibble(
          outcome = outcome_nm, group = grp, variable = var, label = label,
          adj_r2 = r2, signed_r2 = r2 * coef_sign, p_value = pval,
          n = nrow(mod_data)
        )
      }, error = function(e)
        cat(sprintf("  [ERROR] %s × %s: %s\n", outcome_nm, var, e$message)))
    }
  }
  cat(sprintf("  Completed outcome: %s\n", outcome_nm))
}

results_multi <- bind_rows(results_multi_list)
cat(sprintf("\nMulti-outcome models fit: %d rows.\n", nrow(results_multi)))

# ----------------------------------------------------------------
# 2.  BUILD HEATMAP DATA
# ----------------------------------------------------------------
# BH-FDR is applied within each CGM outcome (separately per column).

heat_df <- results_multi %>%
  group_by(outcome) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    sig_star    = case_when(p_adj < 0.001 ~ "***", p_adj < 0.01 ~ "**",
                            p_adj < 0.05  ~ "*",   TRUE ~ ""),
    cell_label  = sig_star,
    text_colour = ifelse(abs(signed_r2 * 100) > 1.2, "white", "grey15"),
    outcome_lab = factor(cgm_outcomes[outcome], levels = cgm_outcomes),
    group       = factor(group, levels = names(domain_colours))
  )

# Restrict to the same predictors shown in Figure 1 (non-zero iAUC adj R²).
keep_labels <- heat_df %>%
  filter(outcome == "meal_55_meal_iauc", adj_r2 > 0) %>% pull(label)
heat_df <- heat_df %>% filter(label %in% keep_labels)

# Row order: domain, then descending iAUC R² (reversed so the first domain /
# largest bar sits at the top, matching Figure 1).
predictor_order <- heat_df %>%
  filter(outcome == "meal_55_meal_iauc") %>%
  arrange(group, desc(adj_r2)) %>% pull(label) %>% rev()
heat_df <- heat_df %>% mutate(label = factor(label, levels = predictor_order))

# ----------------------------------------------------------------
# 3.  FIGURE 4  —  [labels] | [domain strip] | [heatmap]
# ----------------------------------------------------------------

cap_r2 <- 2.5   # colour scale capped at ±2.5% (outliers squished to the limit)

domain_bar_df <- heat_df %>% distinct(label, group) %>% mutate(x = 1)

# Left panel: predictor labels, right-aligned.
p_ylabels <- ggplot(domain_bar_df, aes(x = x, y = label)) +
  geom_blank() +
  scale_x_continuous(expand = expansion(0)) +
  scale_y_discrete(position = "right") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11, base_family = "Helvetica") +
  theme(
    panel.grid        = element_blank(),
    panel.background  = element_blank(),
    axis.text.x       = element_blank(),
    axis.ticks        = element_blank(),
    axis.text.y.right = element_text(size = 12, hjust = 1, colour = "black"),
    plot.margin       = margin(5, 2, 5, 5)
  )

# Middle panel: domain colour strip.
p_domain_bar <- ggplot(domain_bar_df, aes(x = x, y = label, fill = group)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = domain_colours, guide = "none") +
  theme_void()

# Right panel: heatmap cells + domain legend.
fig4_main <- ggplot(heat_df, aes(x = outcome_lab, y = label)) +
  geom_tile(aes(fill = signed_r2 * 100), colour = "white", linewidth = 0.45) +
  geom_text(aes(label = cell_label, colour = text_colour), size = 4) +
  scale_colour_identity() +
  geom_point(aes(shape = group), size = 0, alpha = 0) +   # dummy layer to drive the domain legend
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-cap_r2, cap_r2), oob = scales::squish,
    name = "Signed\nadj R² (%)", labels = function(x) paste0(x, "%"),
    na.value = "grey92"
  ) +
  scale_shape_manual(
    values = rep(22, length(domain_colours)),
    labels = domain_labels_clean, name = "Domain",
    guide  = guide_legend(
      override.aes = list(fill = unname(domain_colours),
                          colour = unname(domain_colours), size = 5, alpha = 1),
      ncol = 1, keywidth = 0.9, keyheight = 0.9)
  ) +
  labs(x = NULL, y = NULL) +
  theme_paper +
  theme(
    panel.grid        = element_blank(),
    axis.text.x       = element_text(size = 12, angle = 90, hjust = 1,
                                     vjust = 0.5, colour = "black"),
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    legend.title      = element_text(size = 13, colour = "black"),
    legend.text       = element_text(size = 12, colour = "black"),
    legend.key.size   = unit(0.40, "cm"),
    legend.background = element_blank()
  )

fig4 <- (p_ylabels | p_domain_bar | fig4_main) +
  plot_layout(widths = c(1.0, 0.05, 0.9))

ggsave("fig4_cgm_multi_outcome_heatmap_final.pdf", fig4, width = 10, height = 14)
cat("Figure 4 (multi-outcome heatmap) saved.\n")

# ----------------------------------------------------------------
# 5.  MULTI-OUTCOME RESULTS TABLE  (written to input/ for make_supp_tables_clean.R)
# ----------------------------------------------------------------

results_multi_out <- results_multi %>%
  group_by(outcome) %>%
  mutate(p_adj_bh = p.adjust(p_value, method = "BH")) %>%   # BH-FDR within each CGM outcome
  ungroup() %>%
  arrange(outcome, group, desc(adj_r2)) %>%
  mutate(
    sig = case_when(
      p_adj_bh < 0.001 ~ "FDR***", p_adj_bh < 0.01 ~ "FDR**",
      p_adj_bh < 0.05  ~ "FDR*",   p_value  < 0.05 ~ "Nominal", TRUE ~ "NS"),
    `Signed R2 pct`   = round(signed_r2 * 100, 3),
    `Adj R2 pct`      = round(adj_r2    * 100, 3),
    `P-value`         = signif(p_value,  3),
    `P-FDR`           = signif(p_adj_bh, 3),
    CGM_Outcome_label = cgm_outcomes[outcome]
  ) %>%
  select(CGM_Outcome = outcome, CGM_Label = CGM_Outcome_label,
         Domain = group, Predictor = label, Variable = variable,
         `Signed R2 pct`, `Adj R2 pct`, `P-value`, `P-FDR`,
         Significance = sig, N = n) %>%
  mutate(Domain = gsub("\n", " ", Domain))   # flatten two-line labels for CSV

dir.create("input", showWarnings = FALSE)
write.csv(results_multi_out, file.path("input", "cgm_multi_outcome_r2_table.csv"),
          row.names = FALSE)
cat("Multi-outcome results table saved to input/cgm_multi_outcome_r2_table.csv\n")
