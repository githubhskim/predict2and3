# ============================================================
# Figure 1 — CGM Exposome adjusted-R² bar chart
# ------------------------------------------------------------
# Bivariate adjusted R² of each predictor against Meal 55 iAUC,
# one bar per predictor, coloured by domain.
#
# Method (following Berry et al., Nat Med 2020): for each predictor,
#   lm(iAUC ~ predictor) -> adjusted R². Separate (non-hierarchical)
#   models, no covariate adjustment ("adjusted" = Wherry penalty for size).
#
# Run:  Rscript fig1_cgm_exposome_r2.R
# Outputs: fig1_cgm_exposome_r2_final.pdf, input/cgm_exposome_r2_table.csv
# ============================================================

if (!exists("df")) source("00_common_setup.R")

# ----------------------------------------------------------------
# 1.  FIT BIVARIATE ADJUSTED-R² MODELS (outcome = iAUC)
# ----------------------------------------------------------------

outcome_var <- "meal_55_meal_iauc"

cat("Fitting models...\n")
results_list <- list()

for (grp in names(predictor_groups)) {
  grp_info <- predictor_groups[[grp]]

  for (i in seq_along(grp_info$vars)) {
    var   <- grp_info$vars[i]
    label <- grp_info$labels[i]

    if (!var %in% names(df)) {
      cat(sprintf("  [MISSING] %s\n", var)); next
    }

    mod_data <- df %>%
      select(outcome = all_of(outcome_var), predictor = all_of(var)) %>%
      mutate(predictor = suppressWarnings(as.numeric(predictor))) %>%
      filter(!is.na(outcome), !is.na(predictor))

    n_obs <- nrow(mod_data)
    if (n_obs < MIN_N) {
      cat(sprintf("  [SKIP n=%d] %s\n", n_obs, var)); next
    }

    tryCatch({
      m    <- lm(outcome ~ predictor, data = mod_data)
      sm   <- summary(m)
      r2   <- pmax(sm$adj.r.squared, 0)              # Wherry-adjusted R²
      fst  <- sm$fstatistic
      pval <- pf(fst[1], fst[2], fst[3], lower.tail = FALSE)

      results_list[[length(results_list) + 1]] <- tibble(
        group = grp, variable = var, label = label,
        adj_r2 = r2, p_value = pval, n = n_obs
      )
    }, error = function(e) cat(sprintf("  [ERROR] %s: %s\n", var, e$message)))
  }
}

results <- bind_rows(results_list)

cat(sprintf("\nSuccessfully modeled %d predictors.\n", nrow(results)))

# ----------------------------------------------------------------
# 2.  PLOT DATA
# ----------------------------------------------------------------
# Drop predictors with 0% explained variance, then order rows by domain and
# descending R² within domain. Levels reversed because ggplot draws level 1 at
# the bottom, so this puts the first domain / largest bar at the top.

plot_df <- results %>%
  filter(adj_r2 > 0) %>%
  mutate(
    group        = factor(group, levels = names(domain_colours)),
    label_with_n = paste0(label, "  (n=", scales::comma(n), ")")
  ) %>%
  arrange(group, desc(adj_r2)) %>%
  mutate(label_ordered = factor(label_with_n, levels = rev(unique(label_with_n))))

# ----------------------------------------------------------------
# 3.  FIGURE 1
# ----------------------------------------------------------------

fig1 <- ggplot(plot_df, aes(x = adj_r2 * 100, y = label_ordered, fill = group)) +
  geom_col(width = 0.7, show.legend = TRUE) +
  geom_text(aes(label = sprintf("%.2f%%", adj_r2 * 100)),
            hjust = -0.1, size = 2.9, colour = "black") +
  scale_fill_manual(values = domain_colours, labels = domain_labels_clean, name = NULL) +
  # sqrt x-scale: starts at 0 (sqrt breaks at negatives); headroom on the right
  # so the value labels are not clipped.
  scale_x_sqrt(expand = expansion(mult = c(0, 0.14)),
               limits = c(0, NA),
               breaks = c(0, 0.5, 1, 2, 3, 5, 7, 9),
               labels = function(x) paste0(x, "%")) +
  guides(fill = guide_legend(ncol = 1, keywidth = 1.0, keyheight = 0.9)) +
  labs(x = "Adjusted R²  (%)", y = NULL) +
  theme_paper +
  theme(
    panel.grid.major.y = element_blank(),   # keep only vertical reference lines
    panel.grid.minor   = element_blank(),
    axis.text.y        = element_text(size = 9),
    axis.text.x        = element_text(size = 9),
    axis.title.x       = element_text(size = 10),
    legend.position    = c(0.66, 0.25),
    legend.text        = element_text(size = 9),
    legend.key.size    = unit(0.40, "cm")
  )

ggsave("fig1_cgm_exposome_r2_final.pdf", fig1, width = 7.5, height = 10)
cat("Figure 1 saved.\n")

# ----------------------------------------------------------------
# 4.  RESULTS TABLE  (written to input/ for make_supp_tables_clean.R)
# ----------------------------------------------------------------

results_out <- results %>%
  arrange(group, desc(adj_r2)) %>%
  mutate(
    p_adj_bh   = p.adjust(p_value, method = "BH"),
    `Adj R² %` = round(adj_r2 * 100, 3),
    `P-value`  = signif(p_value, 3),
    `BH-FDR q` = signif(p_adj_bh, 3)
  ) %>%
  select(Domain = group, Predictor = label, Variable = variable,
         `Adj R² %`, `P-value`, `BH-FDR q`, N = n) %>%
  mutate(Domain = gsub("\n", " ", Domain))   # flatten two-line labels for CSV

dir.create("input", showWarnings = FALSE)
write.csv(results_out, file.path("input", "cgm_exposome_r2_table.csv"),
          row.names = FALSE)
cat("Results table saved to input/cgm_exposome_r2_table.csv\n")
