library(tidyverse)
library(survey)
library(gt)
library(gtsummary)

## Setup

df <- readRDS("data/nhanes_L_analytic.rds")
df$WTDR2D[is.na(df$WTDR2D)] <- 0

## Build variables

df <- df |>
  mutate(
    choline_dens = 1000 * choline_mean / energy_mean,
    log_alt = log(LBXSATSI),
    log_ast = log(LBXSASSI),
    log_ggt = log(LBXSGTSI),
    race3 = fct_collapse(
      factor(RIDRETH3),
      white = "Non-Hispanic White",
      black = "Non-Hispanic Black",
      other_level = "other"
    ),
    chol_q = cut(
      choline_dens,
      breaks = quantile(
        choline_dens[in_aim2],
        probs = seq(0, 1, 0.25),
        na.rm = TRUE
      ),
      include.lowest = TRUE,
      labels = c("Q1 (Lowest)", "Q2", "Q3", "Q4 (Highest)")
    ),
    low_chol_obese = as.numeric(chol_q == "Q1 (Lowest)" & bmi_cat == "Obese"),
    bmi_chol_group = factor(
      paste0(bmi_cat, ", ", 
             if_else(chol_q == "Q1 (Lowest)", "low choline", "higher choline")),
      levels = c("Normal, higher choline", "Normal, low choline",
                 "Overweight, higher choline", "Overweight, low choline",
                 "Obese, higher choline", "Obese, low choline")
    )
  )

## Survey design

des_full <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTDR2D,
  data = df,
  nest = TRUE
)

des2 <- subset(
  des_full,
  in_aim2 %in% TRUE & bmi_cat != "Underweight" & !is.na(choline_dens)
)
des2 <- update(des2, bmi_cat = fct_drop(bmi_cat))

n_total <- nrow(des2$variables)
cat("Design df:", degf(des2), "\n")

## Sample sizes

des2$variables |>
  count(bmi_cat, chol_q) |>
  pivot_wider(names_from = chol_q, values_from = n) |>
  print()

## Table 1: liver enzymes by BMI category and choline quartile
## Cells with SE > 40% of the estimate are suppressed (NCHS guidance)

enzyme_cells <- map_dfr(
  c(ALT = "log_alt", AST = "log_ast", GGT = "log_ggt"),
  function(outcome) {
    svyby(
      reformulate(outcome),
      ~bmi_cat + chol_q,
      des2,
      svymean,
      na.rm = TRUE,
      vartype = c("se", "ci")
    ) |>
      as.data.frame() |>
      rename(est = 3) |>
      transmute(
        bmi_cat,
        chol_q,
        value = if_else(
          se > 0.4 * abs(est),
          "—",
          sprintf("%.1f (%.1f, %.1f)", exp(est), exp(ci_l), exp(ci_u))
        )
      )
  },
  .id = "enzyme"
)

tbl_1 <- enzyme_cells |>
  pivot_wider(names_from = chol_q, values_from = value) |>
  arrange(enzyme, bmi_cat) |>
  gt(groupname_col = "enzyme", rowname_col = "bmi_cat") |>
  tab_header(
    title = "Table 1. Serum liver enzymes by BMI category and choline intake quartile",
    subtitle = paste0(
      "Aim 2 Sample: US adults, NHANES August 2021–August 2023 (N = ",
      scales::comma(n_total), ")"
    )
  ) |>
  tab_spanner(
    label = "Choline intake quartile (mg/1,000 kcal)",
    columns = -bmi_cat
  ) |>
  tab_stubhead(label = "BMI category") |>
  cols_align(align = "center", columns = -bmi_cat) |>
  tab_footnote(
    footnote = "Values are survey-weighted geometric means (95% CI) in U/L.",
    locations = cells_title(groups = "title")
  ) |>
  tab_footnote(
    footnote = "Estimates suppressed (—) where the standard error exceeded 40% of the estimate, per NCHS guidance.",
    locations = cells_title(groups = "title")
  ) |>
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(11),
    heading.align = "left",
    column_labels.font.weight = "bold",
    row_group.font.weight = "bold",
    table.border.top.style = "solid",
    table.border.bottom.style = "solid",
    table_body.hlines.style = "none",
    data_row.padding = px(5)
  )

## Plot

alt_plot_data <- svyby(
  ~log_alt,
  ~bmi_cat + chol_q,
  des2,
  svymean,
  na.rm = TRUE,
  vartype = "ci"
) |>
  as.data.frame() |>
  transmute(
    bmi_cat,
    chol_q,
    alt = exp(log_alt),
    ci_lower = exp(ci_l),
    ci_upper = exp(ci_u)
  )

p_aim2 <- ggplot(alt_plot_data, aes(x = chol_q, y = alt, color = bmi_cat)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  labs(
    x = NULL,
    y = "ALT (U/L)",
    color = "BMI category",
    title = "ALT by choline quartile and BMI category"
  )

## Main test: low choline and obese vs. all
## Running for each liver enzyme

fit_enzyme <- function(outcome, exposure = "low_chol_obese") {
  f <- as.formula(
    paste(outcome, "~", exposure, "+ RIDAGEYR + sex + race3 + INDFMPIR")
  )
  svyglm(f, design = des2)
}

m_alt <- fit_enzyme("log_alt")
m_ast <- fit_enzyme("log_ast")
m_ggt <- fit_enzyme("log_ggt")

summary(m_alt)
summary(m_ast)
summary(m_ggt)

## BMI category by low vs. higher choline

m_alt_groups <- fit_enzyme("log_alt", "bmi_chol_group")
m_ast_groups <- fit_enzyme("log_ast", "bmi_chol_group")
m_ggt_groups <- fit_enzyme("log_ggt", "bmi_chol_group")

## Enzymes differ across the six groups?

regTermTest(m_alt_groups, ~bmi_chol_group)
regTermTest(m_ast_groups, ~bmi_chol_group)
regTermTest(m_ggt_groups, ~bmi_chol_group)

## Table 2: adjusted contrast across the three enzymes
## Bonferroni threshold for 3 outcomes: p < 0.0167

tidy_enzyme <- function(model) {
  tbl_regression(
    model,
    include = low_chol_obese,
    exponentiate = TRUE,
    label = list(low_chol_obese ~ "Low choline + Obesity")
  )
}

tbl_2 <- tbl_merge(
  list(tidy_enzyme(m_alt), tidy_enzyme(m_ast), tidy_enzyme(m_ggt)),
  tab_spanner = c("**ALT**", "**AST**", "**GGT**")
) |>
  modify_caption(
    "**Table 2. Adjusted difference in liver enzymes among adults with both low choline intake and obesity, compared with all other adults**"
  ) |>
  modify_footnote(
    everything() ~ paste(
      "Survey-weighted linear regression on log-transformed enzyme values;",
      "estimates are ratios of geometric means.",
      "Adjusted for age, sex, race and ethnicity, and family income-to-poverty ratio.",
      "Low choline defined as the lowest quartile of choline density; obesity as BMI ≥ 30.",
      "Significance set at a Bonferroni-corrected P < 0.017 for three outcomes."
    )
  )

## Table 3: liver enzymes by six BMI-choline groups

tidy_groups <- function(model) {
  tbl_regression(
    model,
    include = bmi_chol_group,
    exponentiate = TRUE,
    label = list(bmi_chol_group ~ "BMI category and choline intake")
  )
}

tbl_3 <- tbl_merge(
  list(tidy_groups(m_alt_groups), tidy_groups(m_ast_groups), tidy_groups(m_ggt_groups)),
  tab_spanner = c("**ALT**", "**AST**", "**GGT**")
) |>
  modify_caption("**Table 3. Liver enzymes by BMI category and choline intake group**") |>
  modify_footnote(
    everything() ~ paste(
      "Survey-weighted linear regression on log-transformed enzyme values;",
      "estimates are ratios of geometric means.",
      "Reference group: normal weight with choline above the lowest quartile.",
      "Adjusted for age, sex, race and ethnicity, and family income-to-poverty ratio."
    )
  )

## Save

dir.create("output", showWarnings = FALSE)

gtsave(tbl_1, "output/table1_enzymes_by_bmi_choline.html")
gtsave(tbl_1, "output/table1_enzymes_by_bmi_choline.png", vwidth = 1100, zoom = 2)

gtsave(as_gt(tbl_2), "output/table2_adjusted_contrast.html")
gtsave(as_gt(tbl_2), "output/table2_adjusted_contrast.png", vwidth = 900, zoom = 2)
gtsave(as_gt(tbl_3), "output/table3_bmi_choline_groups.png", vwidth = 1100, zoom = 2)

ggsave(
  "output/aim2_plot_alt.png",
  plot = p_aim2,
  width = 8,
  height = 5,
  units = "in",
  dpi = 300
)