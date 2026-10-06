library(tidyverse)
library(survey)
library(gt)
library(gtsummary)

# Setup

df <- readRDS("data/nhanes_L_analytic.rds")

ana <- df |>
  filter(in_aim1, bmi_cat != "Underweight") |>
  mutate(bmi_cat = fct_drop(bmi_cat))

des <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTDR2D,
  data = ana,
  nest = TRUE
)

# Sex and BMI category

# Counts below the AI by sex and BMI category
sex_counts <- ana |>
  filter(!is.na(below_ai), 
         !is.na(sex), 
         !is.na(bmi_cat)) |>
  group_by(sex, bmi_cat) |>
  summarise(n_below_ai = sum(below_ai == 1), .groups = "drop")

# Total sample size by sex
sex_totals <- ana |>
  filter(!is.na(below_ai), 
         !is.na(sex), 
         !is.na(bmi_cat)) |>
  count(sex, name = "N")

Female_N <- sex_totals |> 
  filter(sex == "Female") |> 
  pull(N)

Male_N <- sex_totals |> 
  filter(sex == "Male") |> 
  pull(N)

# Weighted prevalence and confidence intervals
prev_sex_bmi <- svyby(
  ~below_ai,
  ~sex + bmi_cat,
  des,
  svymean,
  na.rm = TRUE,
  vartype = "ci"
) |>
  as.data.frame() |>
  transmute(
    sex,
    bmi_cat,
    prevalence = 100 * as.numeric(below_ai),
    ci_lower = 100 * as.numeric(ci_l),
    ci_upper = 100 * as.numeric(ci_u)
  ) |>
  left_join(sex_counts, by = c("sex", "bmi_cat")) |>
  mutate(
    cell = sprintf(
      "%d, %.1f%% (%.1f–%.1f)",
      n_below_ai,
      prevalence,
      ci_lower,
      ci_upper
    )
  )

# Table data
table_sex_bmi <- prev_sex_bmi |>
  select(bmi_cat, sex, cell) |>
  pivot_wider(names_from = sex, values_from = cell) |>
  select(bmi_cat, Female, Male)

tbl_aim1 <- table_sex_bmi |>
  gt() |>
  tab_header(
    title = "Prevalence of choline intake below the AI",
    subtitle = "Unweighted count below AI, survey-weighted % (95% CI)"
  ) |>
  cols_label(
    bmi_cat = "BMI category",
    Female = md(paste0("Female<br>(N = ", scales::comma(Female_N), ")")),
    Male = md(paste0("Male<br>(N = ", scales::comma(Male_N), ")"))
  ) |>
  cols_align(align = "left", columns = bmi_cat) |>
  cols_align(align = "center", columns = c(Female, Male)) |>
  tab_options(
    table.font.size = px(13),
    heading.title.font.size = px(18),
    heading.subtitle.font.size = px(11),
    column_labels.font.weight = "bold",
    data_row.padding = px(8)
  ) |>
  opt_row_striping()

# Prevalence plot by sex and BMI category
p_aim1_sex_bmi <- ggplot(
  prev_sex_bmi,
  aes(x = bmi_cat, y = prevalence, color = sex)
) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  labs(
    x = NULL,
    y = "% below choline AI",
    color = "Sex",
    title = "Prevalence of choline intake below the AI"
  )

# By BMI category

# Counts below the AI by BMI category
bmi_counts <- ana |>
  filter(!is.na(below_ai), 
         !is.na(bmi_cat)) |>
  group_by(bmi_cat) |>
  summarise(n_below_ai = sum(below_ai == 1), 
            .groups = "drop")

# Total analytic sample size
bmi_N <- ana |>
  filter(!is.na(below_ai), 
         !is.na(bmi_cat)) |>
  nrow()

prev_bmi <- svyby(
  ~below_ai,
  ~bmi_cat,
  des,
  svymean,
  na.rm = TRUE,
  vartype = "ci"
) |>
  as.data.frame() |>
  transmute(
    bmi_cat,
    prevalence = 100 * as.numeric(below_ai),
    ci_lower = 100 * as.numeric(ci_l),
    ci_upper = 100 * as.numeric(ci_u)
  ) |>
  left_join(bmi_counts, by = "bmi_cat") |>
  mutate(
    overall = sprintf(
      "%d, %.1f%% (%.1f–%.1f)",
      n_below_ai,
      prevalence,
      ci_lower,
      ci_upper
    )
  )

tbl_bmi <- prev_bmi |>
  select(bmi_cat, overall) |>
  gt() |>
  tab_header(
    title = "Prevalence of choline intake below the AI by BMI category",
    subtitle = "Unweighted count below AI, survey-weighted % (95% CI)"
  ) |>
  cols_label(
    bmi_cat = "BMI category",
    overall = md(paste0("Overall<br>(N = ", scales::comma(bmi_N), ")"))
  ) |>
  cols_align(align = "left", columns = bmi_cat) |>
  cols_align(align = "center", columns = overall) |>
  tab_options(
    table.font.size = px(13),
    heading.title.font.size = px(18),
    heading.subtitle.font.size = px(11),
    column_labels.font.weight = "bold",
    data_row.padding = px(8)
  ) |>
  opt_row_striping()

# Prevalence plot by BMI category
p_aim1_bmi <- ggplot(prev_bmi, aes(x = bmi_cat, y = prevalence)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  labs(
    x = NULL,
    y = "% below choline AI",
    title = "Prevalence of choline intake below the AI by BMI category"
  ) 

# Aim 1: Women by menopause and BMI

ana_f <- ana |>
  filter(sex == "Female", menopause %in% c("pre", "post"))

des_f <- svydesign(
  ids = ~SDMVPSU,
  strata = ~SDMVSTRA,
  weights = ~WTDR2D,
  data = ana_f,
  nest = TRUE
)

# Counts below the AI by menopause status and BMI category
menopause_counts <- ana_f |>
  filter(!is.na(below_ai), 
         !is.na(menopause), 
         !is.na(bmi_cat)) |>
  group_by(menopause, 
           bmi_cat) |>
  summarise(n_below_ai = sum(below_ai == 1), 
            .groups = "drop")

# Total sample size by menopause status
menopause_totals <- ana_f |>
  filter(!is.na(below_ai), 
         !is.na(menopause), 
         !is.na(bmi_cat)) |>
  count(menopause, name = "N")

pre_N <- menopause_totals |> 
  filter(menopause == "pre") |> 
  pull(N)

post_N <- menopause_totals |> 
  filter(menopause == "post") |> 
  pull(N)

prev_menopause <- svyby(
  ~below_ai,
  ~menopause + bmi_cat,
  des_f,
  svymean,
  na.rm = TRUE,
  vartype = "ci"
) |>
  as.data.frame() |>
  transmute(
    menopause,
    bmi_cat,
    prevalence = 100 * as.numeric(below_ai),
    ci_lower = 100 * as.numeric(ci_l),
    ci_upper = 100 * as.numeric(ci_u)
  ) |>
  left_join(menopause_counts, by = c("menopause", "bmi_cat")) |>
  mutate(
    cell = sprintf(
      "%d, %.1f%% (%.1f–%.1f)",
      n_below_ai,
      prevalence,
      ci_lower,
      ci_upper
    )
  )

table_menopause <- prev_menopause |>
  select(bmi_cat, menopause, cell) |>
  pivot_wider(names_from = menopause, 
              values_from = cell) |>
  select(bmi_cat, pre, post)

tbl_aim1_f <- table_menopause |>
  gt() |>
  tab_header(
    title = "Prevalence of choline intake below the AI among women",
    subtitle = "Unweighted count below AI, survey-weighted % (95% CI)"
  ) |>
  cols_label(
    bmi_cat = "BMI category",
    pre = md(paste0("Pre-menopausal<br>(N = ", scales::comma(pre_N), ")")),
    post = md(paste0("Post-menopausal<br>(N = ", scales::comma(post_N), ")"))
  ) |>
  cols_align(align = "left", columns = bmi_cat) |>
  cols_align(align = "center", columns = c(pre, post)) |>
  tab_options(
    table.font.size = px(13),
    heading.title.font.size = px(18),
    heading.subtitle.font.size = px(11),
    column_labels.font.weight = "bold",
    data_row.padding = px(8)
  ) |>
  opt_row_striping()

# Prevalence plot by menopause status and BMI category
p_aim1_menopause <- ggplot(
  prev_menopause,
  aes(x = bmi_cat, y = prevalence, color = menopause)
) +
  geom_point(position = position_dodge(width = 0.5), size = 2.5) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  labs(
    x = NULL,
    y = "% below choline AI",
    color = "Menopausal status",
    title = "Prevalence of choline intake below the AI among women"
  )

# Table 1: Baseline characteristics by sex

tbl_1 <- des |>
  tbl_svysummary(
    by = sex,
    include = c(
      RIDAGEYR,
      RIDRETH3,
      DMDEDUC2,
      bmi_cat,
      choline_mean,
      below_ai
    ),
    label = list(
      RIDAGEYR ~ "Age, years",
      RIDRETH3 ~ "Race/ethnicity",
      DMDEDUC2 ~ "Education",
      bmi_cat ~ "BMI category",
      choline_mean ~ "Choline intake, mg/day",
      below_ai ~ "Below choline AI"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n_unweighted} ({p}%)"
    ),
    digits = list(
      all_continuous() ~ 1,
      all_categorical() ~ c(0, 1)
    )
  ) |>
  add_overall(last = TRUE) |>
  add_p() |>
  modify_header(
    label ~ "**Characteristic**",
    stat_1 ~ paste0(
      "**Female, n = ",
      sum(ana$sex == "Female", na.rm = TRUE),
      "**"
    ),
    stat_2 ~ paste0(
      "**Male, n = ",
      sum(ana$sex == "Male", na.rm = TRUE),
      "**"
    ),
    stat_0 ~ paste0("**Overall, N = ", nrow(ana), "**")
  ) |>
  modify_spanning_header(c(stat_1, stat_2) ~ "**Sex**") |>
  modify_footnote(
    all_stat_cols() ~
      "Unweighted n (weighted %) for categorical variables; weighted mean (SD) for continuous variables."
  )

# Aim 1 tables
gt::gtsave(
  tbl_aim1,
  "output/aim1_table_sex_bmi.png",
  vwidth = 900,
  vheight = 450,
  zoom = 2
)

gt::gtsave(
  tbl_bmi,
  "output/aim1_table_bmi.png",
  vwidth = 700,
  vheight = 400,
  zoom = 2
)

gt::gtsave(
  tbl_aim1_f,
  "output/aim1_table_menopause_bmi.png",
  vwidth = 900,
  vheight = 450,
  zoom = 2
)

# Aim 1 plots
ggsave(
  "output/aim1_plot_sex_bmi.png",
  plot = p_aim1_sex_bmi,
  width = 8,
  height = 5,
  units = "in",
  dpi = 300
)

ggsave(
  "output/aim1_plot_bmi.png",
  plot = p_aim1_bmi,
  width = 7,
  height = 5,
  units = "in",
  dpi = 300
)

ggsave(
  "output/aim1_plot_menopause_bmi.png",
  plot = p_aim1_menopause,
  width = 8,
  height = 5,
  units = "in",
  dpi = 300
)

# Baseline table
gt::gtsave(
  gtsummary::as_gt(tbl_1),
  "output/table1_baseline_by_sex.png",
  vwidth = 1200,
  vheight = 1000,
  zoom = 2
)