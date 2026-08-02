# load packages
library(tidyverse)

# read merged file (from 01_fetch_nhanes.R)
df <- readRDS("data/nhanes_L_merged_raw.rds")

# sex-specific Adequate Intake (AI) for choline, mg/day
ai_men   <- 550
ai_women <- 425

df <- df |>
  mutate(
    #sex
    sex = case_when(
      RIAGENDR == 1 ~ "male",
      RIAGENDR == 2 ~ "female",
      TRUE          ~ NA_character_
    ),
    
    #BMI category (Aim 1)
    bmi_cat = case_when(
      BMXBMI < 18.5              ~ "underweight",
      BMXBMI >= 18.5 & BMXBMI < 25 ~ "normal",
      BMXBMI >= 25   & BMXBMI < 30 ~ "overweight",
      BMXBMI >= 30              ~ "obese",
      TRUE                       ~ NA_character_
    ),
    bmi_cat = factor(bmi_cat,
                     levels = c("normal", "overweight", "obese", "underweight")),
    # average choline across the two recall days
    choline_mean = rowMeans(cbind(DR1TCHL, DR2TCHL), na.rm = FALSE),
    #sex-specific AI threshold and derived choline measures
    ai_threshold = if_else(sex == "male", ai_men, ai_women),
    below_ai     = choline_mean < ai_threshold,          # Aim 1 outcome
    choline_pct_ai = 100 * choline_mean / ai_threshold,  # Aim 2 exposure
    # average energy across days (for plausibility screen)
    energy_mean = rowMeans(cbind(DR1TKCAL, DR2TKCAL), na.rm = FALSE),
    # menopausal status (women only)
    # RHD043: reason for no period in past 12 mo (7 = natural menopause, 3 = hysterectomy)
    # RHQ031: still having regular periods
    menopause = case_when(
      sex == "male"    ~ NA_character_,
      RHD043 == 7      ~ "post",
      RHD043 == 3      ~ "hysterectomy",   # set aside in primary analysis
      RHQ031 == 1      ~ "pre",
      RHD043 %in% c(1, 2) ~ "pre",
      TRUE             ~ NA_character_
    ),
    
    # exclusion components
    excl_energy = (sex == "male"   & (energy_mean < 500 | energy_mean > 8000)) |
      (sex == "female" & (energy_mean < 500 | energy_mean > 5500)),
    
    # heavy alcohol (ALQ130 = avg drinks/day on drinking days; ALQ121 = frequency)
    heavy_alcohol = case_when(
      sex == "male"   & ALQ130 > 2 ~ TRUE,
      sex == "female" & ALQ130 > 1 ~ TRUE,
      TRUE                          ~ FALSE
    ),
    
    # self-reported liver condition (MCQ160L ever, MCQ170L still)
    liver_disease = MCQ160L == 1,
    # self-reported hepatitis B (HEQ010). NOTE: hep C self-report dropped this cycle;
    # confirm/replace with lab markers when those files are added.
    hepatitis = HEQ010 == 1
  )

# analytic sample flags (no rows dropped)
df <- df |>
  mutate(
    # Aim 1: adults 20+, two valid recalls, measured BMI, valid 2-day weight
    in_aim1 = RIDAGEYR >= 20 &
      !is.na(WTDR2D) & WTDR2D > 0 &
      !is.na(BMXBMI) &
      !is.na(DR1TCHL) & !is.na(DR2TCHL) &
      !excl_energy %in% TRUE,
    #Aim 2: Aim 1 plus non-missing ALT and the liver-related exclusions
    in_aim2 = in_aim1 &
      !is.na(LBXSATSI) &
      !(heavy_alcohol %in% TRUE) &
      !(liver_disease %in% TRUE) &
      !(hepatitis %in% TRUE)
  )

##quick exclusion funnel (for the participant flow diagram)
cat("Full merged sample:      ", nrow(df), "\n")
cat("Adults 20+:              ", sum(df$RIDAGEYR >= 20, na.rm = TRUE), "\n")
cat("In Aim 1 sample:         ", sum(df$in_aim1, na.rm = TRUE), "\n")
cat("In Aim 2 sample:         ", sum(df$in_aim2, na.rm = TRUE), "\n")

# save analytic dataset (still full rows, with flags + derived vars)
saveRDS(df, "data/nhanes_L_analytic.rds")
write_csv(df, "data/nhanes_L_analytic.csv")