## Rows are NOT dropped. Exclusions are stored as flags so the survey design
## can be declared on the full sample

## load packages
library(tidyverse)

## read merged file (from 01_fetch_nhanes.R)
df <- readRDS("data/nhanes_L_merged_raw.rds")

## coerce factor/character columns to plain character for predictable matching
chr <- function(x) as.character(x)

## sex-specific Adequate Intake (AI) for choline, mg/day
ai_men   <- 550
ai_women <- 425

df <- df |>
  mutate(
    
    ## ---- sex ----
    sex = case_when(
      chr(RIAGENDR) == "Male"   ~ "male",
      chr(RIAGENDR) == "Female" ~ "female",
      TRUE                      ~ NA_character_
    ),
    
    ## ---- BMI category ----
    bmi_cat = case_when(
      BMXBMI < 18.5                ~ "underweight",
      BMXBMI >= 18.5 & BMXBMI < 25 ~ "normal",
      BMXBMI >= 25   & BMXBMI < 30 ~ "overweight",
      BMXBMI >= 30                 ~ "obese",
      TRUE                         ~ NA_character_
    ),
    bmi_cat = factor(bmi_cat,
                     levels = c("normal", "overweight", "obese", "underweight")),
    
    ## ---- valid dietary recalls (NHANES reliability flags) ----
    day1_ok = chr(DR1DRSTZ) == "Reliable and met the minimum criteria",
    day2_ok = chr(DR2DRSTZ) == "Reliable and met the minimum criteria",
    
    ## ---- choline: mean of the two recall days (requires both) ----
    choline_mean = rowMeans(cbind(DR1TCHL, DR2TCHL), na.rm = FALSE),
    
    ## ---- AI threshold and derived choline measures ----
    ai_threshold = case_when(
      sex == "male"   ~ ai_men,
      sex == "female" ~ ai_women,
      TRUE            ~ NA_real_
    ),
    below_ai       = as.numeric(choline_mean < ai_threshold),  # 1 = below AI
    choline_pct_ai = 100 * choline_mean / ai_threshold,
    
    ## ---- energy: mean of the two days (plausibility screen) ----
    energy_mean = rowMeans(cbind(DR1TKCAL, DR2TKCAL), na.rm = FALSE),
    
    ## ---- menopausal status (women only) ----
    menopause = case_when(
      sex == "male"                                        ~ NA_character_,
      chr(RHD043) == "Menopause/Change of life"            ~ "post",
      chr(RHD043) == "Hysterectomy"                        ~ "hysterectomy",
      chr(RHQ031) == "Yes"                                 ~ "pre",
      chr(RHD043) %in% c("Pregnancy", "Breast feeding")    ~ "pre",
      TRUE                                                 ~ NA_character_
    ),
    
    ## ---- exclusion components ----
    ## implausible energy intake (NA stays NA, handled in the flags below)
    excl_energy = case_when(
      is.na(energy_mean) | is.na(sex)                                  ~ NA,
      sex == "male"   & (energy_mean < 500 | energy_mean > 8000)       ~ TRUE,
      sex == "female" & (energy_mean < 500 | energy_mean > 5500)       ~ TRUE,
      TRUE                                                             ~ FALSE
    ),
    
    ## heavy alcohol: ALQ130 = avg drinks/day on drinking days.
    ## 777 = Refused, 999 = Don't know -> set to NA before thresholding.
    alq130_num = suppressWarnings(as.numeric(chr(ALQ130))),
    alq130_num = if_else(alq130_num %in% c(777, 999), NA_real_, alq130_num),
    heavy_alcohol = case_when(
      is.na(alq130_num) | is.na(sex)      ~ FALSE,   # no evidence of heavy use
      sex == "male"   & alq130_num > 2    ~ TRUE,
      sex == "female" & alq130_num > 1    ~ TRUE,
      TRUE                                 ~ FALSE
    ),
    
    ## self-reported liver condition (MCQ160L: ever told you had a liver condition)
    liver_disease = chr(MCQ160L) == "Yes",
    
    ## self-reported hepatitis B (HEQ010).
    ## NOTE: hep C self-report (HEQ030) was dropped this cycle; add lab-confirmed
    ## hepatitis markers before finalizing the Aim 2 exclusions.
    hepatitis = chr(HEQ010) == "Yes"
  )

## ---- analytic sample flags (no rows dropped) ----
df <- df |>
  mutate(
    
    ## Aim 1: adults 20+, two reliable recalls, measured BMI, positive 2-day
    ## weight, a derivable AI comparison, and plausible energy intake
    in_aim1 = RIDAGEYR >= 20 &
      day1_ok %in% TRUE & day2_ok %in% TRUE &
      !is.na(WTDR2D) & WTDR2D > 0 &
      !is.na(BMXBMI) &
      !is.na(below_ai) &
      excl_energy %in% FALSE,
    
    ## Aim 2: Aim 1 plus a measured ALT and the liver-related exclusions
    in_aim2 = in_aim1 &
      !is.na(LBXSATSI) &
      !(heavy_alcohol %in% TRUE) &
      !(liver_disease %in% TRUE) &
      !(hepatitis %in% TRUE)
  )

## ---- verification: fail loudly if a recode silently produced all-NA ----
stopifnot(
  "sex is all NA (check RIAGENDR labels)"       = any(!is.na(df$sex)),
  "bmi_cat is all NA"                           = any(!is.na(df$bmi_cat)),
  "ai_threshold is all NA (depends on sex)"     = any(!is.na(df$ai_threshold)),
  "below_ai is all NA"                          = any(!is.na(df$below_ai)),
  "menopause is all NA (check RHD043 labels)"   = any(!is.na(df$menopause)),
  "heavy_alcohol never TRUE (check ALQ130)"     = any(df$heavy_alcohol, na.rm = TRUE),
  "no one in Aim 1 sample"                      = sum(df$in_aim1, na.rm = TRUE) > 0
)

## ---- derivation checks ----
cat("\n--- sex ---\n");           print(table(df$sex, useNA = "ifany"))
cat("\n--- bmi_cat (adults) ---\n")
print(table(df$bmi_cat[df$RIDAGEYR >= 20], useNA = "ifany"))
cat("\n--- below_ai (Aim 1 sample) ---\n")
print(table(df$below_ai[df$in_aim1], useNA = "ifany"))
cat("\n--- menopause ---\n");     print(table(df$menopause, useNA = "ifany"))
cat("\n--- heavy_alcohol ---\n"); print(table(df$heavy_alcohol, useNA = "ifany"))

## ---- exclusion funnel ----
cat("\n--- exclusion funnel ---\n")
cat("Full merged sample:      ", nrow(df), "\n")
cat("Adults 20+:              ", sum(df$RIDAGEYR >= 20, na.rm = TRUE), "\n")
cat("  + two reliable recalls:",
    sum(df$RIDAGEYR >= 20 & df$day1_ok %in% TRUE & df$day2_ok %in% TRUE,
        na.rm = TRUE), "\n")
cat("In Aim 1 sample:         ", sum(df$in_aim1, na.rm = TRUE), "\n")
cat("In Aim 2 sample:         ", sum(df$in_aim2, na.rm = TRUE), "\n")

## ---- save analytic dataset (full rows, with flags + derived vars) ----
saveRDS(df, "data/nhanes_L_analytic.rds")
write_csv(df, "data/nhanes_L_analytic.csv")