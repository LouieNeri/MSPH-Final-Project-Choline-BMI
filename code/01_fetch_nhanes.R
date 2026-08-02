## load packages
library(nhanesA)
library(tidyverse)

## pull each file and keep only needed variables (all "_L" = Aug 2021-Aug 2023 cycle)

# demographics + survey design vars
demo <- nhanes("DEMO_L") |> 
  select(SEQN, RIAGENDR, RIDAGEYR, RIDRETH3, DMDEDUC2,
         INDFMPIR, SDMVPSU, SDMVSTRA, WTINT2YR, WTMEC2YR)

# day 1 nutrients (choline, energy, weights)
diet1 <- nhanes("DR1TOT_L") |>        
  select(SEQN, WTDRD1, WTDR2D, DR1TKCAL, DR1TCHL, DR1DRSTZ, DR1DAY)

# day 2 nutrients (choline)
diet2 <- nhanes("DR2TOT_L") |>        
  select(SEQN, DR2TKCAL, DR2TCHL, DR2DRSTZ)

# body measures (BMI)
bmx <- nhanes("BMX_L") |>
  select(SEQN, BMXBMI)

# biochemistry (ALT, AST, GGT)
biopro <- nhanes("BIOPRO_L") |>       
  select(SEQN, LBXSATSI, LBXSASSI, LBXSGTSI)

# reproductive health (menopause)
rhq <- nhanes("RHQ_L") |>             
  select(SEQN, RHQ031, RHD043)

# alcohol (aim 2 exclusion)
alq <- nhanes("ALQ_L") |>              
  select(SEQN, ALQ121, ALQ130)

# medical conditions (liver disease exclusion)
mcq <- nhanes("MCQ_L") |>              
  select(SEQN, MCQ160L, MCQ170L)

# hepatitis B self-report (aim 2 exclusion)
heq <- nhanes("HEQ_L") |>           
  select(SEQN, HEQ010)

# merge all files on SEQN
df <- list(demo, diet1, diet2, bmx, biopro, rhq, alq, mcq, heq) |>
  reduce(left_join, by = "SEQN")

#checking merged data
nrow(df)                               # full sample
summary(df$DR1TCHL)                    # choline looks reasonable?
summary(df$LBXSATSI)                   # ALT looks reasonable?

## save merged raw data
saveRDS(df, "nhanes_L_merged_raw.rds")
write_csv(df, "nhanes_L_merged_raw.csv")   # csv archive copy