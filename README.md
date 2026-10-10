# Choline Intake Adequacy and Liver Enzymes by Body Weight in US Adults

This respository contains analysis code for a cross-sectional study of dietary choline intake, body mass index (BMI), and liver enzymes in US adults, using NHANES August 2021 to August 2023.

## Study Aims

This study has two aims:

1.  Estimate the prevalence of choline intake less than the AI for each BMI category, based on gender and menopausal status.
2.  Is the association between choline intake and liver enzymes (ALT, AST, and GGT) modified by obesity?

## Data and Variables

All data come from the National Health and Nutrition Examination Survey (NHANES), from the National Center for Health Statistics.

NHANES documentation: <https://wwwn.cdc.gov/nchs/nhanes>

| File               | Contents                                 |
|--------------------|------------------------------------------|
| DEMO_L             | Demographics and survey design variables |
| DR1TOT_L, DR2TOT_L | Day 1 and Day 2 dietary recall totals    |
| BMX_L              | Body measures (BMI)                      |
| BIOPRO_L           | Standard biochemistry profilex           |
| RHQ_L              | Reproductive health (menopausal status)  |
| ALQ_L              | Alcohol use                              |
| MCQ_L              | Medical conditions                       |
| HEQ_L              | Hepatitis status                         |

## Scripts and Reproducibility

Run the scripts in order:

| Script | Function |
|----|----|
| `01_fetch_nhanes.R` | Downloads all the nhanes datasets and combines them by study participant id. |
| `02_clean_data.R` | Creates new variables needed for the analysis |
| `03_aim_1_analysis.R` | Analysis for Aim 1 estimating % less than adequate intake of choline by BMI category and stratified by sex/gender and menopause status. |
| `04_aim2_analysis.R` | Analysis for Aim 2 models for liver enzymes and contrasts for main effect of choline, and secondary analysis. |

## R Version and Packages used

- R 4.3.3
- Packages:

``` r
install.packages(c("nhanesA", "tidyverse", "survey", "gt", "gtsummary"))
```

## Folder Structure

```         
choline-bmi-nhanes/
├── code/        # code used to conduct the analysis      
├── data/        # raw and cleaned data files
├── output/      # figures and tables produced from the code
└── docs/        # supporting files including manuscript
```

## License

This document is placed in the public domain under the terms of [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/). You may copy, adapt, and distribute this document as you wish. See the `LICENSE` file for more information. The NHANES data are in the public domain.
