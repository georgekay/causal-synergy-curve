# NHANES 8-SDOH Analytic File Codebook

File: `nhanes_8sdoh_analytic.csv`

Derived by the NHANES 8-SDOH preprocessing pipeline for the Bayesian PH applied illustration.

Source: `nhanes_harmonized_public_use_extract_source.csv`, a harmonized public-use NHANES 1999-2018 extract with linked mortality follow-up fields, survey-design variables, and pre-derived Bundy-style SDOH indicators.

## Inclusion criteria

- Non-Hispanic Black or non-Hispanic White adults
- Age 20-74 years at NHANES examination
- Linked mortality follow-up time available
- Complete data on all eight SDOH indicators
- Non-missing sex and pooled MEC weight

## Endpoint

- `followup_years`: linked mortality follow-up time in years, based on `PERMTH_EXM` with `PERMTH_INT` fallback
- `death_10yr`: all-cause death within 10 years of follow-up

## SDOH domain coding

All SDOH domains are binary and coded so that:

- `0` = favorable benchmark
- `1` = unfavorable or constrained state

| Manuscript domain | Analytic column | Source column | Favorable benchmark |
|---|---|---|---|
| Employment | `D_employment` | `sdoh_unemployed` | employed or not in the unfavorable employment category |
| Family income | `D_income` | `sdoh_low_income` | not low income |
| Food security | `D_food_security` | `sdoh_food_insecure` | food secure |
| Education | `D_education` | `sdoh_less_than_high_school` | high-school education or higher |
| Access to care | `D_access_to_care` | `sdoh_no_regular_care` | has a regular source of care |
| Health insurance | `D_health_insurance` | `sdoh_no_private_insurance` | has private insurance |
| Housing instability | `D_housing_instability` | `sdoh_rent_or_other` | owns home or lives in the favorable housing category |
| Partnership | `D_partnership` | `sdoh_not_partnered` | married or living with a partner |

## Survey and covariate fields

- `analysis_race`: Black or White analytic race group
- `analysis_age`: age in years
- `sex_female`: female sex indicator
- `cycle`: NHANES survey cycle
- `pooled_mec_weight`: pooled MEC examination weight used for standardization
- `SDMVPSU`, `SDMVSTRA`: NHANES survey design variables retained for design-sensitivity extensions

## Group architecture definitions

The manuscript reports an eight-domain architecture and a three-group architecture. The three-group architecture is a prespecified stability layer:

- Economic: employment, family income, housing instability
- Healthcare/access: access to care, health insurance
- Education/social: food security, education, partnership

The grouping keeps the eight-domain results visible while providing a lower-variance architecture summary for race-stratified comparison.
