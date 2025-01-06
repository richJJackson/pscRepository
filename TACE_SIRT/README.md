Disease: Hepatocellular Carcinoma (HCC)
Treatments: Transarterial Chemoembolisation (TACE) and Selective Internal Radiation Therapy (SIRT)
Aim:  To compare overall survival between TACE and SIRT

---

# Process:

## Clean Datasets
1. Make sure that covariates are defined in the same way across both datasets (e.g. same column names and same units)
2. Define the primary outcome (Overall survival) for each dataset.  This is defined in days in the TACE dataset (turn into months) but needs to be calculated as the difference from the start of treatment until death (or censoring) in the SIRT dataset

## Summary Statistics
1. Summary table summsing data between TACE and SIRT data cohort
2. Kaplan Meier graphs looking at Overall survival by covariates
   
## TACE Model
1. Create on the TACE data using appropriate methodology (stepwise???) to the TACE 2 data
2. Validate using bootstrapping

## PSC comparison   
1. Compare SIRT data to TACE model for overall comparison
2. Perform sub-group analysis to see if treatment effect changes based on covariates
