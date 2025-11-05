# MIRA: Multiple Imputation for Rescue Analgesia

A comprehensive simulation study evaluating different statistical methods for handling rescue medication use in acute pain clinical trials, with focus on bunionectomy as the clinical model. This project introduces a reference-based multiple imputation approach to handling rescue medication.

## Project Overview

This simulation study compares various statistical methods for analyzing pain clinical trials where participants may take rescue medication. The study uses bunionectomy pain trials as the clinical model and implements reference-based multiple imputation techniques.

## Core Simulation Files

- **`data_generation.R`**: Functions for generating counterfactual acute pain trial data

- **`rescue_simulation.R`**: Core rescue medication modeling functions

- **`analysis_functions.R`**: Statistical analysis methods

- **`utils.R`**: Utility and helper functions

- **`run_simulation.R`**: Main simulation execution script

- **`sim_study.R`**: High-level simulation study design

## Output Files

- **Simulation results** (`.rds` files): Generated locally when running simulations
  - Not tracked in Git due to large file sizes

## Project Configuration

- **`MIRA.Rproj`**: RStudio project file
- **`.gitignore`**: Git configuration to exclude temporary files and simulation results due to file size (`.rds` files)

## Rescue Medication Handling Approaches

1. **Multiple Imputation for Rescue Analgesia (MIRA)**: Reference-based imputation 
2. **LOCF (Last Observation Carried Forward)**: Most common method carrying forward last pre-rescue observation
3. **BOCF (Baseline Observation Carried Forward)**: Conservative approach using baseline values
4. **WOCF (Worst Observation Carried Forward)**: Worst-case imputation strategy
5. **ROOT**: Responder Outcome Over Time (composite strategy)
6. **Observed**: No imputation (treatment policy strategy)


## Running the Simulation
The simulation parameters are set in sim_study.R.

## Other Notes
Simulation results are saved as .rds files locally. These files are not tracked in Git due to their large size.
