# libraries for parallelization
library(future)
library(furrr)
library(parallel)
conflicted::conflicts_prefer(dplyr::select, dplyr::filter, dplyr::lag)

# source all required functions
source("data_generation.R")
source("rescue_simulation.R") 
source("analysis_functions.R")
source("utils.R")
source("run_simulation.R")  

#########################################################################################
# Rescue parameters

# Rescue probability table
rescue_probs_bunion <- tibble(
  pain_score = 0:10,
  overall_prob = c(0.000416, 0.000708, 0.00304, 0.00578, 0.0177, 0.0381, 
                   0.0534, 0.0840, 0.122, 0.144, 0.221)
)

# Rescue effect data 
rescue_effect_data <- tribble(
  ~time_hr,      ~gamma_pred,        ~sd_hc,
  0,              0.000000000,      0,
  0.166666667,   -0.038009242,      0.67872,
  0.333333333,   -0.188237897,      1.11441,
  0.5,           -0.437308426,      1.49463,
  0.666666667,   -0.745206819,      1.76421,
  0.833333333,   -1.071561263,      1.94214,
  1,             -1.383915282,      2.05751,
  1.166666667,   -1.659628454,      2.17690,
  1.333333333,   -1.885169398,      2.18667,
  1.5,           -2.054442976,      2.27788,
  1.666666667,   -2.166917955,      2.33746,
  1.833333333,   -2.225902973,      2.31482,
  2,             -2.237107578,      2.30741,
  2.5,           -2.055563311,      2.24939,
  3,             -1.696395688,      2.29637,
  3.5,           -1.299967007,      2.24378,
  4,             -0.943573298,      2.27236,
  4.5,           -0.657087424,      2.15042,
  5,             -0.442869055,      2.13658,
  5.5,           -0.290698573,      2.04421,
  6,             -0.186692703,      2.17562
)

# Scheduled observation points
sched_points <- c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8,
                   10, 12, 14, 18, 22, 24, 26, 28, 30, 32, 34, 36, 38, 42, 46, 48)

##########################################################################################
# detect available cores
n_cores <- max(1, parallel::detectCores() - 2)

# set up parallel processing
plan(multisession, workers = n_cores)
options(future.rng.onMisuse = "ignore")

# full simulation study
full_results <- run_simulation_study(n_simulations = 5000, master_seed = 98765)

# create metadata
simulation_metadata <- list(
  master_seed = 98765,
  n_simulations = 5000,
  rng_kind = RNGkind()[1],
  r_version = R.version.string,
  package_versions = sessionInfo()$otherPkgs,
  n_cores = n_cores,
  simulation_date = Sys.time(),
  conditions = expand_grid(
    n_subjects = c(50, 100, 200),
    placebo_endline = 4,
    active_strength = c("null", "small", "medium", "large")
  ),
  note = "Rescue-dependent methods (observed, MIRA) run with all 3 rescue strengths; others use 'not_applicable'"
)

# save results
results_package <- list(
  contrasts = full_results$contrasts,
  group_estimates = full_results$group_estimates,
  rescue_summary = full_results$rescue_summary,
  metadata = simulation_metadata
)

saveRDS(results_package,
        file = paste0("simulation_results_", Sys.Date(), ".rds"))