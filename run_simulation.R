#########################################################################################################
# run_single_simulation: generates and analyzes a single trial for one set of conditions
#########################################################################################################
run_single_simulation <- function(condition, sim_id, master_seed) {
  
  # set simulation-specific seed
  sim_seed <- master_seed + sim_id * 1000 + 
    hash_condition(condition$n_subjects, condition$placebo_endline, condition$active_strength)
  
  set.seed(sim_seed)
  
  # get treatment parameters
  treatment_params <- get_treatment_params(condition$active_strength, condition$placebo_endline)
  
  # store all analysis results
  all_contrasts <- list()
  all_group_estimates <- list()
  
  # create base metadata for this trial
  base_trial_metadata <- list(
    simulation_id = sim_id,
    n_subjects = condition$n_subjects,
    placebo_endline = condition$placebo_endline,
    active_strength = condition$active_strength,
    scenario_type = ifelse(condition$active_strength == "null", "null", "alternative"),
    initial_rng_seed = sim_seed
  )
  
  # generates trial data
  trial_data_complete <- sim_trial_data(treatment_params, condition$n_subjects)
  
  # extract the different datasets
  rescue_datasets <- trial_data_complete$rescue_datasets
  
  # calculate rescue exposure 
  trial_data_base <- rescue_datasets[["predicted"]]
  rescue_data <- calculate_rescue_exposure(trial_data_base)
  rescue_summary <- rescue_data$by_arm %>%
    mutate(
      simulation_id = base_trial_metadata$simulation_id,
      n_subjects = base_trial_metadata$n_subjects,
      placebo_endline = base_trial_metadata$placebo_endline,
      active_strength = base_trial_metadata$active_strength,
      scenario_type = base_trial_metadata$scenario_type
    )
  
  # run rescue-independent analyses 
  rescue_independent_methods <- c("counterfactual", "locf", "bocf", "wocf", "root")
  
  for (method in rescue_independent_methods) {
    
    trial_metadata <- modifyList(base_trial_metadata, list(rescue_effect_strength = "not_applicable"))
    
    tryCatch({
      
      if (method == "root") {
        # ROOT analysis
        root_data <- calculate_root(trial_data_base, threshold = 0.5, rescue_window = 6)
        root_analysis <- analyze_root_zibeta(root_data, reference_group = "placebo")
        
        # add metadata to contrasts
        if (nrow(root_analysis$contrasts) > 0 && root_analysis$contrasts$converged[1]) {
          contrast_result <- root_analysis$contrasts %>%
            mutate(
              outcome_type = "root",
              analysis_method = "root",
              simulation_id = trial_metadata$simulation_id,
              n_subjects = trial_metadata$n_subjects,
              placebo_endline = trial_metadata$placebo_endline,
              active_strength = trial_metadata$active_strength,
              rescue_effect_strength = trial_metadata$rescue_effect_strength,
              scenario_type = trial_metadata$scenario_type,
              initial_rng_seed = trial_metadata$initial_rng_seed
            )
          all_contrasts[[paste0("root_", method)]] <- contrast_result
        }
        
        # add metadata to group estimates
        if (nrow(root_analysis$group_estimates) > 0 && root_analysis$group_estimates$converged[1]) {
          group_result <- root_analysis$group_estimates %>%
            mutate(
              outcome_type = "root",
              analysis_method = "root",
              simulation_id = trial_metadata$simulation_id,
              n_subjects = trial_metadata$n_subjects,
              placebo_endline = trial_metadata$placebo_endline,
              active_strength = trial_metadata$active_strength,
              rescue_effect_strength = trial_metadata$rescue_effect_strength,
              scenario_type = trial_metadata$scenario_type,
              initial_rng_seed = trial_metadata$initial_rng_seed
            )
          all_group_estimates[[paste0("root_", method)]] <- group_result
        }
        
      } else {
        # SPID analysis (for LOCF, WOCF, BOCF)
        spid_data <- calculate_spid(trial_data_base, analysis_type = method)
        spid_analysis <- analyze_spid_ancova(spid_data, reference_group = "placebo")
        
        # add metadata to contrasts
        if (nrow(spid_analysis$contrasts) > 0 && spid_analysis$contrasts$converged[1]) {
          contrast_result <- spid_analysis$contrasts %>%
            mutate(
              outcome_type = "spid",
              analysis_method = method,
              simulation_id = trial_metadata$simulation_id,
              n_subjects = trial_metadata$n_subjects,
              placebo_endline = trial_metadata$placebo_endline,
              active_strength = trial_metadata$active_strength,
              rescue_effect_strength = trial_metadata$rescue_effect_strength,
              scenario_type = trial_metadata$scenario_type,
              initial_rng_seed = trial_metadata$initial_rng_seed
            )
          all_contrasts[[paste0("spid_", method)]] <- contrast_result
        }
        
        # add metadata to group estimates  
        if (nrow(spid_analysis$group_estimates) > 0 && spid_analysis$group_estimates$converged[1]) {
          group_result <- spid_analysis$group_estimates %>%
            mutate(
              outcome_type = "spid",
              analysis_method = method,
              simulation_id = trial_metadata$simulation_id,
              n_subjects = trial_metadata$n_subjects,
              placebo_endline = trial_metadata$placebo_endline,
              active_strength = trial_metadata$active_strength,
              rescue_effect_strength = trial_metadata$rescue_effect_strength,
              scenario_type = trial_metadata$scenario_type,
              initial_rng_seed = trial_metadata$initial_rng_seed
            )
          all_group_estimates[[paste0("spid_", method)]] <- group_result
        }
      }
      
    }, error = function(e) {
      # record any failed analyses
      failed_contrast <- data.frame(
        contrast = "model_failed",
        estimate = NA_real_,
        se = NA_real_,
        lower_ci = NA_real_,
        upper_ci = NA_real_,
        p_value = NA_real_,
        significant = FALSE,
        analysis_type = method,
        converged = FALSE,
        outcome_type = ifelse(method == "root", "root", "spid"),
        analysis_method = method,
        simulation_id = trial_metadata$simulation_id,
        n_subjects = trial_metadata$n_subjects,
        placebo_endline = trial_metadata$placebo_endline,
        active_strength = trial_metadata$active_strength,
        rescue_effect_strength = trial_metadata$rescue_effect_strength,
        scenario_type = trial_metadata$scenario_type,
        initial_rng_seed = trial_metadata$initial_rng_seed,
        n_imputations = NA
      )
      
      failed_group <- data.frame(
        treatment_group = "model_failed",
        estimate = NA_real_,
        se = NA_real_,
        lower_ci = NA_real_,
        upper_ci = NA_real_,
        analysis_type = method,
        converged = FALSE,
        outcome_type = ifelse(method == "root", "root", "spid"),
        analysis_method = method,
        simulation_id = trial_metadata$simulation_id,
        n_subjects = trial_metadata$n_subjects,
        placebo_endline = trial_metadata$placebo_endline,
        active_strength = trial_metadata$active_strength,
        rescue_effect_strength = trial_metadata$rescue_effect_strength,
        scenario_type = trial_metadata$scenario_type,
        initial_rng_seed = trial_metadata$initial_rng_seed,
        n_imputations = NA
      )
      
      all_contrasts[[paste0("failed_", method)]] <- failed_contrast
      all_group_estimates[[paste0("failed_", method)]] <- failed_group
    })
  }
  
  # run rescue-dependent analyses for each rescue effect version
  rescue_dependent_methods <- c("observed", "mira")
  rescue_strengths <- c("predicted", "underestimated", "overestimated")
  
  for (rescue_strength in rescue_strengths) {
    
    # use the pre-generated dataset for this rescue strength
    trial_data_rescue <- rescue_datasets[[rescue_strength]]
    
    # run rescue-dependent analyses for this rescue strength
    for (method in rescue_dependent_methods) {
      
      trial_metadata <- modifyList(base_trial_metadata, list(rescue_effect_strength = rescue_strength))
      
      tryCatch({
        
        if (method == "mira") {
          # MIRA analysis
          spid_analysis <- analyze_mira(trial_data_rescue, hcapap_rescue_effect, n_imputations = 100)
        } else {
          # Observed analysis
          spid_data <- calculate_spid(trial_data_rescue, analysis_type = method)
          spid_analysis <- analyze_spid_ancova(spid_data, reference_group = "placebo")
        }
        
        # add metadata to contrasts
        if (nrow(spid_analysis$contrasts) > 0 && spid_analysis$contrasts$converged[1]) {
          contrast_result <- spid_analysis$contrasts %>%
            mutate(
              outcome_type = "spid",
              analysis_method = method,
              simulation_id = trial_metadata$simulation_id,
              n_subjects = trial_metadata$n_subjects,
              placebo_endline = trial_metadata$placebo_endline,
              active_strength = trial_metadata$active_strength,
              rescue_effect_strength = trial_metadata$rescue_effect_strength,
              scenario_type = trial_metadata$scenario_type,
              initial_rng_seed = trial_metadata$initial_rng_seed,
              n_imputations = ifelse(method == "mira" && "n_imputations" %in% names(spid_analysis$contrasts), 
                                     spid_analysis$contrasts$n_imputations[1], NA)
            )
          all_contrasts[[paste0("spid_", method, "_", rescue_strength)]] <- contrast_result
        }
        
        # add metadata to group estimates  
        if (nrow(spid_analysis$group_estimates) > 0 && spid_analysis$group_estimates$converged[1]) {
          group_result <- spid_analysis$group_estimates %>%
            mutate(
              outcome_type = "spid",
              analysis_method = method,
              simulation_id = trial_metadata$simulation_id,
              n_subjects = trial_metadata$n_subjects,
              placebo_endline = trial_metadata$placebo_endline,
              active_strength = trial_metadata$active_strength,
              rescue_effect_strength = trial_metadata$rescue_effect_strength,
              scenario_type = trial_metadata$scenario_type,
              initial_rng_seed = trial_metadata$initial_rng_seed,
              n_imputations = ifelse(method == "mira" && "n_imputations" %in% names(spid_analysis$group_estimates), 
                                     spid_analysis$group_estimates$n_imputations[1], NA)
            )
          all_group_estimates[[paste0("spid_", method, "_", rescue_strength)]] <- group_result
        }
        
      }, error = function(e) {
        # record any failed analyses
        failed_contrast <- data.frame(
          contrast = "model_failed",
          estimate = NA_real_,
          se = NA_real_,
          lower_ci = NA_real_,
          upper_ci = NA_real_,
          p_value = NA_real_,
          significant = FALSE,
          analysis_type = method,
          converged = FALSE,
          outcome_type = "spid",
          analysis_method = method,
          simulation_id = trial_metadata$simulation_id,
          n_subjects = trial_metadata$n_subjects,
          placebo_endline = trial_metadata$placebo_endline,
          active_strength = trial_metadata$active_strength,
          rescue_effect_strength = trial_metadata$rescue_effect_strength,
          scenario_type = trial_metadata$scenario_type,
          initial_rng_seed = trial_metadata$initial_rng_seed,
          n_imputations = NA
        )
        
        failed_group <- data.frame(
          treatment_group = "model_failed",
          estimate = NA_real_,
          se = NA_real_,
          lower_ci = NA_real_,
          upper_ci = NA_real_,
          analysis_type = method,
          converged = FALSE,
          outcome_type = "spid",
          analysis_method = method,
          simulation_id = trial_metadata$simulation_id,
          n_subjects = trial_metadata$n_subjects,
          placebo_endline = trial_metadata$placebo_endline,
          active_strength = trial_metadata$active_strength,
          rescue_effect_strength = trial_metadata$rescue_effect_strength,
          scenario_type = trial_metadata$scenario_type,
          initial_rng_seed = trial_metadata$initial_rng_seed,
          n_imputations = NA
        )
        
        all_contrasts[[paste0("failed_", method, "_", rescue_strength)]] <- failed_contrast
        all_group_estimates[[paste0("failed_", method, "_", rescue_strength)]] <- failed_group
      })
    }
  }
  
  # Combine all results
  trial_contrasts <- bind_rows(all_contrasts)
  trial_group_estimates <- bind_rows(all_group_estimates)
  
  return(list(
    contrasts = trial_contrasts,
    group_estimates = trial_group_estimates,
    rescue_summary = rescue_summary,
    metadata = base_trial_metadata
  ))
}

#########################################################################################################
# hash_condition: creates deterministic seed for each simulation
#########################################################################################################
hash_condition <- function(n_subjects, placebo_endline, active_strength) {
  sum(utf8ToInt(paste0(n_subjects, placebo_endline, active_strength))) %% 1000
}

#########################################################################################################
# run_simulation_study: runs simulation study with parallelization
#########################################################################################################
run_simulation_study <- function(n_simulations = 2000, master_seed = 12345) {
  
  # define condition combinations
  conditions <- expand_grid(
    n_subjects = c(50, 100, 200),
    placebo_endline = c(4,6),
    active_strength = c("null", "small", "medium", "large")
  )
  
  n_conditions <- nrow(conditions)
  cat("Running", n_conditions, "conditions with", n_simulations, "simulations each\n")
  cat("Total simulations:", n_conditions * n_simulations, "\n")
  cat("Note: Rescue-dependent methods (observed, MIRA) will be run with 3 rescue strengths each\n")
  
  start_time <- Sys.time()
  
  # set up parallel processing
  plan(multisession, workers = n_cores)
  
  # run simulations in parallel at the condition level
  all_results <- future_map(1:n_conditions, function(condition_idx) {
    
    condition <- conditions[condition_idx, ]
    
    # run simulations for the condition sequentially to maintain seed control
    condition_results <- map(1:n_simulations, function(sim_id) {
      run_single_simulation(condition, sim_id, master_seed)
    })
    
    # combine results for the condition
    condition_contrasts <- map_dfr(condition_results, "contrasts")
    condition_group_estimates <- map_dfr(condition_results, "group_estimates")
    condition_rescue_summary <- map_dfr(condition_results, "rescue_summary")
    
    return(list(
      contrasts = condition_contrasts,
      group_estimates = condition_group_estimates,
      rescue_summary = condition_rescue_summary  
    ))
    
  }, .options = furrr_options(
    seed = TRUE, 
    packages = c("MASS", "gamlss", "tidyverse", "stringr", "Matrix")
  ))
  
  # combine results across all conditions
  combined_contrasts <- map_dfr(all_results, "contrasts")
  combined_group_estimates <- map_dfr(all_results, "group_estimates")
  combined_rescue_summary <- map_dfr(all_results, "rescue_summary")  
  
  total_time <- difftime(Sys.time(), start_time, units = "hours")
  cat("\nSimulation study completed in", round(total_time, 2), "hours\n")
  
  return(list(
    contrasts = combined_contrasts,
    group_estimates = combined_group_estimates,
    rescue_summary = combined_rescue_summary
  ))
}