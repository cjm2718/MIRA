####################################################################################################
# sim_rescue_timing: simulate rescue times based on specified probabilities
####################################################################################################
sim_rescue_timing <- function(counterfactual_data, rescue_prob_table, sched_points, rescue_window = 6, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # helper function to get rescue probability for a given pain score
  get_rescue_prob <- function(pain_score) {
    rescue_prob_table$overall_prob[rescue_prob_table$pain_score == pain_score]
  }
  
  # process each subject to determine rescue timing
  unique_subjects <- unique(counterfactual_data$subject_id)
  
  rescue_timing_data <- map_dfr(unique_subjects, function(subj_id) {
    
    subj_data <- filter(counterfactual_data, subject_id == subj_id) %>%
      arrange(time_hours)
    
    n_timepoints <- nrow(subj_data)
    
    # initialize tracking variables
    rescue_taken <- rep(0, n_timepoints)
    rescue_active <- rep(0, n_timepoints)
    last_rescue_time <- -Inf
    
    # sequential evaluation for rescue decisions (at 15-min intervals)
    for (i in seq_len(n_timepoints)) {
      current_time <- subj_data$time_hours[i]
      current_pain <- subj_data$pain_score[i]
      
      # check if enough time has passed since last rescue and minimum lockout from study start
      if (current_time >= 1.5 && current_time >= (last_rescue_time + rescue_window)) {
        rescue_prob <- get_rescue_prob(current_pain)
        
        if (runif(1) < rescue_prob) {
          rescue_taken[i] <- 1
          last_rescue_time <- current_time
        }
      }
      
      # determine if rescue is currently active 
      rescue_active[i] <- ifelse(current_time > last_rescue_time && 
                                   current_time <= (last_rescue_time + rescue_window), 1, 0)
    }
    
    # return rescue timing data for a subject
    tibble(
      subject_id = subj_id,
      time_hours = subj_data$time_hours,
      counterfactual_pain = subj_data$pain_score,
      rescue_taken = rescue_taken,
      rescue_active = rescue_active,
      scheduled = subj_data$time_hours %in% sched_points
    )
  })
  
  # generate individual rescue response characteristics
  individual_intercepts <- rnorm(length(unique_subjects), mean = 0, sd = 1)
  names(individual_intercepts) <- unique_subjects
  
  return(list(
    rescue_timing_data = rescue_timing_data,
    individual_intercepts = individual_intercepts
  ))
}

####################################################################################################
# sim_base_rescue_effects: simulate base rescue effects for a single instance of rescue-taking
####################################################################################################
sim_base_rescue_effects <- function(rescue_timing_data, rescue_effect_data, 
                                         individual_intercepts,
                                         between_subject_prop = 0.8,
                                         rescue_ar1_rho = 0.6, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  results <- map_dfr(names(individual_intercepts), function(subj_id) {
    
    subj_data <- filter(rescue_timing_data, subject_id == subj_id) %>%
      arrange(time_hours)
    
    n_timepoints <- nrow(subj_data)
    
    # use subject's individual rescue response
    individual_intercept <- individual_intercepts[subj_id]
    
    # initialize variables
    observed_pain <- subj_data$counterfactual_pain
    rescue_residuals <- numeric(n_timepoints)
    rescue_effects <- numeric(n_timepoints)  
    last_rescue_time <- -Inf
    
    # apply rescue effects based on existing timing
    for (i in seq_len(n_timepoints)) {
      current_time <- subj_data$time_hours[i]
      
      # update last rescue time if rescue was taken at this timepoint
      if (subj_data$rescue_taken[i] == 1) {
        last_rescue_time <- current_time
      }
      
      # calculate rescue effects only when rescue is active
      if (subj_data$rescue_active[i] == 1) {
        time_since_rescue_hr <- current_time - last_rescue_time
        
        # population effect from reference data
        population_effect <- approx(rescue_effect_data$time_hr, 
                                    rescue_effect_data$gamma_pred,
                                    xout = time_since_rescue_hr, 
                                    rule = 2)$y
        
        # target variability at this timepoint
        target_sd <- approx(rescue_effect_data$time_hr, 
                            rescue_effect_data$sd_hc,
                            xout = time_since_rescue_hr, 
                            rule = 2)$y
        
        # calculate individual effect
        individual_sd <- sqrt(between_subject_prop) * target_sd
        individual_effect <- individual_intercept * individual_sd
        
        # calculate residual effect with AR(1) correlation
        residual_sd <- sqrt(1 - between_subject_prop) * target_sd
        
        if (i == 1 || subj_data$rescue_active[i-1] == 0) {
          # start of new rescue episode
          rescue_residuals[i] <- rnorm(1, 0, residual_sd)
        } else {
          # continue existing rescue episode with AR(1) correlation
          rescue_residuals[i] <- rescue_ar1_rho * rescue_residuals[i-1] + 
            rnorm(1, 0, residual_sd * sqrt(1 - rescue_ar1_rho^2))
        }
        
        # combine all components for BASE rescue effect
        total_rescue_effect <- population_effect + individual_effect + rescue_residuals[i]
        
        # store the rescue effect (before applying to pain)
        rescue_effects[i] <- total_rescue_effect
        
        # apply to observed pain (rounded)
        observed_pain[i] <- pmax(0, pmin(10, subj_data$counterfactual_pain[i] + 
                                           round(total_rescue_effect)))
      }
    }
    
    # return only scheduled + rescue timepoints
    obs_points <- subj_data$time_hours %in% sched_points | subj_data$rescue_taken == 1
    
    tibble(
      subject_id = subj_id,
      time_hours = subj_data$time_hours[obs_points],
      counterfactual_pain = subj_data$counterfactual_pain[obs_points],
      rescue_taken = subj_data$rescue_taken[obs_points],
      rescue_active = subj_data$rescue_active[obs_points],
      observed_pain = observed_pain[obs_points],
      base_rescue_effect = rescue_effects[obs_points],  
      scheduled = subj_data$time_hours[obs_points] %in% sched_points
    ) %>%
      mutate(treatment = case_when(
        str_detect(subject_id, "placebo") ~ "placebo",
        str_detect(subject_id, "active_low") ~ "active_low", 
        str_detect(subject_id, "active_med") ~ "active_med",
        str_detect(subject_id, "active_high") ~ "active_high",
        str_detect(subject_id, "null") ~ "null"
      ))
  })
  
  return(results)
}

####################################################################################################
# apply_rescue_adjustment: apply adjustment factors to base rescue effects
####################################################################################################
apply_rescue_adjustment <- function(base_data, rescue_effect_strength = "predicted") {
  
  # Define adjustment factors
  adjustment_factor <- switch(rescue_effect_strength,
                              "predicted" = 1.0,
                              "underestimated" = 1.5,  # rescue 50% stronger than expected
                              "overestimated" = 0.5,   # rescue 50% weaker than expected
                              1.0
  )
  
  # Apply adjustment factor to rescue effects and recalculate observed pain
  adjusted_data <- base_data %>%
    mutate(
      # Scale the rescue effect by adjustment factor
      adjusted_rescue_effect = base_rescue_effect * adjustment_factor,
      # Recalculate observed pain using adjusted rescue effect
      observed_pain = pmax(0, pmin(10, counterfactual_pain + round(adjusted_rescue_effect)))
    ) %>%
    select(-base_rescue_effect, -adjusted_rescue_effect)  # Clean up intermediate columns
  
  return(adjusted_data)
}

#########################################################################################################
# sim_trial_data: wrapper function to generate simualted trial data (counterfactual, rescue and observed)
#########################################################################################################
sim_trial_data <- function(treatment_params_list, n_subjects, 
                                          between_subject_prop = 0.8,
                                          rescue_ar1_rho = 0.6) {
  
  # step 1: generate counterfactual trajectories
  counterfactual_data <- map_dfr(treatment_params_list, function(params) {
    do.call(sim_bunion, c(params, list(n_subjects = n_subjects, seed = NULL)))
  })
  
  # step 2: simulate rescue timing and generate individual rescue characteristics
  rescue_timing_result <- sim_rescue_timing(counterfactual_data, rescue_probs_bunion, sched_points)
  
  # step 3: generate base rescue effects as per reference distribution
  base_rescue_data <- sim_base_rescue_effects(rescue_timing_result$rescue_timing_data, 
                                                   rescue_effect_data, 
                                                   rescue_timing_result$individual_intercepts,
                                                   between_subject_prop,
                                                   rescue_ar1_rho)
  
  # step 4: apply different adjustment factors to the base effects
  rescue_strengths <- c("predicted", "underestimated", "overestimated")
  
  rescue_datasets <- map(rescue_strengths, function(strength) {
    apply_rescue_adjustment(base_rescue_data, strength)
  })
  names(rescue_datasets) <- rescue_strengths
  
  return(list(
    counterfactual_data = counterfactual_data,
    rescue_timing = rescue_timing_result$rescue_timing_data,
    individual_intercepts = rescue_timing_result$individual_intercepts,
    rescue_datasets = rescue_datasets
  ))
}