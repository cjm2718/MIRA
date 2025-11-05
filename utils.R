#########################################################################################
# treatment parameter definitions 
#########################################################################################
# Default parameters (placebo endline = 4)
placebo_params_default <- list(
  baseline_mean = 7.0,
  endline_mean = 4.0,
  decay_rate = 0.08,
  dose1_drop = 0.5,      
  dose2_drop = 0.25,      
  rebound_factor = 1,
  treatment_label = "placebo"
)

active_low_params_default <- list(
  baseline_mean = 7.0,
  endline_mean = 3.5,
  decay_rate = 0.08,
  dose1_drop = 1,    
  dose2_drop = 0.5,   
  rebound_factor = 1,
  treatment_label = "active_low"
)

active_med_params_default <- list(
  baseline_mean = 7.0,
  endline_mean = 3.0,
  decay_rate = 0.10,
  dose1_drop = 1.5,    
  dose2_drop = 0.75,    
  rebound_factor = 1, 
  treatment_label = "active_med"
)

active_high_params_default <- list(
  baseline_mean = 7.0,
  endline_mean = 2.5,
  decay_rate = 0.12,
  dose1_drop = 2,    
  dose2_drop = 1,    
  rebound_factor = 1, 
  treatment_label = "active_high"
)

# Null parameters (identical to placebo for Type I error testing)
null_params_default <- list(
  baseline_mean = 7.0,
  endline_mean = 4.0,
  decay_rate = 0.08,
  dose1_drop = 0.5,      
  dose2_drop = 0.25,      
  rebound_factor = 1,
  treatment_label = "null"
)

# adjusted parameters (placebo endline = 6, all shifted up by 2)
placebo_params_adjusted <- modifyList(placebo_params_default, list(endline_mean = 6.0))
active_low_params_adjusted <- modifyList(active_low_params_default, list(endline_mean = 5.5))
active_med_params_adjusted <- modifyList(active_med_params_default, list(endline_mean = 5.0))
active_high_params_adjusted <- modifyList(active_high_params_default, list(endline_mean = 4.5))
null_params_adjusted <- modifyList(null_params_default, list(endline_mean = 6.0))

# function to get treatment parameters based on strength and placebo endline
get_treatment_params <- function(active_strength, placebo_endline) {
  
  if (placebo_endline == 4) {
    placebo_params <- placebo_params_default
    active_low_params <- active_low_params_default
    active_med_params <- active_med_params_default
    active_high_params <- active_high_params_default
    null_params <- null_params_default
  } else {
    placebo_params <- placebo_params_adjusted
    active_low_params <- active_low_params_adjusted
    active_med_params <- active_med_params_adjusted
    active_high_params <- active_high_params_adjusted
    null_params <- null_params_adjusted
  }
  
  treatment_params <- list(placebo_params)
  
  if (active_strength == "null") {
    treatment_params <- c(treatment_params, list(null_params))
  } else if (active_strength == "small") {
    treatment_params <- c(treatment_params, list(active_low_params))
  } else if (active_strength == "medium") {
    treatment_params <- c(treatment_params, list(active_med_params))
  } else if (active_strength == "large") {
    treatment_params <- c(treatment_params, list(active_high_params))
  }
  
  return(treatment_params)
}

#####################################################################################################################
# calculate_rescue_exposure: calculates the exposure to rescue at the subject level and group level
#####################################################################################################################
calculate_rescue_exposure <- function(trial_data) {
  
  # calculate rescue exposure for each subject
  subject_exposure <- trial_data %>%
    arrange(subject_id, time_hours) %>%
    group_by(subject_id, treatment) %>%
    summarise(
      total_study_time = max(time_hours) - min(time_hours),
      rescue_active_time = sum(rescue_active * c(diff(time_hours), 0)),  
      rescue_exposure_pct = rescue_active_time / total_study_time * 100,
      n_rescue_events = sum(rescue_taken),
      .groups = "drop"
    )
  
  # calculate treatment arm summaries
  arm_summaries <- subject_exposure %>%
    group_by(treatment) %>%
    summarise(
      n_subjects = n(),
      mean_rescue_exposure_pct = mean(rescue_exposure_pct),
      sd_rescue_exposure_pct = sd(rescue_exposure_pct),
      median_rescue_exposure_pct = median(rescue_exposure_pct),
      mean_rescue_events = mean(n_rescue_events),
      sd_rescue_events = sd(n_rescue_events),
      prop_subjects_with_rescue = mean(n_rescue_events > 0),
      .groups = "drop"
    )
  
  return(list(
    by_subject = subject_exposure,
    by_arm = arm_summaries
  ))
}
