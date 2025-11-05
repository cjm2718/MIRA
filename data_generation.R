require(tidyverse)
library(purrr)
library(scales)
library(gamlss)
library(Matrix)
library(stringr)
library(scales)
library(ggplot2)
library(zoo)
library(MASS)


###########################################################################################
# sim_bunion: primary simulation function for single arm of a bunionectomy trial
###########################################################################################
sim_bunion <- function(
    
  # number of subjects and time points for analysis 
  n_subjects = 50,
  time_points = seq(0, 48, by = 0.25),  # both scheduled and potential unscheduled observations every 15 min
  
  # trajectory parameters (higher decay rate gives faster initial decline in pain scores)
  baseline_mean = 7.0,
  endline_mean = 4.0,
  decay_rate = 0.08,
  
  # dose effect parameters (q6h dosing pattern)
  dose1_drop = 0.0,      # magnitude of initial drop from first dose (0-6h)
  dose2_drop = 0.0,      # magnitude of drop from second dose (6-12h)  
  dose_decay = 0.8,      # how fast each dose effect fades
  rebound_factor = 0.6,  # proportion of rebound before next dose (0=no rebound, 1=full rebound)
  
  # baseline eligibility
  min_baseline = 5.0,
  
  # between-subject variability parameters (slope SD smaller than observed data to reflect no rescue)
  baseline_sd = 1.2,
  slope_sd = 0.04,
  baseline_slope_cor = -0.3,
  
  # within-subject parameters (residual SD smaller than observed data to reflect no rescue)
  residual_sd = 1.4,
  ar1_rho = 0.6,
  
  # other options
  seed = NULL,
  treatment_label = "Placebo"
) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # function for mean trajectory with dosing effects
  mean_trajectory <- function(t) {
    # base exponential recovery (natural healing)
    base_recovery <- baseline_mean + (endline_mean - baseline_mean) * (1 - exp(-decay_rate * t))
    
    # add dosing effects for first 12 hours (q6h dosing)
    total_dose_effect <- 0
    
    if (t <= 12) {
      # first dose effect (0-6h): drop to maximum at ~1.5h, then rebound
      if (t <= 6) {
        if (t <= 1.5) {
          # drop phase: maximum drop at 1.5 hours
          dose1_effect <- -dose1_drop * (t/1.5)
        } else {
          # rebound phase: gradual recovery from 1.5h to 6h
          # when rebound_factor = 1, should return to 0 (no dose effect)
          # when rebound_factor = 0, should stay at -dose1_drop
          dose1_effect <- -dose1_drop * (1 - rebound_factor) + 
            (-dose1_drop * rebound_factor) * (1 - (t-1.5)/4.5)
        }
        total_dose_effect <- total_dose_effect + dose1_effect
      }
      
      # second dose effect (6-12h): similar pattern shifted by 6 hours
      if (t > 6 && t <= 12) {
        t_since_dose2 <- t - 6
        if (t_since_dose2 <= 1.5) {
          # drop phase: maximum drop at 7.5h (1.5h after dose 2)
          dose2_effect <- -dose2_drop * (t_since_dose2/1.5)
        } else {
          # rebound phase: gradual recovery from 7.5h to 12h  
          dose2_effect <- -dose2_drop * (1 - rebound_factor) + 
            (-dose2_drop * rebound_factor) * (1 - (t_since_dose2-1.5)/4.5)
        }
        total_dose_effect <- total_dose_effect + dose2_effect
      }
    }
    
    return(base_recovery + total_dose_effect)
  }
  
  # population mean trajectory
  mu_pop <- vapply(time_points, mean_trajectory, numeric(1))
  
  # covariance matrix for random effects
  Sigma_re <- matrix(c(baseline_sd^2,
                       baseline_slope_cor * baseline_sd * slope_sd,
                       baseline_slope_cor * baseline_sd * slope_sd,
                       slope_sd^2), 2, 2)
  
  # simulate individual subject
  simulate_subject <- function(subject_id) {
    # generate random effects
    re <- mvrnorm(1, mu = c(0, 0), Sigma = Sigma_re)
    random_intercept <- re[1]
    random_slope <- re[2]
    
    # subject-specific mean trajectory
    mu_subj <- mu_pop + random_intercept + random_slope * time_points
    
    # generate AR(1) correlation structure for errors
    errors <- numeric(length(time_points))
    errors[1] <- rnorm(1, 0, residual_sd)
    
    if (length(time_points) > 1) {
      for (i in 2:length(time_points)) {
        errors[i] <- ar1_rho * errors[i-1] + 
          rnorm(1, 0, residual_sd * sqrt(1 - ar1_rho^2))
      }
    }
    
    # generate pain score trajectory
    pain <- round(pmax(0, pmin(10, mu_subj + errors)))
    
    # apply baseline eligibility shift if needed
    baseline <- pain[1]
    if (baseline < min_baseline) {
      shift_amount <- min_baseline - baseline
      pain_final <- round(pmax(0, pmin(10, pain + shift_amount)))
    } else {
      pain_final <- pain
    }
    
    tibble(
      subject_id = paste0(treatment_label, "_", sprintf("%03d", subject_id)),
      treatment = treatment_label,
      time_hours = time_points,
      pain_score = pain_final
    )
  }
  
  # generate all subjects
  map_dfr(seq_len(n_subjects), simulate_subject)
}
