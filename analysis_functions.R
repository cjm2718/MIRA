############################################################################################################################
# calculate_spid: calculates SPID and baseline pain scores for analysis
############################################################################################################################
calculate_spid <- function(trial_data, analysis_type = "counterfactual", rescue_window = 6) {
  
  if (analysis_type == "counterfactual") {
    # counterfactual uses only scheduled timepoints since no rescue is used
    analysis_data <- trial_data %>% 
      filter(scheduled == TRUE) %>%
      arrange(subject_id, time_hours) %>%
      mutate(pain_for_analysis = counterfactual_pain)
    
  } else {
    # all other methods use both scheduled AND unscheduled rescue timepoints
    analysis_data <- trial_data %>% 
      filter(scheduled == TRUE | rescue_taken == 1) %>% 
      arrange(subject_id, time_hours)
    
    if (analysis_type == "observed_only") {
      # use observed data with no imputation
      analysis_data <- analysis_data %>%
        mutate(pain_for_analysis = observed_pain)
      
    } else {
      # apply imputation methods (locf, bocf, wocf)
      analysis_data <- analysis_data %>%
        group_by(subject_id) %>%
        arrange(time_hours) %>%
        mutate(
          # identify rescue episodes and create imputation windows
          rescue_times = ifelse(rescue_taken == 1, time_hours, NA),
          last_rescue_time = zoo::na.fill(zoo::na.locf(rescue_times, na.rm = FALSE), -Inf),
          in_imputation_window = (time_hours - last_rescue_time) <= rescue_window & 
            last_rescue_time > -Inf
        ) %>%
        mutate(
          baseline_pain = first(observed_pain),
          # get pain score at time of rescue for locf
          rescue_pain = ifelse(rescue_taken == 1, observed_pain, NA),
          last_rescue_pain = zoo::na.fill(zoo::na.locf(rescue_pain, na.rm = FALSE), NA),
          # get highest observed pain up to rescue time for wocf
          max_pain_to_rescue = zoo::rollapply(observed_pain, width = 1:length(observed_pain), 
                                              FUN = max, align = "right", fill = NA, partial = TRUE),
          # apply imputation
          pain_for_analysis = case_when(
            !in_imputation_window ~ observed_pain,
            analysis_type == "locf" ~ last_rescue_pain,
            analysis_type == "bocf" ~ baseline_pain,
            analysis_type == "wocf" ~ max_pain_to_rescue,
            TRUE ~ observed_pain
          )
        ) %>%
        ungroup()
    }
  }
  
  # calculate spid using correct trapezoidal rule formula
  spid_results <- analysis_data %>%
    group_by(subject_id, treatment) %>%
    arrange(time_hours) %>%
    mutate(
      baseline_pain = first(pain_for_analysis),
      pid = baseline_pain - pain_for_analysis  # pain intensity difference from baseline
    ) %>%
    # apply trapezoidal rule: sum of (pid_j + pid_j+1)/2 * (t_j+1 - t_j)
    mutate(
      next_pid = lead(pid),
      time_interval = lead(time_hours) - time_hours,
      trapezoid_area = ifelse(!is.na(next_pid), 
                              (pid + next_pid) / 2 * time_interval, 
                              0)
    ) %>%
    summarise(
      spid_0_48 = sum(trapezoid_area, na.rm = TRUE),
      baseline_pain = first(baseline_pain),
      analysis_type = analysis_type,
      .groups = "drop"
    )
  
  return(spid_results)
}



####################################################################################################################
# analyze_spid_ancova: runs ANCOVA model controlling for baseline pain, returns both contrasts and group estimates #
####################################################################################################################
analyze_spid_ancova <- function(spid_data, reference_group = "placebo") {
  
  # set treatment to factor with reference level
  spid_data <- spid_data %>%
    mutate(treatment = factor(treatment, levels = c(reference_group, 
                                                    setdiff(unique(treatment), reference_group))))
  # track convergence
  converged <- FALSE
  contrasts <- data.frame()
  group_estimates <- data.frame()
  
  tryCatch({
    # fit ANCOVA model
    model <- lm(spid_0_48 ~ treatment + baseline_pain, data = spid_data)
    converged <- TRUE
    
    # extract coefficient estimates for treatment effect
    coef_summary <- summary(model)$coefficients
    exclude_rows <- c("(Intercept)", "baseline_pain")
    treatment_coefs <- coef_summary[!rownames(coef_summary) %in% exclude_rows, , drop = FALSE]
    
    if (nrow(treatment_coefs) > 0) {
      contrasts <- as.data.frame(treatment_coefs) %>%
        tibble::rownames_to_column("contrast") %>%
        mutate(
          significant = `Pr(>|t|)` < 0.05,
          lower_ci = Estimate - qt(0.975, model$df.residual) * `Std. Error`,
          upper_ci = Estimate + qt(0.975, model$df.residual) * `Std. Error`,
          analysis_type = first(spid_data$analysis_type),
          converged = TRUE
        ) %>%
        rename(estimate = Estimate, se = `Std. Error`, p_value = `Pr(>|t|)`) %>%
        select(contrast, estimate, se, lower_ci, upper_ci, p_value, significant, 
               analysis_type, converged)
    }
    
    # calculate adjusted marginal means for each treatment group
    treatment_levels <- levels(spid_data$treatment)
    baseline_mean <- mean(spid_data$baseline_pain, na.rm = TRUE)
    
    newdata <- data.frame(
      treatment = factor(treatment_levels, levels = treatment_levels),
      baseline_pain = baseline_mean
    )
    
    predictions <- predict(model, newdata = newdata, se.fit = TRUE)
    
    group_estimates <- data.frame(
      treatment_group = treatment_levels,
      estimate = predictions$fit,
      se = predictions$se.fit,
      lower_ci = predictions$fit - qt(0.975, model$df.residual) * predictions$se.fit,
      upper_ci = predictions$fit + qt(0.975, model$df.residual) * predictions$se.fit,
      analysis_type = first(spid_data$analysis_type),
      converged = TRUE
    )
    
  }, error = function(e) {
    # for capture if model fails to converge
    contrasts <- data.frame(
      contrast = "model_failed",
      estimate = NA,
      se = NA,
      lower_ci = NA,
      upper_ci = NA,
      p_value = NA,
      significant = FALSE,
      analysis_type = ifelse("analysis_type" %in% names(spid_data), 
                             first(spid_data$analysis_type), "unknown"),
      converged = FALSE
    )
    
    group_estimates <- data.frame(
      treatment_group = "model_failed",
      estimate = NA,
      se = NA,
      lower_ci = NA,
      upper_ci = NA,
      analysis_type = ifelse("analysis_type" %in% names(spid_data), 
                             first(spid_data$analysis_type), "unknown"),
      converged = FALSE
    )
  })
  
  return(list(
    contrasts = contrasts,
    group_estimates = group_estimates
  ))
}

###############################################################################################################################
# analyze_mira: runs ANCOVA models by MIRA approach controlling for baseline pain, returns both contrasts and group estimates #
###############################################################################################################################
analyze_mira <- function(trial_data, rescue_effect_data, n_imputations = 100, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # use all scheduled AND unscheduled timepoints
  analysis_data <- trial_data %>% 
    filter(scheduled == TRUE | rescue_taken == 1) %>% 
    arrange(subject_id, time_hours)
  
  # identify timepoints affected by rescue
  rescue_affected_data <- analysis_data %>%
    group_by(subject_id) %>%
    arrange(time_hours) %>%
    mutate(
      rescue_times = ifelse(rescue_taken == 1, time_hours, NA),
      last_rescue_time = zoo::na.fill(zoo::na.locf(rescue_times, na.rm = FALSE), -Inf),
      time_since_rescue = time_hours - last_rescue_time,
      affected_by_rescue = time_since_rescue <= 6 & time_since_rescue >= 0 & !is.infinite(time_since_rescue)
    ) %>%
    ungroup()
  
  # store results from each imputation
  imputation_contrasts <- list()
  imputation_group_estimates <- list()
  imputation_converged <- logical(n_imputations)
  
  # perform multiple imputations
  for (imp_num in 1:n_imputations) {
    
    # create imputed dataset
    imputed_data <- rescue_affected_data %>%
      mutate(
        # for affected timepoints, subtract expected rescue effect from observed pain to estimate counterfactual
        imputed_pain = ifelse(
          affected_by_rescue,
          {
            # get mean effect from reference data
            expected_rescue_effect <- approx(rescue_effect_data$time_hr, rescue_effect_data$gamma_pred,
                                             xout = time_since_rescue, rule = 2)$y
            
            # get variability from reference data
            rescue_effect_sd <- approx(rescue_effect_data$time_hr, rescue_effect_data$sd_hc,
                                       xout = time_since_rescue, rule = 2)$y
            
            # sample rescue effect from distribution
            sampled_rescue_effect <- rnorm(n(), mean = expected_rescue_effect, sd = rescue_effect_sd)
            
            # impute counterfactual as observed - rescue effect
            pmax(0, pmin(10, observed_pain - sampled_rescue_effect))
          },
          observed_pain  # unchanged for non-rescue-affected timepoints
        )
      )
    
    # calculate spid for this imputed dataset using trapezoidal rule
    spid_imputed <- imputed_data %>%
      group_by(subject_id, treatment) %>%
      arrange(time_hours) %>%
      mutate(
        baseline_pain = first(observed_pain), 
        pid = baseline_pain - imputed_pain  
      ) %>%
      # apply trapezoidal rule: sum of (pid_j + pid_j+1)/2 * (t_j+1 - t_j)
      mutate(
        next_pid = lead(pid),
        time_interval = lead(time_hours) - time_hours,
        trapezoid_area = ifelse(!is.na(next_pid), 
                                (pid + next_pid) / 2 * time_interval, 
                                0)
      ) %>%
      summarise(
        spid_0_48 = sum(trapezoid_area, na.rm = TRUE),
        baseline_pain = first(baseline_pain),
        .groups = "drop"
      )
    
    # analyze this imputed dataset using same ANCOVA approach as other methods
    spid_imputed <- spid_imputed %>%
      mutate(
        treatment = factor(treatment, levels = c("placebo", 
                                                 setdiff(unique(treatment), "placebo"))),
        analysis_type = "mira"
      )
    
    tryCatch({
      analysis_result <- analyze_spid_ancova(spid_imputed, reference_group = "placebo")
      
      if (nrow(analysis_result$contrasts) > 0 && analysis_result$contrasts$converged[1]) {
        imputation_contrasts[[imp_num]] <- analysis_result$contrasts
        imputation_group_estimates[[imp_num]] <- analysis_result$group_estimates
        imputation_converged[imp_num] <- TRUE
      } else {
        imputation_converged[imp_num] <- FALSE
      }
      
    }, error = function(e) {
      imputation_converged[imp_num] <- FALSE
    })
  }
  
  # apply rubin's rules to combine imputations
  if (sum(imputation_converged) >= 2) {
    # use only converged imputations
    valid_contrasts <- imputation_contrasts[imputation_converged]
    valid_group_estimates <- imputation_group_estimates[imputation_converged]
    m <- length(valid_contrasts)
    
    # pool contrasts
    pooled_contrasts <- map_dfr(seq_len(nrow(valid_contrasts[[1]])), function(i) {
      estimates <- map_dbl(valid_contrasts, ~ .x$estimate[i])
      ses <- map_dbl(valid_contrasts, ~ .x$se[i])
      
      # rubin's rules
      pooled_estimate <- mean(estimates)
      within_var <- mean(ses^2)
      between_var <- var(estimates)
      total_var <- within_var + (1 + 1/m) * between_var
      pooled_se <- sqrt(total_var)
      df <- (m - 1) * (1 + within_var / ((1 + 1/m) * between_var))^2
      
      data.frame(
        contrast = valid_contrasts[[1]]$contrast[i],
        estimate = pooled_estimate,
        se = pooled_se,
        lower_ci = pooled_estimate - qt(0.975, df) * pooled_se,
        upper_ci = pooled_estimate + qt(0.975, df) * pooled_se,
        p_value = 2 * (1 - pt(abs(pooled_estimate / pooled_se), df)),
        significant = abs(pooled_estimate / pooled_se) > qt(0.975, df),
        analysis_type = "mira",
        n_imputations = m,
        converged = TRUE
      )
    })
    
    # pool group estimates
    pooled_group_estimates <- map_dfr(seq_len(nrow(valid_group_estimates[[1]])), function(i) {
      estimates <- map_dbl(valid_group_estimates, ~ .x$estimate[i])
      ses <- map_dbl(valid_group_estimates, ~ .x$se[i])
      
      # rubin's rules
      pooled_estimate <- mean(estimates)
      within_var <- mean(ses^2)
      between_var <- var(estimates)
      total_var <- within_var + (1 + 1/m) * between_var
      pooled_se <- sqrt(total_var)
      df <- (m - 1) * (1 + within_var / ((1 + 1/m) * between_var))^2
      
      data.frame(
        treatment_group = valid_group_estimates[[1]]$treatment_group[i],
        estimate = pooled_estimate,
        se = pooled_se,
        lower_ci = pooled_estimate - qt(0.975, df) * pooled_se,
        upper_ci = pooled_estimate + qt(0.975, df) * pooled_se,
        analysis_type = "mira",
        n_imputations = m,
        converged = TRUE
      )
    })
    
  } else {
    # if MI fails, return NAs for contrast and group results
    pooled_contrasts <- data.frame(
      contrast = "active_vs_placebo",
      estimate = NA,
      se = NA, 
      lower_ci = NA,
      upper_ci = NA,
      p_value = NA,
      significant = FALSE,
      analysis_type = "mira",
      n_imputations = 0,
      converged = FALSE
    )
    
    pooled_group_estimates <- data.frame(
      treatment_group = "model_failed",
      estimate = NA,
      se = NA,
      lower_ci = NA,
      upper_ci = NA,
      analysis_type = "mira",
      n_imputations = 0,
      converged = FALSE
    )
  }
  
  return(list(
    contrasts = pooled_contrasts,
    group_estimates = pooled_group_estimates
  ))
}

###################################################################################
# calculate_root: creates subject-level dataset for ROOT analysis           
###################################################################################
calculate_root <- function(trial_data, threshold = 0.5, rescue_window = 6) {
  
  root_data <- trial_data %>% 
    filter(scheduled == TRUE | rescue_taken == 1) %>%
    arrange(subject_id, time_hours)
  
  root_results <- root_data %>%
    group_by(subject_id, treatment) %>%
    arrange(time_hours) %>%
    mutate(
      baseline_pain = first(observed_pain),
      percent_reduction = (baseline_pain - observed_pain) / baseline_pain,
      rescue_times = ifelse(rescue_taken == 1, time_hours, NA),
      last_rescue_time = zoo::na.fill(zoo::na.locf(rescue_times, na.rm = FALSE), -Inf),
      affected_by_rescue = (time_hours - last_rescue_time) <= rescue_window & 
        last_rescue_time > -Inf,
      is_responder = percent_reduction >= threshold & !affected_by_rescue,
      time_interval = c(diff(time_hours), 0)
    ) %>%
    summarise(
      total_time = sum(time_interval),
      responder_time = sum(time_interval[is_responder]),
      root_0_48 = responder_time / total_time,
      baseline_pain = first(baseline_pain),
      .groups = "drop"
    )
  
  return(root_results)
}

#####################################################################################################################
# analyze_root_zibeta: runs ZI beta model controlling for baseline pain, returns both contrasts and group estimates #
#####################################################################################################################
analyze_root_zibeta <- function(root_data, reference_group = "placebo") {
  
  # set treatment to factor with reference level
  root_data <- root_data %>%
    mutate(treatment = factor(treatment, levels = c(reference_group, 
                                                    setdiff(unique(treatment), reference_group))))
  
  outcome <- "root_0_48"
  treatment_var <- "treatment"
  baseline_vars <- "baseline_pain"
  vars <- c(treatment_var, baseline_vars)
  
  formula_mu <- reformulate(vars, response = outcome)
  formula_pi <- reformulate(vars)
  
  tryCatch({
    fit <- gamlss(
      formula = formula_mu,
      sigma.formula = ~ 1,   
      nu.formula = formula_pi,
      family = BEZI,
      data = root_data,
      control = gamlss.control(trace = FALSE)
    )
    
    # levels and covariate means for adjusted marginal means
    trt_levels <- levels(root_data[[treatment_var]])
    baseline_means <- sapply(root_data[baseline_vars], mean, na.rm = TRUE)
    newdata <- as.data.frame(matrix(rep(baseline_means, times = length(trt_levels)),
                                    nrow = length(trt_levels), byrow = TRUE))
    colnames(newdata) <- baseline_vars
    newdata[[treatment_var]] <- factor(trt_levels, levels = trt_levels)
    newdata <- newdata[, c(treatment_var, baseline_vars)]
    
    # design matrices
    X_mu <- model.matrix(reformulate(vars), data = newdata)
    X_pi <- model.matrix(reformulate(vars), data = newdata)
    
    # coefficients 
    beta_mu <- coef(fit, what = "mu")
    beta_pi <- coef(fit, what = "nu")
    names(beta_mu) <- paste0("mu_", names(beta_mu))
    names(beta_pi) <- paste0("pi_", names(beta_pi))
    
    # variance-covariance matrix
    V_all <- vcov(fit)
    mu_names <- paste0("mu_", names(coef(fit, "mu")))
    sigma_names <- paste0("sigma_", names(coef(fit, "sigma")))
    pi_names <- paste0("pi_", names(coef(fit, "nu")))
    rownames(V_all) <- colnames(V_all) <- c(mu_names, sigma_names, pi_names)
    keep_names <- c(mu_names, pi_names)
    V <- V_all[keep_names, keep_names]
    
    # linear predictors with inverse links
    mu_hat <- plogis(as.numeric(X_mu %*% beta_mu))
    pi_hat <- plogis(as.numeric(X_pi %*% beta_pi))
    
    # adjusted marginal means
    f <- mu_hat * (1 - pi_hat)
    
    # gradients
    grad_mu <- (1 - pi_hat) * mu_hat * (1 - mu_hat) * X_mu
    grad_pi <- -mu_hat * pi_hat * (1 - pi_hat) * X_pi
    grad_list <- lapply(seq_along(mu_hat), function(i) c(grad_mu[i, ], grad_pi[i, ]))
    
    # group estimates
    group_estimates <- map_dfr(seq_along(trt_levels), function(i) {
      var_i <- as.numeric(t(grad_list[[i]]) %*% V %*% grad_list[[i]])
      se_i <- sqrt(var_i)
      
      data.frame(
        treatment_group = trt_levels[i],
        estimate = f[i],
        se = se_i,
        lower_ci = f[i] - qnorm(0.975) * se_i,
        upper_ci = f[i] + qnorm(0.975) * se_i,
        converged = TRUE
      )
    })
    
    # contrasts
    ref_idx <- which(trt_levels == reference_group)
    non_ref_indices <- which(trt_levels != reference_group)
    
    if (length(non_ref_indices) > 0) {
      contrasts <- map_dfr(non_ref_indices, function(j) {
        i <- ref_idx
        delta <- f[j] - f[i]
        var_j <- as.numeric(t(grad_list[[j]]) %*% V %*% grad_list[[j]])
        var_i <- as.numeric(t(grad_list[[i]]) %*% V %*% grad_list[[i]])
        cov_ij <- as.numeric(t(grad_list[[j]]) %*% V %*% grad_list[[i]])
        delta_var <- var_j + var_i - 2 * cov_ij
        delta_se <- sqrt(delta_var)
        z <- delta / delta_se
        p <- 2 * (1 - pnorm(abs(z)))
        
        data.frame(
          contrast = paste0("treatment", trt_levels[j]),
          estimate = delta, se = delta_se,
          lower_ci = delta - qnorm(0.975) * delta_se,
          upper_ci = delta + qnorm(0.975) * delta_se,
          p_value = p, significant = p < 0.05, converged = TRUE
        )
      })
    } else {
      contrasts <- data.frame()
    }
    
    return(list(
      contrasts = contrasts,
      group_estimates = group_estimates
    ))
    
  }, error = function(e) {
    return(list(
      contrasts = data.frame(
        contrast = "model_failed", estimate = NA_real_, se = NA_real_, 
        lower_ci = NA_real_, upper_ci = NA_real_, p_value = NA_real_,
        significant = FALSE, converged = FALSE
      ),
      group_estimates = data.frame(
        treatment_group = "model_failed", estimate = NA_real_, se = NA_real_,
        lower_ci = NA_real_, upper_ci = NA_real_, converged = FALSE
      )
    ))
  })
}

