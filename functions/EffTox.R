run_EffTox <- function(data,
                       N,
                       cohort_size,
                       efficacy_hurdle,
                       toxicity_hurdle,
                       p_e,
                       p_t,
                       eff0,
                       tox1,
                       eff_star,
                       tox_star,
                       priors_EffTox,
                       seed = 123,
                       chains = 1,
                       iter = 1000,
                       warmup = 250) {
  
  # Extract dose levels
  doses <- data$profile_info$dose
  ndoses <- length(doses)
  
  # Truncate to N patients
  df_DLT <- data$data_ph1_DLT[1:N, ]
  df_rsp <- data$data_rsp[1:N, ]
  
  # Initialize trial state
  p_id <- 1
  outcome_str <- ""
  dose_history <- c()
  current_dose <- 1
  
  # Iteratively run the EffTox design
  while (p_id <= N) {
    # Assign next cohort
    cohort_ids <- p_id:min(p_id + cohort_size - 1, N)
    
    # Potential outcomes at current dose
    tox <- df_DLT[cohort_ids, current_dose]
    eff <- df_rsp[cohort_ids, current_dose]
    
    # Build outcome string for this cohort
    cohort_str <- make_outcome_str(tox, eff, current_dose)
    
    # Append to outcome history
    if (outcome_str == "") {
      outcome_str <- cohort_str
    } else {
      outcome_str <- paste(outcome_str, cohort_str)
    }
    
    # Update dose history
    dose_history <- c(dose_history, rep(current_dose, length(cohort_ids)))
    
    # Fit EffTox model
    fit_EffTox <- stan_efftox(
      outcome_str = outcome_str,
      real_doses = doses,
      efficacy_hurdle = efficacy_hurdle,
      toxicity_hurdle = toxicity_hurdle,
      p_e = p_e,
      p_t = p_t,
      eff0 = eff0,
      tox1 = tox1,
      eff_star = eff_star,
      tox_star = tox_star,
      priors = priors_EffTox,
      seed = seed,
      refresh = 0,
      chains = chains,
      iter = iter,
      warmup = warmup
    )
    
    # Get next recommended dose
    next_dose <- fit_EffTox$recommended_dose
    
    # EffTox sometimes recommend to stop when there is not enough responses
    # we will overwrite this decision when not all doses have been explored
    # we will recommend the next dose to be the maximum dose for which the 
    # posterior mean probabilities that toxicity at the doses is acceptable is larger than the toxicity hurdle
    if(is.na(next_dose)  & max(fit_EffTox$doses) < ndoses)
    {
      candidate_dose <- current_dose + 1 
      # force to go to next dose,if next dose has not yet been explored
      if(!(candidate_dose %in% dose_history))
      {
        next_dose <- candidate_dose
      } else
      {
        # if next dose has been explored, check if it has acceptable toxicity
        check_tox <- fit_EffTox$prob_acc_tox[candidate_dose] > p_t
        next_dose <- ifelse(check_tox,candidate_dose,NA)
      }
    }
    
    # Check for stopping
    if (is.na(next_dose) | next_dose < 1 | next_dose > ndoses) {
      cat("Trial stopped early at", max(cohort_ids), "patients.\n")
      break
    }
    
    # Update dose and patient ID
    current_dose <- next_dose
    p_id <- max(cohort_ids) + 1
  }
  # prepare output
  out <- list()
  # recommended phase 2 dose (RP2D)
  RP2D <- fit_EffTox$recommended_dose
  out$RP2D <- RP2D
  # profile info for RP2D
  if(is.na(RP2D))
  {
    out$RP2D_info <- NA
  } else
  {
    out$RP2D_info <- list("dose" = doses[RP2D],
                          "DLT_rate" = data$profile_info$DLT_rate[RP2D],
                          "rsp_rate" = data$profile_info$rsp_rate[RP2D],
                          "median_OS" = data$profile_info$median_OS[RP2D])
  }
  # table of dose levels, actual dose values and patient assigment
  summary_tab <- table(factor(fit_EffTox$dat$doses,levels = 1:ndoses)) |>
    as.data.frame()
  colnames(summary_tab) <- c("dose_level","num_pts")
  summary_tab$dose <- doses
  out$summary_tab <- summary_tab
  # dose escalation history
  out$dose_history <- fit_EffTox$dat$doses
  # patient DLT outcomes
  out$DLT_outcome <- fit_EffTox$dat$tox
  # patient rsp outcomes
  out$rsp_outcome <- fit_EffTox$dat$eff
  return(out)
}


# Helper function: convert patient outcomes to outcome string
make_outcome_str <- function(tox, eff=NULL, dose_level) {
  if(is.null(eff))
  {
    letters <- sapply(seq_along(tox), function(i)
      {
      if(tox[i]==1) "T"
      else if(tox[i]==0) "N"
    })
  } else
  {
    letters <- sapply(seq_along(tox), function(i) {
      if (tox[i] == 1 & eff[i] == 0) "T"
      else if (tox[i] == 0 & eff[i] == 1) "E"
      else if (tox[i] == 0 & eff[i] == 0) "N"
      else if (tox[i] == 1 & eff[i] == 1) "B"
    })
  }
  paste0(dose_level, paste0(letters, collapse = ""))
}
