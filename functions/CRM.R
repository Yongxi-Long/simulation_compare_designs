require(trialr)
run_CRM <- function(data,
                    N,
                    cohort_size = 1,
                    target_DLT,
                    tox_bound = 0.9,
                    skeleton=NULL,
                    model="empiric",
                    a0 = NULL,
                    alpha_mean = NULL, alpha_sd = NULL,
                    beta_mean = NULL, beta_sd = sqrt(1.34),
                    beta_shape = NULL, beta_inverse_scale = NULL,
                    seed = 123,
                    chains = 1,
                    iter = 1000,
                    warmup = 250,
                    drug_detail = TRUE)
{
  # Extract dose levels
  doses <- data$profile_info$dose
  ndoses <- length(doses)
  
  # Truncate to N patients
  df_DLT <- data$data_ph1_DLT[1:N, ]
  
  # if user does not supply a skeleton, then use the default
  if(is.null(skeleton))
  {
    skeleton <- dfcrm::getprior(halfwidth = 0.05,target = target_DLT,
                                nu = floor(median(1:ndoses)),
                                nlevel = ndoses,model = model)
  }
  
  # Initialize trial state
  p_id <- 1
  outcome_str <- ""
  dose_history <- c()
  current_dose <- 1
  stopped <- 0
  
  while(stopped==0 & p_id <= N)
  {
    # Assign next cohort
    cohort_ids <- p_id:min(p_id + cohort_size - 1, N)
    
    # Potential outcomes at current dose
    tox <- df_DLT[cohort_ids, current_dose]
    
    # Build outcome string for this cohort
    cohort_str <- make_outcome_str(tox, dose_level=current_dose)
    
    # Append to outcome history
    if (outcome_str == "") {
      outcome_str <- cohort_str
    } else {
      outcome_str <- paste(outcome_str, cohort_str)
    }
    
    # Update dose history
    dose_history <- c(dose_history, rep(current_dose, length(cohort_ids)))
    
    fit_CRM <- trialr::stan_crm(outcome_str =outcome_str, 
                                skeleton = skeleton, 
                                target = target_DLT, 
                                model = model,
                                a0 = a0,
                                alpha_mean = alpha_mean, alpha_sd = alpha_sd,
                                beta_mean = beta_mean, beta_sd= beta_sd,
                                beta_shape = beta_shape, beta_inverse_scale = beta_inverse_scale,
                                seed = seed,
                                refresh = 0,
                                chains = chains,
                                iter = iter,
                                warmup = warmup)
    
    # get the next recommend dose
    next_dose <- fit_CRM$recommended_dose
    #########################
    # Over dose control
    #########################
    # get posterior probability that the current dose is above the MTD
    # P(p_T(d) > target_DLT|data), i.e., Pr(toxicity prob at dose d exceeds target DLT)
    # get a sample from the posterior
    prob_exceeds_target_DLT <- fit_CRM |>
      tidybayes::gather_draws(prob_tox[dose]) |>
      dplyr::group_by(dose) |>
      dplyr::summarise(prob_exceeds_target_DLT = mean(.value > target_DLT)) |>
      dplyr::select(prob_exceeds_target_DLT) |>
      unlist()
    # see if the recommended dose has Pr(toxicity prob at dose d exceeds target DLT) <= tox_bound
    admissible_doses <- which(prob_exceeds_target_DLT <= tox_bound) |> as.numeric()
    if(length(admissible_doses)==0)
    {
      next_dose <- NA
    } else
    {
      if(!(next_dose %in% admissible_doses))
      {
        next_dose <- max(admissible_doses)
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
  } # end while
  
  #prepare output
  out <- list()
  # recommended phase 2 dose (RP2D)
  RP2D <- fit_CRM$recommended_dose
  # check the feasibility of the recommend dose
  if(!(RP2D %in% admissible_doses)) 
  {
    if(length(admissible_doses)==0) 
    {
      RP2D <- NA
    } else
    {
      RP2D <- max(admissible_doses)
    }
  }
  out$RP2D <- RP2D
  # profile info for RP2D
  if(is.na(RP2D) | drug_detail == FALSE)
  {
    out$RP2D_info <- NA
  } else
  {
    out$RP2D_info <- list("dose" = doses[RP2D],
                          "DLT_rate" = data$profile_info$DLT_rate[RP2D],
                          "rsp_rate" = data$profile_info$rsp_rate[RP2D],
                          "median_OS" = data$profile_info$median_OS[RP2D])
  }
  # table of dose levels, actual dose values and patient assignment
  summary_tab <- table(factor(fit_CRM$dat$doses,levels = 1:ndoses)) |>
    as.data.frame()
  colnames(summary_tab) <- c("dose_level","num_pts")
  summary_tab$dose <- doses
  out$summary_tab <- summary_tab
  # dose escalation history
  out$dose_history <- fit_CRM$dat$doses
  # patient DLT outcomes
  out$DLT_outcome <- fit_CRM$dat$tox
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
