run_BOIN12 <- function(data,
                       N,
                       cohort_size = 3,
                       phi_t, phi_e,
                       c_t = 0.95, c_e = 0.9,
                       u1 = 100, u2, u3,u4=0,
                       n_star1 = 6, n_star2 = 8)
{
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
  DLT_outcome <- c()
  rsp_outcome <- c()
  
  # get BOIN boundary if there is dose exploration rule
  if(!is.null(n_star2))
  {
    bound <- BOIN::get.boundary(target = phi_t,
                                ncohort = 100,
                                cohortsize = cohort_size,
                                p.saf = 0.6*phi_t,
                                p.tox = 1.4*phi_t)
    lambda_d <- bound$lambda_d
  }
  
  # Iteratively run the BOIN12 design
  while (p_id <= N) {
    # Assign next cohort
    cohort_ids <- p_id:min(p_id + cohort_size - 1, N)
    
    # Potential outcomes at current dose
    tox <- df_DLT[cohort_ids, current_dose]; DLT_outcome <- c(DLT_outcome,tox)
    eff <- df_rsp[cohort_ids, current_dose]; rsp_outcome <- c(rsp_outcome,eff)
    
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
    
    # fit the BOIN12 model
    design_BOIN12 <- escalation::get_boin12(num_doses = ndoses,
                             phi_t = phi_t,phi_e = phi_e,
                             u1=u1,u2=u2,u3=u3,u4=u4,
                             n_star = n_star1,
                             c_t = c_t,c_e = c_e) 
    fit_BOIN12 <- design_BOIN12 |> escalation::fit(outcome_str)
    
    # Get next recommended dose
    next_dose <- fit_BOIN12$recommended_dose
    
    # Extra dose-exploration rule
    # If the number of patients treated at the current dose d is greater than some value, n_star2 here
    # The observed toxicity rate of this dose is less than the de-escalation boundary lambda_d (from the univariate BOIN)
    # The next higher dose has never been used for treating patients, 
    # we escalate the dose for treating the next cohort of patients.
    
    if(!is.null(n_star2))
    { 
      if(sum(dose_history==current_dose) >= n_star2 & max(dose_history) < ndoses)
      {
        # the observed toxicity rate at current dose
        p_t_current_dose <- summary(fit_BOIN12) |>
          data.frame() |>
          dplyr::filter(dose == current_dose) |>
          dplyr::select(empiric_tox_rate) |>
          unlist() |>
          as.numeric()
        if(p_t_current_dose < lambda_d)
        {
          next_dose <- current_dose + 1
        }
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
  } # end of while
  
  # prepare output
  out <- list()
  # recommended phase 2 dose (RP2D)
  RP2D <- fit_BOIN12$recommended_dose
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
  summary_tab <- table(factor(dose_history,levels = 1:ndoses)) |>
    as.data.frame()
  colnames(summary_tab) <- c("dose_level","num_pts")
  summary_tab$dose <- doses
  out$summary_tab <- summary_tab
  # dose escalation history
  out$dose_history <- dose_history
  # patient DLT outcomes
  out$DLT_outcome <- DLT_outcome
  # patient rsp outcomes
  out$rsp_outcome <- rsp_outcome
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
# res=run_BOIN12(data = data,N=N,u2=u2,u3=u3,
#                phi_t = phi_t , phi_e = phi_e)
# res$dose_history
# res$RP2D
