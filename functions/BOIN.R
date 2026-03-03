run_BOIN <- function(data,
                     N,
                     cohort_size = 3,
                     target_DLT,
                     p_L = NULL, p_U = NULL,
                     tox_bound = 0.95,
                     prior = c(1,1),
                     drug_detail = TRUE
                     )
{
  # Extract dose levels
  doses <- data$profile_info$dose
  ndoses <- length(doses)
  
  # Truncate to N patients
  df_DLT <- data$data_ph1_DLT[1:N, ]
  
  # get the boundaries for escalation
  if(is.null(p_L)) {p_L <- 0.6*target_DLT}
  if(is.null(p_U)) {p_U <- 1.4*target_DLT}
  bound <- BOIN::get.boundary(target = target_DLT,
                              ncohort = 100,
                              cohortsize = cohort_size,
                              p.saf = p_L,
                              p.tox = p_U,
                              cutoff.eli = tox_bound)
  lambda_e <- bound$lambda_e
  lambda_d <- bound$lambda_d
  
  # Initialize trial state
  p_id <- 1
  outcome <- numeric()
  dose_history <- c()
  current_dose <- 1
  stopped <- 0
  # vector to store observed DLT at each dose level
  DLT_per_dose <- rep(0,ndoses)
  # vector to store estimated DLT rate for each dose
  DLT_rate_per_dose <- rep(NA,ndoses)
  # vector to store number of patients assigned to each dose level
  n_per_dose <- rep(0,ndoses)
  # vector to store posterior beta parameters for each dose
  a <- b <- rep(NA,ndoses)
  # vector to store posterior probability that the current dose has DLT > target DLT
  tox_prob_per_dose <- rep(NA,ndoses)
  # vector to store inadmissible doses
  inadmissible_doses  <- numeric(0)
  
  # start dose escalation
  while(stopped==0 & p_id <= N)
  {
    # Assign next cohort
    cohort_ids <- p_id:min(p_id + cohort_size - 1, N)
    
    # Potential outcomes at current dose
    tox <- df_DLT[cohort_ids, current_dose]
    
    # update posterior beta-binomial parameters
    n_per_dose[current_dose] <- n_per_dose[current_dose] + length(cohort_ids)
    DLT_per_dose[current_dose] <- DLT_per_dose[current_dose] + sum(tox)
    a[current_dose] <- prior[1]+DLT_per_dose[current_dose]
    b[current_dose] <- prior[2]+n_per_dose[current_dose]-DLT_per_dose[current_dose]
    DLT_rate_per_dose[current_dose] <- DLT_per_dose[current_dose]/n_per_dose[current_dose]
    
    # Append to outcome history
    outcome <- c(outcome,tox)
    
    # Update dose history
    dose_history <- c(dose_history, rep(current_dose, length(cohort_ids)))
    
    # compare to two boundaries to decide (de)-escalation
    if(DLT_rate_per_dose[current_dose] <= lambda_e)
    {
      # escalate to the next dose
      next_dose <- current_dose + 1
    } else if (DLT_rate_per_dose[current_dose] >= lambda_d)
    {
      # de-escalate to lower dose
      next_dose <- current_dose - 1
    } else
    {
      # stay at the current dose
      next_dose <- current_dose
    }
    
    #########################
    # Over dose control
    #########################
    # get posterior probability that the current dose is above the MTD
    # P(p_T(d) > target_DLT|data), i.e., Pr(toxicity prob at dose d exceeds target DLT)
    tox_prob_per_dose[current_dose] <- pbeta(target_DLT,
                                             a[current_dose],
                                             b[current_dose],
                                             lower.tail = FALSE)
    # check if any dose is inadmissible: P(p_T(d) > target_DLT|data) > tox_bound and num at this dose >=3
    if(any(tox_prob_per_dose > tox_bound,na.rm = TRUE))
    {
      # the lowest dose that is too toxic
      lowest_dose_too_toxic <- min(which(tox_prob_per_dose > tox_bound))
      # this dose and all doses above are inadmissible
      inadmissible_doses <- c(inadmissible_doses,lowest_dose_too_toxic:ndoses)
    }
    # for doses with information, check if they are admissible by the specified toxicity bound
    # doses not explored yet are also admissible
    admissible_doses <- union(which(tox_prob_per_dose <= tox_bound),
                              which(is.na(tox_prob_per_dose))) |> as.numeric()
    admissible_doses <- setdiff(admissible_doses,
                                inadmissible_doses)
    if(length(admissible_doses)==0)
    {
      next_dose <- NA
    } else
    {
      # in case all admissible doses are lower than the next dose
      # or the next dose is 0
      next_dose <- min(next_dose,max(admissible_doses))
    }
    # Check for stopping
    if (is.na(next_dose) | next_dose < 1 | next_dose > ndoses) {
      cat("Trial stopped early at", max(cohort_ids), "patients.\n")
      stopped <- 1
    }
    
    # Update dose and patient ID
    current_dose <- next_dose
    p_id <- max(cohort_ids) + 1
  } # end of while() loop
  if(stopped)
  {
    RP2D <- NA
  } else
  {
    # estimate MTD from trial data using Isonotic regression
    mtd.obj <- BOIN::select.mtd(target = target_DLT,
                                npts = n_per_dose,
                                ntox = DLT_per_dose,
                                cutoff.eli = tox_bound,
                                boundMTD = TRUE,
                                p.tox = p_U)
    # recommended phase 2 dose (RP2D)
    RP2D <- mtd.obj$MTD
  }

  #prepare output
  out <- list()
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
  summary_tab <- table(factor(dose_history,levels = 1:ndoses)) |>
    as.data.frame()
  colnames(summary_tab) <- c("dose_level","num_pts")
  summary_tab$dose <- doses
  out$summary_tab <- summary_tab
  # dose escalation history
  out$dose_history <- dose_history
  # patient DLT outcomes
  out$DLT_outcome <- outcome
  return(out)
}

