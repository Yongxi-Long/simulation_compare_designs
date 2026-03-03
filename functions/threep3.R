run_threep3 <- function(data,
                        N,
                        A = 3, B= 3,
                        C=1, D=1, E=1,
                        quite = TRUE,
                        drug_detail = TRUE)
{
  # Extract dose levels
  doses <- data$profile_info$dose
  ndoses <- length(doses)
  
  # Truncate to N patients
  df_DLT <- data$data_ph1_DLT[1:N, ]
  
  # Initialize trial state
  p_id <- 1
  outcome <- numeric()
  dose_history <- c()
  current_dose <- 1
  stopped <- 0
  # record number of patients and DLTs per dose
  # vector to store observed DLT at each dose level
  DLT_per_dose <- rep(0,ndoses)
  # vector to store estimated DLT rate for each dose
  DLT_rate_per_dose <- rep(NA,ndoses)
  # vector to store number of patients assigned to each dose level
  n_per_dose <- rep(0,ndoses)
  # MTD
  MTD <- NA
  
  # start dose escalation
  while(stopped==0 & p_id <= N)
  {
    # Assign next cohort
    # determine cohort size based on whether or not this dose has been explored
    cohort_size <- ifelse(n_per_dose[current_dose]==0,A,B)
    cohort_ids <- p_id:min(p_id + cohort_size - 1, N)
    
    # Potential outcomes at current dose
    tox <- df_DLT[cohort_ids, current_dose]
    
    n_per_dose[current_dose] <- n_per_dose[current_dose] + length(cohort_ids)
    DLT_per_dose[current_dose] <- DLT_per_dose[current_dose] + sum(tox)
    DLT_rate_per_dose[current_dose] <- DLT_per_dose[current_dose]/n_per_dose[current_dose]
    
    # Append to outcome history
    outcome <- c(outcome,tox)
    
    # Update dose history
    dose_history <- c(dose_history, rep(current_dose, length(cohort_ids)))
    
    # compare to decision rules to decide (de)-escalation
    # if this dose is explored for the first time
    if(n_per_dose[current_dose]==A)
    {
      if(DLT_rate_per_dose[current_dose] < C/A)
      {
        #escalation
        next_dose <- current_dose + 1
        # if there is no more dose levels to escalate, add B more to current dose level
        if(next_dose>ndoses)
        {
          next_dose <- current_dose
        }
      } else if (DLT_rate_per_dose[current_dose] > D/A)
      {
        #de-escalate
        next_dose <- current_dose - 1
      } else
      {
        # stay at the current dose
        next_dose <- current_dose
      }
    } else
    {
      # if the dose is re-visited
      if(DLT_rate_per_dose[current_dose] <= E/(A+B))
      {
        # escalation
        next_dose <- current_dose + 1
        # if there is no more dose levels to escalate, stop and declare current dose to be MTD
        if(next_dose>ndoses)
        {
          stopped <- 1
          MTD <- current_dose
        }
      } else
      {
        # break, previous dose is MTD
        stopped <- 1
        MTD <- ifelse(current_dose - 1 > 0,current_dose - 1,NA)
      }
    }
    
    # Check for stopping
    if (is.na(next_dose) | next_dose < 1 | next_dose > ndoses) {
      if(!quite)
      {
        cat("Trial stopped early at", max(cohort_ids), "patients.\n")
      }
      break
    }
    
    # Update dose and patient ID
    current_dose <- next_dose
    p_id <- max(cohort_ids) + 1
  } # end of while() loop

  #prepare output
  out <- list()
  # recommended phase 2 dose (RP2D)
  RP2D <- MTD
  
  out$RP2D <- RP2D
  # profile info for RP2D
  if(is.na(RP2D)| drug_detail==FALSE)
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
