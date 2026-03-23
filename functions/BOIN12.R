run_BOIN12 <- function(data,
                       N,
                       cohort_size = 3,
                       phi_t, phi_e,
                       c_t = 0.95, c_e = 0.9,
                       u1 = 100, u2, u3,u4=0,
                       u_b = NULL,
                       n_star1 = 6, n_star2 = 8)
{
  
# Extract dose levels
doses <- data$profile_info$dose
ndoses <- length(doses)

# calculate the default value of utility threshold u_b if not supplied
if(is.null(u_b))
{
  # the highest utility that's deemed undesirable
  u_underscore = u1*phi_e*(1-phi_t)+u2*(1-phi_t)*(1-phi_e)+u3*phi_t*phi_e
  # u_b lies in the middle of u_underscore and 100
  u_b <- u_underscore + (100-u_underscore)/2
  # standardize it
  u_b <- u_b/100
}

# Truncate to N patients
df_DLT <- data$data_ph1_DLT[1:N, ]
df_rsp <- data$data_rsp[1:N, ]

# Initialize trial state
p_id <- 1
dose_history <- c()
current_dose <- 1
DLT_outcome <- c()
rsp_outcome <- c()
# vector to store posterior beta parameters
a_T <- b_T <- a_E <- b_E <- a_U <- b_U <- rep(1,ndoses)
# prior for utility
U_prior <- c(1,1)
# vector to store num of patients assigned to each dose
num_pts <- rep(0,ndoses)
# vector to store observed number for 4 outcome combinations
num_tox_noeff <- num_notox_noeff <- num_notox_eff <- num_tox_eff <- rep(0,ndoses)
# vector of observed DLT rate
DLT_rate_obs <- rep(NA,ndoses)
# stop indicator
stopped <- 0

# get BOIN boundary if there is dose exploration rule
bound <- BOIN::get.boundary(target = phi_t,
                              ncohort = 100,
                              cohortsize = cohort_size,
                              p.saf = 0.6*phi_t,
                              p.tox = 1.4*phi_t)
lambda_d <- bound$lambda_d
lambda_e <- bound$lambda_e

# Iteratively run the BOIN12 design
while (p_id <= N) {
  # Assign next cohort
  cohort_ids <- p_id:min(p_id + cohort_size - 1, N)
  
  # Potential outcomes at current dose
  tox <- df_DLT[cohort_ids, current_dose]; DLT_outcome <- c(DLT_outcome,tox)
  eff <- df_rsp[cohort_ids, current_dose]; rsp_outcome <- c(rsp_outcome,eff)

  # Update dose history
  dose_history <- c(dose_history, rep(current_dose, length(cohort_ids)))
  
  # update number of patients per dose
  num_pts[current_dose] <- num_pts[current_dose] + length(cohort_ids)
  
  # update number of patients with the 4 types of outcomes
  num_notox_eff[current_dose] <- num_notox_eff[current_dose] + sum(tox==0 & eff == 1)
  num_tox_eff[current_dose] <- num_tox_eff[current_dose] + sum(tox==1 & eff ==1)
  num_tox_noeff[current_dose] <- num_tox_noeff[current_dose] + sum(tox==1 & eff==0)
  num_notox_noeff[current_dose] <- num_notox_noeff[current_dose] + sum(tox==0 & eff ==0)
  
  # update posterior parameters for the standardized utility
  x_thisdose <- ((u1*num_notox_eff + u2*num_notox_noeff+
    u3 * num_tox_eff + u4*num_tox_noeff)/100)[current_dose]
  a_U[current_dose] <- U_prior[1] + x_thisdose
  b_U[current_dose] <- U_prior[2] + num_pts[current_dose] - x_thisdose
  
  # update desirability for each dose P(u(d) > u_b|data)
  prob_desirable <- pbeta(q=u_b,shape1 = a_U,shape2 = b_U,lower.tail = FALSE)
  
  # update posterior beta parameters
  a_T[current_dose] <- a_T[current_dose]+sum(tox); b_T[current_dose] <- b_T[current_dose]+sum(1-tox)
  a_E[current_dose] <- a_E[current_dose]+sum(eff); b_E[current_dose] <- b_E[current_dose]+sum(1-eff)
  # calculated observed DLT rate
  DLT_rate_obs[current_dose] <- (num_tox_eff+num_tox_noeff)[current_dose]/num_pts[current_dose]
  
  # define admissible doses
  prob_excess_tox <- pbeta(q=phi_t,shape1 = a_T, shape2 = b_T,lower.tail = FALSE)
  prob_futile_eff <- pbeta(q=phi_e,shape1 = a_E, shape2 = b_E)
  admissible_doses <- which(prob_excess_tox < c_t & prob_futile_eff < c_e)
  
  # Get next recommended dose
  if(DLT_rate_obs[current_dose] >= lambda_d)
  {
    next_dose <- current_dose - 1
  } else if ((DLT_rate_obs[current_dose] <= lambda_e) | (DLT_rate_obs[current_dose] > lambda_e & DLT_rate_obs[current_dose] < lambda_d & a_T[current_dose] + b_T[current_dose] - 2 < n_star1))
  {
    # choose a dose from (j-1,j,j+1) with maximum utility and admissible dose pool
    candidate_doses <- c(current_dose - 1, current_dose, current_dose +1)
    # Keep only valid dose levels
    candidate_doses <- candidate_doses[
      candidate_doses >= 1 & candidate_doses <= ndoses
    ]
    # Keep only admissible doses
    candidate_doses <- candidate_doses[which(candidate_doses %in% admissible_doses)]
    next_dose <- ifelse(length(candidate_doses)>0,
                        candidate_doses[which.max(prob_desirable[candidate_doses])],
                        NA)
  } else if (DLT_rate_obs[current_dose] > lambda_e & DLT_rate_obs[current_dose] < lambda_d & a_T[current_dose] + b_T[current_dose] - 2 >= n_star1)
  {
    # choose a dose from (j-1,j) with maximum utility and admissible dose pool
    candidate_doses <- c(current_dose - 1, current_dose)
    # Keep only valid dose levels
    candidate_doses <- candidate_doses[
      candidate_doses >= 1 & candidate_doses <= ndoses
    ]
    # Keep only admissible doses
    candidate_doses <- candidate_doses[which(candidate_doses %in% admissible_doses)]
    next_dose <- ifelse(length(candidate_doses)>0,
                        candidate_doses[which.max(prob_desirable[candidate_doses])],
                        NA)
  }
  
  
  # Extra dose-exploration rule
  # If the number of patients treated at the current dose d is greater than some value, n_star2 here
  # The observed toxicity rate of this dose is less than the de-escalation boundary lambda_d (from the univariate BOIN)
  # The next higher dose has never been used for treating patients, 
  # we escalate the dose for treating the next cohort of patients.
  if(!is.null(n_star2))
  { 
    if(sum(dose_history==current_dose) >= n_star2 & max(dose_history) < ndoses)
    {
      if(DLT_rate_obs[current_dose] < lambda_d)
      {
        next_dose <- current_dose + 1
      }
    }
  }
  
  # Check for stopping
  if ((is.na(next_dose) | next_dose < 1 | next_dose > ndoses) & max(cohort_ids) < N) {
    cat("Trial stopped early at", max(cohort_ids), "patients.\n")
    stopped <- 1
    break
  }
  
  # Update dose and patient ID
  current_dose <- next_dose
  p_id <- max(cohort_ids) + 1
} # end of while

# prepare output
out <- list()
# recommended phase 2 dose (RP2D)
# first identify the MTD
# estimate MTD from trial data using Isotonic regression
MTD <- BOIN::select.mtd(target = phi_t,
                            npts = a_T + b_T  -2,
                            ntox = a_T - 1,
                            cutoff.eli = c_t,
                            boundMTD = TRUE)$MTD
candidate_doses <- admissible_doses[admissible_doses<=MTD]
if(length(candidate_doses) > 0 & stopped == 0)
{
  # posterior mean of the utility
  utility_alldose <- a_U/(a_U+b_U)
  # find the dose that is admissible, below the MTD and with largest mean utility
  RP2D <- candidate_doses[which.max(utility_alldose[candidate_doses])]
} else
{
  RP2D <- NA
}
# RP2D <- ifelse(stopped,NA,which.max(utility_alldose))
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
