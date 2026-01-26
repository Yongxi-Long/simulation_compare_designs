#' @title RCT design evaluation for phase 3
#' @description
#' Evaluate the performance under a randomized controlled trial design
#'  description
#' @param data A list of J data frames for the J doses containing the Phase 3 survival data to be evaluated,
#' must have the following four columns:
#' id: patient id
#' trt: treatment assignment
#' eventtime: time to event
#' status: event indicator
#' @param design_para A named list of design parameters, can be either the output of
#' RCT_design or user-input
#' @export RCT_eval
#' @returns Performance metrics under RCT design
#' reject_H0: A vector of length J indicating whether H0 is rejected for each dose j (j=1,..,J)
#' hazard_ratio: A vector of length J giving the hazard ratio of treatment vs. control
#'  for each dose
#' median_surv: A vector of length J giving the median survival for each dose
#' rmean_surv: A vector of length J giving the restricted mean survival until the maximum followup time (maxt)
#' for each dose
#' num_pts: A vector of length J giving the number of patients treated at each dose
require(survival)
run_RCT <- function(data,
                   design_para,
                   quit_prob = NULL,
                   suspend_threshold = NULL,
                   suspend_window = 30)
{
  # check if design parameters are complete
  check <- sapply(c("alpha","N"), function(i) is.null(design_para[[i]]))
  if(any(check))
  {
    stop("Design parameters for RCT are not complete!")
  }
  # extract number of doses 
  ndoses <- length(data$profile_info$dose)
  # significance level for the test
  alpha <- design_para$alpha
  # number of patients
  N <- design_para$N
  
  # chop data
  data_surv <- sapply(data$data_surv, function(x)
    {
    x[1:N,]
  },simplify = FALSE)
  
  # initialize trial suspension indicator
  suspend_tox <- when_suspend_tox<- rep(0,ndoses)
  # initialize patient crossover proportion
  crossover_prop <- rep(0,ndoses)
  if(!is.null(suspend_threshold) | !is.null(quit_prob))
  {
    if(is.null(data$data_ph2_DLT))
    {
      stop("Need to supply toxicity data if want to implement dose modification!")
    } else if (!is.null(suspend_threshold) & !is.null(quit_prob))
    {
      stop("Can only implement either trial suspension or treatment crossover, not both!")
    } else
    {
      data_tox <- data$data_ph2_DLT[1:N,] 
      data_copy <- data_surv
      # patient with DLT randomly assigned to dose escalation, those on control can't have dose modification
      randomization_mat <- sapply(data_surv, function(data_thisdose) data_thisdose$trt)
      if(!is.null(quit_prob))
      {
        quit_index <- data_tox & matrix((runif(length(data_tox)) < quit_prob),
                                              nrow = nrow(data_tox),
                                              ncol = ncol(data_tox)) & randomization_mat
        # proportion of patients who went through crossover
        crossover_prop <- colMeans(quit_index)
        # exchange efficacy data to de-escalated dose after toxicity has crossed boundary
        if(any(crossover_prop>0))
        {
          for(j in 1:ndose)
          {
              # random sample from controls
              controls <- data_copy[[j]][data_copy[[j]]$trt==0,]
              data_surv[[j]][which(quit_index[,j]==TRUE),
                        c("eventtime","status")] <- controls[sample(nrow(controls),
                                                                    size=sum(quit_index[,j]==TRUE),
                                                                    replace = TRUE),
                                                             c("eventtime","status")]
          }
        }
      }
      if(!is.null(suspend_threshold))
      {
        idx <- lapply(data_tox, function(col)
        {
          cummean <- cumsum(col)/seq_along(col)
          idx <- which(cummean > suspend_threshold)
          idx <- idx[idx>=suspend_window]
          if(length(idx) > 0) idx <- min(idx)
        })
        suspend_tox <- sapply(idx, function(i) length(i) >0)
        # see when the trial is suspended due to excess toxicity for each dose
        when_suspend_tox[suspend_tox] <- idx |>
          unlist() |>
          as.numeric()
      }
    }
  } 
  
  # analysis survival data for each dose level
  res <- sapply(data_surv, function(data_thisdose)
    {
    
    fit_km <- survfit(Surv(eventtime, status) ~ trt, data = data_thisdose)
    mod_cox <- coxph(Surv(eventtime,status) ~ trt, data = data_thisdose,ties = "breslow")
    HR <- exp(mod_cox$coefficients["trt"]) |>
      unname()
    p_value <- summary(mod_cox)$logtest["pvalue"] |>
      unname()
    # significantly superior
    reject_H0 <- p_value < alpha & HR < 1
    rmean_surv <- summary(fit_km)$table[2,"rmean"]
    median_surv <- summary(fit_km)$table[2,"median"]
    return(c(reject_H0=reject_H0, 
             HR=HR, 
             median_surv = median_surv,
             rmean_surv = rmean_surv))
  })
  
  # prepare output
  out <- list()
  # record num of patients per dose
  num_pts_per_dose <- (1-suspend_tox)*N+suspend_tox*when_suspend_tox
  
  out$reject_H0 <- res["reject_H0",]
  out$suspend_tox <- suspend_tox
  out$crossover_prop <- crossover_prop
  out$num_pts_per_dose <- num_pts_per_dose
  # out$HR=res["HR",]
  # out$median_surv = res["median_surv",]
  # out$rmean_surv = res["rmean_surv",]
  return(out)
}


#' @title randomized controlled trial design for phase 3
#' @description Use the two-arm randomized controlled design. By default the outcome
#' is time-to-event 
#' @param alpha,beta type I and type II error rate. Power=1-beta 
#' @param p_trt probability of being randomized to the treatment arm
#' @param hazard0,hazard1 event hazard for the control and the treatment group respectively
#' @param t_FU time duration for study follow-up
#' @param t_accrual time duration for subject accrual
#' @param hazard_dropout drop out hazard (equal) for both groups, default is zero
#' @export RCT_design
#' @returns design parameters for RCT
#' N: total sample size
#' alpha: significance level of the test
require(gsDesign)
RCT_design <- function(alpha = 0.05, beta = 0.2,
                       p_trt = 0.5,
                       hazard0, hazard1,
                       t_FU, t_accrual,
                       hazard_dropout = 0)
{
  library(gsDesign)
  design_para <- nSurvival(alpha = alpha,
                           beta = beta,
                           lambda1 = hazard0,
                           lambda2 = hazard1,
                           Ts = t_FU,
                           Tr = t_accrual,
                           eta = hazard_dropout,
                           ratio = (1-p_trt)/p_trt)
  # get the sample size
  N <- round(design_para$n)
  return(list(N = N, alpha = alpha))
}