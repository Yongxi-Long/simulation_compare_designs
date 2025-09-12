#' @title Generate specific profile
#' @description
#' The function extract corresponding DLT probability/Response probability/Median survival
#' from supplied parametric curves and selected doses
#' @name gen_profile
#' @param drug_profile drug profile, should be a function that takes a dose
#' value as input and outputs the DLT probability/Response probability/Median survival
#' @param dose_levels a vector of dose levels to be investigated
#' @return a vector of DLT probability/Response probability/Median survival for the supplied dose levels
#' @export gen_profile

gen_profile <- function(drug_profile,dose_levels)
{
  # source the corresponding profile function based on drug profile name
  source(paste0("profiles/",drug_profile,".R"))
  sapply(dose_levels,function(dose) do.call(drug_profile,list(dose)))
}

#' @title Simulate Phase I trial data
#' @description
#' The function generate a N*J matrix of DLT indicators based on supplied
#' dose-toxicity profile and sample size. Constraint is that if the patient has 
#' DLT at j-th dose level, he/she must also have DLT at higher dose levels
#' @name gen_DLT
#' @param nsim number of iterations for the simulation
#' @param DLT_rates A vector of length J, the j-th element specifies the 
#' probability of experiencing dose limiting toxicity (DLT) for the j-th dose
#' User input or the output of function gen_profile
#' @param drug_profile have to provide a function of dose-toxicity profile and 
#' a vector of dose levels if do not supply the DLT_rates directly, the function
#' then internally calls gen_profile() if dose_levels is also supplied
#' @param dose_levels a vector of dose levels to be investigated
#' @param N sample size
#' @export gen_DLT
#' @returns data_list
#' A list of data frames of length nsim. Each data frame has N rows and J 
#' columns, where J is the number of dose levels. The n-th row and the j-th entry 
#' specifies whether the n-th patient has a DLT (0/1) at the j-th dose level.
gen_DLT <- function(nsim,
                    DLT_rates=NULL,
                    drug_profile=NULL,dose_levels=NULL,
                    N) 
{
  if(is.null(DLT_rates))
  {
    if(is.null(drug_profile) | is.null(dose_levels))
    {
      stop("Please provide either a vector of DLT rates or a dose-toxicity profile function and selected dose levels")
    } else
    {
      DLT_rates <- gen_profile(drug_profile,dose_levels)
    }
  }
  J <- length(DLT_rates)
  
  data_list <- replicate(nsim,
                         {
                           DLT_mat <- matrix(NA,nrow = N,ncol = J)
                           for(j in J:1)
                           {
                             if(j==J)
                             {
                               DLT_mat[,J] <- rbinom(n=N,size=1,prob =DLT_rates[J])
                             } else
                             {
                               # DLT at this dose is a subset from next higher dose
                               tox_this_dose <- DLT_mat[,j+1]
                               DLT <- DLT_rates[j]/DLT_rates[j+1]
                               tox_this_dose[tox_this_dose==1] <- rbinom(n=sum(tox_this_dose),size=1,prob=DLT)
                               DLT_mat[,j] <- tox_this_dose
                             }
                           }
                           colnames(DLT_mat) <- paste0("dose",1:J)
                           return(data.frame(DLT_mat))
                         },simplify = FALSE)
  return(data_list)
}
# check
# tmp <- gen_DLT(nsim=10,
#                DLT_rates = c(0.1,0.2,0.22,0.3),
#                N=30)

#' @title Simulate Phase II trial data
#' @description
#' The function generate a N*J matrix of objective response indicators based on supplied
#' dose-response profile and sample size. 
#' @name gen_rsp
#' @param nsim number of iterations for the simulation
#' @param rsp_rates The output of function gen_profile:
#' a vector of length J, the j-th element specifies the 
#' probability of having a response for the j-th dose
#' @param drug_profile have to provide a function of dose-toxicity profile and 
#' a vector of dose levels if do not supply the DLT_rates directly
#' @param dose_levels a vector of dose levels to be investigated
#' @param N sample size
#' @export gen_rsp
#' @returns data_list
#' A list of data frames of length nsim. Each data frame has N rows and J 
#' columns, where J is the number of dose levels. The n-th row and the j-th entry 
#' specifies whether the n-th patient has a response (0/1) at the j-th dose level.
gen_rsp <- function(nsim,
                    rsp_rates=NULL,
                    drug_profile=NULL,dose_levels=NULL,
                    N) {
  
  if(is.null(rsp_rates))
  {
    if(is.null(drug_profile) | is.null(dose_levels))
    {
      stop("Please provide either a vector of response rates or a dose-response profile function and selected dose levels")
    } else
    {
      rsp_rates <- gen_profile(drug_profile,dose_levels)
    }
  }
  J <- length(rsp_rates)
  data_list <- replicate(nsim,
                       {
                         rsp_mat <- sapply(rsp_rates, function(p)
                           rbinom(n=N,size = 1,prob = p))
                         colnames(rsp_mat) <- paste0("dose",1:J)
                         return(data.frame(rsp_mat))
                       },simplify = FALSE)
  return(data_list)
}
# check
# tmp <- gen_rsp(nsim = 10,
#                rsp_rates = c(0.2,0.23,0.4,0.44),
#                N=142)

#' @title Simulate Phase III trial data
#' @description
#' The function generate a data frame containing the individual 
#' survival data. The function internally calls simsurv with one extra parameter 
#' that specifies the randomization ratio
#' @name  gen_surv
#' @param nsim number of iterations for the simulation
#' @param N sample size
#' @param median_surv0: median survival for the control group. Then internally
#' the function will convert it to hazard assuming the survival probability follows
#' an exponential distribution
#' @param median_surv1: (a vector of) median survivals for the treatment group. The j-th 
#' element gives the median survival under the j-th dose. 
#' @param use_simsurv: default is FALSE. If users wish to generate more complicated
#' survival data, they can set this to be TRUE and supply parameters for the simsurv()
#' function. 
#' Note: if users wish to simulate survival data for multiple doses using simsurv() function,
#' then all the corresponding parameters should be lists, with the j-th element
#' corresponding to parameters of the j-th dose
#' @param ndose number of doses, only needed if set use_simsurv = TRUE. Otherwise the
#' length of median_surv1 is taken as the number of doses.
#' @param p_trt probability of being randomized to the treated group,
#'  default is 0.5
#' @param maxt maximum follow-up time (will apply administrative censoring for survival
#' times larger than this)
#' @... Other parameters see documentation of simsurv()
#' @export gen_surv
#' @returns data_list
#' A list (length nsim) of a list (length J, J being the number of doses) of data frames.
#' Each data frame has four columns
#' id: patient id
#' trt: treatment assignment
#' eventtime: time to event
#' status: event indicator
require(simsurv)
gen_surv <- function(nsim,
                     N,
                     median_surv0 = NULL,
                     median_surv1 = NULL,
                     p_trt = 0.5,
                     use_simsurv = FALSE,
                     ndose = NULL,
                     maxt,
                     ...
                     )
{
  # get all input arguments
  call <- match.call(expand.dots = TRUE)
  # generate covariates (1:1 randomization)
  covs_df <- replicate(nsim,
                       expr =
                         {
                           data.frame(id=1:N,
                                      trt=rbinom(N,size=1,prob = p_trt))
                         },simplify = FALSE)
  if(use_simsurv)
  {
    default_args <- formals(simsurv)
    call_args <- sapply(call[2:length(call)], function(i)
      eval(i,parent.frame())) # first one is the function name, not needed
    # change default values to user-input ones, if there is any
    # extract supplied parameters for each dose
    call_args_alldoses = sapply(1:ndose, function(j)
    {
      sapply(call_args, function(i) {
        arg <- tryCatch(i[[j]],
                        error = function(e) NA)
        if(is.na(arg))
        {
          arg <- i[[1]]
        }
        return(arg)
      },simplify = FALSE)
    },simplify = FALSE)
  final_args_alldoses <- sapply(1:ndose, function(j)
  {
    final_args <- default_args
    final_args[names(call_args)] <- call_args_alldoses[[j]]
    # keep arguments relevant for simsurv only
    final_args <- final_args[names(default_args)]
    return(final_args)
  },simplify = FALSE)
  data_list <- sapply(covs_df, function(x)
    {
    sapply(final_args_alldoses, function(final_args_thisdose)
    {
      final_args_thisdose$x <- x
      data <- do.call(simsurv::simsurv,final_args_thisdose)
      data <- dplyr::inner_join(data,x,by="id")
      return(data)
    },simplify = FALSE)
  },simplify = FALSE)
# end if(use_simsurv)
} else
{
  # use median survival and exponential survival
  # hazard of control
  hazard0 <- log(2)/median_surv0
  # hazard(s) of treatment
  hazard1 <- log(2)/median_surv1
  # betas: log hazard ratios
  betas <- log(hazard1/hazard0)
  data_list <- sapply(covs_df, function(x)
    {
    sapply(betas, function(i)
    {
      names(i) <- "trt"
      data <- simsurv::simsurv(dist = "exponential",
                      lambdas = hazard0,
                      betas = i,
                      x = x,
                      maxt = maxt)
      data <- dplyr::inner_join(data,x,by="id")
      return(data)
    },simplify = FALSE)
  },simplify = FALSE)
}
return(data_list)
}