#Exposure - Response profile 
exp_rsp <- function(data_exposure,Cmax=4200,k=1.045263,
                    L = 0.7)
{
  dose_levels <- unique(data_exposure[,"dose"])
  data_rsp = data.frame(sapply(dose_levels, function(d)
  {
    Cmin <- data_exposure[data_exposure$dose==d,] |>
      dplyr::select(Cmin) |>
      unlist() |>
      as.vector()
    rsp_rates <- L*(1-exp(-k*Cmin/1000))/(1-exp(-Cmax/1000*k))
    rsp_outcomes <- rbinom(n=length(rsp_rates),size=1,prob = rsp_rates)
    return(rsp_outcomes)
  }))
  colnames(data_rsp) <- paste0("dose",dose_levels)
  return(data_rsp)
}

