#Exposure - Toxicity profile 
exp_tox <- function(data_exposure,Cmax=4200,type="short-term")
{
  dose_levels <- unique(data_exposure[,"dose"])
  if(type=="short-term")
  {
    # L <- 0.5
    # a <- -5.5
    # k <- 2
    L <- 0.6
    a <- -6
    k <- 2.5
  } 
  else if (type == "long-term")
  {
    L <- 0.95
    a <- -4.5
    # # k <- 1.7202
    k <- 1.86
  } else
  {
    stop("Unknown toxicity type, must be either 'short-term' or 'long-term'!")
  }
  Z2 = a + k*(Cmax/1000)
  data_DLT = data.frame(sapply(dose_levels, function(d)
    {
    Cmin <- data_exposure[data_exposure$dose==d,] |>
      dplyr::select(Cmin) |>
      unlist() |>
      as.vector()
    Z1 = a + k*((Cmin)/1000)
    DLT_rates <- L  / ((1 + exp(-Z1))*(1 - exp(-Z2)))
    DLT_outcomes <- rbinom(n=length(DLT_rates),size=1,prob = DLT_rates)
    return(DLT_outcomes)
  }))
  colnames(data_DLT) <- paste0("dose",dose_levels)
  return(data_DLT)
}

