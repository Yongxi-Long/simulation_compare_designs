#Exposure - hazard
exp_surv <- function(data_exposure,h_min=0.02,
                     median_surv0 = 20,
                     p_trt = 0.5,
                     a=0.3274917,b=-0.002728778,Cmax=4200,
                     maxt = 60){
  dose_levels <- unique(data_exposure[,"dose"])
  data_surv = sapply(dose_levels, function(d)
  {
    data_exposure_thisdose <- data_exposure[data_exposure$dose==d,]
    # randomize
    data_exposure_thisdose$trt <- rbinom(n=nrow(data_exposure_thisdose),
                                        size = 1,
                                        prob = p_trt)
    data_exposure_thisdose$hazard <- (data_exposure_thisdose$Cmin < Cmax)*(a*exp(b*data_exposure_thisdose$Cmin) + h_min)+
      (data_exposure_thisdose$Cmin >= Cmax)*h_min
      #a*exp(b*data_exposure_thisdose$Cmin)
      
    # control hazard is constant
    data_exposure_thisdose$hazard[data_exposure_thisdose$trt==0] <- log(2)/median_surv0
    # generate time-to-event from an exponential survival
    data_exposure_thisdose$eventtime <- rexp(n=nrow(data_exposure_thisdose),rate = data_exposure_thisdose$hazard)
    data_exposure_thisdose$status <- as.integer(data_exposure_thisdose$eventtime <= maxt)
    data_exposure_thisdose$eventtime[data_exposure_thisdose$eventtime>=maxt] <- maxt
    # create survival data frame
    df <- data_exposure_thisdose |>
      rename(id = ID) |>
      dplyr::select(id,trt,eventtime,status)
    return(df)
  },simplify = FALSE)
return(data_surv)
}

