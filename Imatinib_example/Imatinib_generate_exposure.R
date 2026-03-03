# Read necessary packages
library(mrgsolve)
library(dplyr)
library(MASS)
#Imatinib doses: 400, 600, 800, 1000
# 2025/4/22 update: modified population characteristics

## Create a population data set ##
imatinib_exp <- function(dose_levels,
                         N,
                         population = "specific",
                         rho = -0.5,
                         time=240){
  if(population=="specific")
  {
    # Create population with individual-specific albumin and WBC levels
    # population_specific <- data.frame(ID = 1:N,  # Adjust to wanted number of IDs
    #                            ALB = rnorm(N, mean = 35, sd = 6.5),  # Generate random albumin level, based on population in article of Demetri et al.
    #                            WBC = rnorm(N, mean = 7.76, sd = 2.91)) # Generate random WBC count (x10^9) level, based on population in article of Demetri et al.
    mu <- c(ALB = 35.5, WBC = 5.59) # Mean levels
    sd_ALB <- 5.0
    sd_WBC <- 1.86
    # rho is Correlation factor, chosen based on plots in article of Demetri et al.
    
    # Create covariance matrix
    Sigma <- matrix(c(sd_ALB^2, rho * sd_ALB * sd_WBC,
                      rho * sd_ALB * sd_WBC, sd_WBC^2), 
                    nrow = 2)
    
    # Generate the population
    data <- mvrnorm(n = N, mu = mu, Sigma = Sigma)
    population <- data.frame(ID = 1:N, ALB = data[,1], WBC = data[,2])
  } else if (population=="basic")
  {
    population <- data.frame(ID = 1:N,  # Adjust to wanted number of IDs
               ALB = 38.3, # Use the standard number
               WBC = 7) # Use the standard value
  } else
  {
    stop("Please supply either population = \"basic\" or \"specific\"! ")
  }
  data_exposure <- sapply(dose_levels, function(dose)
    {
    # Define events for all individuals
    events <- ev_expand(amt = dose, cmt = "GUT", ii = 24, addl = 9, ID = 1:N) # Adjust to wanted dosing regimen
    
    # Merge individual-specific weights and ages with the events
    merged_df <- merge(population, events, by = "ID") # Adjust population if wanted
    
    # Set times for simulation
    times  <- tgrid(start = 0, end = 240, delta = 1) # Adjust if wanted
    
    # Simulation code
    PK_res <- mrgsim_d(mod, data = merged_df, tgrid = times)
    
    # Extract Cmin levels
    PK_res <- as.data.frame(PK_res) |>
      mutate(Cmin = IPRED*1000)
    # option to output Cmin at all times
    if(time == "all")
    {
      Cmin <- PK_res[, c("ID","time","Cmin")] 
    } else
    {
      Cmin <- PK_res[PK_res$time == time, c("ID","time","Cmin")] # Adjust to other Cmin points if wanted
    }
    ##########################
    # need to be changed later
    ##########################
    # replace NAN values with mean
    Cmin$Cmin[is.na(Cmin$Cmin)] <- mean(Cmin$Cmin,na.rm = TRUE)
    
    # Merge with patient characteristics
    out <- merge(merged_df[,c("ID","ALB","WBC","amt")],Cmin,by = "ID") |>
      rename(dose = amt) |>
      dplyr::select(ID,WBC,ALB,dose,time,Cmin)
    return(out)
  },simplify = FALSE)
  data_exposure <- bind_rows(data_exposure) |>
    arrange(ID,dose)
  return(data_exposure)
}
