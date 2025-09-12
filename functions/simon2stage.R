run_simon2stage <- function(data,
                             design_para,
                             suspend_threshold = NULL)
{
  # check if design parameters are complete
  check <- sapply(c("r1","r2","n1","N"), function(i) is.null(design_para[[i]]))
  if(any(check))
  {
    stop("Design parameters for simon2stage are not complete!")
  }
  # extract number of doses 
  ndoses <- length(data$profile_info$dose)
  # extract design parameters
  r1 <- design_para$r1; r2 <- design_para$r2
  n1 <- design_para$n1; N <- design_para$N
  # chop data
  data_rsp <- data$data_rsp[1:N,]
  # evaluate trial under each dose
  # first see if there is trial suspension possibility
  # initialize trial suspension indicator
  suspend_tox <- when_suspend_tox<- rep(0,ndoses)
  if(!is.null(suspend_threshold))
  { 
    if(is.null(data$data_ph2_DLT))
    {
      stop("Need to supply toxicity data if want to implement dose modification!")
    } else
    {
      data_tox <- data$data_ph2_DLT[1:N,] |> as.matrix()
      tmp <- first_crossing(data_tox,suspend_threshold)
      suspend_tox <- !is.na(tmp)
      # see when the trial is suspended due to excess toxicity for each dose
      when_suspend_tox[suspend_tox] <- tmp[suspend_tox]
    }
  }
  # if stopped for futility after 1st stage
  stop_futility <- cumsum(data_rsp)[n1,] <= r1
  # if make it into stage 2 and success at the end
  reject_H0 <- cumsum(data_rsp)[N,] > r2 & !stop_futility
  # check if already suspended due to toxicity
  if(sum(suspend_tox)>0)
  {
    reject_H0[suspend_tox] <- 0
  } 
  # record num of patients per dose
  num_pts_per_dose <- (1-suspend_tox)*(n1*stop_futility + N*(!stop_futility))+suspend_tox*when_suspend_tox

  # prepare output
  out <- list()
  out$reject_H0 <- reject_H0
  out$suspend_tox <- suspend_tox
  out$num_pts_per_dose <- num_pts_per_dose
  return(out)
}


# helper function to design a Simon's two stage trial
simon2stage_design <- function(p0 = NULL,
                               p1 = NULL,
                               alpha = 0.05, beta = 0.2,
                               method = "Minimax",
                               nmax = 200)
{
  # get design parameters
  design_para <- tryCatch(ph2simon(pu=p0,
                                   pa=p1,
                                   ep1=alpha,
                                   ep2=beta,
                                   nmax = nmax),
                          error = function(e) {
                            e$message
                          })
  if(class(design_para)!="ph2simon")
  {
    stop(design_para)
  } else
  {
    #summary(ph2design)
    # use the corresponding minimax or optimal design specified
    # by the method parameter
    r1 <- twostage.admissible(design_para)[method,"r1"]
    n1 <- twostage.admissible(design_para)[method,"n1"]
    r2 <- twostage.admissible(design_para)[method,"r"]
    N <- twostage.admissible(design_para)[method,"n"]
  }
  return(list(r1=r1,n1=n1,r2=r2,N=N))
}

# helper function to track the cumulative mean as you move down rows,
# and find the first row index where it exceeds threshold X.
first_crossing <- function(mat, threshold) {
  N <- nrow(mat)
  p <- ncol(mat)
  cross_at <- rep(NA, p)
  
  # Loop over columns
  for (j in 1:p) {
    cum_means <- cumsum(mat[, j]) / seq_len(N)
    idx <- which(cum_means > threshold)[1]  # first crossing
    if (!is.na(idx)) {
      cross_at[j] <- idx
    }
  }
  
  return(cross_at)
}

# Example
set.seed(123)
mat <- matrix(rbinom(50, 1, 0.3), ncol = 5)  # 10 rows × 5 columns
threshold <- 0.4
first_crossing(mat, threshold)

