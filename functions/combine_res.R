combine_res <- function(res_ph1,
                        res_ph2=NULL,
                        res_ph3=NULL)
{
  # get number of doses
  ndoses <- length(res_ph1[[1]]$summary_tab$dose_level)
  # extract RP2D from res_ph1
  RP2D <- sapply(res_ph1,function(x) x$RP2D)
  # distribution of RP2D from phase I
  RP2D_tab <-  table(factor(RP2D,levels = 1:ndoses), useNA = "ifany") |>
    as.data.frame() 
  colnames(RP2D_tab) <- c("dose_level","freq_selected")
  RP2D_tab <- RP2D_tab |>
    mutate(prob_selected = freq_selected/sum(freq_selected))
  if(!is.null(res_ph2))
  {
    # see for which iteration a dose is recommend (not NA)
    idx <- !is.na(RP2D)
    # Phase II rejection matrix for each iteration (row) and each dose (column)
    ph2_success_per_dose <- t(sapply(res_ph2, function(i) i$reject_H0))
    # conditional (given Phase I) Phase II rejection vector per drug
    rej_ph12 <- rep(NA,length(RP2D))
    rej_ph12[idx] <- ph2_success_per_dose[cbind(which(idx),RP2D[idx])]
    PoS_ph12 <- mean(rej_ph12,na.rm=TRUE)
    
    # probability of trial suspension due to excess toxicity in Phase II
    ph2_suspension_per_dose <- t(sapply(res_ph2, function(i) i$suspend_tox))
    TS_ph12 <- rep(NA,length(RP2D))
    TS_ph12[idx] <- ph2_suspension_per_dose[cbind(which(idx),RP2D[idx])]
    probTS_ph12 <- mean(TS_ph12,na.rm=TRUE)
    if(!is.null(res_ph3))
    {
      # Phase III rejection matrix for each iteration (row) and each dose (column)
      ph3_success_per_dose <- t(sapply(res_ph3, function(i) i$reject_H0))
      # conditional (given Phase I) Phase III rejection vector per drug
      rej_ph13 <- rep(NA,length(RP2D))
      rej_ph13[idx] <- ph3_success_per_dose[cbind(which(idx),RP2D[idx])]
      PoS_ph123 <- mean(rej_ph12 & rej_ph13, na.rm=TRUE)
      
      # probability of trial suspension due to excess toxcity in phase 2&3
      ph3_suspension_per_dose <- t(sapply(res_ph3, function(i) i$suspend_tox))
      TS_ph13 <- rep(NA,length(RP2D))
      TS_ph13[idx] <- ph3_suspension_per_dose[cbind(which(idx),RP2D[idx])]
      probTS_ph123 <- mean(TS_ph12 | TS_ph13,na.rm=TRUE)
      
      out <- list("Table of RP2D" = RP2D_tab,
                  "Phase I-II probability of significance"=PoS_ph12,
                  "Phase I-III probability of significance"=PoS_ph123,
                  "Phase II probability of toxicity suspension" = probTS_ph12,
                  "Phase II&III probability of toxicity suspension" = probTS_ph123)
    } else
    {
      out <- list("Table of RP2D" = RP2D_tab,
                  "Phase I-II probability of significance"=PoS_ph12,
                  "Phase II probability of toxicity suspension" = probTS_ph12)
    }
  } else
  {
    out <- list("Table of RP2D" = RP2D_tab)
  }
  
  return(out)
}