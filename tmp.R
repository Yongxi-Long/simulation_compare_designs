target <- 0.25
skeleton <- c(0.05, 0.15, 0.25, 0.4, 0.6)
outcomes <- '2NN 3NN 4TT'
fit1 <- stan_crm(outcomes, skeleton = skeleton, target = target, 
                 model = 'empiric', beta_sd = sqrt(1.34), seed = 123,
                 beta_mean=0,alpha_mean = 0.3)
fit1$recommended_dose
fit1$prob_mtd
a=as.data.frame(fit1)
# probability of p(dose1) > 0.3
mean(a$`prob_tox[3]`>0.3)
