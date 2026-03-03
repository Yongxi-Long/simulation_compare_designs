tox_fun_shortterm <- function(Cmin)
{
  L <- 0.5
  a <- -5.5
  k <- 2
  Cmax <- 4200
  Z1 = a + k*((Cmin)/1000)
  Z2 = a + k*(Cmax/1000)
  DLT_rate <- L  / ((1 + exp(-Z1))*(1 - exp(-Z2)))
  return(DLT_rate)
}
tox_fun_tmp <- function(Cmin)
{
  L <- 0.6
  a <- -6
  k <- 2.5
  Cmax <- 4200
  Z1 = a + k*((Cmin)/1000)
  Z2 = a + k*(Cmax/1000)
  DLT_rate <- L  / ((1 + exp(-Z1))*(1 - exp(-Z2)))
  return(DLT_rate)
}

ggplot(data.frame(x=c(0,6000)),aes(x=x))+
  geom_function(fun=tox_fun_shortterm,linetype=1,color="#49548A")+
  geom_function(fun=tox_fun_tmp,linetype=2,color="#B14E56")+
  theme_bw()
