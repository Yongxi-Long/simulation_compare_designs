# to get the exposure-toxicity curve, we need two points of (Cmin,DLT rate)
# we make the assumption that at the maximum Cmin of 4200, the DLT rate is 80%
# another point is the DLT rate at the toxicity threshold Cmin = 3180,
# we can try different values to tun the steepness of the curve
DLT_3180 <- 0.6
DLT_4200 <- 0.85
tox_system <- function(k)
{
  a <- -5.3035
  Ctox <- 3180
  Cmax <- 4200
  Ztox <- a + k*(Ctox/1000)
  Zmax <- a + k*(Cmax/1000)
  DLT_tox <- DLT_4200  / ((1 + exp(-Ztox))*(1 - exp(-Zmax)))
  return(c(
           DLT_tox - DLT_3180))
}
# solve for a,k using optim
k <- uniroot(f=tox_system,interval = c(-5,5))$root
k
