# The calculations assume that the system is at equilibrium. Mean at-age parameters are provided.

rm(list=ls())
library(ggplot2)
source("functions.R") # contains functions "survivorship_F" and "MSYcalc"

# Read in the data sets

D <- read.csv("Pre_workshop/inputs.csv")

# at-age variables from age 1 to 9+ group

w_age <- D$Weight # weight-at-age in kg; ages 1 to 9+
m_age <- D$Maturity # maturity; ages 1 to 9+
v_age <- D$Vulnerability # vulnerability or selectivity; ages 1 to 9+

# other Parameters

M <- D$M[1] # assumed natural mortality rate
h <- D$h[1] # assumed steepness of the Beverton-Holt stock recruitment relationship
R0 <- D$R0[1] # model estimated unfished equilibrium recruitment (units = millions of age 1 recruits)

########################################################################################################################
# Move h/R0/phi0 parameterization of B-H SRR to a/b parameterization 
########################################################################################################################

# survivorship_F: function from functions.R (calculates survivorship-at-age from f = F (default is F=0), M = M, n_ages = number of age classes (from age of recruitment to plus group), sel = vulnerability-at-age when F > 0)
# survivorship_F <- function(f=0,M,n_ages,sel)
l_age <- survivorship_F(M=M,n_ages=length(w_age)) # unfished survivorship-at-age (F=0)
phi0 <- sum(l_age*w_age*m_age) # unfished spawning biomass per recruit

# B-H a and b
BHa <- 4*h/(phi0*(1-h)) # estimated Beverton-Holt a
BHb <- 1/R0*(BHa-1/phi0) # estimated Beverton-Holt b

########################################################################################################################
# Calculate equilibrium SSBmsy
########################################################################################################################

# MSYcalc: 
# function from functions.R returns a list with Fmsy, msy, SSBmsy from:
# M=M, waa=weight-at-age, mat=maturity-at-age, sel=vulnerability-at-age, Beverton-Holt a and b
# MSYcalc <- function(M,waa,mat,sel,a,b)

calc <- MSYcalc(M=M,waa=w_age,mat=m_age,sel=v_age,a=BHa,b=BHb) #biomass units (kg X million) = kt

#####################################################################################################################
# Data Frame for SPR
########################################################################################################################

# set up a data frame with various F, calculate spawning biomass per recruit (phiF) and SPR for various F
DF <- data.frame(f=seq(0.01,3,0.01),phiF=NA,SPR=NA)

for(i in 1:nrow(DF)){
  l_f_age <- survivorship_F(f=DF$f[i],M=M,n_ages=length(v_age),sel=v_age,message=F)
  DF$phiF[i] <- sum(l_f_age*w_age*m_age)
}

DF$SPR <- DF$phiF/phi0

########################################################################################################################
# Solutions
########################################################################################################################

#2. phi0 (spawning biomass per recruit)

phi0

#4.SPR at F = 0.2 and F = 0.5

DF$SPR[DF$f==0.2]
DF$SPR[DF$f==0.5]

#5. F at SPR40% and 35%

DF$abs_SPR_minus40 <- abs(DF$SPR - 0.4)
DF$abs_SPR_minus35 <- abs(DF$SPR - 0.35)

f_at_SPR_40 <- DF$f[which(DF$abs_SPR_minus40 == min(DF$abs_SPR_minus40))]
f_at_SPR_35 <- DF$f[which(DF$abs_SPR_minus35 == min(DF$abs_SPR_minus35))]

#6. For M = 0.2 and assuming no stock recruitment relationship

Eq_SSB_at_SPR_40 = DF$phiF[DF$f==f_at_SPR_40] * R0
Eq_SSB_at_SPR_35 = DF$phiF[DF$f==f_at_SPR_35] * R0

#6. For M = 0.3 and assuming no stock recruitment relationship

l_age_M.3 <- survivorship_F(M=0.3,n_ages=length(w_age)) # unfished survivorship-at-age (F=0)
phi0_M.3  <- sum(l_age_M.3 *w_age*m_age)

DF_M.3 <- data.frame(f=seq(0.01,3,0.01),phiF=NA,SPR=NA)

for(i in 1:nrow(DF_M.3)){
  l_f_age <- survivorship_F(f=DF_M.3$f[i],M=0.3,n_ages=length(v_age),sel=v_age,message=F)
  DF_M.3$phiF[i] <- sum(l_f_age*w_age*m_age)
}

DF_M.3$SPR <- DF_M.3$phiF/(phi0_M.3)


DF_M.3$abs_SPR_minus40 <- abs(DF_M.3$SPR - 0.4)
DF_M.3$abs_SPR_minus35 <- abs(DF_M.3$SPR - 0.35)

f_at_SPR_40_M.3 <- DF_M.3$f[which(DF_M.3$abs_SPR_minus40 == min(DF_M.3$abs_SPR_minus40))]
f_at_SPR_35_M.3 <- DF_M.3$f[which(DF_M.3$abs_SPR_minus35 == min(DF_M.3$abs_SPR_minus35))]

Eq_SSB_at_SPR_40_M.3 <- DF_M.3$phiF[DF_M.3$f==f_at_SPR_40_M.3] * R0
Eq_SSB_at_SPR_35_M.3 <- DF_M.3$phiF[DF_M.3$f==f_at_SPR_35_M.3] * R0


#7. SSB0

SSB0 = R0/(BHa-BHb*R0)

BHa_h.9 <- 4*0.9/(phi0*(1-0.9)) # estimated Beverton-Holt a
BHb_h.9 <- 1/R0*(BHa_h.9-1/phi0) # estimated Beverton-Holt b

SSB0_h.9 = R0/(BHa_h.9-BHb_h.9*R0) #assume R0 stays constant

BHa_M.3 <- 4*h/(phi0_M.3*(1-h)) # estimated Beverton-Holt a
BHb_M.3 <- 1/R0*(BHa_M.3-1/phi0_M.3) # estimated Beverton-Holt b

SSB0_M.3 <- R0/(BHa_M.3-BHb_M.3*R0) #assume R0 stays constant

#8. SSB at 50% Rmax

Rmax <- BHa/BHb
SSB_50rmax <- 0.5*Rmax/(BHa-BHb*0.5*Rmax)

Rmax_h.9 <- BHa_h.9 /BHb_h.9 
SSB_50rmax_h.9  <- 0.5*Rmax_h.9/(BHa_h.9-BHb_h.9*0.5*Rmax_h.9)

Rmax_M.3 <- BHa_M.3 /BHb_M.3 
SSB_50rmax_M.3  <- 0.5*Rmax_M.3/(BHa_M.3-BHb_M.3*0.5*Rmax_M.3)

#9. SSBmsy

SSBmsy <- calc$SSBmsy
SSBmsy_h.9 <- MSYcalc(M=M,waa=w_age,mat=m_age,sel=v_age,a=BHa_h.9,b=BHb_h.9)$SSBmsy
SSBmsy_M.3 <- MSYcalc(M=0.3,waa=w_age,mat=m_age,sel=v_age,a=BHa_M.3,b=BHb_M.3)$SSBmsy

#10. SSBmsy/SSB0

SSBmsy/SSB0
SSBmsy_h.9/SSB0_h.9
SSBmsy_M.3/SSB0_M.3
