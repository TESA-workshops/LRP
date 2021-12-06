########################################################################################################################
# Exercise 2: SSB0 and SSBmsy calculations, plots
########################################################################################################################

rm(list=ls())
library(ggplot2)
source(paste0(getwd(),"/functions.R")) # contains functions "survivorship_F" and "MSYcalc"

# Read in the two data sets
AA <- read.csv(paste0(getwd(),"/Exercise 2/ex2_at_age_data.csv"))
D <- read.csv(paste0(getwd(),"/Exercise 2/ex2_data.csv"))

########################################################################################################################
# Move h/R0/phi0 parameterization of B-H SRR to a/b parameterization and calculate "unfished equilibrium SSB"
########################################################################################################################

# Model assumes a Beverton-Holt SRR with steepness h = 0.75
h <- 0.75 # model assumed steepness
R0 <- 2.395766 # model estimated equilibrium unfished recruitment [provided]
M <- 0.3 # model assumed natural mortality rate

# survivorship_F: function from functions.R (calculates survivorship-at-age from f=F (default is F=0), M=M, n_ages = number of age classes (including 0 and plus group), sel=vulnerability-at-age when F !=0)
# survivorship_F <- function(f=0,M,n_ages,sel)
l_age <- survivorship_F(M=M,n_ages=length(AA$w_age)) # survivorship-at-age (unfished, F=0)
phi0 <- sum(l_age*AA$w_age*AA$m_age) # unfished spawning biomass  per recruit

# B-H a and b
BHa <- 4*h/(phi0*(1-h)) # estimated Beverton-Holt a
BHb <- 1/R0*(BHa-1/phi0) # estimated Beverton-Holt b

# Unfished equilibrum SSB0
SSB0 <- R0*phi0

########################################################################################################################
# Calculate equilibrium SSBmsy and ratio of 
########################################################################################################################

# MSYcalc: function from functions.R returns a list with Fmsy, msy, SSBmsy from M=M, waa=weight-at-age, mat=maturity-at-age, sel=vulnerability-at-age, Beverton-Holt a and b
# MSYcalc <- function(M,waa,mat,sel,a,b)

calc <- MSYcalc(M=M,waa=AA$w_age,mat=AA$m_age,sel=AA$v_age,a=BHa,b=BHb)

Fmsy <- calc$Fmsy
msy <- calc$msy
SSBmsy <- calc$SSBmsy

#Ratio of SSBmsy to SSB0

SSBmsy/SSB0

########################################################################################################################
# Plots from data frame AA
########################################################################################################################

# set up a data frame with various F, calculate spawning biomass per recruit (phiF) and SPR for various F
f=seq(0.01,3,0.01)
DF <- data.frame(f=f,phiF=NA,SPR=NA,YPR=NA)
YPR_a <- SURV_a <- matrix(NA,ncol=length(AA$w_age),nrow=nrow(DF))

for(i in 1:nrow(DF)){
  SURV_a[i,] <- l_f_age <- survivorship_F(f=DF$f[i],M=M,n_ages=length(AA$v_age),sel=AA$v_age,message=F)
  DF$phiF[i] <- sum(l_f_age*AA$w_age*AA$m_age)
  for(j in 1:ncol(YPR_a)){
    YPR_a[i,j] <- SURV_a[i,j]*AA$w_age[j]*(1-exp(-(0.3+f[i]*AA$v_age[j])))*f[i]*AA$v_age[j]/(0.3+f[i]*AA$v_age[j])
  }
}

DF$YPR <- rowSums(YPR_a) 
DF$eq_rec_f <- (1/BHb*(BHa-1/DF$phiF))
DF$eq_rec_f[DF$eq_rec_f<0] <- 0
DF$Yield <- DF$YPR*DF$eq_rec_f 
DF$SPR <- DF$phiF/phi0
DF$eq_SSB <- DF$eq_rec_f/(BHa-DF$eq_rec_f*BHb)

DF$abs_SPR_minus40 <- abs(DF$SPR - 0.4)
DF$abs_SPR_minus50 <- abs(DF$SPR - 0.5)

f_at_SPR_40 <- DF$f[which(DF$abs_SPR_minus40 == min(DF$abs_SPR_minus40))]
f_at_SPR_50 <- DF$f[which(DF$abs_SPR_minus50 == min(DF$abs_SPR_minus50))]

Eq_SSB_at_SPR_40 = DF$phiF[DF$f==f_at_SPR_40] * R0
Eq_SSB_at_SPR_50 = DF$phiF[DF$f==f_at_SPR_50] * R0



LRP1 = 146
LRP21 = 233
LRP22 = 297
LRP3 = mean(c(min(D$SSB),median(D$SSB[26:50]),Eq_SSB_at_SPR_40,0.2*SSB0,0.4*calc$SSBmsy))
LRP4 = min(D$SSB)
LRP5 = 185

Rmax <- BHa/BHb
SSB_50rmax <- 0.5*Rmax/(BHa-BHb*0.5*Rmax)

LRP6 = mean(c(min(D$SSB),SSB_50rmax,0.2*SSB0,0.4*calc$SSBmsy))


#Plot
ggplot(D) + geom_path(mapping=aes(y=SSB,x=Year)) + theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0) +
  geom_hline(yintercept=LRP1, color = "purple") + 
  geom_hline(yintercept=LRP4, color = "purple") + 
  geom_hline(yintercept=LRP5, color = "purple") + 
  
  geom_hline(yintercept=LRP21, color = "blue") + 
  geom_hline(yintercept=LRP22, color = "blue") + 
  
  geom_hline(yintercept=LRP6, color = "green") +
  geom_hline(yintercept=LRP3, color = "green") +
  
  geom_hline(yintercept=0.2*SSB0, color = "red",linetype=2) +
  geom_hline(yintercept=0.4*calc$SSBmsy, color = "darkred",linetype=2) 
  


ggplot() + geom_point(D[!is.na(D$Rec),],mapping=aes(y=Rec,x=SSB),colour="blue") + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=c(0,12.5)) + expand_limits(x=0) + geom_function(fun=function(x) BHa*x/(1+BHb*x)) +
  geom_vline(xintercept=LRP1, color = "purple") + 
  geom_vline(xintercept=LRP4, color = "purple") + 
  geom_vline(xintercept=LRP5, color = "purple") + 
  
  geom_vline(xintercept=LRP21, color = "blue") + 
  geom_vline(xintercept=LRP22, color = "blue") + 
  
  geom_vline(xintercept=LRP6, color = "green") +
  geom_vline(xintercept=LRP3, color = "green") +
  
  geom_vline(xintercept=0.2*SSB0, color = "red",linetype=2) +
  geom_vline(xintercept=0.4*calc$SSBmsy, color = "darkred",linetype=2) 


mDF <- rbind(c(0,phi0,1,0,R0,0,SSB0),DF)

ggplot() + geom_path(mDF,mapping=aes(y=Yield,x=eq_SSB)) + theme_classic() + labs(x="SSB (kt)", y="SSB (kt") + expand_limits(y=0) + expand_limits(x=0) +
  geom_vline(xintercept=LRP1, color = "purple") + 
  geom_vline(xintercept=LRP4, color = "purple") + 
  geom_vline(xintercept=LRP5, color = "purple") + 
  
  geom_vline(xintercept=LRP21, color = "blue") + 
  geom_vline(xintercept=LRP22, color = "blue") + 
  
  geom_vline(xintercept=LRP6, color = "green") +
  geom_vline(xintercept=LRP3, color = "green") +
  
  geom_vline(xintercept=0.2*SSB0, color = "red",linetype=2) +
  geom_vline(xintercept=0.4*calc$SSBmsy, color = "darkred",linetype=2) 
                     

maxyield = max(DF$Yield)

ggplot() + geom_path(mDF,mapping=aes(y=Yield/maxyield,x=eq_SSB)) + theme_classic() + labs(x="SSB (kt)", y="Relative Yield") + expand_limits(y=0) + expand_limits(x=0) +
  geom_vline(xintercept=LRP1, color = "purple") + 
  geom_vline(xintercept=LRP4, color = "purple") + 
  geom_vline(xintercept=LRP5, color = "purple") + 
  
  geom_vline(xintercept=LRP21, color = "blue") + 
  geom_vline(xintercept=LRP22, color = "blue") + 
  
  geom_vline(xintercept=LRP6, color = "green") +
  geom_vline(xintercept=LRP3, color = "green") +
  
  geom_vline(xintercept=calc$SSBmsy, color = "orange",size=1) +
  
  geom_vline(xintercept=0.2*SSB0, color = "red",linetype=2) +
  geom_vline(xintercept=0.4*calc$SSBmsy, color = "darkred",linetype=2) 
