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

# Plots for Weight-at-age, Maturity-at-age, and vulnerability-at-age
ggplot(AA,aes(y=w_age/1000,x=age)) + geom_path() + theme_classic() + labs(x="Age", y="Weight (kg)") + expand_limits(y=0) + expand_limits(x=0) + scale_x_continuous(breaks=0:11)
ggplot(AA,aes(y=m_age,x=age)) + geom_path() + theme_classic() + labs(x="Age", y="Maturity") + expand_limits(y=0) + expand_limits(x=0) + scale_x_continuous(breaks=0:11)
ggplot(AA,aes(y=v_age,x=age)) + geom_path() + theme_classic() + labs(x="Age", y="Vulnerability") + expand_limits(y=0) + expand_limits(x=0) + scale_x_continuous(breaks=0:11)

########################################################################################################################
# Plots from data frame D
########################################################################################################################
#Stock recruitment pairs (labels are years)
ggplot(D[!is.na(D$Rec),],aes(y=Rec,x=SSB,label=Year)) + geom_point(colour="blue") + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=c(0,12.5)) + expand_limits(x=0) + geom_text(mapping=aes(y=Rec,x=SSB,label=Year),nudge_y = 0.5,size=3) 

#Stock recruitment pairs and model fit
ggplot() + geom_point(D[!is.na(D$Rec),],mapping=aes(y=Rec,x=SSB),colour="blue") + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=c(0,12.5)) + expand_limits(x=0) + geom_function(fun=function(x) BHa*x/(1+BHb*x)) 

#Historical Recruitment
ggplot(D[!is.na(D$Rec),],aes(y=Rec,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="Recruitment (10^9)") + expand_limits(y=0) + expand_limits(x=50) 

#Historical SSB 
ggplot(D,aes(y=SSB,x=Year)) + geom_path() +theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0)

#Historical Catch
ggplot(D,aes(y=Catch,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="Catch (kt)") + expand_limits(y=0)

#Historical F
ggplot(D,aes(y=f,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="F") + expand_limits(y=0) 

#Plot Equilibrium SSB0 (purple), Equilibrium SSBmsy (green)
ggplot(D) + geom_path(mapping=aes(y=SSB,x=Year)) + theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0) +
  geom_hline(yintercept=SSB0, color = "purple") + 
  geom_hline(yintercept=SSBmsy, color = "green") 

#Plot Historical F and Equilibrium Fmsy (green)
ggplot(D,aes(y=f,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="F") + expand_limits(y=0) +
  geom_hline(yintercept=Fmsy, color = "green")  
  

#Plot Acoustic index years 25-50
  #Add x yr moving average (example 3 years)
  D$MA_Index <- NA
  D$MA_Index[D$Year%in%28:50] <- apply(cbind(D$Acoustic_Index[D$Year%in%28:50],
                                             D$Acoustic_Index[D$Year%in%27:49],
                                             D$Acoustic_Index[D$Year%in%26:48]),1,mean)

  #Add loess smoother (example span = 0.5)
  lsmooth <- loess(Acoustic_Index ~ Year,data=D,span=0.5)  
  D$loess_Index <- NA
  D$loess_Index[D$Year%in%26:50] <- predict(lsmooth)

ggplot(D[!is.na(D$Acoustic_Index),]) + geom_path(mapping=aes(y=Acoustic_Index,x=Year),size=1) + theme_classic() + labs(x="Year", y="Acoustic SSB (kt)") + expand_limits(y=0,x=0) +
  geom_path(data=D[!is.na(D$MA_Index),],mapping=aes(y=MA_Index,x=Year),color="red") +
  geom_path(data=D[!is.na(D$loess_Index),],mapping=aes(y=loess_Index,x=Year),color="blue") 


#Plot relative exploitation rate [example: catch/(smoothed index)] years 25-50
ggplot(D[!is.na(D$Acoustic_Index),]) + geom_path(mapping=aes(y=Catch/loess_Index,x=Year)) + theme_classic() + labs(x="Year", y="Relative Exploitation Rate") + expand_limits(y=0,x=0) 


#Plot stock status indicator (example: suppose LRP = 30% SSB0)
ggplot(D) + 
  geom_path(mapping=aes(y=SSB/(0.3*SSB0),x=Year),color="blue") +
  theme_classic() + labs(x="Year", y="Ratio SSB/LRP") + expand_limits(y=0) +
  geom_hline(yintercept=1, linetype="dashed") 


#Plot variability in stock status indicator in year 50 (assume ratio of SSB/SSB0 is normally distributed with CV 20% in each year, example: suppose LRP = 30% SSB0)
myindicator <- D$SSB
mylrp <- 0.3*SSB0

  D$ratio_mean <- myindicator/(mylrp)
  D$ratio_p25 <- myindicator/(mylrp) - 0.675*0.2*myindicator/(mylrp) #25th percentile, Z score = -0.675
  D$ratio_p75 <- myindicator/(mylrp) + 0.675*0.2*myindicator/(mylrp) #75th percentile, Z score = 0.675
  D$ratio_p025 <- myindicator/(mylrp) - 1.96*0.2*myindicator/(mylrp) #2.5th percentile, Z score = -1.96
  D$ratio_p095 <- myindicator/(mylrp) + 1.96*0.2*myindicator/(mylrp) #97.5th percentile, Z score = 1.96
  
ggplot(D) + 
  geom_ribbon(aes(x=Year,ymin=ratio_p025,ymax=ratio_p095),fill='#e6e6e6') +
  geom_ribbon(aes(x=Year,ymin=ratio_p25,ymax=ratio_p75),fill='#ababab') +
  geom_path(mapping=aes(y=ratio_mean,x=Year),color="blue") +
  theme_classic() + labs(x="Year", y="Ratio Indicator/LRP") + expand_limits(y=0) +
  geom_hline(yintercept=1, linetype="dashed") 