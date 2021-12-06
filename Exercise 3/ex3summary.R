########################################################################################################################
# Exercise 3: SSB0 and SSBmsy calculations, plots
########################################################################################################################

rm(list=ls())
library(ggplot2)
source(paste0(getwd(),"/functions.R")) # contains functions "survivorship_F" and "MSYcalc"

# Read in data

Data <- readRDS(paste0(getwd(),"/Exercise 3/ex3_data.rda"))

WAA <- Data[[1]] # weight-at-age in g years 1-50
MAT <- Data[[2]] # maturity-at-age years 1-50
VUL <- Data[[3]] # vulnerability-at-age years 1-50

D <- Data[[4]] # Model estimated parameters by year 

#Calculate Surplus Production P(t) = B(t+1) - B(t) + CATCH(t)
D$P<-NA
D$P[1:49] <- D$B[2:50]-D$B[1:49]+D$Catch[1:49]

########################################################################################################################
# Move h/R0/phi0 parameterization of B-H SRR to a/b parameterization and calculate "unfished equilibrium SSB"
########################################################################################################################

# Model assumes a Beverton-Holt SRR with steepness h = 0.7
h <- 0.75 # model assumed steepness
R0 <- 2.395766 # model estimated equilibrium unfished recruitment [provided]
M <- 0.3 # model assumed natural mortality rate

# mean unfished spawning biomass per recruit (phi0) over fist 5 years (5 years = mean generation time)
phi0_5<-mean(D$phi0[1:5])

# B-H a and b. Here we assume that the stock recruitment relationship is estimated using phi0 from the first 5 years
BHa <- 4*h/(phi0_5*(1-h)) # estimated Beverton-Holt a
BHb <- 1/R0*(BHa-1/phi0_5) # estimated Beverton-Holt b

# Some SSB0 options
 # equilibrium unfished SSB0 (using mean phi0 over fist 5 years)
 SSB0 <- R0*phi0_5
 SSB0_l10 <- R0*mean(D$phi0[41:50])
 
 # equilibrium unfished SSB0 (using mean phi0 over entire time series)
 SSB0x <- R0*mean(D$phi0)

########################################################################################################################
# Calculate SSBmsy using WAA, MAT, VUL from a specific time period using the "MSYcalc" function
# Note: examples below are for the first 5 years and the last 10 years (years 41-50) of the historical time series 
########################################################################################################################

# MSYcalc: function from source.R returns a list with Fmsy, msy, SSBmsy from M=M, waa=weight-at-age, mat=maturity-at-age, sel=vulnerability-at-age, Beverton-Holt a and b
# MSYcalc <- function(M,waa,mat,sel,a,b)

calc_years1_5 <- MSYcalc(M=0.3,waa=apply(WAA[1:5,],2,mean), mat=apply(MAT[1:5,],2,mean), sel=apply(VUL[1:5,],2,mean),a=BHa, b=BHb)
calc_last_10_years <- MSYcalc(M=0.3,waa=apply(WAA[41:50,],2,mean), mat=apply(MAT[41:50,],2,mean), sel=apply(VUL[41:50,],2,mean),a=BHa, b=BHb)
calc_last_20_years <- MSYcalc(M=0.3,waa=apply(WAA[31:50,],2,mean), mat=apply(MAT[31:50,],2,mean), sel=apply(VUL[31:50,],2,mean),a=BHa, b=BHb)
calc_first_20_years <- MSYcalc(M=0.3,waa=apply(WAA[1:20,],2,mean), mat=apply(MAT[1:20,],2,mean), sel=apply(VUL[1:20,],2,mean),a=BHa, b=BHb)

calc_years1_5$Fmsy
calc_years1_5$SSBmsy
calc_years1_5$msy

# Ratio of SSBmsy to SSB0 using average biological parameters over the first 5 years

SSBmsy_SSB0 <- calc_years1_5$SSBmsy / SSB0

########################################################################################################################
# Plots from data frame D
########################################################################################################################

LRP1 = 0.4*calc_last_20_years$SSBmsy

D$SSB[20]/LRP1
D$SSB[50]/LRP1


D$LRP1 = c(rep(NA,30),rep(0.4*calc_last_20_years$SSBmsy,20))
LRP2 = D$SSB[22] #Brec
D$LRP3 = c(rep(NA,40),rep(0.2*SSB0_l10,10))
D$LRP32 = c(rep(NA,40),rep(0.2*SSB0,10))
D$LRP4 = c(rep(NA,25),rep(121,25)) #0.4SSNfspr40%
D$LRP5 = c(rep(0.4*calc_first_20_years$SSBmsy,20),rep(NA,20),rep(0.4*calc_last_10_years$SSBmsy,10))
LRP6 = 0.2*SSB0 
LRP62 = D$SSB[20]  #Brec

#Stock recruitment pairs  (labels are years)
ggplot(D[!is.na(D$Rec),],aes(y=Rec,x=SSB,label=Year)) + geom_point(colour="blue") + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=c(0,12.5)) + expand_limits(x=0) + geom_text(mapping=aes(y=Rec,x=SSB,label=Year),nudge_y = 0.5,size=3) 

#Stock recruitment pairs and model fit
ggplot() + geom_point(D[!is.na(D$Rec),],mapping=aes(y=Rec,x=SSB),colour="blue") + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=c(0,12.5)) + expand_limits(x=0) + geom_function(fun=function(x) BHa*x/(1+BHb*x)) 

#Historical Recruitment
ggplot(D[!is.na(D$Rec),],aes(y=Rec,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="Recruitment (10^9)") + expand_limits(y=0,x=50) 

#Historical SSB 
ggplot(D,aes(y=SSB,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0)

#Historical Catch
ggplot(D,aes(y=Catch,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="Catch (kt)") + expand_limits(y=0)

#Historical F
ggplot(D,aes(y=f,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="F") + expand_limits(y=0) 

#Surplus Production
ggplot(D[!is.na(D$P),],aes(y=P,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="Surplus Production (kt)") + expand_limits(y=0,x=50) +   geom_hline(yintercept=0)   

ggplot(D[!is.na(D$P),],aes(y=P,x=SSB)) + geom_point() + theme_classic() + labs(x="SSB (kt)", y="Surplus Production (kt)") + expand_limits(y=0,x=0) +   geom_hline(yintercept=0)   

ggplot(D[!is.na(D$P),],aes(y=P,x=SSB)) + geom_path() + theme_classic() + labs(x="SSB (kt)", y="Surplus Production (kt)") + expand_limits(y=0,x=0) +   geom_hline(yintercept=0)  +  geom_text(mapping=aes(y=P,x=SSB,label=Year),nudge_y = 0.5,size=4,colour="blue") 

#Plot Historical SSB with some metrics that could be used for LRPs
ggplot(D) + 
  geom_path(mapping=aes(y=SSB,x=Year)) +
  theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0) +
  geom_path(mapping=aes(y=LRP1,x=Year), linetype="dashed", color = "red",size=1) +  
  geom_hline(yintercept=LRP2, linetype="dashed", color = "purple",size=1) +  
  geom_path(mapping=aes(y=LRP3,x=Year), linetype="dashed", color = "green",size=1) +  
  #geom_hline(yintercept=LRP3, color = "red",size=1) +  
  geom_path(mapping=aes(y=LRP4,x=Year),  color = "blue",size=1) +  
  #geom_hline(yintercept=LRP4, color = "blue",size=1) +  
  geom_path(mapping=aes(y=LRP5+5,x=Year),  color = "red",linetype=1, size=1) +  
  geom_hline(yintercept=LRP6, color = "green",size=1) +  
  geom_hline(yintercept=LRP62, linetype="dashed", color = "purple",size=1)   
  















#Plot stock status indicator based on SSB (examples below)
 mylrp <- (0.2*D$dSSB0a) # example 20% dynamic SSB0
 mylrp <- (0.5*calc_last_10_years$SSBmsy) # 50% equilibrium SSBmsy using data from last 10 years

ggplot(D) + geom_path(mapping=aes(y=SSB/mylrp,x=Year),color="blue") + theme_classic() + labs(x="Year", y="Ratio SSB/LRP") + expand_limits(y=0) +
  geom_hline(yintercept=1, linetype="dashed")

########################################################################################################################
# F indicators
########################################################################################################################

#Equilibrium Fmsy using annual data for WAA, MAT, VUL

D$F_msy_annual <- NA

for(i in 1:nrow(D)){
  D$F_msy_annual[i] <- MSYcalc(M=0.3,waa=apply(WAA[i,],2,mean), mat=apply(MAT[i,],2,mean), sel=apply(VUL[i,],2,mean),a=BHa, b=BHb)$Fmsy
}

ggplot(D,aes(y=f,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="F") + expand_limits(y=0) +
  geom_hline(yintercept=calc_last_10_years$Fmsy, color = "blue") +   
  geom_path(mapping=aes(y=F_msy_annual,x=Year),color="purple")        

########################################################################################################################
# Empirical Indicators Below (if needed)
########################################################################################################################

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
