########################################################################################################################
# Exercise 4: LRPs for Model 1 with low equilibrium catch and Model 2 with high equilibrium catch
########################################################################################################################

rm(list=ls())
library(ggplot2)
source(paste0(getwd(),"/functions.R")) # contains functions "survivorship_F" and "MSYcalc"

# Read in data

Data <- readRDS(paste0(getwd(),"/Exercise 4/ex4_data.rda"))

WAA <- Data$WAA # weight-at-age for years 1:50
MAT <- Data$MAT # maturity-at-age for years 1:50
VUL1 <- as.data.frame(Data$VUL1) # vulnerability-at-age for years 1:50 for model 1
VUL2<- as.data.frame(Data$VUL2) # vulnerability-at-age for years 1:50 for model 2
D1<- Data$D1 # Model estimated parameters by year for model 1
D2 <- Data$D2 # Model estimated parameters by year for model 2

########################################################################################################################
# Move h/R0/phi0 parameterization of B-H SRR to a/b parameterization and calculate "unfished equilibrium SSB"
########################################################################################################################
# Model assumes a Beverton-Holt SRR with steepness h = 0.75 and M = 0.3
h <- 0.75 # model assumed steepness
R0_1 <- 2.395766 # model estimated equilibrium unfished recruitment [provided]
M_1 <- 0.3 # model assumed natural mortality rate

# mean unfished spawning biomass per recruit (phi0) over fist 5 years (5 years = mean generation time)
phi0_5y_1<-mean(D1$phi0[1:5])

# B-H a and b. Here we assume that the stock recruitment relationship is estimated using phi0 from the first 5 years
BHa_1 <- 4*h/(phi0_5y_1*(1-h)) # estimated Beverton-Holt a
BHb_1 <- 1/R0_1*(BHa_1-1/phi0_5y_1) # estimated Beverton-Holt b

# equilibrium unfished SSB0 (using phi0 over entire time series)
SSB0_1 <- R0_1*mean(D1$phi0)

# annual "equilibrum" SSB0
D1$annualSSB0 <- R0_1*D1$phi0

########################################################################################################################
# Move h/R0/phi0 parameterization of B-H SRR to a/b parameterization and calculate "unfished equilibrium SSB"
########################################################################################################################
# Model assumes a Beverton-Holt SRR with steepness h = 0.75 and M = 0.3
h <- 0.75 # model assumed steepness
R0_2 <- 3.15979# model estimated equilibrium unfished recruitment [provided]
M_2 <- 0.3 # model assumed natural mortality rate

# mean unfished spawning biomass per recruit (phi0) over fist 5 years (5 years = mean generation time)
phi0_5y_2<-mean(D2$phi0[1:5])

# B-H a and b. Here we assume that the stock recruitment relationship is estimated using phi0 from the first 5 years
BHa_2 <- 4*h/(phi0_5y_2*(1-h)) # estimated Beverton-Holt a
BHb_2 <- 1/R0_2*(BHa_2-1/phi0_5y_2) # estimated Beverton-Holt b

# equilibrium unfished SSB0 (using phi0 over entire time series)
SSB0_2 <- R0_2*mean(D2$phi0)

# annual "equilibrum" SSB0
D2$annualSSB0 <- R0_2*D2$phi0

########################################################################################################################
# Calculate SSBmsy using WAA, MAT, VUL from first 5 years and last 10 years [Note: 50 years of data]
# Note: MSYcalc function can be used to estimate MSY reference points for various time periods. 
########################################################################################################################

# MSYcalc: function from source.R returns a list with Fmsy, msy, SSBmsy from M=M, waa=weight-at-age, mat=maturity-at-age, sel=vulnerability-at-age, Beverton-Holt a and b
# MSYcalc <- function(M,waa,mat,sel,a,b)

Model1_calc_years1_5 <- MSYcalc(M=M_1,waa=apply(WAA[1:5,],2,mean), mat=apply(MAT[1:5,],2,mean), sel=apply(VUL1[1:5,],2,mean),a=BHa_1, b=BHb_1)
Model2_calc_years1_5 <- MSYcalc(M=M_2,waa=apply(WAA[1:5,],2,mean), mat=apply(MAT[1:5,],2,mean), sel=apply(VUL2[1:5,],2,mean),a=BHa_2, b=BHb_2)

Model1_calc_last_10_years <- MSYcalc(M=M_1,waa=apply(WAA[41:50,],2,mean), mat=apply(MAT[41:50,],2,mean), sel=apply(VUL1[41:50,],2,mean),a=BHa_1, b=BHb_1)
Model2_calc_last_10_years <- MSYcalc(M=M_2,waa=apply(WAA[41:50,],2,mean), mat=apply(MAT[41:50,],2,mean), sel=apply(VUL2[41:50,],2,mean),a=BHa_2, b=BHb_2)

########################################################################################################################
# Plots from data frame D1 and D2. Blue plots are model1 and Purple plots are model2
########################################################################################################################

#Stock recruitment pairs  (labels are years)
ggplot(D1[!is.na(D1$Rec),],aes(y=Rec,x=SSB,label=Year)) + geom_point(colour="blue") + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=0) + expand_limits(x=0) + geom_text(mapping=aes(y=Rec,x=SSB,label=Year),nudge_y = 0.5,size=3) 
ggplot(D2[!is.na(D2$Rec),],aes(y=Rec,x=SSB,label=Year)) + geom_point(colour="purple") + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=0) + expand_limits(x=0) + geom_text(mapping=aes(y=Rec,x=SSB,label=Year),nudge_y = 0.5,size=3) 

#Stock recruitment pairs and model fit
ggplot() + geom_point(D1[!is.na(D1$Rec),],mapping=aes(y=Rec,x=SSB),colour="blue") + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=c(0,12.5)) + expand_limits(x=0) + geom_function(fun=function(x) BHa_1*x/(1+BHb_1*x)) 
ggplot() + geom_point(D2[!is.na(D2$Rec),],mapping=aes(y=Rec,x=SSB),colour="purple") + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=c(0,12.5)) + expand_limits(x=0) + geom_function(fun=function(x) BHa_2*x/(1+BHb_2*x)) 

#Historical Recruitment
ggplot() + geom_path(D1[!is.na(D1$Rec),],mapping=aes(y=Rec,x=Year),colour="blue") + theme_classic() + labs(x="Year", y="Recruitment (10^9)") + expand_limits(y=0,x=50) 
ggplot() + geom_path(D2[!is.na(D2$Rec),],mapping=aes(y=Rec,x=Year),colour="purple") + theme_classic() + labs(x="Year", y="Recruitment (10^9)") + expand_limits(y=0,x=50) 

#Historical SSB 
ggplot() + geom_path(D1,mapping=aes(y=SSB,x=Year),colour="blue") + theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0)
ggplot() + geom_path(D2,mapping=aes(y=SSB,x=Year),colour="purple") + theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0)

#Historical Catch
  M1_hist_wq_catch = data.frame(Year = -20:0, Catch = 4)
  M2_hist_wq_catch = data.frame(Year = -20:0, Catch = 80)

ggplot() + geom_path(D1,mapping=aes(y=Catch,x=Year),colour="blue") + theme_classic() + labs(x="Year", y="Catch (kt)") + expand_limits(y=0) +
  geom_path(M1_hist_wq_catch,mapping=aes(y=Catch,x=Year),colour="blue",linetype=2)
ggplot() + geom_path(D2,mapping=aes(y=Catch,x=Year),colour="purple") + theme_classic() + labs(x="Year", y="Catch (kt)") + expand_limits(y=0) +
  geom_path(M2_hist_wq_catch,mapping=aes(y=Catch,x=Year),colour="purple",linetype=2)

#Historical F
ggplot() + geom_path(D1,mapping=aes(y=f,x=Year),colour="blue") + theme_classic() + labs(x="Year", y="F") + expand_limits(y=0) +
  geom_hline(yintercept=Model1_calc_last_10_years$Fmsy, linetype="dashed", color = "green")  

ggplot() + geom_path(D2,mapping=aes(y=f,x=Year),colour="purple") + theme_classic() + labs(x="Year", y="F") + expand_limits(y=0) +
  geom_hline(yintercept=Model2_calc_last_10_years$Fmsy, linetype="dashed", color = "green")  

#Plot Historical SSB with some metrics that could be used for LRPs
#See code from Exercise 3 to generate additional plots.
ggplot(D1) + 
  geom_path(mapping=aes(y=SSB,x=Year),colour="blue") +
  theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0) +
  geom_hline(yintercept=0.2*SSB0_1, linetype="dashed", color = "green")  

ggplot(D2) + 
  geom_path(mapping=aes(y=SSB,x=Year),colour="purple") +
  theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0) +
  geom_hline(yintercept=0.2*SSB0_2, linetype="dashed", color = "green")  

########################################################################################################################
# Plot of Stock Status over time using SSB
########################################################################################################################

# Example 1: mean over both models
my_mean_limit <- (0.2*SSB0_1 + 0.2*SSB0_2)/2 # example 20% equilibrium unfished SSB0 (using phi0 over fist 5 years)
MEAN_SSB <- data.frame(mean_SSB=(D1$SSB+D2$SSB)/2,Year=1:50)

ggplot() + geom_path(MEAN_SSB, mapping=aes(y=mean_SSB/my_mean_limit,x=Year),color="blue") + theme_classic() + labs(x="Year", y="Ratio SSB/LRP") + expand_limits(y=0) +
  geom_hline(yintercept=1, linetype="dashed")

# Example 2: one model 
mylrp <- (0.2*SSB0_1) # example 20%  SSB0

ggplot(D1) + geom_path(mapping=aes(y=SSB/mylrp,x=Year),color="blue") + theme_classic() + labs(x="Year", y="Ratio SSB/LRP") + expand_limits(y=0) +
  geom_hline(yintercept=1, linetype="dashed")


mylrp2 <- (0.2*SSB0_2) # example 20%  SSB0

ggplot(D2) + geom_path(mapping=aes(y=SSB/mylrp2,x=Year),color="purple") + theme_classic() + labs(x="Year", y="Ratio SSB/LRP") + expand_limits(y=0) +
  geom_hline(yintercept=1, linetype="dashed")