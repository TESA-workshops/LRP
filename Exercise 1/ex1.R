########################################################################################################################
# Exercise 1: Plot of Catch by MU and plots of indices with example moving averages and loess smoothers
########################################################################################################################

rm(list=ls())
library(ggplot2)

#Read in data inputs
C <- read.csv(paste0(getwd(),"/Exercise 1/ex1_catch.csv"))
D <- read.csv(paste0(getwd(),"/Exercise 1/ex1_indices.csv"))

head(C)
tail(D)

#Calculate CPUE for purse seine fleet in MU1
D$PS_CPUE <- (D$PS_Catch_MU1)/D$PS_Effort

########################################################################################################################
# Plots (entire stock area)
########################################################################################################################

# 1. Catch by MU
ggplot() + geom_area(data=C, mapping=aes(y=Catch_kt,x=Year,fill=MU)) + theme_classic() + labs(x='Year', y="Total Catch (kt)") + scale_x_continuous(limits = c(0,50), expand = c(0, 0)) + scale_y_continuous(limits = c(0,300), expand = c(0, 0)) + scale_fill_manual(values=c("lightskyblue2","plum3","lightgoldenrod"))
 

# 2. MU1 and MU2 bottom trawl survey with an example moving average and loess smoother if required

  #Add moving average for the index (example 3 years)
  D$MA_BT1_2  <- NA
  D$MA_BT1_2[11:50] <- apply(cbind(D$BT_Index_MU1_2[D$Year%in%11:50],
                             D$BT_Index_MU1_2[D$Year%in%10:49],
                             D$BT_Index_MU1_2[D$Year%in%9:48]),1,mean)

  #Add loess smoother for the index (example span = 0.5)
  lsmooth1 <- loess(BT_Index_MU1_2 ~ Year,data=D,span=0.5)  
  D$loess_BT1_2 <- NA
  D$loess_BT1_2[D$Year%in%9:50] <- predict(lsmooth1)

  #Edit plot as necessary, including horizontal line for LRP
ggplot(D[!is.na(D$BT_Index_MU1_2),]) + geom_path(mapping=aes(y=BT_Index_MU1_2,x=Year)) + theme_classic() + labs(x="Year", y="MU1 and MU2 Bottom Trawl Survey Index (kt)") + expand_limits(y=0,x=50) +
  scale_x_continuous(limits = c(0,50), expand = c(0, 0)) +
  geom_path(data=D[!is.na(D$MA_BT1_2),],mapping=aes(y=MA_BT1_2,x=Year),color="red") +
  geom_path(data=D[!is.na(D$loess_BT1_2),],mapping=aes(y=loess_BT1_2,x=Year),color="blue") +
  geom_hline(yintercept=1000, linetype="dashed", color = "orange",size=1)

# 3. MU1 and MU3 acoustic survey with an example moving average and loess smoother if required
 
  #Add moving average for the index (example 3 years)
  D$MA_As1_3 <- NA
  D$MA_As1_3[28:50] <- apply(cbind(D$Ac_Index_MU1_3[D$Year%in%28:50],
                             D$Ac_Index_MU1_3[D$Year%in%27:49],
                             D$Ac_Index_MU1_3[D$Year%in%26:48]),1,mean)

  #Add loess smoother for the index (example span = 0.5)
  lsmooth2 <- loess(Ac_Index_MU1_3 ~ Year,data=D,span=0.5)  
  D$loess_As1_3 <- NA
  D$loess_As1_3[D$Year%in%26:50] <- predict(lsmooth2)

  #Edit plot as necessary, including horizontal line for LRP
ggplot(D[!is.na(D$Ac_Index_MU1_3),]) + geom_path(mapping=aes(y=Ac_Index_MU1_3,x=Year)) + theme_classic() + labs(x="Year", y="MU1 and MU3 Acoustic SSB Index (kt)") + expand_limits(y=0,x=50) +
  scale_x_continuous(limits = c(0,50), expand = c(0, 0)) +
  geom_path(data=D[!is.na(D$MA_As1_3),],mapping=aes(y=MA_As1_3,x=Year),color="red") +
  geom_path(data=D[!is.na(D$loess_As1_3),],mapping=aes(y=loess_As1_3,x=Year),color="blue") +
  geom_hline(yintercept=750, linetype="dashed", color = "orange",size=1)

########################################################################################################################
# Plots (MU1)
########################################################################################################################

# 1. MU1 PS Catch
ggplot(D[!(is.na(D$PS_Catch_MU1)),],aes(y=PS_Catch_MU1,x=Year)) + geom_path() + theme_classic() + labs(x='Year', y="MU1 Purse Seine Catch (kt)") + expand_limits(y=0)

# 2. MU1 PS Effort
ggplot(D[!(is.na(D$PS_Effort_MU1)),],aes(y=PS_Effort_MU1,x=Year)) + geom_path() + theme_classic() + labs(x='Year', y="MU1 Purse Seine Effort (trips)") + expand_limits(y=0)

# 3. MU1 PS CPUE
ggplot(D[!(is.na(D$PS_CPUE)),],aes(y=PS_CPUE,x=Year)) + geom_path() + theme_classic() + labs(x='Year', y="MU1 Purse Seine CPUE (kt/trip)") + expand_limits(y=0)


# 4. MU1 bottom trawl survey with an example moving average and loess smoother if required

  #Add moving average for the index (example 3 years)
  D$MA_BT1 <- NA
  D$MA_BT1[11:50] <- apply(cbind(D$BT_Index_MU1[D$Year%in%11:50],
                             D$BT_Index_MU1[D$Year%in%10:49],
                             D$BT_Index_MU1[D$Year%in%9:48]),1,mean)

  #Add loess smoother (example span = 0.5)
  lsmooth3 <- loess(BT_Index_MU1 ~ Year,data=D,span=0.5)  
  D$loess_As1 <- NA
  D$loess_As[D$Year%in%9:50] <- predict(lsmooth3)

  #Edit plot as necessary, including horizontal line for LRP
ggplot(D[!is.na(D$BT_Index_MU1),]) + geom_path(mapping=aes(y=BT_Index_MU1,x=Year)) + theme_classic() + labs(x="Year", y="MU1 Bottom Trawl Survey Index (kt)") + expand_limits(y=0,x=50) +
  scale_x_continuous(limits = c(0,50), expand = c(0, 0)) +
  geom_path(data=D[!is.na(D$MA_BT1),],mapping=aes(y=MA_BT1,x=Year),color="red") +
  geom_path(data=D[!is.na(D$loess_As),],mapping=aes(y=loess_As,x=Year),color="blue") +
  geom_hline(yintercept=250, linetype="dashed", color = "orange",size=1)


# 5. MU1 Acoustic Survey with an example moving average and loess smoother if required

  #Add moving average for the index (example 3 years)
  D$MA_As1 <- NA
  D$MA_As1[28:50] <- apply(cbind(D$Ac_Index_MU1[D$Year%in%28:50],
                             D$Ac_Index_MU1[D$Year%in%27:49],
                             D$Ac_Index_MU1[D$Year%in%26:48]),1,mean)

  #Add loess smoother (example span = 0.5)
  lsmooth4 <- loess(Ac_Index_MU1 ~ Year,data=D,span=0.5)  
  D$loess_As1 <- NA
  D$loess_As1[D$Year%in%26:50] <- predict(lsmooth4)
  
  #Edit plot as necessary, including horizontal line for LRP
ggplot(D[!is.na(D$Ac_Index_MU1),]) + geom_path(mapping=aes(y=Ac_Index_MU1,x=Year)) + theme_classic() + labs(x="Year", y="MU1 Acoustic Index of SSB (kt)") + expand_limits(y=0,x=50) +
  scale_x_continuous(limits = c(0,50), expand = c(0, 0)) +
  geom_path(data=D[!is.na(D$MA_As1),],mapping=aes(y=MA_As1,x=Year),color="red") +
  geom_path(data=D[!is.na(D$loess_As1),],mapping=aes(y=loess_As1,x=Year),color="blue") +
  geom_hline(yintercept=750, linetype="dashed", color = "orange",size=1)


