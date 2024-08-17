library(nlme)
library(MASS)
# setwd("\\\\nask.man.ac.uk\\home$\\Desktop\\Longitudial Data\\Course Work")
setwd("C:\\Users\\mz17\\Desktop\\Stat App II Final Project")
cd4=read.table("cd4_new.txt",header=T)
cd4data<-data.frame(cd4)
# adding two items to the data.frame
cd4data$time2<-(cd4data$time)^2
cd4data$time3<-(cd4data$time)^3
# group the data
cd4data<-groupedData(cd4~time|subject,data=cd4data,FUN=mean,outer=~age)
head(cd4data)
####################################################################################
# Model Selection
#-----------------------------------------------------------------------
# Response Transformation
cd4lm<-lm(cd4~time+time2+time3+cig+drug+cesd+time:cig,data=cd4data)
summary(cd4lm)
boxcox(cd4lm)
# Y=sqrt(cd4) sqrt transformation
#-----------------------------------------------------------------------
# scatterplot matrix
pairs(~time+cd4+age+cig+drug+partner+cesd,data=cd4)
#-----------------------------------------------------------------------
# Comparison for mean functions
cd4lm1<-lm(sqrt(cd4)~time+cig+drug+cesd,data=cd4data)
summary(cd4lm1)
cd4lm2<-lm(sqrt(cd4)~time+time2+time3+cig+drug+cesd+time:cig,data=cd4data)
summary(cd4lm2)
anova(cd4lm1,cd4lm2)
plot(cd4lm2,which=1)
acf(resid(cd4lm2)) # expected to see autocorrelation for the same patient
#-----------------------------------------------------------------------
# Linear Mixed model with AR1 error covariance
cd4lmm.1.I<-lme(sqrt(cd4)~time+time2+time3+cig+drug+cesd+time:cig,data=cd4data,
                random=~1|subject,correlation=corAR1(form=~time|subject),method="ML")
summary(cd4lmm.1.I)
cd4lmm.un.ar1<-lme(sqrt(cd4)~time+time2+time3+cig+drug+cesd+time:cig,data=cd4data,
                   random=~time|subject,correlation=corAR1(form=~time|subject),method="ML")
summary(cd4lmm.un.ar1)
anova(cd4lmm.1.I,cd4lmm.un.ar1)
AIC(cd4lm1,cd4lm2,cd4lmm.1.I,cd4lmm.un.ar1)
# comparison for the covariance matrix for random effects
cd4lmm.un.ar1<-lme(sqrt(cd4)~time+time2+time3+cig+drug+cesd+time:cig,data=cd4data,
                   random=~time|subject,correlation=corAR1(form=~time|subject),method="ML")
summary(cd4lmm.un.ar1)
cd4lmm.diag.ar1<-update(cd4lmm.un.ar1,random=pdDiag(~time))
summary(cd4lmm.diag.ar1)
cd4lmm.inde.ar1<-update(cd4lmm.un.ar1,random=pdIdent(~time))
summary(cd4lmm.inde.ar1)
anova(cd4lmm.un.ar1,cd4lmm.diag.ar1,cd4lmm.inde.ar1)
# From the LR-test and AIC&BIC cd4lmm.diag.ar1 is the best among them.
cd4lmm.diag.ar4<-update(cd4lmm.un.ar1,correlation=corARMA(form=~time|subject,p=4,q=0))
summary(cd4lmm.diag.ar4)
acf(resid(cd4lmm.diag.ar4))
pacf(resid(cd4lmm.diag.ar4))
anova(cd4lmm.diag.ar1,cd4lmm.diag.ar4)
AIC(cd4lm1,cd4lm2,cd4lmm.1.I,cd4lmm.un.ar1,cd4lmm.diag.ar1,cd4lmm.inde.ar1)
plot(cd4lmm.inde.ar1)
summary(cd4lmm.diag.ar1)
# Based on the model selected in a
# cd4lmm.diag.ar1
# (1)
summary(cd4lmm.diag.ar1)
# Value Std.Error DF t-value p-value
# time -2.196960 0.1561324 2000 -14.07113 0.0000
# cig 0.641313 0.1291044 2000 4.96740 0.0000
# anova.lme(cd4lmm.diag.ar1)
# CD4+ cell decrease over the time.
# CD4+ dell increase with the number of packs of cigarettes smoked per day
# cd4data$cd4
plot(cd4data$time,cd4data$cd4)
f<-function(x){(28.527470-2.196960*x-0.440518*x^2+0.118350*x^3)^2}
x<-seq(-4,6,by=.1)
y<-f(x)
lines(x,y,lwd=2,col=2)
abline(v=0,lwd=2,col=4)
reduce0<-lme(sqrt(cd4)~time2+time3+cig+drug+cesd+time:cig,data=cd4data,
             random=pdDiag(~time),correlation=corAR1(form=~time|subject),method="ML")
anova(reduce0,cd4lmm.diag.ar1)
reduce1<-lme(sqrt(cd4)~time+time2+time3+cig+drug+cesd,data=cd4data,
             random=pdDiag(~time),correlation=corAR1(form=~time|subject),method="ML")
anova(reduce1,cd4lmm.diag.ar1)
# (2)
# Value Std.Error DF t-value p-value
# drug 0.473593 0.3132456 2000 1.51189 0.1307
# cesd -0.043725 0.0136695 2000 -3.19876 0.0014
# At the 10% significant level,
# the drug using dosen't significantly affect the change of the CD4+ cells
# At the 5% significant level,
# (score of center for epidemiological studies of depression)
# the cesd significantly affects the amount of CD4+ cells.
# However this effect is quite small and negative.
#####################################################################################
cd4lmm.un.ar1<-lme(sqrt(cd4)~time+time2+time3+cig+drug+cesd+time:cig,data=cd4data,
                   random=~time|subject,correlation=corAR1(form=~time|subject),method="ML")
summary(cd4lmm.un.ar1)
# AIC BIC logLik
# 14228.91 14303.96 -7101.456
# Phi1
# 0.4679066
intervals(cd4lmm.un.ar1)
# lower est. upper
# Phi1 0.3994602 0.4679066 0.5311654
# same as
cd4lmm.un.ar1<-lme(sqrt(cd4)~time+time2+time3+cig+drug+cesd+time:cig,data=cd4data,
                   random=~time|subject,correlation=corAR1(form=~time),method="ML")
summary(cd4lmm.un.ar1)
# Phi1
# 0.4679066
str(cd4lmm.un.ar1)
arima(resid(cd4lmm.un.ar1)[2371:2376],order=c(1,0,0))