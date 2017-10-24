########### SARIMA MODEL OF MONTHLY ADR ##############

# Import Data
myd=read.table("Monthly.csv",header=T, sep=',')

#Load Libraries
library(tseries)
library(fBasics)
library(forecast)
library(lmtest) 
library(fUnitRoots)

# Data Exploration
ADR=myd$ADR
hist(myd$ADR, xlab="Monthly ADR", prob=TRUE,main="Histogram")
xfit<-seq(min(myd$ADR),max(myd$ADR),length=40)
yfit<-dnorm(xfit,mean=mean(myd$ADR),sd=sd(myd$ADR)) 
lines(xfit, yfit, col="blue", lwd=2) 
qqnorm(myd$ADR)
qqline(myd$ADR, col = 2) 
normalTest(ADR,method=c("jb")) 

#Create TS
ADRts =ts(myd[,2], start=c(2005,1), freq=12)
plot(ADRts,main="Fig 1_1 Monthly ADR TS Plot",ylab="ADR")

par(mfcol=c(1,2))
hist(ADRts, xlab="ADR", prob=TRUE, main="Fig 2_1 Histogram of ADR")
qqnorm(ADRts)
qqline(ADRts, col = 2)
normalTest(ADRts,method=c("jb"))

#create ACF
acf(as.vector(ADRts),lag.max=30, main="Fig3_1(a) ACF of ADR TS")
pacf(as.vector(ADRts),lag.max=30, main="PACF of time series data")

#Apply first difference
dx=diff(ADRts) 
acf(as.vector(dx),lag.max=50, main="Fig3_1(b) ACF of 1st diff ADR")
pacf(as.vector(dx),lag.max=50, main="PACF of time series data")

#Apply seasonal difference
sdx=diff(dx,12)
acf(as.vector(sdx),lag.max=50, main="Fig3_1(c) ACF of DSDX")
Pacf(as.vector(sdx),lag.max=50, main="PACF of DSDX starts")

#Select Model
m0=auto.arima(ADRts, ic="bic", seasonal=T, stationary=F, trace=T)
m1=Arima(ADRts,order=c(1,0,1),seasonal = list(order=c(0,1,1),period=12),method="ML")
m1
coeftest(m1)
normalTest(m1$residuals,method = c("jb"))
par(mfcol=c(1,1))
qqnorm(m1$residuals)
qqline(m1$residuals,col=2)
Acf(m1$residuals,lag.max=50, main="ACF_ADR M1 Residuals")
Box.test(m1$residuals, 4, "Ljung-Box", fitdf=length(coef(m1)))
Box.test(m1$residuals, 6, "Ljung-Box", fitdf=length(coef(m1)))

m2=Arima(ADRts,order=c(0,1,1),seasonal=list(order=c(0,1,1),period=12), method="ML")
m2
coeftest(m2)
normalTest(m2$residuals,method = c("jb"))
qqnorm(m2$residuals)
qqline(m2$residuals,col=2)
Acf(m2$residuals)
Box.test(m2$residuals, 4, "Ljung-Box", fitdf=length(coef(m2)))
Box.test(m2$residuals, 6, "Ljung-Box", fitdf=length(coef(m2)))

m3=Arima(ADRts,order=c(0,1,1),seasonal=list(order=c(1,0,1),period=12), method="ML")
m3
coeftest(m3)
normalTest(m3$residuals,method = c("jb"))
qqnorm(m3$residuals)
qqline(m3$residuals,col=2)
Acf(m3$residuals)
Box.test(m3$residuals, 4, "Ljung-Box", fitdf=length(coef(m3)))
Box.test(m3$residuals, 6, "Ljung-Box", fitdf=length(coef(m3)))

#6 Forecast
f1=forecast.Arima(m1,h=12)
f1
plot(f1,include=50, main="ADR 12 monthes Forecasts")
plot(m1$x)
lines(fitted(m1),col=2)

f2=forecast.Arima(m2,h=12)
f2
plot(f2,include=50)
plot(m2$x)
lines(fitted(m2),col=2)

f3=forecast.Arima(m3,h=12)
f3
plot(f3,include=50)
plot(m3$x)
lines(fitted(m3),col=2)

#7 Backtesting
source("backtest.R")
pm1=backtest(m1,ADRts,120,1)
pm2=backtest(m2,ADRts,120,1)
pm3=backtest(m3,ADRts,120,1)

#### SARIMA MONTHLY Occupied Rooms FORECAST#####
#1 Normal?
OR=myd[,5]
normalTest(OR,method = c("jb"))
qqnorm(OR)
qqline(OR,col=2)

#2 Plot TS
ORts=ts(OR,start=c(2005,1),frequency = 12)
plot(ORts,main = "Fig 1_2 Monthly Occupired Rms TS Plot")

par(mfcol=c(1,2))
hist(ORts, xlab="Occulied Rms", prob=TRUE, main="Fig 2_2 Histogram Occ")
qqnorm(ORts)
qqline(ORts, col = 2)
normalTest(ORts,method=c("jb"))

#3 Stationary? Seasonality?
acf(as.vector(ORts),lag.max=30, main="Fig3_2(a) ACF of OR TS")
pacf(as.vector(ORts),lag.max=30, main="PACF of time series data")

#Apply first difference
dx=diff(ORts) 
acf(as.vector(dx),lag.max=50, main="Fig3_2(b) ACF of 1st diff OR")
pacf(as.vector(dx),lag.max=50, main="PACF of time series data")

#Apply seasonal difference
sdx=diff(dx,12)
acf(as.vector(sdx),lag.max=50, main="Fig3_2(c) ACF of DSDX")
Pacf(as.vector(sdx),lag.max=50, main="PACF of DSDX starts")

#4 Identify p q d D
m5=auto.arima(ORts,stationary = F,seasonal = T,trace = T)

#5 Fit Model and Residuals Analysis
m5=Arima(ORts, order=c(1,1,1), seasonal=list(order=c(0,1,2), period=12), method="ML")
coeftest(m5)
plot(m5$residuals)
qqnorm(m5$residuals)
qqline(m5$residuals,col=2)
normalTest(m5$residuals,method=c("jb"))
Box.test(m5$residuals,lag=6,type=c("Ljung"),fitdf=length(coef(m5)))
Box.test(m5$residuals,lag=8,type=c("Ljung"),fitdf=length(coef(m5)))
Acf(m5$residuals)
Pacf(m5$residuals)
m6=Arima(ORts, order=c(0,1,1), seasonal=list(order=c(0,1,2), period=12), method="ML")
coeftest(m6)
plot(m6$residuals)
qqnorm(m6$residuals)
qqline(m6$residuals,col=2)
normalTest(m6$residuals,method=c("jb"))
Box.test(m6$residuals,lag=6,type=c("Ljung"),fitdf=length(coef(m6)))
Box.test(m6$residuals,lag=8,type=c("Ljung"),fitdf=length(coef(m6)))
Acf(m6$residuals)
Pacf(m6$residuals)
m7=Arima(ORts, order=c(0,1,6), seasonal=list(order=c(0,1,2), period=12), method="ML")
coeftest(m7)
m7=Arima(ORts, order=c(0,1,6), seasonal=list(order=c(0,1,2), period=12), method="ML",fixed=c(NA,0,NA,0,0,NA,NA,NA))
coeftest(m7)
m7=Arima(ORts, order=c(0,1,6), seasonal=list(order=c(0,1,2), period=12), method="ML",fixed=c(NA,0,0,0,0,NA,NA,NA))
coeftest(m7)
plot(m7$residuals)
qqnorm(m7$residuals)
qqline(m7$residuals,col=2)
normalTest(m7$residuals,method=c("jb"))
Box.test(m7$residuals,lag=8,type=c("Ljung"),fitdf=4)
Box.test(m7$residuals,lag=10,type=c("Ljung"),fitdf=4)
Acf(m7$residuals)
Pacf(m7$residuals)

#6 Forecast
f5=forecast.Arima(m5,h=12)
f5
plot(f5,include=50)
plot(m5$x)
lines(fitted(m5),col=2)

f6=forecast.Arima(m6,h=12)
f6
plot(f6,include=50)
plot(m6$x)
lines(fitted(m6),col=2)

f7=forecast.Arima(m7,h=12)
f7
plot(f7,include=50, main="Occupied Rm 12mth Forecasts")
plot(m7$x)
lines(fitted(m7),col=2)

#7 Backtesting
source("backtest.R")
pm5=backtest(m5,ORts,120,1)
pm6=backtest(m6,ORts,120,1)
pm7=backtest(m7,ORts,120,1)

#### SARIMA MONTHLY RevPAR FORECAST#####
#1 Normal?
RP=myd[,4]
normalTest(RP,method = c("jb"))
qqnorm(RP)
qqline(RP,col=2)

#2 Plot TS
RPts=ts(RP,start=c(2005,1),frequency = 12)
plot(RPts,main = "Fig 1_3 Monthly RevPAR TS Plot")

par(mfcol=c(1,2))
hist(RPts, xlab="RevPAR", prob=TRUE, main="Fig 2_3 Histogram RevPAR")
qqnorm(RPts)
qqline(RPts, col = 2)
normalTest(RPts,method=c("jb"))

#ACF
acf(as.vector(RPts),lag.max=30, main="Fig3_3(a) ACF of RevPAR TS")
pacf(as.vector(RPts),lag.max=30, main="PACF of time series data")

#Apply first difference
dx=diff(RPts) 
acf(as.vector(dx),lag.max=50, main="Fig3_3(b) ACF of 1st diff RevPAR")
pacf(as.vector(dx),lag.max=50, main="PACF of time series data")

#Apply seasonal difference
sdx=diff(dx,12)
acf(as.vector(sdx),lag.max=50, main="Fig3_3(c) ACF of DSDX")
Pacf(as.vector(sdx),lag.max=50, main="PACF of DSDX starts")

#4 Identify p q d D
m8=auto.arima(RPts,stationary = F,seasonal = T,trace = T)

#5 Fit Model and Residuals Analysis
m8=Arima(RPts, order=c(2,0,0), seasonal=list(order=c(2,1,1), period=12), method="ML")
coeftest(m8)
plot(m8$residuals)
qqnorm(m8$residuals)
qqline(m8$residuals,col=2)
normalTest(m8$residuals,method=c("jb"))
Box.test(m8$residuals,lag=6,type=c("Ljung"),fitdf=length(coef(m8)))
Box.test(m8$residuals,lag=8,type=c("Ljung"),fitdf=length(coef(m8)))
m9=Arima(RPts, order=c(2,0,0), seasonal=list(order=c(0,1,1), period=12), method="ML")
coeftest(m9)
plot(m9$residuals)
qqnorm(m9$residuals)
qqline(m9$residuals,col=2)
normalTest(m9$residuals,method=c("jb"))
Box.test(m9$residuals,lag=5,type=c("Ljung"),fitdf=length(coef(m9)))
Box.test(m9$residuals,lag=7,type=c("Ljung"),fitdf=length(coef(m9)))
Acf(m9$residuals)
Pacf(m9$residuals)
########################### BASED ON ACF AND PACF##################
m10=Arima(RPts, order=c(2,0,1), seasonal=list(order=c(0,1,1), period=12), method="ML")
coeftest(m10)
plot(10$residuals)
qqnorm(m10$residuals)
qqline(m10$residuals,col=2)
normalTest(m10$residuals,method=c("jb"))
Box.test(m10$residuals,lag=5,type=c("Ljung"),fitdf=length(coef(m10)))
Box.test(m10$residuals,lag=7,type=c("Ljung"),fitdf=length(coef(m10)))
Acf(m10$residuals)
Pacf(m10$residuals)
################### BEST MODEL SO FAR######################
m11=Arima(RPts, order=c(1,0,1), seasonal=list(order=c(1,1,0), period=12), method="ML")
coeftest(m11)
plot(m11$residuals)
qqnorm(m11$residuals)
qqline(m11$residuals,col=2)
normalTest(m11$residuals,method=c("jb"))
Box.test(m11$residuals,lag=5,type=c("Ljung"),fitdf=length(coef(m11)))
Box.test(m11$residuals,lag=7,type=c("Ljung"),fitdf=length(coef(m11)))
Acf(m11$residuals)
Pacf(m11$residuals)

#6 Forecast
f8=forecast.Arima(m8,h=12)
f8
plot(f8,include=50)
plot(m8$x)
lines(fitted(m8),col=2)

f9=forecast.Arima(m9,h=12)
f9
plot(f9,include=50)
plot(m9$x)
lines(fitted(m9),col=2)

f10=forecast.Arima(10,h=12)
f10
plot(f10,include=50)
plot(m10$x)
lines(fitted(m10),col=2)

f11=forecast.Arima(m11,h=12)
f11
plot(f11,include=50, main="RevPAR 12 monthes Forecastes")
plot(m11$x)
lines(fitted(m11),col=2)

#7 Backtesting
source("backtest.R")
pm8=backtest(m8,RPts,120,1)
pm9=backtest(m9,RPts,120,1)
pm10=backtest(m10,RPts,120,1)
pm11=backtest(m11,RPts,120,1)

###########Conclusion##############
myd$Month.Dummy.Variables[1:132]= c("Jan","Feb","March","April","May","June","July","Aug","Sep","Oct","Nov","Dec")
myd$Month.Dummy.Variables[133:143] = c("Jan","Feb","March","April","May","June","July","Aug","Sep","Oct","Nov")
myd$Jan <- (myd$Month.Dummy.Variables=="Jan")*1
myd$Feb <- (myd$Month.Dummy.Variables=="Feb")*1
myd$March <- (myd$Month.Dummy.Variables=="March")*1
myd$April <- (myd$Month.Dummy.Variables=="April")*1
myd$May <- (myd$Month.Dummy.Variables=="May")*1
myd$June <- (myd$Month.Dummy.Variables=="June")*1
myd$July <- (myd$Month.Dummy.Variables=="July")*1
myd$Aug <- (myd$Month.Dummy.Variables=="Aug")*1
myd$Sep <- (myd$Month.Dummy.Variables=="Sep")*1
myd$Oct <- (myd$Month.Dummy.Variables=="Oct")*1
myd$Nov <- (myd$Month.Dummy.Variables=="Nov")*1
myd$Dec <- (myd$Month.Dummy.Variables=="Dec")*1
boxplot(myd$ADR~myd$Month.Dummy.Variables,las=1,main="BoxPlot of Monthly ADR",col=c(3,2),ylab="ADR in Dollars")
box()
axis(2,at=seq(100,300,10),seq(100,300,10),las=1)
par(mfcol=c(2,1))
acf(as.vector(diff(ADRts,12)),lag.max=50, main="ADR Model ACF Plot")
Pacf(as.vector(diff(ADRts,12)),lag.max=50, main="ADR Model PACF Plot")
par(mfcol=c(1,1))
boxplot(myd$Occupied.Rooms.in.thousands.~myd$Month.Dummy.Variables,las=1,main="BoxPlot of Occupied Rooms",col=c(3,2),ylab="Occupied Rooms In Thousands")
box()
axis(2,at=seq(100,300,10),seq(100,300,10),las=1)
acf(as.vector(sdx),lag.max=50, main="Occupied Rooms Model ACF Plot")
Pacf(as.vector(sdx),lag.max=50, main="Occupied Rooms Model PACF Plot")
boxplot(myd$RevPAR~myd$Month.Dummy.Variables,las=1,main="BoxPlot of Monthly RevPAR",col=c(3,2),ylab="RevPAR in Dollars")
box()
axis(2,at=seq(40,300,10),seq(40,300,10),las=1)
acf(as.vector(sdx),lag.max=50, main="RevPAR Model ACF Plot")
Pacf(as.vector(sdx),lag.max=50, main="RevPAR Model PACF Plot")

