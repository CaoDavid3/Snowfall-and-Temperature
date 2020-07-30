#------------------------------------
# Stat 370
# Name:David Cao and Nirudika Velupillai
# Instructor: Cristina Anton 
# Project
#------------------------------------
library(readr)
require(TSA)

lag.plot1=function(data1,max.lag=1,corr=TRUE,smooth=FALSE){ 
  name1=paste(deparse(substitute(data1)),"(t-",sep="")
  name2=paste(deparse(substitute(data1)),"(t)",sep="")
  data1=as.ts(data1)
  max.lag=as.integer(max.lag)
  prow=ceiling(sqrt(max.lag))
  pcol=ceiling(max.lag/prow)
  a=acf(data1,max.lag,plot=FALSE)$acf[-1]
  par(mfrow=c(prow,pcol), mar=c(2.5, 4, 2.5, 1), cex.main=1.1, font.main=1)
  for(h in 1:max.lag){                       
    plot(lag(data1,-h), data1, xy.labels=FALSE, main=paste(name1,h,")",sep=""), ylab=name2, xlab="") 
    if (smooth==TRUE) 
      lines(lowess(ts.intersect(lag(data1,-h),data1)[,1],
                   ts.intersect(lag(data1,-h),data1)[,2]), col="red")
    if (corr==TRUE)
      legend("topright", legend=round(a[h], digits=2), text.col ="blue", bg="white", x.intersp=0)
  }
}

diag1=function(xdata, fitit)
{
  n = length(xdata)
  k = length(fitit$coef)
  BIC = -2*fitit$loglik+k*log(n)
  
  AIC = -2*fitit$loglik+2*k
  AICc = AIC+2*(1+k)*(k+2)/(n-k-2)
  list(AIC=AIC, AICc=AICc, BIC=BIC)
}

erro=function(xdata, preddata)
{
  #xdata are the real values
  # preddata are the predicted values
  # xdata and preddata should have the same dimension
  n = length(preddata)
  m=length(xdata)
  e=xdata-preddata# the error
  MSE=mean(e*e)
  MPE=mean(e/xdata)
  MAE=mean(abs(e))
  MAPE=mean(abs(e/xdata))   
  list(MPE=MPE, MSE=MSE, MAE=MAE, MAPE=MAPE)
}
Edmonton <- read_csv("EMonthData.csv")
EdmontonFull <- Edmonton[,4] # full data set
ESnowFull <- ts(EdmontonFull, start = 1967, frequency = 12)

ESnow <- ts(EdmontonFull, start=1967,frequency=12)
win.graph(width=10, height=10,pointsize=12)
plot.ts(ESnow, type="o", pch = 8, cex = 0.75, main="Snow fall (cm) annually in 1967-2017")

summary(ESnow)
win.graph(width=10, height=10,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
hist(ESnow, prob=TRUE, 12)   # histogram    
lines(density(ESnow))     # smooth it - ?density for details 
qqnorm(ESnow)             # normal Q-Q plot  
qqline(ESnow)             # add a line    

win.graph(width=10, height=10,pointsize=12)
par(mfrow=c(1,2))        # set up the graphics 
acf(ESnow, lag.max = 50)
pacf(ESnow, lag.max = 50)
eacf(ESnow)

win.graph(width=15, height=15,pointsize=12)
lag.plot1(ESnow,12,corr=TRUE,smooth=TRUE)

win.graph(width=15, height=15,pointsize=12)
bct<-BoxCox.ar(ESnow + 1,method='yw')
bct$ci
bct$mle

logESnow = log(ESnow + 1)
win.graph(width=10, height=10,pointsize=12)
plot.ts(logESnow, type="o", pch = 8, cex = 0.75, main="Snow fall(cm) annually in 1937 to 2015")

#------------------------------------------
# Take the difference of the time series by 12
#------------------------------------------
DESnow <- diff(logESnow,12)
win.graph(width=10, height=10,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
acf(DESnow, lag.max = 50)
pacf(DESnow, lag.max = 50)
eacf(DESnow)

win.graph(width=10, height=10,pointsize=12)
plot.ts(DESnow, type="o", pch = 8, cex = 0.75, main="Snow fall(cm) annually in 1937")

summary(DESnow)
win.graph(width=10, height=10,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
hist(DESnow, prob=TRUE, 12)   # histogram    
lines(density(DESnow))     # smooth it - ?density for details 
qqnorm(DESnow)             # normal Q-Q plot  
qqline(DESnow)             # add a line    

#-------------------------
#fit1
#-------------------------
fit1 = arima(logESnow, order = c(0,0,0), seasonal = list(order = c(2,1,2), period = 12))
fit1$coef #(parameter estimates)
sqrt(diag(fit1$var.coef))# standard errors
fit1$sigma2 # noise variance
res1=residuals(fit1)

win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
hist(res1, prob=TRUE, 12)   # histogram    
lines(density(res1))     # smooth it - ?density for details 
qqnorm(res1)             # normal Q-Q plot  
qqline(res1)             # add a line 

win.graph(width=9.5, height=12,pointsize=12)
par(mfrow=c(2,1))
acf(res1, lag.max =80)
pacf(res1, lag.max=80)

win.graph(width=9.5, height=7,pointsize=12)
tsdiag(fit1,gof=20,omit.initial=T)
detectAO(fit1) # detect outliers 
detectAO(fit1, robust=F)
detectIO(fit1)

#----------------
# fit 2
#----------------
fit2=arima(logESnow,order=c(0,0,0),seasonal=list(order=c(2,1,2), period=12),xreg=as.numeric (seq (logESnow) == 245,409))
fit2$coef #(parameter estimates)
sqrt(diag(fit2$var.coef))# standard errors
fit2$sigma2 # noise variance
res2=residuals(fit2)

win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
hist(res2, prob=TRUE, 12)   # histogram    
lines(density(res2))     # smooth it - ?density for details 
qqnorm(res2)             # normal Q-Q plot  
qqline(res2)             # add a line 

win.graph(width=9.5, height=12,pointsize=12)
par(mfrow=c(2,1))
acf(res2, lag=40)
pacf(res2, lag=40)

win.graph(width=9.5, height=7,pointsize=12)
tsdiag(fit2,gof=20,omit.initial=T)
detectAO(fit2) # detect outliers 
detectAO(fit2, robust=F)
detectIO(fit2)

#-----------------------------
# Prediction
#-----------------------------
###################  Predictions
win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(2,1)) 
plot(fit2,n.ahead=10,type='b', pch=16, newxreg = 0)
plot(fit2,n.ahead=10,type='b', pch=16, xlim=c(2015,2019), newxreg = 0)
points(time(logESnow), logESnow, col = 2, pch = 16, type = 'b')

#Prediction errors 
xESnow=logESnow[603:612]
prem1=predict(fit2, n.ahead = 10, newxreg = 0)
preddata1=(prem1$pred)
preddata1
erro(xESnow, preddata1)

diag1(logESnow,fit2)

#Prediction errors
xESnow=logESnow[603:612]
xESnow <- exp(xESnow) +1

win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(1,1)) 
#plot(fit2,n.ahead=10,type='b', pch=16, newxreg = 0)
plot(fit2,n.ahead=10,type='b', pch=16, xlim=c(2015,2019), newxreg = 0, main="Prediction for ARIMA(0,0,0) X (2,1,2)12")
points(time(logESnow), logESnow, col = c("#006666"), pch = 20, type = 'b')
abline(v = c(2015,2016,2017,2018,2019), lty = 2)

prem1=predict(fit2, n.ahead = 10, newxreg = 0)
preddata1=(prem1$pred)
preddata1
erro(xESnow, preddata1)

diag1(logESnow,fit2)

#-----------------------------------------------
# Comparision between the two 
#-----------------------------------------------
logESnowFull = log(ESnowFull + 1)
ESnowD <- diff(logESnowFull, 12)

edm1 <- read_csv("update.csv")
edmt = edm1[,3]
edm_ts = ts(edmt,start = 1967, frequency=12) #complete set
win.graph(width=9.5, height=7,pointsize=12)
plot.ts(edm_ts, type="o", main="Edmonton Stony Plain Temperature from 1967 to 2017")

edm_tsd <- diff(edm_ts,12)

##################################################
win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(3,1))        # set up the graphics 
acf(edm_tsd, lag=40)
acf(ESnowD, lag=40)

#not good because we don't have periodic data and it is rather a continuous model
win.graph(width=9.5, height=7,pointsize=12)
periodogram(edm_tsd)
win.graph(width=9.5, height=7,pointsize=12)
periodogram(ESnowD)
#raw periodogram-better choice but not consistent
win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(2,1))  
spec.pgram(edm_tsd,taper=0,detrend="False",log="no")
spec.pgram(ESnowD,taper=0,detrend="False",log="no")

#####################
win.graph(width=9.5, height=7,pointsize=12)
#try different parameters for smoothing L=9
par(mfrow=c(2,1))
k=kernel("daniell",4)  
plot(spec.pgram(edm_tsd,k,taper=0,plot=FALSE),plot.type="marginal",log="no")
plot(spec.pgram(ESnowD,k,taper=0,plot=FALSE),plot.type="marginal",log="no")
#####################
# log_10, dB=decibels
win.graph(width=9.5, height=7,pointsize=12)
#try different parameters for smoothing L=9
par(mfrow=c(2,1))
k=kernel("daniell",4)  
plot(spec.pgram(edm_tsd,k,taper=0,plot=FALSE),plot.type="marginal",log="dB")
plot(spec.pgram(ESnowD,k,taper=0,plot=FALSE),plot.type="marginal",log="dB")
##########################
#log 
win.graph(width=9.5, height=7,pointsize=12)
#try different parameters for smoothing L=9
par(mfrow=c(2,1))
k=kernel("daniell",4)  
plot(spec.pgram(edm_tsd,k,taper=0,plot=FALSE),plot.type="marginal",log="yes")
plot(spec.pgram(ESnowD,k,taper=0,plot=FALSE),plot.type="marginal",log="yes")
##########################
win.graph(width=9.5, height=7,pointsize=12)
#now try l=15-maybe too smooth
par(mfrow=c(2,1))
k=kernel("daniell",7)   
plot(spec.pgram(edm_tsd,k,taper=0,plot=FALSE),plot.type="marginal",log="no", main="Nasdaq-Smoothered periodogram, L=15")
plot(spec.pgram(ESnowD,k,taper=0,plot=FALSE),plot.type="marginal",log="no")
#####################


###### plot on the same graph
x=ts(cbind(ESnowD,edm_tsd))
k=kernel("daniell",7)  
s=spec.pgram(x,k,taper=0)#graph of log spectrum including a 95% confidence iterval
win.graph(width=9.5, height=7,pointsize=12)
plot(s, plot.type="marginal", log="no", main="Smoothered periodograms, L=15, ESnowD(red), edm_tsd(Black)" )
#the bi-coherence
win.graph(width=9.5, height=7,pointsize=12)
plot(s, plot.type="coh", ci.lty=2)
#more smoothing seems to be needed
k=kernel("daniell",25)  
s=spec.pgram(x,k,taper=0)
win.graph(width=9.5, height=7,pointsize=12)
plot(s, plot.type="coh", ci.lty=2)



