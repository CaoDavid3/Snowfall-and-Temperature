
# Working with Complete Data
#----------------------------------------------------------------------------------------------------------------------------------------------------------
library(readxl)
(edm1 <- read.csv("~/Downloads/update.csv"))
(edmt = edm1[,2])
#Training set
(edmt.1 = edmt[1:602])



# Plot the Data
require(TSA)
(edm_ts = ts((edmt.1),start = 1967, frequency=12)) #complete set
plot.ts(edm_ts, type="o", main="Edmonton Stony Plain Temperature from 1967 to 2017")

# Lag Plots
lag.plot=function(data1,max.lag=1,corr=TRUE,smooth=FALSE){ 
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


lag.plot(edm_ts,60,corr=TRUE, smooth=TRUE)

# Checking for Normality
hist(edm_ts, probability=TRUE, 12, main="Edmonton Mean Temperature") #normality does not hold
qqnorm(edm_ts) #same goes for the qq plots
qqline(edm_ts)
summary(edm_ts)

# Analyze the ACF and PACF of the Original Data
acf(edm_ts, lag.max=100, main="Original Temperature in Edmonton")
acf(edm_ts,lag.max=36)
pacf(edm_ts,lag.max=100, main= "Original Temperature in Edmonton ")
pacf(edm_ts, lag.max=30)

# BoxCox to find type of transformation
edm_ts1 = edm_ts+27 #Because we have negative values
var=BoxCox.ar(edm_ts1,method="yw") 
var$ci
var$mle 
# Based on the MLE, we were suggested to take the square of the data

sqrtran =(edm_ts)^2
trnfm=plot.ts(sqrtran, type="o", main= "Square Transformation of Mean Temperature")
# However, after doing so, it seems like it didn't do too much to the variance,
# as our original data did not look like a transformation was needed, but clearly
# the effects of outliers, as seen in both cases

# Trying to come up with a model
#Since we have a pattern every 12 months, let's try seasonal difference

seasondiff=diff(edm_ts, 12)
par(mfrow=c(2,1))
acf(seasondiff,lag.max=100, main="Seasonal Difference (12)") #cut off after lag 1?
pacf(seasondiff,lag.max=100, main="Seasonal Difference (12)") 
# After taking a seasonal difference, we realized that it was essential that
# we also considered a seasonal MA, as the ACF cuts off after lag 1
# and tails off after every lag

# To check whether we have the presence of two seasons, we checked the
# periodogram, which has confirmed that we only have 1. 
b=periodogram(edm_ts)
b$freq
b$spec
# when we counted to see at which frequency the peak occurs, 
# we found that it was at 52,which gave us approx 0.083 = 1/12 
b1=b$freq[52]
b1


#Based on the seasonal differenced time series, our immediate thought
# was to try a seasonal MA with difference model
model = arima(edm_ts, order=c(0,0,0),seasonal=list(order=c(0,1,1), period=12))
model
modelres=model$residuals
hist(modelres,probability=TRUE)
lines(density(modelres))
qqnorm(modelres)
qqline(modelres)
acf(modelres,lag.max=100)
pacf(modelres,lag.max=100)
tsdiag(model,gof=30)


#Try 1. 
model1=arima(edm_ts,order=c(3,0,3),seasonal=list(order=c(0,1,1), period=12))
model1
res1=model1$residuals
summary(res1)
hist(res1,probability=TRUE, 12, main="Residuals of (3,0,3)X(0,1,1)12")
lines(density(res1)) 
qqnorm(res1, main="Normal Q-Q Plot of Residuals Under (3,0,3)X(0,1,1)12")
qqline(res1)
par(mfrow=c(2,1))
acf(res1,lag.max=100, main="ACF of Residuals Under (3,0,3)X(0,1,1)12")
pacf(res1,lag.max=100, main="PACF of Residuals Under (3,0,3)X(0,1,1)12")
tsdiag(model1,gof=30)

tabl=armasubsets(y=res1,nar=3,nma=3,y.name='test',ar.method='ols')
plot(tabl)
detectAO(model1,robust=FALSE)
detectAO(model1)
detectIO(model1)
outlier1=arima(edm_ts,order=c(3,0,3),seasonal=list(order=c(0,1,1), period=12), xreg=as.numeric(seq(edm_ts)==c(73,146)))
outlier1
outres = outlier1$residuals
hist(outres,probability = TRUE,12, main="After Considering Outliers (3,0,3)X(0,1,1)12")
lines(density(outres))
summary(outres)
qqnorm(outres, main="Outliers + Residuals Under (3,0,3)X(0,1,1)12")
qqline(outres)
acf(outres,lag.max=100, main= "ACF After Outliers")
pacf(outres,lag.max=100, main= "PACF After Outliers")
tsdiag(outlier1, gof=30)
#-------------------------------------------------------------------------------------------


# Trying code from online sources (Google)
library(tsoutliers)
library(TSA)
library(lmtest)
library(astsa)
(outliers_excess_ts <- tso(edm_ts))
help("tso")

par(mfrow=c(1,1))
ao <- filter(edm_ts, filter = 0, method = "recursive")
plot(ao, main = "Additive Outlier - TC delta = 0", type = "s")
outliers_excess_ts <- tso(edm_ts, types = c("TC", "AO", "LS", "IO", "SLS"))
outliers_excess_ts # transient change outlier occuring at 1981 was identified
plot(outliers_excess_ts)
outliers_excess_ts$outliers
outliers_excess_ts$outliers$ind
outliers_excess_ts$outliers$time
(n = length(edm_ts))
(mo_tc = outliers("TC",outliers_excess_ts$outliers$ind))
(tc = outliers.effects(mo_tc,n))
(coefhat = as.numeric(outliers_excess_ts$outliers["coefhat"]))
(tc_effect = coefhat*tc)
tc_effect_ts = ts(tc_effect, frequency = frequency(edm_ts), start=start(edm_ts))
excess_wo_ts = edm_ts-tc_effect_ts
plot(cbind(edm_ts, excess_wo_ts, tc_effect_ts))

sarima(excess_wo_ts, p=1,d=0,q=0,P=0, Q=1, D=1, S=12)

arimax{TSA}
arimax_model <- arimax(edm_ts,
                       order = c(1,0,0),
                       seasonal = list(order = c(0,1,1), period = 12),
                       xtransf = data.frame(I1 = (1*(seq(edm_ts) == outliers_excess_ts$outliers$ind))),                                                                                      
                       transfer = list(c(1,0)),
                       method='ML')

summary(arimax_model)


mo <- outliers("AO", 10)
ao <- outliers.effects(mo, n)
plot(ao, type = "h", main = "AO: additive outlier")



#--------------------------------------

# Predicting with ARIMA(3,0,3)X(0,1,1)12
nrow(edmt)
trained = edmt[1:602,]
(trained)
# Time series of the training set
(edm_train1 = ts(na.omit(trained),start = 1967, frequency=12))
plot.ts(edm_train1, type="o", col="blue", main="Working With Training Set")

#ACF and PACF of the training set
acf(edm_train1,lag.max=100, main="ACF of the Training Set")
pacf(edm_train1,lag.max=100, main="PACF of the Training Set")
acf(edm_train1,lag.max=100)
pacf(edm_train1,lag.max=100)

#Fitting a model using the training set
# It wouldnt let me consider outliers, why?
fittrain=arima(edm_train1,order=c(3,0,3),seasonal=list(order=c(0,1,1), period=12))
plot(fittrain,n.ahead=10,type="b",pch=16, main="Prediction Using Training Set 1967-2017")
plot(fittrain,n.ahead=10,type="b",pch=16,xlim=c(2017,2018), main="Prediction Using Training Set 2017-2018")
grid(NULL,NA,lwd=2)
points(time(edm_ts),(edm_ts),type="b", col="blue", pch=16)


erro=function(xdata, preddata)
{
  #xdata are the real values
  # preddata are the predicted values
  # xdata and preddata should have the same dimension
  n = length(preddata)
  m=length(xdata)
  e=t(xdata)-t(preddata) # the error
  MSE=mean(e*e)
  MPE=mean(e/t(xdata))
  MAE=mean(abs(e))
  MAPE=mean(abs(e/t(xdata)))   
  list(MPE=MPE, MSE=MSE, MAE=MAE, MAPE=MAPE)
}

xdata1=edmt[603:612,] # Test Set
prem=predict(fittrain,n.ahead=10) #Predicting based on training set under fit01
prem
preddata=prem$pred
erro(xdata1,preddata)


diag1=function(xdata, fitit)
{
  n = length(xdata)
  k = length(fitit$coef)
  BIC = -2*fitit$loglik+k*log(n)
  
  AIC = -2*fitit$loglik+2*k
  AICc = AIC+2*(1+k)*(k+2)/(n-k-2)
  list(AIC=AIC, AICc=AICc, BIC=BIC)
}

diag1(edm_ts,fittrain)

#------------------------------------------------------------------------------#
#Prediction Using the seasonally differenced data
(edmts = ts(na.omit(trained),start = 1967, frequency=12))
seasonal = diff(edmts,12)
fit01=arima(seasonal,order=c(3,0,3),seasonal=list(order=c(0,0,1), period=12))
plot(fit01,n.ahead=10,type="b",pch=16, main="Applying Seasonal Difference 1967-2017")
plot(fit01,n.ahead=10,type="b",pch=16,xlim=c(2014,2018), main="Applying Seasonal Difference 2014-2018")
grid(NULL,NA,lwd=2)
season1=diff(edm_ts,12) #seasonal difference of the original
points(time((season1)),(season1),type="b", col="blue", pch=16)

erro=function(xdata, preddata)
{
  #xdata are the real values
  # preddata are the predicted values
  # xdata and preddata should have the same dimension
  n = length(preddata)
  m=length(xdata)
  e=t(xdata)-t(preddata) # the error
  MSE=mean(e*e)
  MPE=mean(e/t(xdata))
  MAE=mean(abs(e))
  MAPE=mean(abs(e/t(xdata)))   
  list(MPE=MPE, MSE=MSE, MAE=MAE, MAPE=MAPE)
}

xdata=edmt[603:612,] # Test Set
prem1=predict(fit01,n.ahead=10) #Predicting based on training set under fit01
prem1
preddata=prem1$pred
erro(xdata1,preddata)


diag1=function(xdata, fitit)
{
  n = length(xdata)
  k = length(fitit$coef)
  BIC = -2*fitit$loglik+k*log(n)
  
  AIC = -2*fitit$loglik+2*k
  AICc = AIC+2*(1+k)*(k+2)/(n-k-2)
  list(AIC=AIC, AICc=AICc, BIC=BIC)
}

diag1(season1,fit01)
#------------------------------------------------------------------------------------

#Try 2. 
model2=arima(edm_ts,order=c(2,1,3),seasonal=list(order=c(0,1,1), period=12))
model2
resid2=model2$residuals
summary(resid2)
hist(resid2,probability=TRUE, 12)
lines(density(resid2)) 
qqnorm(resid2)
qqline(resid2)
par(mfrow=c(2,1))
acf(resid2,lag.max=100)
pacf(resid2,lag.max=100)
tsdiag(model2,gof=30)

tabl=armasubsets(y=resid2,nar=5,nma=5,y.name='test',ar.method='ols')
plot(tabl)
detectAO(model2,robust=FALSE)
detectAO(model2)
detectIO(model2)
outlierfit2=arima(edm_ts,order=c(2,1,3),seasonal=list(order=c(0,1,1), period=12), xreg=as.numeric(seq(edm_ts)==c(146)))
outlierfit2
resid22 = outlierfit2$residuals
hist(resid22,probability = TRUE,12)
lines(density(resid22))
summary(resid22)
qqnorm(resid22)
qqline(resid22)
acf(resid22,lag.max=100)
pacf(resid22,lag.max=100)
tsdiag(outlierfit2, gof=30)

# Predictions Using Training Set
fit02=arima(edm_ts1,order=c(2,1,3),seasonal=list(order=c(0,1,1), period=12))
plot(fit02,n.ahead=10,type="b",pch=16)
plot(fit02,n.ahead=10,type="b",pch=16,xlim=c(2017,2018))
grid(NULL,NA,lwd=2)
points(time((edm_ts)),(edm_ts),type="b", col="purple", pch=16)


erro=function(xdata, preddata)
{
  #xdata are the real values
  # preddata are the predicted values
  # xdata and preddata should have the same dimension
  n = length(preddata)
  m=length(xdata)
  e=t(xdata)-t(preddata) # the error
  MSE=mean(e*e)
  MPE=mean(e/t(xdata))
  MAE=mean(abs(e))
  MAPE=mean(abs(e/t(xdata)))   
  list(MPE=MPE, MSE=MSE, MAE=MAE, MAPE=MAPE)
}

(xdata2=edmt[603:611]) # Test Set
prem2=predict(model2,n.ahead=10) #Predicting based on training set under fit01
prem2
preddata2=prem2$pred
erro(xdata2,preddata2)


diag1=function(xdata, fitit)
{
  n = length(xdata)
  k = length(fitit$coef)
  BIC = -2*fitit$loglik+k*log(n)
  
  AIC = -2*fitit$loglik+2*k
  AICc = AIC+2*(1+k)*(k+2)/(n-k-2)
  list(AIC=AIC, AICc=AICc, BIC=BIC)
}

diag1(edm_ts1,fit02)
#----------------------------------------------------------------------
# What if we worked on the seasonal data from the start
seasondiff=diff(edm_ts1,12)
fit022=arima(seasondiff,order=c(2,1,3),seasonal=list(order=c(0,0,1), period=12))
plot(fit022,n.ahead=10,type="b",pch=16)
plot(fit022,n.ahead=10,type="b",pch=16,xlim=c(2017,2018))
grid(NULL,NA,lwd=2)
season=diff(edm_ts,12)
points(time((season)),(season),type="b", col="purple", pch=16)


xdata2=edmt[603:612,] # Test Set
prem22=predict(fit022,n.ahead=10) #Predicting based on training set under fit01
prem22
preddata22=prem22$pred
erro(xdata2,preddata22)


diag1=function(xdata, fitit)
{
  n = length(xdata)
  k = length(fitit$coef)
  BIC = -2*fitit$loglik+k*log(n)
  
  AIC = -2*fitit$loglik+2*k
  AICc = AIC+2*(1+k)*(k+2)/(n-k-2)
  list(AIC=AIC, AICc=AICc, BIC=BIC)
}

diag1(season,fit022)
#--------------------------------------------------------------------------------------------

