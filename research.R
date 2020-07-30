h<-read.csv("nas_big.csv")
require(TSA)
dim(h)
hh=rev(h[,7])
h_small<-hh[1:105]
ts_hs<-ts(h_small, start=1, frequency=1)
ts_hl<-ts(hh, start=c(2009,4,1), frequency=52)
win.graph(width=9.5, height=7,pointsize=12)
plot.ts(ts_hs, type="o", main="Nasdaq-from April 1st 2009 to April 1st 2011", xlab="week")

win.graph(width=9.5, height=7,pointsize=12)
plot.ts(ts_hl, type="o", main="Nasdaq-136")

ts_hsd=diff(ts_hs)
win.graph(width=9.5, height=7,pointsize=12)
plot.ts(ts_hsd, type="o", main="DiffNasdaq-105")
win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
acf(ts_hsd, lag=40)
pacf(ts_hsd, lag=40)

########################
g<-read.csv("tsx_big.csv")
dim(g)
gg=rev(g[,7])
g_small<-gg[1:105]
ts_gs<-ts(g_small,start=1, frequency=1)
ts_gl<-ts(gg,start=c(2009,4,1), frequency=52)
win.graph(width=9.5, height=7,pointsize=12)
plot.ts(ts_gs, type="o", main="TSX-from April 1st 2009 to April 1st 2011", xlab="week")

win.graph(width=9.5, height=7,pointsize=12)
plot.ts(ts_gl, type="o", main="TSX-136")
ts_gsd=diff(ts_gs)
win.graph(width=9.5, height=7,pointsize=12)
plot.ts(ts_gsd, type="o", main="DiffTSX-105")
win.graph(width=9.5, height=7,pointsize=12)
par(mfrow=c(2,1))        # set up the graphics 
acf(ts_gsd, lag=40)
pacf(ts_gsd, lag=40)
#########################
# if you want both on one graph
its=ts(cbind(ts_hl,ts_gl))
win.graph(width=9.5, height=7,pointsize=12)
plot.ts(its, type="o", plot.type="single",col=1:2, lty=1:2,main="NASdaq (Black -), TSX (Red --) adjusted closing price from April 1st 2009 to November 1st 2011", xlab="week", ylab="Nasdaq, TSX adjusted closing price")
abline(v=106,col="blue")

####################
par(mfrow=c(3,1))        # set up the graphics 
acf(ts_gsd, lag=40)
acf(ts_hsd, lag=40)
win.graph(width=9.5, height=7,pointsize=12)#diff series
ccf(ts_gsd, ts_hsd,lag=60) # cross correlation function

win.graph(width=9.5, height=7,pointsize=12)#initial series

ccf(ts_gs, ts_hs,lag=60)#problems because not stationary

###########################
#not good because we don't have periodic data and it is rather a continuous model
win.graph(width=9.5, height=7,pointsize=12)
periodogram(ts_gsd)
win.graph(width=9.5, height=7,pointsize=12)
periodogram(ts_hsd)
win.graph(width=9.5, height=7,pointsize=12)
#raw periodogram-better choice but not consistent
par(mfrow=c(2,1))  
spec.pgram(ts_gsd,taper=0,detrend="False",log="no")
spec.pgram(ts_hsd,taper=0,detrend="False",log="no")
#####################
#####################
win.graph(width=9.5, height=7,pointsize=12)
#try different parameters for smoothing L=9
par(mfrow=c(2,1))
k=kernel("daniell",4)  
plot(spec.pgram(ts_gsd,k,taper=0,plot=FALSE),plot.type="marginal",log="no")
plot(spec.pgram(ts_hsd,k,taper=0,plot=FALSE),plot.type="marginal",log="no")
#####################
# log_10, dB=decibels
win.graph(width=9.5, height=7,pointsize=12)
#try different parameters for smoothing L=9
par(mfrow=c(2,1))
k=kernel("daniell",4)  
plot(spec.pgram(ts_gsd,k,taper=0,plot=FALSE),plot.type="marginal",log="dB")
plot(spec.pgram(ts_hsd,k,taper=0,plot=FALSE),plot.type="marginal",log="dB")
##########################
#log 
win.graph(width=9.5, height=7,pointsize=12)
#try different parameters for smoothing L=9
par(mfrow=c(2,1))
k=kernel("daniell",4)  
plot(spec.pgram(ts_gsd,k,taper=0,plot=FALSE),plot.type="marginal",log="yes")
plot(spec.pgram(ts_hsd,k,taper=0,plot=FALSE),plot.type="marginal",log="yes")

##########################
win.graph(width=9.5, height=7,pointsize=12)
#now try l=15-maybe too smooth
par(mfrow=c(2,1))
k=kernel("daniell",7)   
plot(spec.pgram(ts_gsd,k,taper=0,plot=FALSE),plot.type="marginal",log="no")
plot(spec.pgram(ts_hsd,k,taper=0,plot=FALSE),plot.type="marginal",log="no", main="Nasdaq-Smoothered periodogram, L=15")

#####################
###### plot on the same graph
x=ts(cbind(ts_gsd,ts_hsd))
k=kernel("daniell",4)  
s=spec.pgram(x,k,taper=0)#graph of log spectrum including a 95% confidence iterval
win.graph(width=9.5, height=7,pointsize=12)
plot(s, plot.type="marginal", log="no", main="Smoothered periodograms, L=9, Nasdaq(red), TSX(Black)" )
#the bi-coherence
win.graph(width=9.5, height=7,pointsize=12)
plot(s, plot.type="coh", ci.lty=2)
#more smoothing seems to be needed
k=kernel("daniell",7)  
s=spec.pgram(x,k,taper=0)
win.graph(width=9.5, height=7,pointsize=12)
plot(s, plot.type="coh", ci.lty=2)
