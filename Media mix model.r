

##source(url("http://www.stat.pitt.edu/stoffer/tsa2/Rcode/itall.R"))
##z<- as.matrix(read.table("C:/Users/mccj07/Documents/Media mix modeling/MMM_data.csv", sep=",", header=T), 104, 4, byrow = TRUE)

### Lets first generate some Random DAta
x1=2000*10^arima.sim(n = 104, list(ar = c(0.8897, -0.178)),sd = sqrt(0.0796))
x2=5000*10^arima.sim(n = 104, list(ar = c(0.9497, -0.1758)),sd = sqrt(0.1296))
x3=1000*10^arima.sim(n = 104, list(ar = c(0.9197, -0.1258)),sd = sqrt(0.0696))
y=10000*10^arima.sim(n = 104, list(ar = c(0.9097, -0.1458), ma = c(-0.279, 0.188)),sd = sqrt(0.1596))

y=y+.75*x2+1.2*x1+.9*x3

xxxy=cbind(x1,x2,x3,y) 
plot(xxxy,type='l')

xxxy 


##  Lets go ahead and log the data
y1.l=log(y)
x1.l=log(x1)
x2.l=log(x2)
x3.l=log(x3)

## lets remove the mean for each series. 
y1.m=mean(y1.l)
x1.m=mean(x1.l)
x2.m=mean(x2.l)
x3.m=mean(x3.l)


y1a=y1.l-y1.m
x1a=x1.l-x1.m
x2a=x2.l-x2.m
x3a=x3.l-x3.m

##   The two series X and Y must be made stationary before starting interpreting 
##   the crosscorrelation function. Means of X and Y and variances of X and Y must 
##   be constant over time. So the Autocorrelation Function of X and Y must be looked at 
##   first of all and appropriate pre-transformations and differencing transformations 
##   must be made until the ACF dies down quickly or has a cut off kind of behavior.

par(mfcol=c(2,2))
acf(y1a,lag=20)
acf(x1a,lag=20)
pacf(y1a,lag=20)
pacf(x1a,lag=20)

par(mfcol=c(2,2))
acf(x2a,lag=20)
acf(x3a,lag=20)
pacf(x2a,lag=20)
pacf(x3a,lag=20)

## in this case we generated staionary time series . so all are stationary



xxxya=cbind(x1a,x2a,x3a,y1a) 
lag.plot(xxxya, 4, do.lines=FALSE) 

x1_f=arima(x1a,order=c(2,0,0))
x2_f=arima(x2a,order=c(2,0,0))
x3_f=arima(x3a,order=c(2,0,0))

x1_f
x2_f
x3_f


###y1=arima(z[,1],order=c(2,0,0))
### choos x2 as a filter 
tsdiag(x1_f, gof.lag=20)
tsdiag(x2_f, gof.lag=20)
tsdiag(x3_f, gof.lag=20)


f1=c(1,-x2_f$coef[1:2])


y1fil=filter(y1a, sides=1,f1)
y1fil=y1fil[3:103]

x1fil=filter(x1a, sides=1,f1)
x1fil=x1fil[3:103]
x2fil=filter(x1a, sides=1,f1)
x2fil=x2fil[3:104]
x3fil=filter(x1a, sides=1,f1)
x3fil=x3fil[3:104]

##   Identify an ARIMA model for the input series X that you made stationary and apply this model to Y.

##   Get the residuals from the model for X and the residuals of the model for Y. This is called Prewhitening the series.


ccf(x1fil,y1fil, ylab = "cross-correlation x1 and y")
ccf(x2fil,y1fil, ylab = "cross-correlation x2 and y")
ccf(x3fil,y1fil, ylab = "cross-correlation x3 and y")

x1f=x1_f$residuals[3:104]
x2f=x2_f$residuals[3:104]
x3f=x3_f$residuals[3:104]



##   Find the cross-correlation function between the residuals. This crosscorrelation function allows 
##   us to find the impulse response function.

###   Use the estimates of the impulse response function to make guesses of the orders of the actual Y and X models.

source("C:/Users/mccj07/Documents/R/source/ccm.R")

xy=cbind(x1f,y1a[3:104],x2f,x3f) 
##<== Combine filtered series.
ccm(xy,20)

##   With the latter get initial estimates of the parameters. 
###  Estimate the parameters for the Y and the X part of the model


(fit2.gls = arima(y1a[3:104], order=c(2,0,2), xreg=cbind(x1f, x2f, x3f))) 
acf(fit2.gls$residuals)
pacf(fit2.gls$residuals)


Box.test(resid(fit2.gls), 12, type="Ljung")
pchisq(12.377, 10, lower=FALSE) 


x1.p=predict(x1_f, n.ahead=13, level=c(80,95))
x2.p=predict(x2_f, n.ahead=13, level=c(80,95))
x3.p=predict(x3_f, n.ahead=13, level=c(80,95))

newxreg=cbind(x1.p$pred,x2.p$pred,x3.p$pred)
y.fore = predict(fit2.gls, n.ahead=13, newxreg=newxreg)  

exp(y.fore$pred)
y.fore$pred+y1.m


U = exp((y.fore$pred + 2*y.fore$se)+y1.m)
L = exp((y.fore$pred - 2*y.fore$se)+y1.m)
miny=min(y,L)
maxy=max(y,U)
ts.plot(y,exp(y.fore$pred+y1.m),col=1:2, ylim=c(miny,maxy))
lines(U, col="blue", lty="dashed")
lines(L, col="blue", lty="dashed")

exy=ts(exp(y.fore$pred+y1.m),start=105)

pred.xxxy=cbind(exp(x1.p$pred+x1.m),exp(x2.p$pred+x2.m),exp(x3.p$pred+x3.m),exy)

pred.xxxy