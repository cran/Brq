plot.Brq <-
function(x, plottype=c("hist", "trace", "ACF", "traceACF", "histACF","tracehist",
 "traceACFhist"),Dim=1,breaks=30,lwd=1,col1=0,col2=1,col3=1,col4=1, ...){
k=ncol(as.matrix(x$beta[,Dim]))

if(k==2) par(mfrow=c(1,2))
if(k==3) par(mfrow=c(1,3))
if(k==4) par(mfrow=c(2,2))
if(k>4 & k<=12) par(mfrow=c(ceiling(k/3),3))
if(k>12) par(mfrow=c(ceiling(k/3),3))

plottype <- match.arg(plottype)

switch(plottype,
# Trace plot
trace=for(i in 1:k){
ts.plot(x$beta[,i],xlab="iterations",ylab="",main= noquote(names(coef(x)))[i],col=col4)
},

# Autocorrelation plot
ACF=for(i in 1:k){
acf(x$beta[,i],main= noquote(names(coef(x)))[i], col=col3)
},

# Trace and Autocorrelation plots 
traceACF={
par(mfrow=c(k,2))
for(i in 1:k){
ts.plot(x$beta[,i],xlab="iterations",ylab="",main= noquote(names(coef(x)))[i],col=col4)
acf(x$beta[,i],main= noquote(names(coef(x)))[i], col=col3)
}},

# Histogram and Autocorrelation plots 
histACF= {
par(mfrow=c(k,2))
for(i in 1:k){
hist(x$beta[,i],breaks=breaks,prob=TRUE, main="",xlab=noquote(names(coef(x)))[i], col=col1)
lines(density(x$beta[,i], adjust=2), lty="dotted", col=col2, lwd=lwd)
acf(x$beta[,i],main= noquote(names(coef(x)))[i], col=col3)
}},

# Trace and Histogram plots
tracehist= {
par(mfrow=c(k,2))
for(i in 1:k){
ts.plot(x$beta[,i],xlab="iterations",ylab="",main= noquote(names(coef(x)))[i],col=col4)
hist(x$beta[,i],breaks=breaks,prob=TRUE, main="",xlab=noquote(names(coef(x)))[i], col=col1)
lines(density(x$beta[,i], adjust=2), lty="dotted", col=col2, lwd=lwd)
}},

# Histogram, Autocorrelation and Trace plots 
traceACFhist={
par(mfrow=c(k,3))
for(i in 1:k){
ts.plot(x$beta[,i],xlab="iterations",ylab="",main= noquote(names(coef(x)))[i],col=col4)
hist(x$beta[,i],breaks=breaks,prob=TRUE, main="",xlab=noquote(names(coef(x)))[i], col=col1)
lines(density(x$beta[,i], adjust=2), lty="dotted", col=col2, lwd=lwd)
acf(x$beta[,i],main= noquote(names(coef(x)))[i], col=col3)
}},

# Histogram plot 
hist=for(i in 1:k){
hist(x$beta[,i],breaks=breaks,prob=TRUE, main="",xlab=noquote(names(coef(x)))[i], col=col1)
lines(density(x$beta[,i], adjust=2), lty="dotted", col=col2, lwd=lwd)
})

}