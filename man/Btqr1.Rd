\name{Btqr}
\alias{Btqr}
\title{
Bayesian tobit quantile regression 
}
\description{
 This function implements the idea of Bayesian tobit quantile regression  employing
 a likelihood function that is based on the asymmetric Laplace distribution (Yu and  Stander, 2007). The asymmetric
 Laplace  error distribution is written as  scale mixtures of normal distributions as in Reed and Yu (2009). 
}
\usage{
Btqr(x,y, tau = 0.5, left = 0,  runs = 11000, burn = 1000, thin=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
\code{Matrix of predictors.} 
}
\item{y}{
\code{Vector of dependent variable. } 
}
  \item{tau}{
 \code{The quantile of interest. Must be between 0 and 1.} 
}
  \item{left}{
 \code{ Left censored point.}
}
  \item{runs}{
 \code{Length of desired Gibbs sampler output.} 
}
  \item{burn}{\code{Number of Gibbs sampler iterations before output is saved.} }
 \item{thin}{\code{thinning parameter of MCMC draws.}}
}

\examples{
# Example 
set.seed(12345)
x <- abs(rnorm(100))
y <- -0.5 + x +(.25 + .25*x)*rnorm(100)
plot(x,y, type="n")
h <-(y > 0)
points(x[h],y[h],cex=.9,pch=16)
points(x[!h],y[!h],cex=.9,pch=1)
y <- pmax(0,y)
for(tau in (2:8)/9){
fit=Brq(y~x,tau=tau, method="Btqr", left=0, runs=1000, burn=500)$coef
# Note: runs =11000 and burn =1000
Xs=sort(x)
Xc=cbind(1,sort(x))
Xcf=Xc\%*\%c(fit)
Xcfp=pmax(0,Xcf)
lines(Xs,Xcfp,col="red")}
}
\author{
  Rahim Alhamzawi
}
