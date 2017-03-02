\name{BTrq}
\alias{BTrq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian Tobit Quantile Regression}
\description{This function implements the idea of Bayesian  tobit quantile regression using
 a likelihood function that is based on the asymmetric Laplace distribution. The asymmetric
 Laplace  error distribution is written as  scale mixtures of normal distributions.
 In this function,   a two-level hierarchical Bayesian model is used. Specifically, I put   zero mean Gaussian priors on the regression coefficients with 
 non informative Jeffreys prior distributions for the unknown variances. }
\usage{BTrq(formula, tau =0.5, runs =11000, burn =1000)}
\arguments{
   \item{formula}{Model formula.}
  \item{tau}{The quantile of interest. Must be between 0 and 1.}
    \item{runs}{Length of desired Gibbs sampler output.}
  \item{burn}{Number of Gibbs sampler iterations before output is saved.}
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
fit=BTrq(y~x,tau=tau, runs=1000, burn=500)$coef
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
\references{

[1] Alhamzawi, R. (2014). Model selection in quantile regression models. Journal of Applied Statistics, 42, 445-458.

}



