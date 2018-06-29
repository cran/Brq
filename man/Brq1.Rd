\name{Bqr}
\alias{Bqr}
\alias{Brq}
\alias{Brq.default}
\alias{Brq.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian Quantile Regression}
\description{This function implements the idea of Bayesian quantile regression employing
 a likelihood function that is based on the asymmetric Laplace distribution (Yu and  Moyeed, 2001). The asymmetric
 Laplace  error distribution is written as  scale mixtures of normal distributions as in Reed and Yu (2009).}
\usage{Bqr(x,y, tau =0.5, runs =11000, burn =1000, thin=1)}
\arguments{
   \item{x}{Matrix of predictors.}
   \item{y}{Vector of dependent variable. }
  \item{tau}{The quantile of interest. Must be between 0 and 1.}
    \item{runs}{Length of desired Gibbs sampler output.}
  \item{burn}{Number of Gibbs sampler iterations before output is saved.}
\item{thin}{thinning parameter of MCMC draws.}
}
\examples{
# Example 1
n <- 100
x <- runif(n=n,min=0,max=5)
y <- 1 + 1.5*x + .5*x*rnorm(n)
Brq(y~x,tau=0.5,runs=2000, burn=500)
Brq(y~x,tau=0.95,runs=2000, burn=500)


# Example 2
n <- 100
x <- runif(n=n,min=0,max=5)
y <- 1 + 1.5*x+ .5*x*rnorm(n)
plot(x,y, main="Scatterplot and Quantile Regression Fit", xlab="x", cex=.5, col="gray")
for (i in 1:5) {
if (i==1) p = .05
if (i==2) p = .25
if (i==3) p = .50
if (i==4) p = .75
if (i==5) p = .95
fit = Brq(y~x,tau=p,runs=1500, burn=500)
# Note: runs =11000 and burn =1000
abline(a=mean(fit$coef[1]),b=mean(fit$coef[2]),lty=i,col=i)
}
abline( lm(y~x),lty=1,lwd=2,col=6)
legend(x=-0.30,y=max(y)+0.5,legend=c(.05,.25,.50,.75,.95,"OLS"),lty=c(1,2,3,4,5,1),
lwd=c(1,1,1,1,1,2),col=c(1:6),title="Quantile")


}
\author{
  Rahim Alhamzawi
}