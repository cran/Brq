\name{BLtqr}
\alias{BLtqr}
\title{
Bayesian tobit quantile regression 
}
\description{
 This function implements the idea of Bayesian  Lasso tobit quantile regression  using a likelihood function that is based on the asymmetric Laplace distribution
  (Rahim, 2016). The asymmetric Laplace  error distribution is written as a scale mixture of normal distributions as in Reed and Yu (2009).  
 This function implements the Bayesian lasso for linear tobit quantile regression models by assigning scale mixture of normal (SMN)
priors on the parameters and independent exponential priors on their variances.
 A Gibbs sampling algorithm for the Bayesian Lasso tobit quantile regression is constructed by sampling the parameters from their full
conditional distributions.
}
\usage{
BLtqr(x,y, tau = 0.5, left = 0,  runs = 11000, burn = 1000, thin=1)
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
  \item{burn}{
 \code{Number of Gibbs sampler iterations before output is saved.} }
 \item{thin}{\code{thinning parameter of MCMC draws.}}
}



\author{
  Rahim Alhamzawi
}


\examples{
# Example 
n <- 150
p=8
Beta=c(5, 0, 0, 0, 0, 0, 0, 0)
x <- matrix(rnorm(n=p*n),n)
x=scale(x)
y <-x\%*\%Beta+rnorm(n)
y=y-mean(y)
y=pmax(0, y)

fit = Brq(y~0+x,tau=0.5, method="BLtqr",runs=5000, burn=1000)
summary(fit)
model(fit)
}
