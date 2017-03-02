\name{BBrq}
\alias{BBrq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian binary Quantile Regression}
\description{This function implements the idea of Bayesian  binary quantile regression using
 a likelihood function that is based on the asymmetric Laplace distribution. The asymmetric
 Laplace  error distribution is written as  scale mixtures of normal distributions.
 In this function,   a two-level hierarchical Bayesian model is used. Specifically, I put   zero mean Gaussian priors on the regression coefficients with 
 non informative Jeffreys prior distributions for the unknown variances. }
\usage{BBrq(formula, tau =0.5, runs =11000, burn =1000)}
\arguments{
   \item{formula}{Model formula.}
  \item{tau}{The quantile of interest. Must be between 0 and 1.}
    \item{runs}{Length of desired Gibbs sampler output.}
  \item{burn}{Number of Gibbs sampler iterations before output is saved.}
}
\examples{
#set.seed(1234)
n <- 200
x=rnorm(n)
ystar <-  1 + x + rnorm(n=n, mean=0, sd=1)
y=as.numeric(ystar>0)
fit = BBrq(y~x,tau=0.5,runs=3000, burn=1000)
# Note: runs =11000 and burn =1000
fit$coef
}
\author{
  Rahim Alhamzawi
}
\references{

[1] Alhamzawi, R. (2014). Model selection in quantile regression models. Journal of Applied Statistics, 42, 445-458.

}



