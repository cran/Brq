\name{BLqr}
\alias{BLqr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian  Lasso quantile regression}
\description{This function implements the idea of Bayesian  Lasso quantile regression  using a likelihood function that is based on the asymmetric Laplace distribution (Rahim, 2016). The asymmetric
 Laplace  error distribution is written as  a scale mixture of normal distributions as in Reed and Yu (2009).  
 This function implements the Bayesian lasso for linear quantile regression models by assigning scale mixture of normal (SMN)
priors on the parameters and independent exponential priors on their variances. We propose an alternative
Bayesian analysis of the Bayesian lasso problem reported in Li, et al. (2010).   A Gibbs sampling algorithm for the Bayesian Lasso quantile regression is constructed by sampling the parameters from their full
conditional distributions.
 
 
 
  }
\usage{BLqr(x,y, tau = 0.5, runs = 11000, burn = 1000, thin=1)}
\arguments{
   \item{x}{Matrix of predictors.}
   \item{y}{Vector of dependent variable. }
  \item{tau}{The quantile of interest. Must be between 0 and 1.}
    \item{runs}{Length of desired Gibbs sampler output.}
  \item{burn}{Number of Gibbs sampler iterations before output is saved.}
  \item{thin}{thinning parameter of MCMC draws.}
}

\author{
  Rahim Alhamzawi
}

\references{

[1]  Alhamzawi, R. (2016). Bayesian variable selection in  quantile regression using asymmetric Laplace distribution. Working paper.

[2]  Reed, C. and  Yu, K. (2009).  A partially collapsed Gibbs sampler for Bayesian quantile regression. \emph{Technical Report}. Department of Mathematical Sciences, Brunel
University. URL: \url{http://bura.brunel.ac.uk/bitstream/2438/3593/1/fulltext.pdf}.

[3]  Li, Q., Xi, R. and Lin, N. (2010). Bayesian regularized quantile regression. Bayesian Analysis,
5(3): 533-56.
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

fit = Brq(y~0+x,tau=0.5, method="BLqr",runs=5000, burn=1000)
summary(fit)
model(fit)
}
