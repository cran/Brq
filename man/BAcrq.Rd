\name{BAcrq}
\alias{BAcrq}
\title{Bayesian adaptive lasso Tobit quantile regression}
\description{This function implements the Bayesian  Tobit quantile  regression with adaptive lasso penalty. The function  support  Tobit quantile regression with left or  right censored point. }
\usage{BAcrq(formula, tau=0.5, method="Lcrq", Ce=0, runs=15000, burn=1000)}
\arguments{
  \item{formula}{Model formula.}
  \item{tau}{The quantile of interest. Must be between 0 and 1.}
  \item{method}{The fitting method: either "Lcrq" for  Tobit quantile regression with left censored point or  "Rcrq" for  Tobit quantile regression 
    with right censored point.}
  \item{Ce}{Censored point for Tobit quantile regression.}
  \item{runs}{Length of desired Gibbs sampler output.}
  \item{burn}{Number of Gibbs sampler iterations before output is saved.}
}
\references{
[1] Alhamzawi, R. and Yu, K. and Benoit,  D. (2012). Bayesian adaptive Lasso quantile regression.  Statistical Modelling 2012; 12(3):  279-297.}
\examples{


# Example  (Tobit quantile regression with left censored point).
n=100
x <- rnorm(n)
y <- 0.5 + x + 0.5*x*rnorm(n)
y <- pmax(0,y)
fit=BAcrq(y~x, tau= 0.50,  method= "Lcrq", Ce=0, runs=2000, burn=1000)
fit$coeff
summary.brq(fit)

}




