
\name{Btopm}
\alias{Btopm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Btopm}
\description{This function uses Bayesian information criteria to select promising quantile regression models.}
\usage{Btopm(formula, tau =0.5, family ="rq", Ce =0, runs =15000, burn =1000)}
\arguments{
   \item{formula}{Model formula.}
  \item{tau}{The quantile of interest. Must be between 0 and 1.}
  \item{Ce}{Censored point for Tobit quantile regression.}
  \item{family}{The fitting method: either "rq" for  quantile regression with continuous response, "bin" for binary 
    quantile regression,  "Lcrq" for  Tobit quantile regression with left censored point or,  "Rcrq" for  Tobit quantile regression 
    with right censored point.}
  \item{runs}{Length of desired Gibbs sampler output.}
  \item{burn}{Number of Gibbs sampler iterations before output is saved.}
}
%\details{For further details see  Alhamzawi and  Yu  (2012).  }
\examples{
# Example 1  (Quantile regression).
data(stackloss)
fit=Btopm(stack.loss ~ stack.x, tau=0.05, family="rq", runs=1000, burn=500)
fit$coeff
fit=Btopm(stack.loss ~ stack.x, tau=0.50, family="rq", runs=1000, burn=500)
fit$coeff
fit=Btopm(stack.loss ~ stack.x, tau=0.95, family="rq", runs=1000, burn=500)
fit$coeff



}




