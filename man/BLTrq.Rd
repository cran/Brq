\name{BLTrq}
\alias{BLTrq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian Lasso Tobit quantile regression}
\description{This function implements the idea of Bayesian  Lasso tobit quantile regression  using a likelihood function that is based on the asymmetric Laplace distribution (Rahim, 2016). The asymmetric
 Laplace  error distribution is written as  scale mixtures of normal distributions as in Reed and Yu (2009). We specify  zero mean Gaussian priors to the 
 regression coefficients with an independent prior distributions for the unknown variances.}
\usage{BLTrq(formula, tau =0.5, runs =11000, burn =1000)}
\arguments{
   \item{formula}{Model formula.}
  \item{tau}{The quantile of interest. Must be between 0 and 1.}
    \item{runs}{Length of desired Gibbs sampler output.}
  \item{burn}{Number of Gibbs sampler iterations before output is saved.}
}
\examples{
set.seed(12345)
n <- 100
p=8
Beta=c(5, 0, 0, 0, 0, 0, 0, 0)
x <- matrix(rnorm(n=p*n),n)
y <-x\%*\%Beta+ rnorm(n)
y=pmax(y,0)
fit=BLTrq(y~0+x,tau=0.5, runs=5000, burn=1000)
# Note: runs =11000 and burn =1000
fit$coef
}

\author{
  Rahim Alhamzawi
}
\references{

[1]  Alhamzawi, R. (2016). Bayesian variable selection in  quantile regression using asymmetric Laplace distribution. Working paper.

[2]  Reed, C. and  Yu, K. (2009).  A partially collapsed Gibbs sampler for Bayesian quantile regression. \emph{Technical Report}. Department of Mathematical Sciences, Brunel
University. URL: \url{http://bura.brunel.ac.uk/bitstream/2438/3593/1/fulltext.pdf}.

}