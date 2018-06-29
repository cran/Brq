\name{BALqr}
\alias{BALqr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian adaptive Lasso quantile regression}
\description{This function implements the idea of Bayesian adaptive Lasso quantile regression employing a likelihood function that is based on
 the asymmetric Laplace distribution. The asymmetric Laplace  error distribution is written as a scale mixture of normals
  as in Reed and Yu (2009). The proposed method  (\code{BALqr}) extends the Bayesian Lasso quantile regression   by allowing different penalization parameters for different regression
   coeffficients (Alhamzawi et al., 2012). }
\usage{BALqr(x,y, tau = 0.5, runs = 11000, burn = 1000, thin=1)}
\arguments{
   \item{x}{Matrix of predictors.}
   \item{y}{Vector of dependent variable. }
  \item{tau}{The quantile of interest. Must be between 0 and 1.}
    \item{runs}{Length of desired Gibbs sampler output.}
  \item{burn}{Number of Gibbs sampler iterations before output is saved.}
  \item{thin}{thinning parameter of MCMC draws.}
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
fit = Brq(y~0+x,tau=0.5, method="BALqr",runs=5000, burn=1000)
summary(fit)

}
\author{
  Rahim Alhamzawi
}

\references{
[1] Alhamzawi, Rahim, Keming Yu, and Dries F. Benoit. (2012). Bayesian adaptive Lasso quantile regression. Statistical Modelling 12.3: 279-297.

[2]  Reed, C. and  Yu, K. (2009).  A partially collapsed Gibbs sampler for Bayesian quantile regression. \emph{Technical Report}. Department of Mathematical Sciences, Brunel
University. URL: \url{http://bura.brunel.ac.uk/bitstream/2438/3593/1/fulltext.pdf}.

}




