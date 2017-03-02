\name{BALrq}
\alias{BALrq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian adaptive Lasso quantile regression}
\description{This function implements the idea of Bayesian adaptive Lasso quantile regression employing a likelihood function that is based on
 the asymmetric Laplace distribution (Yu and  Moyeed, 2001). The asymmetric Laplace  error distribution is written as  scale mixtures of normal distributions
  as in Reed and Yu (2009). The method extends the Bayesian Lasso quantile regression  by allowing different penalization parameters for different regression
   coeffficients (Alhamzawi et al., 2012).}
\usage{BALrq(formula, tau =0.5, runs =11000, burn =1000)}
\arguments{
   \item{formula}{Model formula.}
  \item{tau}{The quantile of interest. Must be between 0 and 1.}
    \item{runs}{Length of desired Gibbs sampler output.}
  \item{burn}{Number of Gibbs sampler iterations before output is saved.}
}
\examples{
# Example 
n <- 150
p=8
Beta=c(5, 0, 0, 0, 0, 0, 0, 0)
x <- matrix(rnorm(n=p*n),n)
x=scale(x)
y <-x\%*\%Beta+ rnorm(n)
y=y-mean(y)
fit = BALrq(y~0+x,tau=0.5,runs=5000, burn=1000)
# Note: runs =11000 and burn =1000
fit$coef

}
\author{
  Rahim Alhamzawi
}

\references{
[1] Alhamzawi, Rahim, Keming Yu, and Dries F. Benoit. (2012). Bayesian adaptive Lasso quantile regression. Statistical Modelling 12.3: 279-297.

[2]  Reed, C. and  Yu, K. (2009).  A partially collapsed Gibbs sampler for Bayesian quantile regression. \emph{Technical Report}. Department of Mathematical Sciences, Brunel
University. URL: \url{http://bura.brunel.ac.uk/bitstream/2438/3593/1/fulltext.pdf}.

[3]  Yu, K. and Moyeed, R.A. (2001).  Bayesian Quantile Regression.  \emph{Statistics  &  Probability Letters}, 54, 437--447.
URL: \url{http://www.sciencedirect.com/science/article/pii/S0167715201001249}.


}

\seealso{
[1]  Li, Q., Xi, R. and Lin, N. (2010). Bayesian regularized quantile regression. Bayesian Analysis,
5(3): 533-56.
}



