\name{BALtqr}
\alias{BALtqr}
\title{Bayesian adaptive Lasso tobit quantile regression}
\description{This function implements the idea of Bayesian adaptive Lasso tobit quantile regression employing a likelihood function that is based on
 the asymmetric Laplace distribution. The asymmetric Laplace  error distribution is written as  a scale mixture of normal distributions
  as in Reed and Yu (2009). The proposed method (\code{BALtqr}) extends the Bayesian Lasso tobit quantile regression  by allowing different penalization parameters for different regression
   coeffficients (Alhamzawi et al., 2013).}
\usage{
BALtqr(x,y, tau = 0.5, left = 0,  runs = 11000, burn = 1000, thin=1)
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
 \code{Number of Gibbs sampler iterations before output is saved.} 
}
\item{thin}{\code{thinning parameter of MCMC draws.}}

}



\author{
  Rahim Alhamzawi
}

\references{
[1] Alhamzawi, Rahim. (2013). Tobit Quantile Regression with the adaptive Lasso penalty. The Fourth International Arab Conference of Statistics, 450 ISSN (1681 6870). 

[2]  Reed, C. and  Yu, K. (2009).  A partially collapsed Gibbs sampler for Bayesian quantile regression. \emph{Technical Report}. Department of Mathematical Sciences, Brunel
University. URL: \url{http://bura.brunel.ac.uk/bitstream/2438/3593/1/fulltext.pdf}.

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

fit = Brq(y~0+x,tau=0.5, method="BALtqr",runs=5000, burn=1000)
summary(fit)
model(fit)
}


