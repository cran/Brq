\name{brq}
\alias{brq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{brq}
\description{This function implements the idea of Bayesian  estimation in  quantile regression models employing
 a likelihood function that is based on the asymmetric Laplace distribution (Yu and  Moyeed, 2009). The asymmetric
 Laplace  error distribution is written as  scale mixtures of normal distributions as in Reed and Yu (2009) and 
 Kozumi and Kobayashi (2009). We specify  zero mean Gaussian priors to the regression coefficients with an independent
 prior distributions for the unknown variances. Then, the variance of the prior is motivated by the power of the correlation between the
 response variable  and the predictors. The combination between variance of the prior and the correlation between the response variable  and the 
predictors  motivates the amount of shrinkage in the prior. This two-level hierarchical Bayesian model is similar to the Bayesian
 Lasso model  in Li et al. (2010). However, this prior is preferable as there is no need to the penalty parameter in 
the Lasso model (Li et al., 2010) neither need to specify any hyper-parameter or any  threshold values.
 We may look to this prior as  adaptively shrinks the quantile regression estimate to zero with shrinkage increasing as long 
as  the regression coefficient  not active in the model. This function also support Tobit and  binary dependent variables.}
\usage{brq(formula, tau =0.5, family ="rq", Ce =0, runs =15000, burn =1000)}
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
fit=brq(stack.loss ~ stack.x, tau=0.05, family="rq", runs=1000, burn=500)
fit$coeff
fit=brq(stack.loss ~ stack.x, tau=0.25, family="rq", runs=1000, burn=500)
fit$coeff
fit=brq(stack.loss ~ stack.x, tau=0.50, family="rq", runs=1000, burn=500)
fit$coeff
fit=brq(stack.loss ~ stack.x, tau=0.75, family="rq", runs=1000, burn=500)
fit$coeff
fit=brq(stack.loss ~ stack.x, tau=0.95, family="rq", runs=1000, burn=500)
fit$coeff


# Example 2
n <- 150
x <- runif(n=n,min=0,max=5)
y <- 1 + 1.5*x+ .5*x*rnorm(n)
plot(x,y, main="Scatterplot and Quantile Regression Fit", xlab="x", cex=.5, col="gray")
for (i in 1:5) {
if (i==1) p = .05
if (i==2) p = .25
if (i==3) p = .50
if (i==4) p = .75
if (i==5) p = .95
fit = brq(y~x,tau=p,family="rq",runs=2000, burn=1000)
abline(a=mean(fit$c[1]),b=mean(fit$c[2]),lty=i,col=i)
}
abline( lm(y~x),lty=1,lwd=2,col=6)
legend(x=-0.30,y=max(y)+0.5,legend=c(.05,.25,.50,.75,.95,"OLS"),lty=c(1,2,3,4,5,1),
lwd=c(1,1,1,1,1,2),col=c(1:6),title="Quantile")


#Example (Tobit quantile regression with left censored point).
n=100
x <- rnorm(n)
y <- 0.5 + x + 0.5*x*rnorm(n)
y <- pmax(0,y)
fit=brq(y~x, tau= 0.50, family="Lcrq", runs=1000, burn=500)
fit$coeff



# Example (Binary quantile regression).
n = 200
x = matrix(runif(n=n*5, min=0, max=10),n)
ystar = -2*x[,1] + 2*x[,2] +rnorm(n)
y <- as.numeric(ystar>0)
fit=brq(y ~ x, tau=0.5, family="bin", runs=1000, burn=500)
fit$coeff

}


\references{

[1]Alhamzawi, R. (2012). Bayesian variable selection in binary and Tobit quantile regression models.

[2] Alhamzawi, R. and  Yu, K.  (2012). Conjugate priors and variable selection for Bayesian quantile regression. \emph{Computational Statistics and Data Analysis}. In press.
URL: \url{http://www.sciencedirect.com/science/article/pii/S0167947312000345}

[3]  Kozumi, H. and  Kobayashi, G. (2011). Gibbs sampling methods for Bayesian quantile regression.  \emph{Journal of Statistical Computation and
Simulation}. {bf  81},  1565-1578. URL: \url{http://www.tandfonline.com/doi/abs/10.1080/00949655.2010.496117}.

[4] Li, Q.,  Xi, R.  and Lin, N.  (2010). Bayesian regularized quantile regression. Bayesian Analysis, 5(3): 533-56.
\url{ba.stat.cmu.edu/journal/2010/vol05/issue03/lixin.pdf}

[5]  Reed, C. and  Yu, K. (2009).  A partially collapsed Gibbs sampler for Bayesian quantile regression. \emph{Technical Report}. Department of Mathematical Sciences, Brunel
University. URL: \url{http://bura.brunel.ac.uk/bitstream/2438/3593/1/fulltext.pdf}.

[6]  Yu, K. and Moyeed, R.A. (2001).  Bayesian Quantile Regression.  \emph{Statistics  &  Probability Letters}, 54, 437--447.
URL: \url{http://www.sciencedirect.com/science/article/pii/S0167715201001249}.


}



