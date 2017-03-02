\name{Brq}
\alias{Brq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian Quantile Regression}
\description{This function implements the idea of Bayesian  estimation in  quantile regression models employing
 a likelihood function that is based on the asymmetric Laplace distribution (Yu and  Moyeed, 2001). The asymmetric
 Laplace  error distribution is written as  scale mixtures of normal distributions as in Kozumi and Kobayashi (2011).
 In this function,   a two-level hierarchical Bayesian model is used. Specifically, I put   zero mean Gaussian priors on the regression coefficients with 
 non informative Jeffreys prior distributions for the unknown variances. }
\usage{Brq(formula, tau =0.5, runs =11000, burn =1000)}
\arguments{
   \item{formula}{Model formula.}
  \item{tau}{The quantile of interest. Must be between 0 and 1.}
    \item{runs}{Length of desired Gibbs sampler output.}
  \item{burn}{Number of Gibbs sampler iterations before output is saved.}
}
\examples{
# Example 1
n <- 100
x <- runif(n=n,min=0,max=5)
y <- 1 + 1.5*x+ .5*x*rnorm(n)
plot(x,y, main="Scatterplot and Quantile Regression Fit", xlab="x", cex=.5, col="gray")
for (i in 1:5) {
if (i==1) p = .05
if (i==2) p = .25
if (i==3) p = .50
if (i==4) p = .75
if (i==5) p = .95
fit = Brq(y~x,tau=p,runs=2000, burn=500)
# Note: runs =11000 and burn =1000
abline(a=mean(fit$c[1]),b=mean(fit$c[2]),lty=i,col=i)
}
abline( lm(y~x),lty=1,lwd=2,col=6)
legend(x=-0.30,y=max(y)+0.5,legend=c(.05,.25,.50,.75,.95,"OLS"),lty=c(1,2,3,4,5,1),
lwd=c(1,1,1,1,1,2),col=c(1:6),title="Quantile")


# Example 2
#data(salary)
#attach(salary)
# y=b0+b1x1+b2x2
#plot(Years,  Salaries,cex=0.5, ylab="Salaries", xlab="Years",main="",col="gray")
#for(p in c(0.05, 0.25, 0.50, 0.75, 0.95)){
#xseq  <-  seq(min(Years), max(Years), len=5000)
#fit <- Brq( Salaries~Years+I(Years^2),tau=p, runs=5000, burn=1000)
#lines(xseq, fit$coef \%*\% rbind(1, xseq,xseq^2),lty=1, lwd=1.5, col="blue")
#}
}
\author{
  Rahim Alhamzawi
}
\references{

[1]Alhamzawi, R. (2014). Model selection in quantile regression models. Journal of Applied Statistics, 42, 445-458.

[2]  Kozumi, H. and  Kobayashi, G. (2011). Gibbs sampling methods for Bayesian quantile regression.  \emph{Journal of Statistical Computation and
Simulation}. {bf  81},  1565-1578. URL: \url{http://www.tandfonline.com/doi/abs/10.1080/00949655.2010.496117}.

[6]  Yu, K. and Moyeed, R.A. (2001).  Bayesian Quantile Regression.  \emph{Statistics  &  Probability Letters}, 54, 437--447.
URL: \url{http://www.sciencedirect.com/science/article/pii/S0167715201001249}.


}



