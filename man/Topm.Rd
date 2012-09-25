\name{Topm}
\alias{Topm}
\title{Model selection}
\description{Variable selection in classical and quantile regression models using the usual classical estimator and a Bayesian information criterion.}
\usage{Topm(formula,tau=0.5, family= "rq",max.steps =1000, C=0.5)}
\arguments{
  \item{formula}{\code{Model formula.} }
  \item{tau}{ \code{Quantile level.} }
  \item{family}{ \code{The fitting method: either "rq"  for quantile regression with continuous response, "Lm" for classical linear regression, or  "Lcrq" for Tobit quantile regression with left censored point.} }
  \item{max.steps}{ \code{Number of the iteration to update the threshold value.} }
  \item{C}{ \code{Critical point.} }}
%\details{}
%\references{}
%\author{}
%\note{}
%\seealso{}
\examples{
# Example 1
set.seed(123456)
Beta <- c(3, 1.5, 0, 0, 2, 0, 0, 0)
#Beta <- c(5, 0, 0, 0, 0, 0, 0, 0)
n <- 100
p <- length(Beta)
S <- matrix(NA,p,p)

for (i in 1:p)
  for (j in 1:p)
     S[i,j] <- 0.5^(abs(i-j))

x <- matrix(rnorm(n*p),n,p) \%*\% t(chol(S))
y <- x \%*\% Beta + rnorm(n)


# Mean regression
Topm(y~0+x, family ="Lm",max.steps =1000)$coeff

# Quantile regression
Topm(y~0+x,tau=0.5, family ="rq",max.steps =1000)$coeff

# Tobit quantile regression
y=pmax(0,y)
yc=rep(0,n)
Topm(y~0+x,tau=0.5, family ="Lcrq",max.steps =1000)$coeff

# Example 2

data(Body)
attach(Body)
Topm(Weight~BiacSk+BiilSk+BitrSk+CheDeSk+CheDiSk+ElbowSk+WristSk+KneeSk+AnkleSk 
      ,tau=0.5, family ="rq",max.steps= 5000)$coeff
}
\references{
[1] Alhamzawi, R.  (2012).  Bayesian variable selection in quantile regression models.

[2] Koenker, R. (2012).  quantreg: Quantile regression. ver 4.79 .
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }
