\name{print.summary.Brq}
\alias{print.summary.Brq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print the summary of Bayesian quantile regression methods.
}
\description{
  \code{print.Brq.summary} is an S3 method that prints the content of an S3 object of class \code{Brq.summary} to the console. 
}
\usage{
%print.summary.Brq(x, ...)
\method{print}{summary.Brq}(x, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
 \code{object of class Brq.} 
}
  \item{\dots}{
 \code{not used.} 
}
}

\value{
Print the summary of Bayesian quantile regression methods.
}

\author{
  Rahim Alhamzawi
}
\examples{
# Simulate data from heteroskedastic regression
set.seed(666)
n <- 100
x <- runif(n=n,min=0,max=10)
y <- 2 + 2*x + 0.6 *x*rnorm(n=n, mean=0, sd=1)

fit <- Brq(y ~ x, tau=0.5, method="Bqr") 

# Return Bayes estimates and credible intervals 
Results<- summary(fit)

# Print the Brq.summary object
Results
}
