\name{summary.Brq}
\alias{summary.Brq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize the output of the \code{Brq} function}
\description{
\code{summary.Brq} is an S3 method that summarizes the output of the \code{Brq} function in an object of class summary.Brq.
  This function produces the Bayesian estimates and credible intervals for the regression coefficients. 
}
\usage{
%summary.Brq(object, ...)
\method{summary}{Brq}(object, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
 \code{an object of class Brq.} 
}
  \item{\dots}{
 \code{not used.} 
}
}

\examples{
# Example 

fit = Brq(stack.loss~stack.x,tau=0.5, method= "Bqr", runs=5000, burn=1000)
summary(fit)
}

\author{
  Rahim Alhamzawi
}

\references{

[1]  Alhamzawi, R. (2018). Brq: An R package for Bayesian Quantile Regression. Working paper.

}