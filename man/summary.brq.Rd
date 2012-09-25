\name{summary.brq}
\alias{summary.brq}
\title{Posterior summary}
\description{This function displays the posterior mean,  standard error  and credible intervals for the regression coefficients.  }
\usage{\method{summary}{brq}(object, ...)}
\arguments{\item{object}{objects of class \code{brq}, originating from a call to \code{brq}}
\item{\dots}{additional arguments for the fitting the model}}
\examples{
# Example1
fit=brq(stack.loss ~ stack.x, tau=0.95, family="rq", runs=1000, burn=500)
summary.brq(fit)
}




