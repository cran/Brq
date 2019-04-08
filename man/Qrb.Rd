\name{BBqr}
\alias{BBqr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Bayesian Binary Quantile Regression}
\description{This function implements the idea of Bayesian Binary quantile regression employing
 a likelihood function that is based on the asymmetric Laplace distribution. }
\usage{BBqr(x,y, tau =0.5, runs =11000, burn =1000, thin=1)}
\arguments{
   \item{x}{Matrix of predictors.}
   \item{y}{Vector of dependent variable. }
  \item{tau}{The quantile of interest. Must be between 0 and 1.}
    \item{runs}{Length of desired Gibbs sampler output.}
  \item{burn}{Number of Gibbs sampler iterations before output is saved.}
\item{thin}{thinning parameter of MCMC draws.}
}

\author{
  Rahim Alhamzawi
}