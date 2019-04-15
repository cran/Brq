DIC=function(object){
# Estimate Deviance Information Criterion (DIC)
#
# References:
#   Bayesian Data Analysis.
#   Gelman, A., Carlin, J., Stern, H., and Rubin D.
#   Second Edition, 2003
llSum = 0
y=object$y
N=dim(object$MuY)[1]
PostMu=apply(object$MuY, 2, mean)
PostVar=apply(object$VarY, 2, mean)

L=sum( dnorm(y, PostMu, PostVar, log=TRUE) )

  for (i in 1:N) {
       m=object$MuY[i, ]
       s=sqrt(object$VarY[i, ])
       llSum = llSum + sum( dnorm(y, m, s, log=TRUE) )
       }

  P   = 2 * (L - (1 / N * llSum))
  dic = -2 * (L - P)
return(dic)
      class(DIC) <- "Brq"

}

