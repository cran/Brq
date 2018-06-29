summary.Brq <-
function(object, ...){
  
 CredInt=apply(object$beta,2,quantile,c(0.025,0.975))

  TAB <- cbind( Coefficient = coef(object),
                L.CredIntv =  CredInt[1,],
                U.CredIntv = CredInt[2,])
  
  colnames(TAB) <- c("Estimate",  "L.CredIntv",  "U.CredIntv")
  
  result <- list(call=object$call,tau=object$tau,coefficients=TAB)  
  class(result) <- "summary.Brq"
  result
}
