summary.Brq <-
function(object, ...){
if(length(object$tau)>1){
p=dim(coef(object))[1]
result=array(,dim = c(p, 3, length(object$tau)))
estim=NULL
for (i in 1:length(object$tau)){
CredInt=apply(object$beta[,,i],2,quantile,c(0.025,0.975))

  TAB <- cbind( Coefficient = coef(object)[,i],
                L.CredIntv =  CredInt[1,],
                U.CredIntv = CredInt[2,])

  colnames(TAB) <- c("Estimate",  "L.CredIntv",  "U.CredIntv")
    colnames(result) <- c("Estimate",  "L.CredIntv",  "U.CredIntv")

  result[,,i] <- TAB
}
for(i in 1:length(object$tau)){
estim$tau=object$tau[i]
estim$coefficients=result[,,i]
print(estim)
}
}else{
  
 CredInt=apply(object$beta,2,quantile,c(0.025,0.975))

  TAB <- cbind( Coefficient = coef(object),
                L.CredIntv =  CredInt[1,],
                U.CredIntv = CredInt[2,])
  
  colnames(TAB) <- c("Estimate",  "L.CredIntv",  "U.CredIntv")
  
  result <- list(call=object$call,tau=object$tau,coefficients=TAB)  
  class(result) <- "summary.Brq"
  result
}
}