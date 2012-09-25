summary.brq <-
function(object,...)
 { 
   tau=object$tau
   beta <- object$beta
   Mean=format(round(apply(beta, 2, mean), 8), quote = FALSE)
   Median=format(round(apply(beta, 2, median), 8), quote = FALSE)
   Sd=format(round(apply(beta, 2, sd), 8), quote = FALSE)
   Quantile=format(round(apply(beta,2,quantile,c(0.025,0.5,0.975)), 8), quote = FALSE)
   result <- data.frame(cbind(names(object$c), Mean,Sd,Median,Quantile[1,], Quantile[3,]))
   names(result)=c("Variable","Mean","Sd", "Median"," 2.5% C.I.","97.5% C.I.")

result
}