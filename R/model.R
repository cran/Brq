model <-
function(object){
welcome<-function(){
    cat("=====  Model selection based on credible intervals ======")
    cat("\n")
    cat("#                                                       #")
    cat("\n")
    cat("#               Author: Rahim Alhamzawi                 #")
    cat("\n")
    cat("#               Contact: rahim.alhamzawi@qu.edu.iq      #")
    cat("\n")
    cat("#                      July, 2018                       #")
    cat("\n")
    cat("#                                                       #")
    cat("\n")
    cat("=========================================================")
    cat("\n")
}
##############################################################
result=NULL
if(length(object$tau)>1){
for(ii in 1:length(object$tau)){
CredInt = apply(object$beta[,,ii], 2, quantile, c(0.025, 0.975))
#Estimate= apply(object$beta[,,ii], 2, mean)
Estimate= object$coefficients[,ii]


for(i in 1:length(CredInt [1,])){
if (sign(CredInt [1,i])==-1 & sign (CredInt [2,i])==1)  Estimate [i]=0 
}
result= cbind(result,Estimate)}
}else{
CredInt = apply(object$beta, 2, quantile, c(0.025, 0.975))
Estimate= coef(object)

for(i in 1:length(CredInt [1,])){
if (sign(CredInt [1,i])==-1 & sign (CredInt [2,i])==1)  Estimate [i]=0 
}
result= cbind(Estimate)

}

welcome()
taulabs <- paste("tau=", format(round(object$tau, 3)))
dimnames(result) <- list(dimnames(object$beta)[[2]], taulabs)
rownames(result)=rownames(coef(object))
 result
}