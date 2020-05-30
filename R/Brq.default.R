Brq.default <-
function(x, y, tau=0.5, method=c("Bqr","BBqr","BLqr","BLBqr","BALqr","Btqr","BLtqr","BALtqr"), left=0, runs= 5000, burn= 1000, thin=1, ...)
{   
set.seed(123456) 
x <- as.matrix(x)
y <- as.numeric(y)
n=dim(x)[1]
p=dim(x)[2]
Coeff=NULL
esti=NULL
method <- match.arg(method,c("Bqr","BBqr","BLqr","BLBqr","BALqr","Btqr","BLtqr","BALtqr"))
Coeff=NULL
Betas <- array(,dim = c(length(seq(burn, runs, thin)), p, length(tau)))
if(length(tau)>1){
for (i in 1:length(tau)){
est= switch(method,
         Bqr   = Bqr(x,y,tau=tau[i], runs=runs, burn=burn, thin=thin),
         BBqr   = BBqr(x,y,tau=tau[i], runs=runs, burn=burn, thin=thin),
         BLqr  = BLqr(x,y,tau=tau[i], runs=runs, burn=burn, thin=thin),
         BLBqr  = BLqr(x,y,tau=tau[i], runs=runs, burn=burn, thin=thin),
         BALqr = BALqr(x,y,tau=tau[i], runs=runs, burn=burn, thin=thin),
         Btqr  = Btqr(x,y,tau=tau[i],left = 0,runs=runs, burn=burn, thin=thin),
         BLtqr = BLtqr(x,y,tau=tau[i], left = 0, runs=runs, burn=burn, thin=thin),
         BALtqr= BALtqr(x,y,tau=tau[i], left = 0, runs=runs, burn=burn, thin=thin))
Betas[,,i]=est$beta
result=est$coefficients
Coeff= cbind(Coeff,result)


}
paste("tau=", format(round(tau, 3)))
taulabs <- paste("tau=", format(round(tau, 3)))
dimnames(Coeff) <- list(dimnames(x)[[2]], taulabs)
esti$beta=Betas
esti$tau <- tau
esti$coefficients=Coeff
esti$call <- match.call()
class(esti) <- "Brq"
esti
} else {
  est= switch(method,
         Bqr   = Bqr(x,y,tau=tau, runs=runs, burn=burn, thin=thin),
         BBqr   = BBqr(x,y,tau=tau, runs=runs, burn=burn, thin=thin),
         BLqr  = BLqr(x,y,tau=tau, runs=runs, burn=burn, thin=thin),
         BLBqr  = BLqr(x,y,tau=tau, runs=runs, burn=burn, thin=thin),
         BALqr = BALqr(x,y,tau=tau, runs=runs, burn=burn, thin=thin),
         Btqr  = Btqr(x,y,tau=tau,left = 0,runs=runs, burn=burn, thin=thin),
         BLtqr = BLtqr(x,y,tau=tau, left = 0, runs=runs, burn=burn, thin=thin),
         BALtqr= BALtqr(x,y,tau=tau, left = 0, runs=runs, burn=burn, thin=thin))

est$tau <- tau
est$fitted.values <- as.vector(x %*% est$coefficients)
est$residuals <- y - est$fitted.values
est$call <- match.call()
class(est) <- "Brq"
est
}
}