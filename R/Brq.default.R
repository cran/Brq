Brq.default <-
function(x, y, tau=0.5, method=c("Bqr","BLqr","BALqr","Btqr","BLtqr","BALtqr"), left=0, runs=2000, burn=1000, thin=1, ...)
{
x <- as.matrix(x)
y <- as.numeric(y)
 method <- match.arg(method)
 est= switch(method,
         Bqr   = Bqr(x,y,tau=tau, runs=runs, burn=burn, thin=thin),
         BLqr  = BLqr(x,y,tau=tau, runs=runs, burn=burn, thin=thin),
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

