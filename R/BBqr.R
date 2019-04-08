BBqr <-
function(x,y,tau=0.5, runs=11000, burn=1000, thin=1) {
    #x:    matrix of predictors.
    #y:    vector of dependent variable. 
    #tau:  quantile level.
    #runs: the length of the Markov chain.
    #burn: the length of burn-in.
    #thin: thinning parameter of MCMC draws
    x <- as.matrix(x)  
    if(ncol(x)==1) {x=x} else {
    x=x
    if (all(x[,2]==1)) x=x[,-2] }

      # Calculate some useful quantities
        n <- nrow(x)
        p <- ncol(x)
   
      # check input
        if (tau<=0 || tau>=1) stop ("invalid tau:  tau should be >= 0 and <= 1. 
               \nPlease respecify tau and call again.\n")
        if(n != length(y)) stop("length(y) not equal to nrow(x)")
        if(n == 0) return(list(coefficients=numeric(0),fitted.values=numeric(0),
              deviance=numeric(0)))
        if(!(all(is.finite(y)) || all(is.finite(x)))) stop(" All values must  be 
              finite and non-missing") 

      # Saving output matrices 
        betadraw  = matrix(nrow=runs, ncol=p)
        sigmadraw = matrix(nrow=runs, ncol=1)

      # Calculate some useful quantities
        xi     = (1 - 2*tau) 
        zeta   = tau*(1-tau)

      # Initial valus
        beta   = rep(0.99, p)
        v      = rep(1, n)
        sigma  = 1
        lambda= 0.05
        Lambda=diag(lambda, p)
        ystar<-y-0.5

      # low and upp
        low<-ifelse(y==1,0,-Inf)
        upp<-ifelse(y==1,Inf,0)

      # Draw from inverse Gaussian distribution
        rInvgauss <- function(n, mu, lambda = 1){
        un <- runif(n)
        Xi <- rchisq(n,1)
        f <- mu/(2*lambda)*(2*lambda+mu*Xi+sqrt(4*lambda*mu*Xi+mu^2*Xi^2))
        s <- mu^2/f
        ifelse(un < mu/(mu+s), s, f)}

      # Draw from a truncated normal distribution
        rtnorm<-function(n,mean=0,sd=1,lower.bound=-Inf,upper.bound=Inf){
        lower<-pnorm(lower.bound,mean,sd)
        upper<-pnorm(upper.bound,mean,sd)
        qnorm(runif(n,lower,upper),mean,sd)}

      # Start the algorithm
        for (iter in 1: runs) {

      # Draw the latent variable v from inverse Gaussian distribution.
        lambda = 1/(2*sigma)
        mu     = 1/(abs(ystar - x%*%beta))
        v      = c(1/rInvgauss(n, mu = mu, lambda = lambda))
      
      # Draw sigma
        shape =   3/2*n 
        rate  = sum((ystar - x%*%beta - xi*v)^2/(4*v))+zeta*sum(v) 
        sigma = 1/rgamma(1, shape= shape, rate= rate)

      # Draw beta
        V=diag(1/(2*sigma*v))
        varcov <- chol2inv(chol(t(x)%*%V%*%x + Lambda))
        betam  <- varcov %*% (t(x)%*%(V %*% (ystar-xi*v)))
        beta   <-betam+t(chol(varcov))%*%rnorm(p)

      # Draw ystar
        mu<-x%*%beta+xi*v
        sd<- sqrt(2*sigma*v)
        mu[which(y==1 & mu<0)]=0
        mu[which(y==0 & mu>0)]=0
        ystar<-rtnorm(n,mean=mu,sd=sd,lower.bound=low,upper.bound=upp)
        ystar= ystar/sd

      # Sort beta and sigma
        betadraw[iter,]  = beta
        sigmadraw[iter,] = sigma
}
        coefficients =apply(as.matrix(betadraw[-(1:burn), ]),2,mean)
        names(coefficients)=colnames(x)
        if (all(x[,1]==1))  names(coefficients)[1]= "Intercept"  

        result <- list(beta = betadraw[seq(burn, runs, thin),],
        sigma  = sigmadraw[seq(burn, runs, thin),],
        coefficients=coefficients)
    
      return(result)
      class(result) <- "BBqr"
      result
}


