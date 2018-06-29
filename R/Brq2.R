BLqr <-
function(x, y, tau=0.5, runs=11000, burn=1000, thin=1) {

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
        Lambdadraw= matrix(nrow=runs, ncol=1)
        sigmadraw = matrix(nrow=runs, ncol=1)

      # Calculate some useful quantities
        xi     = (1 - 2*tau) 
        zeta   = tau*(1-tau)

      # Initial valus
        beta   = rep(1, p)
        s      = rep(1, p)
        v      = rep(1, n)
        Lambda2 = 1
        sigma  = 1

      # Hyperparameters
        a = 0.1
        b = 0.1

      # Draw from inverse Gaussian distribution
        rInvgauss <- function(n, mu, lambda = 1){
        un <- runif(n)
        Xi <- rchisq(n,1)
        f <- mu/(2*lambda)*(2*lambda+mu*Xi+sqrt(4*lambda*mu*Xi+mu^2*Xi^2))
        s <- mu^2/f
        ifelse(un < mu/(mu+s), s, f)}

      # Start the algorithm
        for (iter in 1: runs) {

      # Draw the latent variable v from inverse Gaussian distribution.
        lambda = 1/(2*sigma)
        mu     = 1/(abs(y - x%*%beta))
        v      = c(1/rInvgauss(n, mu = mu, lambda = lambda))
      
      # Draw the latent variable s from inverse Gaussian distribution.
        lambda= Lambda2
        mu    = sqrt(lambda/(beta^2/sigma) )
        s     =c(1/rInvgauss(p, mu = mu, lambda = lambda))

      # Draw sigma
        shape =   p/2 + 3/2*n 
        rate  = sum((y - x%*%beta - xi*v)^2 / (4*v) )+zeta*sum(v) + sum(beta^2/(2*s)) 
        sigma = 1/rgamma(1, shape= shape, rate= rate)
        
      # Draw beta
        V=diag(1/(2*v))
        invA <- chol2inv(chol(t(x)%*%V%*%x + diag(1/s)) )
        betam  <- invA%*%(t(x)%*%(V %*% (y-xi*v)))
        varcov=sigma*invA
        beta   <-betam+t(chol(varcov))%*%rnorm(p) 
               
      # Draw Lambda2
        tshape  = p + a 
        trate   = sum(s)/2 + b
        Lambda2 = rgamma(1, shape=tshape, rate=trate)

      # Sort beta and sigma
        betadraw[iter,]  = beta
        Lambdadraw[iter,]= Lambda2
        sigmadraw[iter,] = sigma
}
        coefficients =apply(as.matrix(betadraw[-(1:burn), ]),2,mean)
        names(coefficients)=colnames(x)
        if (all(x[,1]==1))  names(coefficients)[1]= "Intercept"  

        result <- list(beta = betadraw[seq(burn, runs, thin),],
        lambda = Lambdadraw[seq(burn, runs, thin),],
        sigma  <- sigmadraw[seq(burn, runs, thin),],
        coefficients=coefficients)
    
      return(result)
      class(result) <- "BLqr"
      result
}


