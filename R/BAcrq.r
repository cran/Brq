

BAcrq <-
function(formula, tau = 0.5, method = "Lcrq", Ce = 0, runs = 15000, burn = 1000) {

    #x:    matrix of predictors.
    #y:    vector of dependent variable. 
    #tau:  quantile level.
    #runs: the length of the Markov chain.
    #burn: the length of burn-in.
  
    x=formula[3][[1]]
    y=formula[2][[1]]
    call <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")

    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf, contrasts)
    x <- as.matrix(x)  
    if(ncol(x)==1) {x=x} else {
    x=x
    if (all(x[,2]==1)) x=x[,-2] }

      # Calculate some useful quantities
       n <- nrow(x)
       p <- ncol(x)
       if (method == "Lcrq") ym <- pmax(Ce, y)
       if (method == "Rcrq") ym <- pmin(Ce, y)
       ind1 = which(ym != Ce)
       ind0 = which(ym == Ce)

      # check input
        if (tau<=0 || tau>=1) stop ("invalid tau:  tau should be >= 0 and <= 1. 
               \nPlease respecify tau and call again.\n")
        if(n != length(y)) stop("length(y) not equal to nrow(x)")
        if(n == 0) return(list(coefficients=numeric(0),fitted.values=numeric(0),
              deviance=numeric(0)))
        if(!(all(is.finite(y)) || all(is.finite(x)))) stop(" All values must  be 
              finite and non-missing") 
      # Saving output matrices 
        betadraw  <- matrix(nrow=runs, ncol=p)
        sigmadraw <- matrix(nrow=runs, ncol=1)
    
        if(all(x[,1]==1))  Xc=x[,-1] else{Xc=x}
        if(all(x[,1]!=1))  Cor.yx=as.vector(cor(y,Xc)) else{Cor.yx=c(1,cor(y,Xc))}

      # parameters
        xi   <- (1 - 2*tau) 
        zeta <- tau*(1-tau)

      # Initial valus
        beta  <- rep(1,p)
        s     <- rep(1,p)
        v     <- rep(1,n)
        sigma <- 1
        err   <- 1e-16

        # Start the algorithm
        for (iter in 1: runs) {

        # Draw the latent y
         XB=x%*%beta
         mean.y <- XB+xi*v
         sd.y   <- sqrt(as.vector(sigma*2*v))
         if(method=="Rcrq"){
         if (!(length(ind0))==FALSE)
         y[ind0]=rtruncnorm(length(ind0),   a= Ce,  b=Inf, mean=mean.y[ind0],   sd=sd.y[ind0])}
         if(method=="Lcrq"){
         if (!(length(ind0))==FALSE)
         y[ind0]=rtruncnorm(length(ind0),   a=-Inf, b=Ce,  mean=mean.y[ind0],   sd=sd.y[ind0])}

      # Draw the latent v from inverse Gaussian distribution.
      #Chhikara, R. and Folks, J. (1989). "The inverse Gaussian distribution".
        mu     <- abs(1 /((y-XB)+ err))
        Ai     <- mu*rchisq(n,1)
        SAi    <- sqrt(2/sigma*Ai+Ai^2)
        L      <- 1/(mu*sigma*(1/sigma+Ai+SAi))
        S      <- mu^2*L
        v      <- c(ifelse((mu+S)*runif(n) > mu, L, 1/S))


       # Draw lambda2
        tshape = 1 + 1e-3
        trate = s/(2 * sigma) +  1e-3
        lambda2 <- 1/rgamma(p, shape = tshape, rate = trate)

        # Draw s
        lambda <- 1/(sigma*lambda2)
        mu <- sqrt(1/(sigma*lambda2*beta^2 + 1e-12)) 
        Ai     <- mu*rchisq(p,1)
        SAi    <- sqrt(4*lambda*Ai+Ai^2)
        Lambda <- 2*lambda
        L      <- 1/(mu/Lambda*(Lambda + Ai +SAi))
        S      <- mu^2*L
        s      <- c(ifelse((mu+S)*runif(p) > mu, L, 1/S))


      # Draw sigma
      #  shape <-  3/2*n
      #  rate  <- (1/4)*sum( (y - XB - xi*v)^2 / (v) )+ zeta*sum(v) + sum(s/(2*lambda2))
      #  sigma <- 1/rgamma(1, shape=shape, rate=1/rate)
        sigma=1

      # Draw  beta
        S=diag(1/s)  
        if(ncol(x)==1) S=1/s
        v=as.vector(v)
        V=diag(1/(2*sigma*v))
        varcov <- chol2inv(chol(t(x)%*%V%*%x + S) )
        betam  <- varcov %*% t(x)%*%V %*% (y-xi*v)
        beta   <-betam+t(chol(varcov))%*%rnorm(p)

        coefficients=beta
        id=which(abs(beta*Cor.yx)/sqrt(1/s)<=   0.1) 
        coefficients[id]=0

      # Sort beta and sigma
          betadraw[iter,] = coefficients
          sigmadraw[iter,]= sigma
      }
      if(ncol(x)==1) coefficients =mean(betadraw) else{
      coefficients =apply( betadraw[-(1:burn), ],2,median)}
      names(coefficients)=colnames(x)
      result <- list(beta = betadraw[-(1:burn),],
      sigma  = sigmadraw[-(1:burn),], coefficients=coefficients,tau=tau)

      return(result)
      class(result) <- "BAcrq"
      result
}