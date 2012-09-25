brq<-
function(formula, tau=0.5, family= "rq", Ce=0, runs=15000, burn=1000) {
   
    #x:    matrix of predictors.
    #y:    vector of dependent variable. 
    #tau:  quantile level.
    #runs: the length of the Markov chain.
    #burn: the length of burn-in.
    #Ce:   censoring  point
   
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

      # Censored data
      if(family=="Lcrq")     ym <- pmax(Ce,y)
      if(family=="Rcrq")     ym <- pmin(Ce,y)
      if(family=="bin")      ym <- y
      if(family=="bin")      Ce <- 0
      if(family=="rq")       ym <- y
      ind1  <- which(ym != Ce)
      ind0  <- which(ym ==Ce)
   
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
        Lambdadraw <- matrix(nrow=runs, ncol=p)

        if(all(x[,1]==1))  Xc=x[,-1] else{Xc=x}
        if(all(x[,1]!=1))  Cor.yx=as.vector(cor(y,Xc)) else{Cor.yx=c(1,cor(y,Xc))}
        Cor.yx= Cor.yx^2
        
      # parameters
        xi     <- (1 - 2*tau) 
        zeta   <- tau*(1-tau)

      # Initial valus
        beta   <- rep(1,p)
        v      <- rep(1,n)
        Lamb   <- rep(0.01,p)
        invLamb<- diag(Lamb)
        sigma  <- 1
        err    <- 1e-16
     
      # Start the algorithm
        for (iter in 1: runs) {

       # Draw the latent y
         XB     <- x%*%beta
         mean.y <- XB + xi*v
         sd.y   <- sqrt(as.vector(sigma*2*v))
         if(family=="bin"){ 
         if (!(length(ind1))==FALSE)
         y[ind1]=rtruncnorm(length(ind1), a=0, b=Inf, mean=mean.y[ind1],sd=sd.y[ind1])
         if (!(length(ind0))==FALSE)
         y[ind0]=rtruncnorm(length(ind0), a=-Inf, b=0, mean=mean.y[ind0],sd=sd.y[ind0])
         }
         if(family=="Rcrq"){
         if (!(length(ind0))==FALSE)
         y[ind0]=rtruncnorm(length(ind0),   a= Ce,  b=Inf, mean=mean.y[ind0],   sd=sd.y[ind0])}
         if(family=="Lcrq"){
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

      # sample Lamb from inverse Gaussian distribution.
        mu     <- 1/ abs(beta^2/Cor.yx + 1e-16)
        Ai     <- mu*rchisq(p,1)
        SAi    <- sqrt(2*Ai+Ai^2)
        L      <- 1/(mu*(1+Ai+SAi))
        S      <- mu^2*L
        Lamb   <- c(ifelse((mu+S)*runif(p) > mu, L, 1/S))

      # Draw sigma
        shape  <-  3*n/2 
        rate   <- sum( (y - XB - xi*v)^2 / (4*v) )+ zeta*sum(v)  + sum(Lamb)/4
        sigma  <- 1/rgamma(1, shape=shape, rate=1/rate)
        if(family=="bin") sigma=1
        if(family=="Lcrq") sigma=1
        if(family=="Rcrq") sigma=1
        
      # Draw  beta 
        invLamb          <- diag(1/(2*Cor.yx*Lamb))
        if(p==1) invLamb <- 1/(2*Cor.yx*Lamb)
        v      <- as.vector(v)
        V      <- diag(1/(2*sigma*v))
        varcov <- chol2inv(chol(t(x)%*%V%*%x + invLamb) )
        betam  <- varcov %*% (t(x)%*%(V %*% (y-xi*v)))
        beta   <- betam+t(chol(varcov))%*%rnorm(p)

        coefficients=beta
        id=which(abs(beta)/sqrt(Lamb)<=   0.1)
        coefficients[id]=0

      # Sort beta and  lambdad
          betadraw[iter,] = coefficients
          Lambdadraw[iter,]= Lamb

      }
      if(length(Cor.yx)==1)   coefficients =mean(betadraw) else{
      coefficients =apply( betadraw[-(1:burn), ],2,median)}
      names(coefficients)=colnames(x)
      result <- list(beta= betadraw[-(1:burn),],
      Lambda  = Lambdadraw[-(1:burn),], coefficients=coefficients, tau=tau)

      return(result)
  class(result) <- "brq"
    result
}

