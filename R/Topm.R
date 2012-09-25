Topm <-
function(formula,tau=0.5, family= "rq", max.steps = 1000, C=0.5){

   x=formula[3][[1]]
   y=formula[2][[1]]
   #family <- match.arg(family)
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
   if(max.steps>=5000) max.steps= 5000
   n=nrow(x)
   p=ncol(x)
   BetaOld=NULL
   
   BetaOld=rq(formula,tau=tau)$coeff  
   if (family=="Lcrq") BetaOld=crq(Curv(y,yc) ~ 0+x, taus = tau, method = "Pow")$coeff    
   if (family=="Lm") BetaOld=lm(formula)$coeff
   beta2=BetaOld
   err= 1e-6
   sigma=1e-1
   Index=rep(0,p)

# The correlation coefficient of x and y
   if(all(x[,1]==1))  Xc=x[,-1] else{Xc=x}
   Cor.yx=abs(cor(y,Xc))

# sd(xj)/sd(y)
   ST=apply(Xc,2,sd)/sd(c(y))
   if(all(x[,1]==1)) Cor.yx=c(1,Cor.yx)
   if(all(x[,1]==1)) ST=c(1,ST)
                        
# Saving output matrices 
   betadraw    <- matrix(nrow=max.steps, ncol=p)

# Start the algorithm
for(iter in 1:max.steps){

   if (iter/30 == as.integer(iter/30)) {beta2=BetaOld}
    
   mu=sqrt((1/(2*sigma))/(BetaOld^2 + err)/2)
   Ai     <- mu*rchisq(p,1)
   SAi    <- sqrt(2/sigma*Ai+Ai^2)
   L      <- 1/(mu*(1+ sigma*(Ai + SAi)))
   S      <- mu^2*L
   Lamb   <- c(ifelse((mu+S)*runif(p) > mu, L, 1/S))
 
   Index =abs(BetaOld*Cor.yx)/sqrt(2*Lamb) + BetaOld^2* ST + abs(BetaOld)/sqrt(2*Lamb)
       
   Zero=which(Index<= C)
   beta2[Zero]=0
   if(sum(beta2)==0) beta2=BetaOld
   betadraw[iter,]=beta2

}
 
    Beta2 = apply(betadraw, 2, median)
    NewIndex0 = which(Beta2 == 0)
    NewIndex1 = which(Beta2 != 0)
    beta2[NewIndex0] = 0
    x = x[, NewIndex1]
    betanew = rq(y ~ 0 + x, tau = tau)$coeff
    if (family=="Lcrq") betanew=crq(Curv(y,yc) ~ 0+ x, taus = tau, method = "Pow")$coeff
    if (family == "Lm")  betanew = lm(y ~ 0 + x)$coeff
    beta2[NewIndex1] = betanew
    names(beta2) = names(BetaOld)
    result <- list(coeff = beta2)
    return(result)

}

