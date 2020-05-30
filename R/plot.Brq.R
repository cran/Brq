plot.Brq <-
function (x, plottype = c("hist", "trace", "ACF", "traceACF", 
    "histACF", "tracehist", "traceACFhist"), Coefficients = 1, breaks = 30, 
    lwd = 1, col1 = 0, col2 = 1, col3 = 1, col4 = 1, ...) 
{   
      call <- match.call()
       mf <- match.call(expand.dots = FALSE)
       mf$drop.unused.levels <- FALSE
   
     Betas=as.matrix(x$beta[, Coefficients])
    k = ncol(as.matrix( Betas))
    if (k == 2) 
        par(mfrow = c(1, 2))
    if (k == 3) 
        par(mfrow = c(1, 3))
    if (k == 4) 
        par(mfrow = c(2, 2))
    if (k > 4 & k <= 12) 
        par(mfrow = c(ceiling(k/3), 3))
    if (k > 12) 
        par(mfrow = c(ceiling(k/3), 3))
    plottype <- match.arg(plottype)
    switch(plottype, trace = for (i in 1:k) {
        ts.plot(Betas[, i], xlab = "iterations", ylab = "", 
            main = noquote(names(coef(x)))[i], col = col4)
    }, ACF = for (i in 1:k) {
        acf(Betas[, i], main = noquote(names(coef(x)))[i], col = col3)
    }, traceACF = {
        par(mfrow = c(k, 2))
        for (i in 1:k) {
            ts.plot(Betas[, i], xlab = "iterations", ylab = "", 
                main = noquote(names(coef(x)))[i], col = col4)
            acf(Betas[, i], main = noquote(names(coef(x)))[i], 
                col = col3)
        }
    }, histACF = {
        par(mfrow = c(k, 2))
        for (i in 1:k) {
            hist(Betas[, i], breaks = breaks, prob = TRUE, main = "", 
                xlab = noquote(names(coef(x)))[i], col = col1)
            lines(density(Betas[, i], adjust = 2), lty = "dotted", 
                col = col2, lwd = lwd)
            acf(Betas[, i], main = noquote(names(coef(x)))[i], 
                col = col3)
        }
    }, tracehist = {
        par(mfrow = c(k, 2))
        for (i in 1:k) {
            ts.plot(Betas[, i], xlab = "iterations", ylab = "", 
                main = noquote(names(coef(x)))[i], col = col4)
            hist(Betas[, i], breaks = breaks, prob = TRUE, main = "", 
                xlab = noquote(names(coef(x)))[i], col = col1)
            lines(density(Betas[, i], adjust = 2), lty = "dotted", 
                col = col2, lwd = lwd)
        }
    }, traceACFhist = {
        par(mfrow = c(k, 3))
        for (i in 1:k) {
            ts.plot(Betas[, i], xlab = "iterations", ylab = "", 
                main = noquote(names(coef(x)))[i], col = col4)
            hist(Betas[, i], breaks = breaks, prob = TRUE, main = "", 
                xlab = noquote(names(coef(x)))[i], col = col1)
            lines(density(Betas[, i], adjust = 2), lty = "dotted", 
                col = col2, lwd = lwd)
            acf(Betas[, i], main = noquote(names(coef(x)))[i], 
                col = col3)
        }
    }, hist = for (i in 1:k) {
        hist(Betas[, i], breaks = breaks, prob = TRUE, main = "", 
            xlab = noquote(names(coef(x)))[i], col = col1)
        lines(density(Betas[, i], adjust = 2), lty = "dotted", 
            col = col2, lwd = lwd)
    })
}