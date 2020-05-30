\name{ImmunogG}
\alias{ImmunogG}
\docType{data}
\title{Immunoglobulin G Data}
\description{ This data set comprises the serum concentration (grams per litre) of immunoglobulin G in 298 children aged from 6 months to 6 years and is discussed in detail by Isaacs 
et al. (1983) and used in Keming et al. (2003).
}
\usage{data("ImmunogG")}
\format{
  A data frame with 298 observations on the following 2 variables.
  \describe{
    \item{\code{IgG}}{The serum concentration (grams per litre) of immunoglobulin G}
    \item{\code{Age}}{Age of children in years.}
  }
}
\source{
[1] Isaacs, D., Altman, D. G., Tidmarsh, C. E., Valman, H. B. and Webster, A. D. B. (1983) Serum Immunoglobin concentrations
 in preschool children measured by laser nephelometry: reference ranges for IgG, IgA and IgM.  J. Clin. Path., 36, 1193-1196.
}
\references{
[1] Yu, K., Lu, Z., & Stander, J. (2003). Quantile regression: applications and current research areas. Journal of the Royal Statistical Society:Series D (The Statistician), 52(3), 331-350.
}
\examples{
data(ImmunogG)
attach(ImmunogG)
y=IgG
x=Age
X=cbind(x, x^2)
plot(Age,IgG)
plot(Age,IgG, main="Quantile Regression Fit", cex=.5, col="gray")
for (i in 1:5) {
if (i==1) p = .05
if (i==2) p = .25
if (i==3) p = .50
if (i==4) p = .75
if (i==5) p = .95
fit = Brq(y~X,tau=p,runs=400, burn=100)
# Note: runs =11000 and burn =1000
curve(fit$coef[1]+fit$coef[2]*x+fit$coef[3]*x^2,col=i,add=TRUE)
}


}

\keyword{ImmunogG}