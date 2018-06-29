\name{Prostate}
\docType{data}
\alias{Prostate}
\title{Prostate Cancer Data}
\description{The data  is taken from a prostate cancer study (Stamey et al., 1989).
  It has 97 rows and 9 columns. The response variable is the logarithm of prostate-specific
antigen. The covariates  are eight clinical measures: the logarithm
of cancer volume (lcavol), the logarithm of prostate
weight (lweight), age, the logarithm of the amount of benign
prostatic hyperplasia (lbph), seminal vesicle invasion
(svi), the logarithm of capsular penetration (lcp), the Gleason
score (gleason) and the percentage Gleason score 4 or 5 (pgg45).
}
\usage{data("Prostate")}
\format{
  The data frame has the following components:
  \itemize{
    \item \bold{lcavol} : The logarithm of cancer volume.
    \item \bold{lweight} : The logarithm of prostate weight.
    \item \bold{age} : Age
    \item \bold{lbph} : The logarithm of the amount of benign prostatic hyperplasia.
    \item \bold{svi} : Seminal vesicle invasion.
    \item \bold{lcp} : The logarithm of capsular penetration.
    \item \bold{gleason} : Gleason score.
    \item \bold{pgg45} : Percentage Gleason scores 4 or 5.
    \item \bold{lpsa} : The logarithm of prostate-specific antigen.
  }
}
\source{
  Stamey, T. A., Kabalin, J. N., McNeal, J. E., Johnstone, I. M., Freiha, F., Redwine, E. A., and  Yang, N. (1989). Prostate specific antigen in the diagnosis and treatment of adenocarcinoma of the prostate. II. Radical prostatectomy treated patients. The Journal of urology, 141(5), 1076-1083.
}
\references{

  Stamey, T. A., Kabalin, J. N., McNeal, J. E., Johnstone, I. M., Freiha, F., Redwine, E. A., & Yang, N. (1989). Prostate specific antigen in the diagnosis and treatment of adenocarcinoma of the prostate. II. Radical prostatectomy treated patients. The Journal of urology, 141(5), 1076-1083.
}

\examples{
data(Prostate)
y=Prostate[,9]
x=as.matrix(Prostate[, -9])
y=y-mean(y)
x=scale(x)
fit=Brq(y~0+x, method="BLqr")
plot(fit, plottype="hist",D=1:8)
plot(fit, plottype="hist",D=1:8, col2=4,lwd=2)
plot(fit, plottype="trace",D=1:8, col4=4)
}

\keyword{datasets}



