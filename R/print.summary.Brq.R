print.summary.Brq <-
function(x, ...)
{
cat("Call:\n")
print(x$call)
cat("\n")
cat("tau:")
print(x$tau)
cat("\n")
print(x$coefficients)
}
