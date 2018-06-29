Brq.formula <-
function(formula, data=list(), ...)
{
mf <- model.frame(formula=formula, data=data)
x <- model.matrix(attr(mf, "terms"), data=mf)
y <- model.response(mf)
est <- Brq.default(x, y, ...)
est$call <- match.call()
est$formula <- formula
est
}
