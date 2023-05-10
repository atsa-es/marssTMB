dat <- t(harborSealWA)[2:4, ]
fit <- MARSS(dat, model=list(Q="unconstrained"))
fitchol <- var_to_cholvar(fit)
Qchol = coef(fitchol, type="matrix")$Q
Q = coef(fit, type="matrix")$Q
Q
crossprod(Qchol)
Q = coef(fit, type="matrix", what="start")$Q
Qchol = coef(fitchol, type="matrix", what="start")$Q
Q
crossprod(Qchol)

fit <- MARSS(dat, model=list(Q=diag(0.1,3)+0.01), control=list(maxit=15))
fitchol <- var_to_cholvar(fit)
Qchol = coef(fitchol, type="matrix")$Q
Q = coef(fit, type="matrix")$Q
Q
crossprod(Qchol)
Q = coef(fit, type="matrix", what="start")$Q
Qchol = coef(fitchol, type="matrix", what="start")$Q
Q
crossprod(Qchol)