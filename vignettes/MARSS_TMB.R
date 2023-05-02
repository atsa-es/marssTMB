## ----setup, include = FALSE---------------------------------------------------
# knitr options
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
run.comparisons <- FALSE
set.seed(1234)

## -----------------------------------------------------------------------------
library(MARSS)
data(lakeWAplankton, package = "MARSS")
phytoplankton <- c("Cryptomonas", "Diatoms", "Greens", "Unicells", "Other.algae")
dat <- as.data.frame(lakeWAplanktonTrans) |>
  subset(Year >= 1980 & Year <= 1989) |>
  subset(select=phytoplankton) |>
  t() |>
  MARSS::zscore()

## -----------------------------------------------------------------------------
mod.list <- list(R='unconstrained', m=1, tinitx=1)

## -----------------------------------------------------------------------------
m1 <- MARSS(dat, model=mod.list, form='dfa', z.score=FALSE, silent = TRUE)
m2 <- MARSS(dat, model=mod.list, form='dfa', z.score=FALSE, silent = TRUE, method="BFGS")

## -----------------------------------------------------------------------------
library(marssTMB)
m3 <- dfaTMB(dat, model=list(m=1, R='unconstrained'))
m4 <- MARSS_tmb(dat, model=mod.list)
m5 <- MARSS_tmb(dat, model=mod.list, control=list(fun.opt="optim"))

## ----echo=FALSE---------------------------------------------------------------
df <- data.frame(name=c("MARSS-EM", "MARSS-BFGS", "dfaTMB-nlminb", "MARSS_tmb-nlminb", "MARSS_tmb-optim-BFGS"),
                 logLik=c(m1$logLik, m2$logLik, m3$logLik, m4$logLik, m5$logLik))
knitr::kable(df)

## ----echo=FALSE---------------------------------------------------------------
library(tidyr)
library(ggplot2)

pars <- data.frame(
  name = c(paste0("R", rownames(coefficients(m1)$R)), paste0("Z", rownames(coefficients(m1)$Z))),
  EM = c(coefficients(m1)$R, coefficients(m1)$Z),
  BFGS = c(coefficients(m2)$R, coefficients(m2)$Z),
  TMB1 = c(as.vector(m3$Estimates$R[lower.tri(m3$Estimates$R,diag=TRUE)]), as.vector(m3$Estimates$Z)),
  TMB2 = c(coefficients(m4)$R, coefficients(m4)$Z))
pars <- pars %>% tidyr::pivot_longer(!name, names_to = "model")

dodge <- position_dodge(width=0.5)
ggplot(pars, aes(x=name, y=value, col=model)) +
  geom_point(position=dodge) +
  ggtitle("same estimates")

## ----include=FALSE------------------------------------------------------------
if(!run.comparisons){
  load("dfa-time-comparisons.rda")
}else{
df <- c()
mods <- list()
for(R in c("diagonal and equal", "diagonal and unequal", "unconstrained")){
  for(m in 1:3){
    mod <- list(m = m, R = R, tinitx=1)
    tfit <- system.time(fit <- MARSS(dat, model=mod, form="dfa", control=list(maxit=10000)))
    df <- rbind(df, data.frame(fun="MARSS", opt.function="EM", m=m, R=R, ncovar = 0, time=tfit[1], logLik=fit$logLik, convergence = fit$convergence))
    mods <- c(mods, list(fit))
    tfit <- system.time(fit <- MARSS(dat, model=mod, form="dfa", method="BFGS"))
    df <- rbind(df, data.frame(fun="MARSS", opt.function="BFGS", m=m, R=R, ncovar = 0, time=tfit[1], logLik=fit$logLik, convergence = fit$convergence))
    mods <- c(mods, list(fit))
    tfit <- system.time(fit <- MARSS_tmb(dat, model=mod, form="dfa"))
    df <- rbind(df, data.frame(fun="TMB", opt.function="nlminb", m=m, R=R, ncovar = 0, time=tfit[1], logLik=fit$logLik, convergence = fit$convergence))
    mods <- c(mods, list(fit))
  }
}
}

## -----------------------------------------------------------------------------
library(dplyr)
df2 <- df |> mutate(mod = paste0(ncovar, "-", fun, "-", opt.function))
ggplot(df2, aes(fill=mod, y=time, x=m)) + 
    geom_bar(position="dodge", stat="identity") +
  facet_wrap(~R, scales = "free_y") +
  scale_y_continuous() +
  ggtitle("TMB is faster esp for R unconstrained")

## -----------------------------------------------------------------------------
# use a simpler R
mod.list2 <- list(m=1, R='diagonal and unequal', tinitx=1)
# add a temperature covariate
temp <- as.data.frame(lakeWAplanktonTrans) |>
    subset(Year >= 1980 & Year <= 1989) |>
    subset(select=Temp)
covar <- t(temp)
t6 <- system.time(m6 <- MARSS_tmb(dat, model=mod.list2, form="dfa", covariates=covar, silent = TRUE))
t7 <- system.time(m7 <- MARSS(dat, model=mod.list2, form="dfa", covariates=covar, silent = TRUE))

## -----------------------------------------------------------------------------
TP <- as.data.frame(lakeWAplanktonTrans) |>
    subset(Year >= 1980 & Year <= 1989) |>
    subset(select=TP)
covar <- rbind(covar, t(TP))
t8 <- system.time(m8 <- MARSS_tmb(dat, model=mod.list2, form="dfa", covariates=covar, silent = TRUE))
t9 <- system.time(m9 <- MARSS(dat, model=mod.list2, form="dfa", covariates=covar, silent = TRUE))

## ----echo=FALSE---------------------------------------------------------------
df <- data.frame(name=rep(c("MARSS-EM", "MARSS_tmb-nlminb"), 2),
                 num_covar = rep(1:2, each = 2),
                 time=c(t6[1], t7[1], t8[1], t9[1]),
                 logLik=c(m6$logLik, m7$logLik, m8$logLik, m9$logLik))
knitr::kable(df)

## ----echo=FALSE---------------------------------------------------------------
library(tidyr)
n1 <- length(coef(m6, type="vector"))
n2 <- length(coef(m8, type="vector"))
pars <- data.frame(
  name = c(rep(names(coef(m6, type="vector")), 2), rep(names(coef(m8, type="vector")), 2)),
  ncovar = as.factor(c(rep(c(1, 1), each=n1), rep(c(2, 2), each = n2))),
  model = c(rep(c("EM", "TMB"), each=n1), rep(c("EM", "TMB"), each = n2)),
  value = c(coef(m6, type="vector"), coef(m7, type="vector"), coef(m8, type="vector"), coef(m9, type="vector")))

## ----echo=FALSE---------------------------------------------------------------
library(ggplot2)
dodge <- position_dodge(width=1)
ggplot(pars, aes(x=name, y=value, col=model)) +
  geom_point(position=dodge) +
  coord_flip() + 
  facet_wrap(~ncovar) +
  ggtitle("estimates should be the same (roughly)")

## ----include=FALSE------------------------------------------------------------
if(!run.comparisons){
  load("dfa-time-comparisons.rda")
}else{
for(R in c("diagonal and equal", "diagonal and unequal", "unconstrained")){
  for(m in 1:3){
    mod <- list(m = m, R = R, tinitx=1)
    tfit <- system.time(fit <- MARSS(dat, model=mod, form="dfa", covariates=covar, control=list(maxit=10000)))
    df <- rbind(df, data.frame(fun="MARSS", opt.function="EM", m=m, R=R, ncovar=2, time=tfit[1], logLik=fit$logLik, convergence = fit$convergence))
    mods <- c(mods, list(fit))
    tfit <- system.time(fit <- MARSS(dat, model=mod, form="dfa", covariates=covar, method="BFGS"))
    df <- rbind(df, data.frame(fun="MARSS", opt.function="BFGS", m=m, R=R, ncovar=2, time=tfit[1], logLik=fit$logLik, convergence = fit$convergence))
    mods <- c(mods, list(fit))
    tfit <- system.time(fit <- MARSS_tmb(dat, model=mod, form="dfa", covariates=covar))
    df <- rbind(df, data.frame(fun="TMB", opt.function="nlminb", m=m, R=R, ncovar=2, time=tfit[1], logLik=fit$logLik, convergence = fit$convergence))
    mods <- c(mods, list(fit))
  }
}
}

## -----------------------------------------------------------------------------
library(dplyr)
df2 <- df |> mutate(mod = paste0(ncovar, "-", fun, "-", opt.function))
ggplot(df2, aes(fill=mod, y=time, x=m)) + 
    geom_bar(position="dodge", stat="identity") +
  facet_wrap(~R, scales = "free_y") +
  scale_y_continuous() +
  ggtitle("TMB is faster than MARSS EM")

