## ----setup, include = FALSE---------------------------------------------------
# knitr options
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
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
m1.em <- MARSS(dat, model=list(R='unconstrained', m=1, tinitx=1), form='dfa', z.score=FALSE, silent = TRUE)

## -----------------------------------------------------------------------------
library(marssTMB)
m1.tmb <- dfaTMB(dat, model=list(m=1, R='unconstrained'))

## -----------------------------------------------------------------------------
library(tidyr)
pars <- data.frame(name = c(paste0("R", rownames(coefficients(m1.em)$R)),
                            paste0("Z", rownames(coefficients(m1.em)$Z))),
                   EM = c(coefficients(m1.em)$R, coefficients(m1.em)$Z),
  TMB = c(as.vector(m1.tmb$Estimates$R[lower.tri(m1.tmb$Estimates$R,diag=TRUE)]),
  as.vector(m1.tmb$Estimates$Z)))
pars <- pars %>% tidyr::pivot_longer(2:3, names_to = "model")

## -----------------------------------------------------------------------------
library(ggplot2)
dodge <- position_dodge(width=0.5)
ggplot(pars, aes(x=name, y=value, col=model)) +
  geom_point(position=dodge) +
  ggtitle("same estimates")

## -----------------------------------------------------------------------------
m1.bfgs <- MARSS(dat, model=list(R='unconstrained', m=2, tinitx=1), form='dfa', z.score=FALSE, silent = TRUE,
               method="BFGS")
m1.tmb <- dfaTMB(dat, model=list(m=2, R='unconstrained'))
c(m1.bfgs$logLik, m1.tmb$logLik)

## ----echo=FALSE---------------------------------------------------------------
ZmatGen <- function(n, m) {
  tZ <- matrix(0.5, nrow = n, ncol = m)
  if (m > 1) {
    for (i in 1:(m - 1)) {
      tZ[i, (m - (m - 2) + (i - 1)):m] <- rep(0, m - 1 - (i - 1))
    }
    return(tZ)
  } else {
    return(tZ)
  }
}
pars <- data.frame(name = c(paste0("R", rownames(coefficients(m1.bfgs)$R)),
                            paste0("Z", rownames(coefficients(m1.bfgs)$Z))),
                   EM = c(coefficients(m1.bfgs)$R, coefficients(m1.bfgs)$Z),
  TMB = c(as.vector(m1.tmb$Estimates$R[lower.tri(m1.tmb$Estimates$R,diag=TRUE)]),
  as.vector(m1.tmb$Estimates$Z[ZmatGen(nrow(dat), 2)!=0])))
pars <- pars %>% tidyr::pivot_longer(2:3, names_to = "model")
dodge <- position_dodge(width=0.5)
ggplot(pars, aes(x=name, y=value, col=model)) +
  geom_point(position=dodge) +
  ggtitle("same estimates")

## -----------------------------------------------------------------------------
mod.list <- list(R='diagonal and equal', m=2, tinitx=1)
m1.bfgs <- MARSS(dat, model=mod.list, form='dfa', z.score=FALSE, silent = TRUE,
               method="BFGS")
m1.tmb <- dfaTMB(dat, model=mod.list)
c(m1.bfgs$logLik, m1.tmb$logLik)

## -----------------------------------------------------------------------------
data(lakeWAplankton, package = "MARSS")
phytoplankton <- c("Cryptomonas", "Diatoms", "Greens", "Unicells", "Other.algae")
dat <- as.data.frame(lakeWAplanktonTrans) |>
  subset(Year >= 1980 & Year <= 1989) |>
  subset(select=phytoplankton) |>
  t() |>
  MARSS::zscore()

# add a temperature covariate
temp <- as.data.frame(lakeWAplanktonTrans) |>
    subset(Year >= 1980 & Year <= 1989) |>
    subset(select=Temp)
covar <- t(temp)
m_cov_tmb <- dfaTMB(dat, model=list(m=1, R='diagonal and unequal'), 
                    EstCovar = TRUE, Covars = covar)
m_cov_tmb$Estimates$D

# add a 2nd covariate
TP <- as.data.frame(lakeWAplanktonTrans) |>
    subset(Year >= 1980 & Year <= 1989) |>
    subset(select=TP)
covar <- rbind(covar, t(TP))
m_cov2_tmb <- dfaTMB(dat, model=list(m=1, R='diagonal and unequal'), 
                    EstCovar = TRUE, Covars = covar)
m_cov2_tmb$Estimates$D

## -----------------------------------------------------------------------------
dat2 <- rbind(dat, dat+matrix(rnorm(nrow(dat)*ncol(dat),0,0.2), nrow(dat), ncol(dat)))
dat2 <- MARSS::zscore(dat2)
mod.list <- list(R='diagonal and unequal', m=3, tinitx=1)
# Control set to match setting in MARSS.tmp()
t1.bfgs <- system.time(m1.bfgs <- MARSS(dat2, 
                                        model=mod.list, form='dfa', 
                                        z.score=FALSE, silent = TRUE, method="BFGS", 
                                        control = list(reltol = 1e-12, maxit=2000)))[1]
t1.em <- system.time(m1.em <- MARSS(dat2, model=mod.list, form='dfa', z.score=FALSE, silent = TRUE))[1]
t1.tmb <- system.time(m1.tmb <- dfaTMB(dat2, model=mod.list))[1]
t1.tmb2 <- system.time(m1.tmb2 <- dfaTMB(dat2, model=mod.list, fun.opt="optim"))[1]

## -----------------------------------------------------------------------------
# Log-likelihood
c(em=m1.em$logLik, kfoptim=m1.bfgs$logLik, TMBnlminb=m1.tmb$logLik, TMBoptim=m1.tmb2$logLik)
# Time
c(em=t1.em, kfoptim=t1.bfgs, TMBnlminb=t1.tmb, TMBoptim=t1.tmb2)

## -----------------------------------------------------------------------------
m1 <- MARSS(dat2, model=mod.list, form='dfa', z.score=FALSE, silent = TRUE, control = list(minit=1, maxit=10))

t2.bfgs <- system.time(m2.bfgs <- MARSS(dat2, model=mod.list, form='dfa', z.score=FALSE, silent = TRUE, method="BFGS", inits = m1))[1]
c(t2.bfgs, m2.bfgs$logLik)

## -----------------------------------------------------------------------------
df <- data.frame(t=1:37, bfgs=tidy(m1.bfgs)$estimate, em=tidy(m1.em)$estimate) |>
  pivot_longer(2:3)
ggplot(df, aes(x=t, y=value, col=name)) + geom_point()

