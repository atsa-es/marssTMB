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
mod.list <- list(R='unconstrained', m=1, tinitx=1)

## -----------------------------------------------------------------------------
t1 <- system.time(m1 <- MARSS(dat, model=mod.list, form='dfa', z.score=FALSE, silent = TRUE))
t2<- system.time(m2 <- MARSS(dat, model=mod.list, form='dfa', z.score=FALSE, silent = TRUE, method="BFGS"))

## -----------------------------------------------------------------------------
library(marssTMB)
t3<- system.time(m3 <- dfaTMB(dat, model=list(m=1, R='unconstrained')))
t4 <- system.time(m4 <- MARSS_tmb(dat, model=mod.list))
t5 <- system.time(m5 <- MARSS_tmb(dat, model=mod.list, control=list(fun.opt="optim")))

## ----echo=FALSE---------------------------------------------------------------
df <- data.frame(name=c("MARSS-EM", "MARSS-BFGS", "dfaTMB-nlminb", "MARSS_tmb-nlminb", "MARSS_tmb-optim-BFGS"),
                 time=c(t1[1], t2[1], t3[1], t4[1], t5[1]),
                 logLik=c(m1$logLik, m2$logLik, m3$logLik, m4$logLik, m5$logLik))
knitr::kable(df)

## ----echo=FALSE---------------------------------------------------------------
library(tidyr)
pars <- data.frame(
  name = c(paste0("R", rownames(coefficients(m1)$R)), paste0("Z", rownames(coefficients(m1)$Z))),
  EM = c(coefficients(m1)$R, coefficients(m1)$Z),
  BFGS = c(coefficients(m2)$R, coefficients(m2)$Z),
  TMB1 = c(as.vector(m3$Estimates$R[lower.tri(m3$Estimates$R,diag=TRUE)]), as.vector(m3$Estimates$Z)),
  TMB2 = c(coefficients(m4)$R, coefficients(m4)$Z))
pars <- pars %>% tidyr::pivot_longer(!name, names_to = "model")

## -----------------------------------------------------------------------------
library(ggplot2)
dodge <- position_dodge(width=0.5)
ggplot(pars, aes(x=name, y=value, col=model)) +
  geom_point(position=dodge) +
  ggtitle("same estimates")

## ----eval=FALSE---------------------------------------------------------------
#  # add a temperature covariate
#  temp <- as.data.frame(lakeWAplanktonTrans) |>
#      subset(Year >= 1980 & Year <= 1989) |>
#      subset(select=Temp)
#  covar <- t(temp)
#  m6 <- MARSS_tmb(dat, model=list(m=1, R='diagonal and unequal'),
#                      EstCovar = TRUE, Covars = covar)
#  m6$Estimates$D
#  
#  # add a 2nd covariate
#  TP <- as.data.frame(lakeWAplanktonTrans) |>
#      subset(Year >= 1980 & Year <= 1989) |>
#      subset(select=TP)
#  covar <- rbind(covar, t(TP))
#  m_cov2_tmb <- dfaTMB(dat, model=list(m=1, R='diagonal and unequal'),
#                      EstCovar = TRUE, Covars = covar)
#  m_cov2_tmb$Estimates$D

