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
system.time(m1.em <- MARSS(dat, model=list(R='unconstrained', m=1, tinitx=1), form='dfa', z.score=FALSE, silent = TRUE))
system.time(m1.bfgs <- MARSS(dat, model=list(R='unconstrained', m=1, tinitx=1), form='dfa', z.score=FALSE, silent = TRUE, method="BFGS"))

## -----------------------------------------------------------------------------
library(marssTMB)
system.time(m1.tmb1 <- dfaTMB(dat, model=list(m=1, R='unconstrained')))

## -----------------------------------------------------------------------------
system.time(m1.tmb2 <- MARSS_tmb(dat, model=list(m=1, R='unconstrained', tinitx=1), control=list(fun.opt="nlminb")))

## -----------------------------------------------------------------------------
library(tidyr)
pars <- data.frame(
  name = c(paste0("R", rownames(coefficients(m1.em)$R)), paste0("Z", rownames(coefficients(m1.em)$Z))),
  EM = c(coefficients(m1.em)$R, coefficients(m1.em)$Z),
  TMB1 = c(as.vector(m1.tmb1$Estimates$R[lower.tri(m1.tmb1$Estimates$R,diag=TRUE)]), as.vector(m1.tmb1$Estimates$Z)),
  TMB2 = c(coefficients(m1.tmb2)$R, coefficients(m1.tmb2)$Z))
pars <- pars %>% tidyr::pivot_longer(2:3, names_to = "model")

## -----------------------------------------------------------------------------
library(ggplot2)
dodge <- position_dodge(width=0.5)
ggplot(pars, aes(x=name, y=value, col=model)) +
  geom_point(position=dodge) +
  ggtitle("same estimates")

