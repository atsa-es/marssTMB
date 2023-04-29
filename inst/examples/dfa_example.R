library(MARSS)
data(lakeWAplankton, package = "MARSS")
phytoplankton <- c("Cryptomonas", "Diatoms", "Greens", "Unicells", "Other.algae")
dat <- as.data.frame(lakeWAplanktonTrans) |>
  subset(Year >= 1980 & Year <= 1989) |>
  subset(select=phytoplankton) |>
  t() |>
  MARSS::zscore()

#fit with MARSS
m1.em <- MARSS(dat, model=list(R='unconstrained', m=1, tinitx=1), form='dfa', z.score=FALSE)
mym1 <- dfaTMB(obs = dat, NumStates = 1, ErrStruc='UNC')
cbind(coefficients(m1.em)$R,as.vector(mym1$Estimates$R[lower.tri(mym1$Estimates$R,diag=T)]))
cbind(coefficients(m1.em)$Z,as.vector(mym1$Estimates$Z)) 