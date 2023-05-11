library(MARSS)
data(lakeWAplankton, package = "MARSS")
phytoplankton <- c("Cryptomonas", "Diatoms", "Greens", "Unicells", "Other.algae")
dat <- as.data.frame(lakeWAplanktonTrans) |>
  subset(Year >= 1980 & Year <= 1989) |>
  subset(select=phytoplankton) |>
  t() |>
  MARSS::zscore()

# fit model
m1 <- MARSS(dat, model=list(m=3, tinitx=1), form="dfa", method="TMB")
m2 <- MARSS(dat, model=list(m=3, tinitx=1), form="dfa")

# fit model with TMB from another fit
# This is not normally how one would do this
m3 <- marssTMB:::MARSSfit.TMB(m2)

# This would be how you would normally do this
m4 <- MARSS(dat, model=list(m=3, tinitx=1), form="dfa", method="TMB", inits=m2)

