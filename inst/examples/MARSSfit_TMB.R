library(MARSS)
data(lakeWAplankton, package = "MARSS")
phytoplankton <- c("Cryptomonas", "Diatoms", "Greens", "Unicells", "Other.algae")
dat <- as.data.frame(lakeWAplanktonTrans) |>
  subset(Year >= 1980 & Year <= 1989) |>
  subset(select=phytoplankton) |>
  t() |>
  MARSS::zscore()

# set-up the model
m1 <- MARSS(dat, model=list(m=3, tinitx=1), form="dfa", method="TMB")
m2 <- MARSS(dat, model=list(m=3, tinitx=1), form="dfa")

