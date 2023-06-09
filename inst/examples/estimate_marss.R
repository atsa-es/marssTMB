library(MARSS)
data(lakeWAplankton, package = "MARSS")
phytoplankton <- c("Cryptomonas", "Diatoms", "Greens", "Unicells", "Other.algae")
dat <- as.data.frame(lakeWAplanktonTrans) |>
  subset(Year >= 1980 & Year <= 1989) |>
  subset(select=phytoplankton) |>
  t() |>
  MARSS::zscore()

# set-up the model
mod <- MARSS(dat, model=list(m=3, tinitx=1), form="dfa", fit=FALSE, silent=TRUE)
# fit
fit <- estimate_marss(mod)

