context("DFA estimation")
library(MARSS)

test_that("1 trend DFA model works", {
  set.seed(1)
  
  data(lakeWAplankton, package = "MARSS")
  phytoplankton <- c("Cryptomonas", "Diatoms", "Greens", "Unicells", "Other.algae")
  dat <- as.data.frame(lakeWAplanktonTrans) |>
    subset(Year >= 1980 & Year <= 1989) |>
    subset(select=phytoplankton) |>
    t() |>
    MARSS::zscore()
  
  m1.tmb <- dfaTMB(dat, model=list(m=1, R='diagonal and equal'))
  expect_equal(m1.tmb$Estimates$logsdObs[1], -0.1125247, tolerance = 0.01)
  expect_equal(c(m1.tmb$Estimates$Z), c(0.247, 0.227, 0.127, 0.368, -0.091), tolerance = 0.01)
  
  m2.tmb <- dfaTMB(dat, model=list(m=1, R='diagonal and unequal'))
  expect_equal(m2.tmb$Estimates$logsdObs, 
               c(-0.143, -0.060, -0.025, -0.887, -0.005), 
               tolerance = 0.01)
  expect_equal(c(m2.tmb$Estimates$Z), c(0.332, 0.220, 0.144, 0.615, -0.033), tolerance = 0.01)

  m3.tmb <- dfaTMB(dat, model=list(m=1, R='unconstrained'))
  expect_equal(m3.tmb$Estimates$logsdObs, 
               c(-0.033,-0.037,-0.018, -0.510, -0.064), 
               tolerance = 0.01)
  expect_equal(c(m3.tmb$Estimates$Z), c(0.115,0.123, 0.089, 0.389, -0.164), tolerance = 0.01)
  
})

test_that("2 trend DFA model works", {
  set.seed(1)
  
  data(lakeWAplankton, package = "MARSS")
  phytoplankton <- c("Cryptomonas", "Diatoms", "Greens", "Unicells", "Other.algae")
  dat <- as.data.frame(lakeWAplanktonTrans) |>
    subset(Year >= 1980 & Year <= 1989) |>
    subset(select=phytoplankton) |>
    t() |>
    MARSS::zscore()
  
  m1.tmb <- dfaTMB(dat, model=list(m=2, R='diagonal and equal'))
  expect_equal(m1.tmb$Estimates$logsdObs[1], -0.2424557, tolerance = 0.01)
  expect_equal(c(m1.tmb$Estimates$Z)[1:5], c(0.394, 0.238, 0.392, 0.432, 0.188), tolerance = 0.01)
  
  m2.tmb <- dfaTMB(dat, model=list(m=2, R='diagonal and unequal'))
  expect_equal(m2.tmb$Estimates$logsdObs, 
               c(-0.173,-0.063,-0.108,-0.848,-0.431), 
               tolerance = 0.01)
  expect_equal(c(m2.tmb$Estimates$Z)[1:5], c(0.416, 0.158, 0.357, 0.567, 0.392), tolerance = 0.01)
  
  m3.tmb <- dfaTMB(dat, model=list(m=1, R='unconstrained'))
  expect_equal(m3.tmb$Estimates$logsdObs, 
               c(-0.033,-0.037,-0.018, -0.510, -0.064), 
               tolerance = 0.01)
  expect_equal(c(m3.tmb$Estimates$Z)[1:5], c(0.115, 0.123, 0.089, 0.389, -0.164), tolerance = 0.01)
  
})
