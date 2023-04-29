library(ggplot2)

set.seed(123)
x = cumsum(rnorm(30))
y = x + rnorm(length(x), 0, 0.01)
estimate_drift = TRUE # U in MARSS
estimate_rho = FALSE # AR(1) parameter, b in MARSS

res <- uniTMB(y)

res$coef

ggplot(res$df, aes(t, pred)) + 
  geom_ribbon(aes(ymin=pred-2*se, ymax = pred+2*se),alpha=0.5) + 
  geom_line() + 
  geom_point(aes(t,y),col="red",alpha=0.5) + 
  xlab("Time") + ylab("Data") + 
  theme_bw()
