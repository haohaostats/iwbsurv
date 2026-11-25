library(iwbsurv)

d <- rotterdam_example()
km <- survival::survfit(survival::Surv(time, status) ~ 1, data = d)
anchor <- stats::stepfun(km$time, c(1, km$surv))(15)

fit <- iwb_fit(
  time = d$time,
  status = d$status,
  t_obs = 5,
  t_max = 15,
  Sa_mean = anchor,
  ess = 50,
  nsim = 2000,
  seed = 2026
)

summary(fit)
predict(fit, times = c(5, 10, 15))
predict(fit, times = c(5, 10, 15), type = "weight")
restricted_survival(fit)
plot(fit, xlab = "Time (years)")
