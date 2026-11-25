library(iwbsurv)

near <- function(x, y, tolerance = 1e-10)
  stopifnot(isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tolerance)))

fails <- function(expr) stopifnot(inherits(try(force(expr), silent = TRUE), "try-error"))

d <- rotterdam_example()
time <- d$time
status <- d$status
set.seed(91)
rng <- .Random.seed
obs <- obs_fit_rp(time, status, 5, nsim = 300, seed = 2026)
stopifnot(identical(rng, .Random.seed))
fit <- iwb_fit(time, status, 5, 15, Sa_mean = 0.39, obs = obs,
               step = 0.1, seed = 2026, keep_draws = TRUE)
stopifnot(identical(rng, .Random.seed))

at <- c(0.5, 2, 4, 7, 10)
ev <- iwbsurv:::obs_eval_rp(obs, at, matrix(obs$theta_hat, nrow = 1))
ref_s <- summary(obs$fit, t = at, type = "survival", ci = FALSE)[[1]]$est
ref_h <- summary(obs$fit, t = at, type = "hazard", ci = FALSE)[[1]]$est
near(as.numeric(exp(-ev$H)), ref_s)
near(as.numeric(ev$h), ref_h)

ext <- iwbsurv:::draw_ext_cgomp(300, rep(0.8, 300), 4, 10, 0.4, seed = 2)
ok <- ext$ok
near(exp(-iwbsurv:::ext_dH(ext$rho[ok], ext$gamma[ok], 6)) * 0.8, ext$Sa[ok])
kap <- iwbsurv:::kappa_star(matrix(1), matrix(0.2), matrix(2),
                           matrix(1), 0.5, 2)
near(kap, 1 / 2.4)

pre <- fit$time <= fit$t_obs
for (nm in c("H_ext", "H_iwb", "H_beta"))
  near(fit$draws[[nm]][, pre], fit$draws$H_obs[, pre])
stopifnot(all(fit$iwb >= 0 & fit$iwb <= 1),
          all(fit$iwb[, "lo"] <= fit$iwb[, "hi"]),
          all(diff(fit$weight$pi) >= -1e-10),
          all(fit$weight$pi[pre] == 0),
          fit$kappa >= 0, fit$kappa <= 1,
          fit$diag$min_h_iwb >= -1e-10)
near(fit$iwb[, "mean"], colMeans(exp(-fit$draws$H_iwb)))
near(fit$iwb[, "lo"],
     apply(exp(-fit$draws$H_iwb), 2, quantile, probs = 0.025, names = FALSE))
near(predict(fit, type = "weight")$effective, fit$kappa * fit$weight$pi)
near(predict(fit, times = at)$mean,
     approx(fit$time, fit$iwb[, "mean"], at)$y)
stopifnot(nrow(as.data.frame(fit)) == 4 * length(fit$time))

toy <- fit
toy$time <- c(0, 4, 10)
toy$iwb <- cbind(mean = 1 - toy$time / 10, lo = 0, hi = 1)
toy$t_max <- 10
near(restricted_survival(toy, 0, 10), 5)
near(restricted_survival(toy, 4, 10), 1.8)

fit2 <- iwb_fit(time, status, 5, 15, t_a = 12, Sa_mean = 0.39,
                obs = obs, step = 0.13, seed = 2026)
stopifnot(15 %in% fit2$time, all(c(5, 12) %in% fit2$time))
near(utils::tail(predict(fit2)$time, 1), 15)

fixture <- read.csv("reference-rotterdam.csv")
actual <- predict(fit, times = fixture$time)
reference_error <- c(
  mean = max(abs(actual$mean - fixture$mean)),
  lower = max(abs(actual$lower - fixture$lower)),
  upper = max(abs(actual$upper - fixture$upper)),
  weight = max(abs(predict(fit, times = fixture$time, type = "weight")$pi -
                   fixture$pi)),
  kappa = abs(fit$kappa - fixture$kappa[1])
)
print(reference_error)
stopifnot(all(reference_error < 1e-6))

fails(iwb_fit(time, status + 1, 5, 15, Sa_mean = 0.39))
fails(iwb_fit(time, status, 5, 15, Sa_mean = 1))
fails(iwb_fit(time, status, 5, 15, t_a = 16, Sa_mean = 0.39))
fails(iwb_fit(time, status, 5, 15, Sa_mean = 0.39, step = 5))
fails(iwb_fit(time, status, 5, 15, Sa_mean = 0.39, obs = list()))
fails(iwb_fit(time + 1, status, 5, 15, Sa_mean = 0.39, obs = obs))
fails(predict(fit, times = 16))
fails(restricted_survival(fit, 10, 4))

fig <- tempfile(fileext = ".pdf")
grDevices::pdf(fig, width = 7, height = 5.5)
plot(fit)
plot(fit, type = "weight")
grDevices::dev.off()
stopifnot(file.info(fig)$size > 1000)
unlink(fig)
cat("All numerical and interface checks passed.\n")
