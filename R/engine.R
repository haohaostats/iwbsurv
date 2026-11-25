rp_basis <- function(x, knots) {
  kmin <- knots[1]; kmax <- knots[length(knots)]
  B <- cbind(1, x)
  if (length(knots) > 2) {
    for (kj in knots[-c(1, length(knots))]) {
      lam <- (kmax - kj) / (kmax - kmin)
      B <- cbind(B, pmax(x - kj, 0)^3 - lam * pmax(x - kmin, 0)^3 -
                    (1 - lam) * pmax(x - kmax, 0)^3)
    }
  }
  B
}

rp_dbasis <- function(x, knots) {
  kmin <- knots[1]; kmax <- knots[length(knots)]
  D <- cbind(0, rep(1, length(x)))
  if (length(knots) > 2) {
    for (kj in knots[-c(1, length(knots))]) {
      lam <- (kmax - kj) / (kmax - kmin)
      D <- cbind(D, 3 * pmax(x - kj, 0)^2 - 3 * lam * pmax(x - kmin, 0)^2 -
                    3 * (1 - lam) * pmax(x - kmax, 0)^2)
    }
  }
  D
}

.obs_fit_rp_core <- function(time, status, t_obs, k = 2, nsim = 2000, seed = 1L) {
  st <- as.integer(status)
  tt <- as.numeric(time)
  st[tt > t_obs] <- 0L
  tt <- pmax(pmin(tt, t_obs), 1e-6)
  if (sum(st) < k + 3) stop("Too few events before t_obs for the spline model.")
  train <- data.frame(time = tt, status = st)

  fit <- flexsurv::flexsurvspline(survival::Surv(time, status) ~ 1, data = train,
                                  k = k, scale = "hazard")
  est <- fit$res.t[, "est"]
  Sig <- fit$cov
  set.seed(seed)
  Z <- matrix(stats::rnorm(nsim * length(est)), nsim, length(est))
  theta <- sweep(Z %*% chol(Sig), 2, est, "+")
  colnames(theta) <- names(est)
  list(fit = fit, knots = as.numeric(fit$knots), theta = theta,
       theta_hat = est, train = train)
}

obs_eval_rp <- function(obs, tt, theta = obs$theta) {
  pos <- tt > 0
  x   <- log(tt[pos])
  B   <- rp_basis(x, obs$knots)
  dB  <- rp_dbasis(x, obs$knots)
  Hp  <- exp(theta %*% t(B))
  hp  <- sweep((theta %*% t(dB)) * Hp, 2, tt[pos], "/")
  M <- nrow(theta)
  H <- matrix(0, M, length(tt)); H[, pos] <- Hp
  h <- matrix(0, M, length(tt)); h[, pos] <- hp
  if (any(!pos)) h[, !pos] <- hp[, 1]
  valid <- rowSums(!is.finite(hp) | hp <= 0) == 0
  list(H = H, h = h, valid = valid)
}

draw_ext_cgomp <- function(M, S_c, t_obs, t_a, Sa_mean, ess = 50,
                           loggamma_mean = log(0.08), loggamma_sd = 0.5,
                           seed = 1L) {
  stopifnot(length(S_c) == M, t_a > t_obs, Sa_mean > 0, Sa_mean < 1)
  set.seed(seed)
  Sa    <- stats::rbeta(M, ess * Sa_mean, ess * (1 - Sa_mean))
  gamma <- stats::rlnorm(M, loggamma_mean, loggamma_sd)
  q     <- Sa / S_c
  ok    <- is.finite(q) & q > 0 & q < 1
  rho   <- rep(NA_real_, M)
  rho[ok] <- gamma[ok] * (-log(q[ok])) / expm1(gamma[ok] * (t_a - t_obs))
  list(Sa = Sa, gamma = gamma, rho = rho, ok = ok, p_inc = mean(!ok))
}

ext_dH <- function(rho, gamma, u) (rho / gamma) * expm1(outer(gamma, u))

ext_h <- function(rho, gamma, u) rho * exp(outer(gamma, u))

weight_iwb <- function(t_obs, tg, dHo, dHe, lh0, lrho, k = 10) {
  if (!requireNamespace("scam", quietly = TRUE)) stop("Package 'scam' is required.")
  Vo <- c(stats::var(lh0),  apply(log(dHo), 2, stats::var))
  Ve <- c(stats::var(lrho), apply(log(dHe), 2, stats::var))
  tt <- c(t_obs, tg)
  r_raw <- Vo / (Vo + Ve)
  if (any(!is.finite(r_raw) | r_raw <= 0 | r_raw >= 1))
    stop("Non-finite or degenerate dispersion ratio; review the model and anchor.")
  dat <- data.frame(t = tt, y = stats::qlogis(r_raw))
  k_eff <- max(5L, min(as.integer(k), length(tt) - 2L))
  fit <- scam::scam(y ~ s(t, k = k_eff, bs = "mpi"), data = dat)
  ell <- as.numeric(stats::predict(fit, newdata = data.frame(t = tt)))

  r_s <- stats::plogis(ell)
  r0  <- r_s[1]
  pi  <- (r_s - r0) / (1 - r0)
  pi[1] <- 0
  pi  <- pmin(pmax(pi, 0), 1 - 1e-12)

  G <- length(tt)
  dell <- numeric(G)
  dell[1] <- (ell[2] - ell[1]) / (tt[2] - tt[1])
  dell[G] <- (ell[G] - ell[G - 1]) / (tt[G] - tt[G - 1])
  if (G > 2) dell[2:(G - 1)] <- (ell[3:G] - ell[1:(G - 2)]) / (tt[3:G] - tt[1:(G - 2)])
  dpi <- pmax(stats::dlogis(ell) * dell / (1 - r0), 0)

  list(t = tt, pi = pi, dpi = dpi, r_raw = r_raw, r_smooth = r_s,
       V_obs = Vo, V_ext = Ve, edf = sum(fit$edf), scam = fit)
}

kappa_star <- function(h_o, h_e, H_o, H_e, pi, dpi) {
  A <- sweep(h_e - h_o, 2, pi, "*") + sweep(H_e - H_o, 2, dpi, "*")
  neg <- A < 0
  if (!any(neg)) return(1)
  min(1, min(h_o[neg] / (-A[neg])))
}

summ_draws <- function(S, level = 0.95) {
  a <- (1 - level) / 2
  cbind(mean = colMeans(S),
        lo   = apply(S, 2, stats::quantile, probs = a,     names = FALSE),
        hi   = apply(S, 2, stats::quantile, probs = 1 - a, names = FALSE))
}

.iwb_fit_core <- function(time, status, t_obs, t_max, t_a = t_max, Sa_mean, ess = 50,
                    step = (t_max - t_obs) / 200, nsim = 2000, k_rp = 2, k_w = 10,
                    loggamma_mean = log(0.08), loggamma_sd = 0.5,
                    level = 0.95, seed = 1L, obs = NULL, keep_draws = FALSE) {
  if (is.null(obs))
    obs <- obs_fit_rp(time, status, t_obs, k = k_rp, nsim = nsim, seed = seed)

  tt <- sort(unique(round(c(seq(0, t_max, by = step), t_obs, t_a, t_max), 10)))
  ev <- obs_eval_rp(obs, tt)
  p_neg <- mean(!ev$valid)
  Ho <- ev$H[ev$valid, , drop = FALSE]
  ho <- ev$h[ev$valid, , drop = FALSE]
  i0 <- which.min(abs(tt - t_obs))
  ie <- which(tt >= t_obs)
  S_c <- exp(-Ho[, i0])

  ext <- draw_ext_cgomp(nrow(Ho), S_c, t_obs, t_a, Sa_mean, ess,
                        loggamma_mean, loggamma_sd, seed = seed + 7919L)
  ok <- ext$ok
  if (sum(ok) < 100) stop("Fewer than 100 compatible draws; increase nsim or review the anchor.")
  Ho <- Ho[ok, , drop = FALSE]
  ho <- ho[ok, , drop = FALSE]
  rho <- ext$rho[ok]
  gam <- ext$gamma[ok]

  u   <- tt[ie] - t_obs
  He  <- Ho
  he  <- ho
  He[, ie] <- Ho[, i0] + ext_dH(rho, gam, u)
  he[, ie] <- ext_h(rho, gam, u)
  if (any(!is.finite(He)) || any(!is.finite(he)))
    stop("Non-finite external tail; review the time units and shape parameters.")
  dHo <- Ho[, ie, drop = FALSE] - Ho[, i0]
  dHe <- He[, ie, drop = FALSE] - He[, i0]

  w <- weight_iwb(t_obs, tt[ie][-1], dHo[, -1, drop = FALSE], dHe[, -1, drop = FALSE],
                  lh0 = log(ho[, i0]), lrho = log(rho), k = k_w)
  pi_full  <- numeric(length(tt)); pi_full[ie]  <- w$pi
  dpi_full <- numeric(length(tt)); dpi_full[ie] <- w$dpi

  kap <- kappa_star(ho[, ie, drop = FALSE], he[, ie, drop = FALSE],
                    Ho[, ie, drop = FALSE], He[, ie, drop = FALSE], w$pi, w$dpi)

  Hb <- Ho + sweep(He - Ho, 2, kap * pi_full, "*")
  hb <- ho + kap * (sweep(he - ho, 2, pi_full, "*") + sweep(He - Ho, 2, dpi_full, "*"))

  v <- (tt[ie] - t_obs) / (t_max - t_obs)
  pib  <- numeric(length(tt)); pib[ie]  <- stats::pbeta(v, 3, 3)
  dpib <- numeric(length(tt)); dpib[ie] <- stats::dbeta(v, 3, 3) / (t_max - t_obs)
  Hbeta <- Ho + sweep(He - Ho, 2, pib, "*")
  hbeta <- ho + sweep(he - ho, 2, pib, "*") + sweep(He - Ho, 2, dpib, "*")

  out <- list(
    time   = tt,
    t_obs  = t_obs, t_max = t_max, t_a = t_a,
    obs    = summ_draws(exp(-Ho), level),
    ext    = summ_draws(exp(-He), level),
    iwb    = summ_draws(exp(-Hb), level),
    beta   = summ_draws(exp(-Hbeta), level),
    weight = list(pi = pi_full, r_raw = w$r_raw, r_smooth = w$r_smooth,
                  t = w$t, edf = w$edf),
    kappa  = kap,
    p_inc  = ext$p_inc,
    diag   = list(n_draws = sum(ok), p_neg = p_neg,
                  min_h_iwb  = min(hb[, ie]),
                  min_h_beta = min(hbeta[, ie]),
                  pi_end     = utils::tail(w$pi, 1)),
    train  = obs$train,
    obs_fit = obs$fit,
    anchor = list(t_a = t_a, Sa_mean = Sa_mean, ess = ess)
  )
  if (keep_draws) out$draws <- list(H_obs = Ho, H_ext = He, H_iwb = Hb, H_beta = Hbeta)
  class(out) <- "iwb_fit"
  out
}
