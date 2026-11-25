.scalar <- function(x, name, lower = -Inf, upper = Inf, integer = FALSE,
                    open_lower = FALSE, open_upper = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x)
  if (ok) ok <- if (open_lower) x > lower else x >= lower
  if (ok) ok <- if (open_upper) x < upper else x <= upper
  if (ok && integer) ok <- x == floor(x)
  if (!ok) stop("Invalid value for '", name, "'.", call. = FALSE)
  invisible(x)
}

.data_check <- function(time, status) {
  if (!is.numeric(time) || !length(time) || any(!is.finite(time)) ||
      any(time < 0))
    stop("'time' must contain finite, non-negative numbers.", call. = FALSE)
  if (!(is.numeric(status) || is.logical(status)) ||
      length(status) != length(time) || anyNA(status) ||
      any(!status %in% c(0, 1)))
    stop("'status' must contain 0 (censored) or 1 (event), one per time.",
         call. = FALSE)
}

.preserve_rng <- function(expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  })
  force(expr)
}

obs_fit_rp <- function(time, status, t_obs, k = 2, nsim = 2000, seed = 1L) {
  .data_check(time, status)
  .scalar(t_obs, "t_obs", 0, open_lower = TRUE)
  .scalar(k, "k", 0, integer = TRUE)
  .scalar(nsim, "nsim", 100, integer = TRUE)
  .scalar(seed, "seed", 0, .Machine$integer.max - 7919, integer = TRUE)
  out <- .preserve_rng(.obs_fit_rp_core(time, status, t_obs, k, nsim, seed))
  out$t_obs <- t_obs
  out$k <- k
  out$nsim <- nsim
  class(out) <- "iwb_observation"
  out
}

iwb_fit <- function(time, status, t_obs, t_max, t_a = t_max, Sa_mean, ess = 50,
                    step = (t_max - t_obs) / 200, nsim = 2000, k_rp = 2, k_w = 10,
                    loggamma_mean = log(0.08), loggamma_sd = 0.5,
                    level = 0.95, seed = 1L, obs = NULL, keep_draws = FALSE) {
  .data_check(time, status)
  .scalar(t_obs, "t_obs", 0, open_lower = TRUE)
  .scalar(t_max, "t_max", t_obs, open_lower = TRUE)
  .scalar(t_a, "t_a", t_obs, t_max, open_lower = TRUE)
  .scalar(Sa_mean, "Sa_mean", 0, 1, open_lower = TRUE, open_upper = TRUE)
  .scalar(ess, "ess", 0, open_lower = TRUE)
  .scalar(step, "step", 1e-10, open_lower = TRUE)
  .scalar(nsim, "nsim", 100, integer = TRUE)
  .scalar(k_rp, "k_rp", 0, integer = TRUE)
  .scalar(k_w, "k_w", 5, integer = TRUE)
  .scalar(loggamma_mean, "loggamma_mean")
  .scalar(loggamma_sd, "loggamma_sd", 0)
  .scalar(level, "level", 0, 1, open_lower = TRUE, open_upper = TRUE)
  .scalar(seed, "seed", 0, .Machine$integer.max - 7919, integer = TRUE)
  if (!is.logical(keep_draws) || length(keep_draws) != 1L || is.na(keep_draws))
    stop("'keep_draws' must be TRUE or FALSE.", call. = FALSE)
  if (any(c(t_obs, t_a, t_max) != round(c(t_obs, t_a, t_max), 10)))
    stop("Cutoff, anchor and horizon require at most 10 decimal places.",
         call. = FALSE)
  if (t_max / step > 1e6)
    stop("The requested grid exceeds one million points.", call. = FALSE)
  grid <- sort(unique(round(c(seq(0, t_max, by = step), t_obs, t_a, t_max), 10)))
  if (sum(grid >= t_obs) < 7L)
    stop("At least seven grid points at or after 't_obs' are required.",
         call. = FALSE)
  if (!is.null(obs)) {
    if (!inherits(obs, "iwb_observation") || obs$t_obs != t_obs)
      stop("'obs' must be an obs_fit_rp() result with the same cutoff.",
           call. = FALSE)
    st <- as.integer(status)
    st[time > t_obs] <- 0L
    train <- data.frame(time = pmax(pmin(as.numeric(time), t_obs), 1e-6),
                        status = st)
    if (!identical(obs$train, train))
      stop("'obs' was fitted to different observations.", call. = FALSE)
  }
  out <- .preserve_rng(.iwb_fit_core(
    time, status, t_obs, t_max, t_a, Sa_mean, ess, step, nsim, k_rp, k_w,
    loggamma_mean, loggamma_sd, level, seed, obs, keep_draws))
  out$level <- level
  out$call <- match.call()
  out$settings <- list(step = step, seed = seed,
                       nsim = if (is.null(obs)) nsim else obs$nsim,
                       k_rp = if (is.null(obs)) k_rp else obs$k, k_w = k_w,
                       loggamma_mean = loggamma_mean, loggamma_sd = loggamma_sd)
  out
}
