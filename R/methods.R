print.iwb_fit <- function(x, ...) {
  cat("Information-weighted survival extrapolation\n")
  cat("Cutoff:", x$t_obs, "  Horizon:", x$t_max, "\n")
  cat("Anchor:", x$anchor$Sa_mean, "at", x$t_a,
      "  Concentration:", x$anchor$ess, "\n")
  cat("Retained draws:", x$diag$n_draws, "  Kappa:", format(x$kappa), "\n")
  cat("Terminal effective borrowing:",
      format(x$kappa * utils::tail(x$weight$pi, 1)), "\n")
  invisible(x)
}

summary.iwb_fit <- function(object, ...) {
  j <- length(object$time)
  data.frame(component = c("iwb", "obs", "ext", "beta"),
             time = object$time[j],
             mean = vapply(object[c("iwb", "obs", "ext", "beta")],
                           function(x) x[j, "mean"], numeric(1)),
             lower = vapply(object[c("iwb", "obs", "ext", "beta")],
                            function(x) x[j, "lo"], numeric(1)),
             upper = vapply(object[c("iwb", "obs", "ext", "beta")],
                            function(x) x[j, "hi"], numeric(1)),
             row.names = NULL)
}

predict.iwb_fit <- function(object, times = object$time,
                            component = c("iwb", "obs", "ext", "beta"),
                            type = c("survival", "weight"), ...) {
  component <- match.arg(component)
  type <- match.arg(type)
  if (!is.numeric(times) || !length(times) || any(!is.finite(times)) ||
      any(times < 0 | times > object$t_max))
    stop("'times' must lie between zero and the fitted horizon.", call. = FALSE)
  interp <- function(y) stats::approx(object$time, y, xout = times)$y
  if (type == "weight")
    return(data.frame(time = times, pi = interp(object$weight$pi),
                      effective = object$kappa * interp(object$weight$pi)))
  y <- object[[component]]
  data.frame(time = times, mean = interp(y[, "mean"]),
             lower = interp(y[, "lo"]), upper = interp(y[, "hi"]))
}

as.data.frame.iwb_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  out <- do.call(rbind, lapply(c("iwb", "obs", "ext", "beta"), function(m) {
    data.frame(component = m, stats::predict(x, component = m))
  }))
  rownames(out) <- row.names
  out
}

restricted_survival <- function(object, from = object$t_obs, to = object$t_max,
                                 component = c("iwb", "obs", "ext", "beta"),
                                 n = 401L) {
  if (!inherits(object, "iwb_fit"))
    stop("'object' must be an iwb_fit result.", call. = FALSE)
  component <- match.arg(component)
  .scalar(from, "from", 0, object$t_max)
  .scalar(to, "to", from, object$t_max, open_lower = TRUE)
  .scalar(n, "n", 2, integer = TRUE)
  tt <- seq(from, to, length.out = n)
  yy <- stats::predict(object, times = tt, component = component)$mean
  sum(diff(tt) * (utils::head(yy, -1) + utils::tail(yy, -1)) / 2)
}

plot.iwb_fit <- function(x, type = c("survival", "weight"),
                         components = c("iwb", "obs", "ext", "beta"),
                         interval = TRUE, xlab = "Time",
                         ylab = NULL, ...) {
  type <- match.arg(type)
  components <- match.arg(components, c("iwb", "obs", "ext", "beta"),
                          several.ok = TRUE)
  colors <- c(iwb = "#1764AB", obs = "#263746", ext = "#7952A3",
              beta = "#CB661B")
  labels <- c(iwb = "IWB", obs = "RP continuation",
              ext = "External tail", beta = "Beta blending")
  ltys <- c(iwb = 1, obs = 3, ext = 2, beta = 4)
  op <- graphics::par(mar = c(7, 4.4, 1.8, 1), las = 1, bty = "l",
                       mgp = c(2.4, 0.7, 0))
  on.exit(graphics::par(op))
  if (is.null(ylab))
    ylab <- if (type == "survival") "Survival probability" else "Borrowing weight"
  graphics::plot(x$time, rep(NA_real_, length(x$time)), type = "n",
                 xlim = c(0, x$t_max), ylim = c(0, 1),
                 xlab = xlab, ylab = ylab, ...)
  if (type == "survival") {
    if (isTRUE(interval) && "iwb" %in% components)
      graphics::polygon(c(x$time, rev(x$time)),
                         c(x$iwb[, "lo"], rev(x$iwb[, "hi"])),
                         col = grDevices::adjustcolor(colors["iwb"], 0.17),
                         border = NA)
    for (m in c(setdiff(components, "iwb"), intersect("iwb", components)))
      graphics::lines(x$time, x[[m]][, "mean"], col = colors[m],
                       lty = ltys[m], lwd = 2.2)
    labs <- labels[components]
    if (isTRUE(interval) && "iwb" %in% components)
      labs[components == "iwb"] <- paste0("IWB (", 100 * x$level, "% interval)")
    graphics::legend("bottom", inset = c(0, -0.30), xpd = NA, bty = "n",
                      legend = labs, col = colors[components],
                      lty = ltys[components], lwd = 2.2, ncol = 2, cex = 0.85)
  } else {
    graphics::lines(x$time, x$weight$pi, col = colors["obs"], lty = 2, lwd = 2.2)
    graphics::lines(x$time, x$kappa * x$weight$pi, col = colors["iwb"], lwd = 2.2)
    graphics::legend("bottom", inset = c(0, -0.30), xpd = NA, bty = "n",
                      legend = c("Unscaled weight", "Effective borrowing"),
                      col = colors[c("obs", "iwb")], lty = c(2, 1),
                      lwd = 2.2, ncol = 2, cex = 0.85)
  }
  invisible(x)
}
