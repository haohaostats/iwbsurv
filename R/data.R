rotterdam_example <- function() {
  d <- survival::rotterdam
  data.frame(time = pmax(d$dtime / 365.25, 1 / 365.25),
             status = as.integer(d$death))
}
