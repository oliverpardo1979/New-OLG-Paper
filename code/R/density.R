type_density <- function(i) {
  a <- 2.4450 / 1.006341403957648
  b <- 0.4655
  c <- 0.2329
  a * exp(-((i - b) / c)^2)
}

make_type_grid <- function(n = 4001L) {
  if (n < 101L) {
    stop_model("La cuadratura requiere al menos 101 nodos.")
  }

  i <- seq(0, 1, length.out = n)
  h <- i[2L] - i[1L]
  quadrature <- rep(h, n)
  quadrature[c(1L, n)] <- h / 2
  raw_weights <- quadrature * type_density(i)
  mass <- sum(raw_weights)

  list(
    i = i,
    density = type_density(i),
    weights = raw_weights / mass,
    raw_mass = mass
  )
}

type_cdf <- function(cutoff, grid = make_type_grid()) {
  cutoff <- min(max(cutoff, 0), 1)
  keep <- grid$i <= cutoff
  if (!any(keep)) {
    return(0)
  }
  sum(grid$weights[keep])
}
