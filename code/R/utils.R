stop_model <- function(..., call. = FALSE) {
  stop(sprintf(...), call. = call.)
}

assert_scalar <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
    stop_model("`%s` debe ser un escalar numérico finito.", name)
  }
  if (x < lower || x > upper) {
    stop_model("`%s` debe pertenecer a [%s, %s].", name, lower, upper)
  }
  invisible(TRUE)
}

max_relative_gap <- function(new, old, floor = 1e-8) {
  max(abs(new - old) / pmax(abs(old), floor))
}

weighted_sum <- function(x, weights) {
  if (length(x) != length(weights)) {
    stop_model("La variable y los pesos de integración tienen longitudes distintas.")
  }
  sum(x * weights)
}

crra_utility <- function(consumption, sigma) {
  value <- rep(-Inf, length(consumption))
  valid <- is.finite(consumption) & consumption > 0
  if (sigma == 1) {
    value[valid] <- log(consumption[valid])
  } else {
    value[valid] <- consumption[valid]^(1 - sigma) / (1 - sigma)
  }
  value
}

copy_list <- function(x) {
  unserialize(serialize(x, NULL))
}
