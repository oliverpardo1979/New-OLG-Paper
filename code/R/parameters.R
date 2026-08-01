default_parameters <- function(regime = c("compete", "pillars", "funded")) {
  regime <- match.arg(regime)

  list(
    # Demografía: un periodo equivale aproximadamente a 40 años.
    m = 0.40,
    n = 1.013^40 - 1,

    # Preferencias.
    sigma = 2.50,
    beta = 0.97^40,
    delta = 1.00^40,

    # Tecnología.
    g = 1.004^40 - 1,
    alpha = 0.40,
    depreciation = 0.00,
    A_informal = 0.10,

    # Perfiles de productividad.
    rho_informal = 2.90,
    rho_formal = 4.64,
    base_informal = 1.00,
    base_formal = 1.00,

    # Impuestos.
    tau_consumption = 0.127,
    tau_labor = 0.125,
    tau_payroll = 0.125,
    tau_capital = 0.320,

    # Transferencias públicas.
    transfer_young = 0.000,
    universal_basic_pension = 0.000,
    universal_basic_saving = 0.000,
    conditional_basic_pension = 0.000,
    conditional_basic_saving = 0.000,
    government_spending = 0.330,

    # Seguridad social.
    tau_pension = 0.125,
    replacement_rate = 0.640,
    regime = regime,
    payg_type_threshold = 0.7358,
    payg_income_threshold = 0.000
  )
}

validate_parameters <- function(param) {
  required <- c(
    "m", "n", "sigma", "beta", "delta", "g", "alpha",
    "depreciation", "A_informal", "rho_informal", "rho_formal",
    "base_informal", "base_formal", "tau_consumption", "tau_labor",
    "tau_payroll", "tau_capital", "transfer_young",
    "universal_basic_pension", "universal_basic_saving",
    "conditional_basic_pension", "conditional_basic_saving",
    "government_spending", "tau_pension", "replacement_rate",
    "regime", "payg_type_threshold", "payg_income_threshold"
  )

  missing <- setdiff(required, names(param))
  if (length(missing) > 0L) {
    stop_model("Faltan parámetros: %s.", paste(missing, collapse = ", "))
  }

  if (!param$regime %in% c("compete", "pillars", "funded")) {
    stop_model("`regime` debe ser `compete`, `pillars` o `funded`.")
  }

  assert_scalar(param$m, "m", 0, 1)
  assert_scalar(param$n, "n", -0.999999, Inf)
  assert_scalar(param$sigma, "sigma", .Machine$double.eps, Inf)
  assert_scalar(param$beta, "beta", 0, Inf)
  assert_scalar(param$delta, "delta", 0, Inf)
  assert_scalar(param$g, "g", -0.999999, Inf)
  assert_scalar(param$alpha, "alpha", .Machine$double.eps, 1 - .Machine$double.eps)
  assert_scalar(param$depreciation, "depreciation", 0, Inf)
  assert_scalar(param$A_informal, "A_informal", 0, Inf)
  assert_scalar(param$tau_consumption, "tau_consumption", -0.999999, Inf)
  assert_scalar(param$tau_labor, "tau_labor", 0, Inf)
  assert_scalar(param$tau_payroll, "tau_payroll", 0, Inf)
  assert_scalar(param$tau_capital, "tau_capital", 0, 1 - .Machine$double.eps)
  assert_scalar(param$tau_pension, "tau_pension", 0, Inf)
  assert_scalar(param$replacement_rate, "replacement_rate", 0, Inf)
  assert_scalar(param$payg_type_threshold, "payg_type_threshold", 0, 1)
  assert_scalar(param$payg_income_threshold, "payg_income_threshold", 0, Inf)

  invisible(TRUE)
}

matlab_parameter_crosswalk <- function() {
  data.frame(
    matlab = c(
      "m", "n", "sig", "beta", "del", "g", "alp", "d", "As",
      "rhos", "rhof", "bases", "basef", "tauc", "taul", "taun",
      "tauk", "tra", "ubp", "ubs", "cbp", "cbs", "G", "taup",
      "b", "thr", "bary", "pill"
    ),
    r = c(
      "m", "n", "sigma", "beta", "delta", "g", "alpha",
      "depreciation", "A_informal", "rho_informal", "rho_formal",
      "base_informal", "base_formal", "tau_consumption", "tau_labor",
      "tau_payroll", "tau_capital", "transfer_young",
      "universal_basic_pension", "universal_basic_saving",
      "conditional_basic_pension", "conditional_basic_saving",
      "government_spending", "tau_pension", "replacement_rate",
      "payg_type_threshold", "payg_income_threshold", "regime"
    ),
    note = c(
      rep("", 27),
      "compete, pillars or funded"
    ),
    stringsAsFactors = FALSE
  )
}
