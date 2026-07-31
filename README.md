# From Competing Pension Regimes to a Pillar System

Paper and R implementation by Oliver Pardo.

This repository contains a two-period OLG model for Colombia with endogenous
informality and pension participation. It compares:

1. Law 100 of 1993: workers choose informality, formal capitalization, or
   formal PAYG;
2. Law 2381 of 2024: formal earnings up to 2.3 SMLMV enter PAYG and only the
   excess enters the funded component.

The reform is simulated as a counterfactual because the integral entry into
force of Law 2381 remains suspended. One model period is approximately 40
years, so transition paths are generational rather than annual forecasts.

## Main outputs

- `main.tex`: paper and algorithm appendix.
- `code/R/policy_transition.R`: pillar mapping, endogenous transition,
  two-cycle validation, microdata, and distributional statistics.
- `code/scripts/run_policy_reform_analysis.R`: reproducible policy analysis.
- `code/tests/run_policy_reform_tests.R`: stationary, transition, fiscal, and
  distributional tests.
- `results/policy_transition_path.csv`: endogenous paths and reported window.
- `results/policy_distribution_summary.csv`: wage and welfare distributions.
- `results/policy_scenario_summary.csv`: Law 100, central pillars, and stress.
- `results/policy_parameter_crosswalk.csv`: source and status of policy inputs.
- `figures/policy_*.png`: paper figures.
- `output/pdf/main.pdf`: compiled paper.

The older benchmark files remain available as `endogenous_*.csv`. Results from
the exogenous predecessor are archived in `results/legacy_exogenous/`.

## Run

From `code/`:

```powershell
Rscript tests/run_endogenous_tests.R
Rscript tests/run_policy_reform_tests.R
Rscript scripts/run_policy_reform_analysis.R
```

The policy script solves 201-type stationary equilibria and a 101-type
transition with exact insertion of endogenous choice roots. The last three
cohorts form a terminal buffer and are excluded from reported transition
figures. Validation accepts either convergence to the pillar steady state or a
nondegenerate two-generation cycle; the central calibration produces the
latter.

## Compile

From the repository root:

```powershell
tectonic --keep-logs --outdir build main.tex
```
