# Pension Regime Choice under Informality

Paper and R implementation by Oliver Pardo.

The model is a two-period OLG economy in which heterogeneous workers choose
simultaneously among:

1. informal employment;
2. formal employment with a funded account;
3. formal employment with a PAYG pension.

Funded contributions are portable. PAYG eligibility increases with worker
type, non-eligible contributors receive a partial refund, and eligible benefits
have a floor and a ceiling. The model obtains positive masses in all three
alternatives without quotas or random-utility shocks.

## Repository structure

- `main.tex`: paper and algorithm appendix.
- `references.bib`: bibliography.
- `code/R/`: model functions.
- `code/scripts/run_endogenous_analysis.R`: reproducible analysis.
- `code/tests/run_endogenous_tests.R`: endogenous-choice validation.
- `results/endogenous_steady_state.csv`: benchmark equilibrium.
- `results/endogenous_sensitivity.csv`: eligibility sensitivity.
- `results/endogenous_micro_choices.csv`: type-level choices and utilities.
- `results/legacy_exogenous/`: archived results from the earlier exogenous
  architecture draft.
- `output/pdf/`: compiled paper.

## Run the model

From `code/`:

```powershell
Rscript tests/run_endogenous_tests.R
Rscript scripts/run_endogenous_analysis.R
```

The analysis first solves a 101-type grid and uses that equilibrium to
initialize the reported 501-type solution.

## Compile the paper

From the repository root:

```powershell
tectonic --keep-logs --outdir build main.tex
```

The endogenous-choice module currently reports stationary equilibria. The
legacy transition code is retained for the exogenous benchmarks and is not
used for the endogenous results in the paper.
