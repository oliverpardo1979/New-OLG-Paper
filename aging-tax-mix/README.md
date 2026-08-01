# Financing Population Aging: tax mix, debt, and informality

This folder contains a separate paper and a self-contained R replication.
It does not modify the pension-regime paper or the public-spending paper in
the repository.

## Reproduce the results

From this directory:

```powershell
Rscript code/tests/test_model.R
Rscript code/scripts/run_paper.R
..\tmp\tectonic\bin\tectonic.exe --keep-logs --outdir build main.tex
```

The R code uses base R only. It calibrates the old-age outlay and the
formal-work fixed cost, solves conditional stationary equilibria, searches
over the tax simplex, generates debt-anchor counterfactuals, and writes every
table and figure used by the paper.

## Interpretation

The annual 2024--2070 series is a sequence of stationary equilibria conditional
on DANE's dependency ratio. It is not a perfect-foresight transition. The
benchmark has a representative household, so it does not measure
intragenerational consumption variance or distributional incidence. These
limitations are stated in the paper.
