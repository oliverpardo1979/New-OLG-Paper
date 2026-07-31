# Pension Pillars under Informality

Paper by Oliver Pardo: “Pension Pillars under Informality: General-Equilibrium Effects of Reallocating Contributions.”

The model is a two-period OLG economy with endogenous formal–informal sector choice,
present-biased saving, funded accounts, PAYG pensions, and an endogenous consumption tax.

## Files

- main.tex: manuscript and numerical-algorithm appendix.
- references.bib: bibliography.
- results/steady_state_results.csv: steady-state outcomes.
- results/transition_diagnostic.csv: selected transition dates.
- output/pdf/: compiled paper.

Compile from the repository root with:

    tectonic --keep-logs --outdir build main.tex

The main comparison uses 501 productivity types and approximately matches PAYG receipts.
The transition uses a coarser grid and is a numerical diagnostic, not a forecast.
