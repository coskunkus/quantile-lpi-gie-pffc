# A General Quantile-Based Lifetime Performance Index — reproduction code

`R` code reproducing every table, figure and reported number in

> Alkhuffash, O. and Kuş, C. *A General Quantile-Based Lifetime Performance
> Index and its Application to the GIE Distribution under Progressive
> First-Failure Censoring.*

The paper proposes a lifetime performance index built from the median and the
interquartile range,

$$C_L^{\xi} = \frac{\xi(0.5) - L}{\xi(0.75) - \xi(0.25)},$$

as an alternative to the classical moment-based index $C_L^{M} = (\mu - L)/\sigma$,
which is undefined whenever the lifetime distribution has no finite variance.
Inference is developed for the generalized inverted exponential (GIE)
distribution under progressive first-failure censoring (PFFC).

This repository holds the code only. The manuscript is not part of it; the
scripts write the LaTeX table fragments, the `\newcommand` macro files and the
figures into `tables/` and `figures/`, which start empty.

## Quick start

```r
install.packages(c("numDeriv", "statmod", "ggplot2", "tidyr", "dplyr",
                   "patchwork", "future", "furrr", "progressr", "progress"))
```

```
cd code
Rscript run_paper.R
```

About ten minutes. That runs the seven validation scripts, Table 1, Figures 1–5,
both data analyses, and every number the running text quotes.

## The two entry points

```
Rscript run_paper.R         # everything except the simulations   (~10 min)
Rscript run_simulations.R   # the four Monte Carlo studies        (hours)
Rscript run_all.R           # both, in the right order            (hours)
```

The Monte Carlo studies of Sections 5.1–5.4 take several hours. Their raw output
is tracked here, in `code/results/`, so you do not have to run them: `run_paper.R`
reads those csv files to build Figures 4–5 and the coverage averages and
replication counts the prose quotes, and rebuilds every table and macro file the
manuscript needs except Tables 2–8 and 11. Delete `code/results/` and run
`run_simulations.R` if you want the simulations rebuilt from nothing; everything
is seeded, so the same numbers should come back.

One dependency runs from the applications into the simulations rather than the
other way. The small-sample block of Section 5.4 mirrors the first application
exactly, reading its shape, scale, limit and group size from the macro file that
`10_realdata1_ballbearings.R` writes. `run_paper.R` checks after each run that
the shipped Section 5.4 results still match the application, and says so if they
have drifted apart. `run_simulations.R` stops with an explanation if that macro
file does not exist yet.

## What is here

```
code/
  00_gie_pffc.R         the model: density, quantiles, moments, likelihood,
                        MLE, delta method, bootstrap, PFFC generator
  01,10,11,12_*.R       figures and the two data analyses
  04,05,06,07_sim_*.R   the four Monte Carlo studies
  08,09_*.R             figures and text macros built from results/
  run_paper.R           entry point: everything but the simulations
  run_simulations.R     entry point: the simulations
  run_all.R             both
  validation/           seven scripts checking the analytical results in the
                        paper against independent numerical calculations
  results/              raw simulation output, as csv
tables/                 written by the scripts: LaTeX table fragments and
figures/                \newcommand macro files, and the figures
```

`code/README.md` documents every script, what it produces and how long it takes.

## Nothing in the manuscript is typed by hand

The claim is meant literally. Every number the paper reports — in a table, in a
caption, or in the running prose — comes from a `\newcommand` that one of these
scripts writes into `tables/values_*.tex`. The manuscript refers to them by name
(`\bbAlphaHat`, `\ctCPACIa`, `\figPeakM`, …), so a sentence three sections away
from a table cannot describe an older version of it. The censored sample of
Section 6.1 is constructed by `10_realdata1_ballbearings.R` and written into the
manuscript the same way, because the construction has to obey the censoring
scheme it illustrates and a hand-typed sample did not.

## Reproducibility

- `MASTER_SEED` is set in `00_gie_pffc.R`; every script derives its own seed
  from it. The two parallel simulations assign one seed per configuration, so
  their results do not depend on the order of traversal and are identical
  whether run serially or in parallel.
- `run_paper.R` and `run_simulations.R` both print `sessionInfo()` when they
  finish.
- The seven validation scripts run as part of `run_paper.R` and stop with an
  error if any check fails. Between them they check the analytical derivatives,
  the location–scale invariance of both indexes, the closed-form moment against
  two independent numerical integrals, the finiteness bounds of Appendix B, the
  closed-form Hessian and the fixed-point iteration, the decomposition of the
  bootstrap coverage, and — against a direct simulation of the life test itself
  — the PFFC sample generator on which every simulated number depends.

## License

MIT, see `LICENSE`.
