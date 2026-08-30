# Reproducibility material

**A General Quantile-Based Lifetime Performance Index and its Application to the GIE Distribution under Progressive First-Failure Censoring**

This directory contains all `R` code needed to reproduce every table, figure and
numerical statement in the paper. The same code is available publicly at
<https://github.com/coskunkus/quantile-lpi-gie-pffc>, which is the address given
in the concluding section of the manuscript.

## Requirements

```r
install.packages(c("numDeriv", "statmod", "ggplot2", "tidyr", "dplyr",
                   "patchwork", "future", "furrr", "progressr", "progress"))
```

`05_sim_ci.R` and `06_sim_heavytail.R` run in parallel by default and need
`future` and `furrr`; set `PARALLEL <- FALSE` at the top of either to run on a
single core instead. All four simulations show a progress bar with an estimate
of the time remaining, which needs `progressr` (`progress` gives the nicer bar
when running interactively). Every script installs anything missing on first
use, and falls back gracefully if the bar cannot be drawn.

R >= 4.1 is assumed. No other software is required.

The results reported in the paper were produced under:

```
R 4.5.2 (2025-10-31 ucrt), x86_64-w64-mingw32/x64, Windows 11
numDeriv  2016.8-1.1
statmod   1.5.1
ggplot2   4.0.3
tidyr     1.3.2
patchwork 1.3.2
progressr 0.18.0
dplyr     1.2.1
```

Nothing here depends on those exact versions, and none of them is pinned.
`run_paper.R`, `run_simulations.R` and `run_all.R` each end by printing
`sessionInfo()`, so any run records its own versions in its own output and the
block above can be checked against it.

## How to run

There are two entry points, and for almost every purpose the first is the one
you want.

```
cd code
Rscript run_paper.R         # everything except the simulations   (~10 min)
Rscript run_simulations.R   # the four Monte Carlo studies        (hours)
```

`run_paper.R` runs the eight validation scripts, Tables 1-4, Figures 1-5, both
applications, and every number the running text quotes. It does not re-run the
Monte Carlo studies of Sections 5.1-5.4. It does not need to: the raw output of
those studies ships with this package in `results/`, and three scripts read it
rather than recompute it. `08_figures_simulation.R` draws Figures 4-5 from it,
`09_values_text.R` takes the coverage averages and replication counts the prose
quotes from it, and `04_sim_bias_mse.R` runs in a reuse mode that reads
`results/sim_bias_mse.csv` and goes straight to Tables 2-4, so those tables and
their captions are rewritten on every run without the study behind them being
repeated. Run `04_sim_bias_mse.R` on its own and it simulates, as before.
Tables 5-8 and 11 are left as they are, since the scripts that write them have
no reuse mode. So `run_paper.R` on its own reproduces the manuscript, and any
table can be checked against the raw output it came from.

`run_simulations.R` rebuilds `results/` from nothing. That is the only reason to
run it. It takes several hours; `05_sim_ci.R` and `06_sim_heavytail.R` dominate
and both run in parallel by default.

One dependency runs from the applications into the simulations, not the other
way: the small-sample block of Section 5.4 mirrors the first application exactly
and reads its shape, scale, limit and group size from the macro file
`10_realdata1_ballbearings.R` writes. So if the application is changed,
`07_sim_smallm.R` has to be re-run. `run_paper.R` checks this for you: it
compares the true index the shipped Section 5.4 block was built at with the one
the current application implies, and says so if they have drifted apart.

`run_all.R` does both, in the right order, from nothing. Several hours.

Individual scripts can also be run on their own; each sources `00_gie_pffc.R`
and sets its own seed, so results do not depend on execution order, with the one
exception noted above (`07` reads what `10` writes).

The scripts write their LaTeX output into `../tables/` and `../figures/`.
Those directories start empty in this repository: the manuscript is not part of
it, and the fragments are regenerated rather than tracked. Copy them into the
manuscript's directory and compile it there:

```
pdflatex paper_QREI_revised.tex     # three times: the table floats need a
                                    # third pass to settle
```

## Files

Every table and every figure in the manuscript is produced by one of these
scripts, and so is every number quoted in the running text. Nothing in the paper
is typed in by hand.

| File | Tables | Figures | Text macros | Runtime |
|---|---|---|---|---|
| `00_gie_pffc.R` | - | - | - | shared definitions only |
| `12_table1_translation.R` | 1 | - | `values_translation.tex` | seconds |
| `01_figures.R` | - | 1, 2, 3 | `values_figures.tex` | seconds |
| `10_realdata1_ballbearings.R` | 9 | - | `values_realdata1.tex` | < 1 min |
| `11_realdata2_guineapig.R` | 10 | - | `values_realdata2.tex` | < 1 min |
| `04_sim_bias_mse.R` | 2, 3, 4 | - | `values_sim.tex` | 20-40 min (seconds in reuse mode) |
| `05_sim_ci.R` | 5, 11 (App. C) | - | - | ~1-3 h, parallel |
| `06_sim_heavytail.R` | 6, 7 | - | `values_heavytail.tex` | ~1-2 h, parallel |
| `07_sim_smallm.R` | 8 | - | `values_smallm.tex` | 30-60 min (needs 10) |
| `08_figures_simulation.R` | - | 4, 5 | - | seconds (needs 04, 05) |
| `09_values_text.R` | - | - | `values_text.tex` | seconds (reads `results/`) |

`run_paper.R` runs every row except `05`, `06` and `07`, and runs `04` in its
reuse mode, which reads `results/sim_bias_mse.csv` instead of simulating.
`run_simulations.R` runs `04`, `05`, `06` and `07`, all four simulating.

### Validation scripts

These do not produce output for the paper; they check that the analytical
results in the paper agree with independent numerical calculations. Each stops
with an error if a check fails.

| File | Checks |
|---|---|
| `validation/v1_check_derivatives.R` | the analytical derivatives in the Remark of Section 3, against `numDeriv` |
| `validation/v2_check_invariance.R` | location-scale invariance of both indexes |
| `validation/v3_check_moment_formula.R` | the closed-form moment of Eq. (3), against quadrature and `integrate()` |
| `validation/v4_check_information.R` | the inequalities and bounds of Appendix B, including heavy-tailed cases with `alpha <= 2` |
| `validation/v5_check_information_and_mle.R` | the closed-form Hessian of Section 3 against `numDeriv`, the explicit MLE of `alpha` given `lambda`, and the fixed-point iteration of Eq. (15) against direct maximisation |
| `validation/v6_check_bootstrap.R` | that the undercoverage of the percentile interval in Table 5 is a property of that interval and not a fault in the code |
| `validation/v7_check_pffc_generator.R` | that `generate_pffc()` reproduces the life test it stands for, against a direct simulation of the experiment |
| `validation/v8_check_score_moments.R` | the score-moment argument of Section 3: the identity `log(1 - e^{-lambda/x}) = -E/(alpha k)` on samples from `generate_pffc()`, the exponential spacings behind it, the two closed forms of the score against `loglik_score()`, and the finiteness of the fourth moments of both components at `alpha` = 0.5, 1 and 2, where the lifetime has no variance |

`v6` answers the one thing in Table 5 that looks like a bug. The percentile
interval covers around 0.85-0.93 while the asymptotic interval sits at nominal,
which is the reverse of the usual expectation. The script reproduces one cell of
the design and decomposes the result: the Monte Carlo bias of the estimate, the
bootstrap's own estimate of that bias, which side each interval misses on, where
each interval is centred, and how often a bootstrap refit fails. It runs through
`bootstrap_index()` and the three `ci_*()` functions themselves rather than
re-implementing them. What comes out is that the estimator is biased at
`m = 25`, that the bootstrap measures the bias accurately, and that the
percentile interval, built from the quantiles of a distribution centred at
`Chat + bias` and hence at roughly `C + 2 * bias`, carries the bias twice, while
`NB` subtracts it and the ACI carries it once. The percentile interval is
*longer* than the ACI and still covers less, so the deficit is a displacement
and not a width. It takes about five minutes and runs as part of `run_paper.R`.

`v7` checks the one function every simulated number depends on. `generate_pffc()`
does not simulate the life test: it uses the equivalence noted by Wu and Kuş
(2009, Sec. 2), by which progressively first-failure-censored order statistics from `F`
are distributed as a progressively type-II censored sample from
`F_k(x) = 1 - (1 - F(x))^k`, which for the GIE is again GIE with the shape
multiplied by `k`, and then generates the type-II sample by the
exponential-spacings form of Balakrishnan and Sandhu (1995). Both steps are
standard, but together they replace the experiment entirely, and a mistake in
either would still produce an ordered sample of the right length. So `v7`
simulates the experiment itself, with `n` groups of `k` units and each failure's group
withdrawn along with `R_i` further groups drawn at random, and compares the two
coordinate by coordinate with a Kolmogorov-Smirnov test, over five designs
covering all four censoring schemes, `alpha <= 2`, and the design of the first
application. It also checks the closed-form marginal of the first failure, that
the group counts reach exactly zero at the `m`-th failure, and that the four
special cases Wu and Kuş list, namely complete sample, first-failure censoring,
progressive type-II censoring and ordinary type-II censoring, all fall out of the
general code. One to two minutes; also part of `run_paper.R`.

## Output

- `results/`: raw numbers as `.csv`. **Tracked in this repository**, because
  the four Monte Carlo studies that produce them take hours and `run_paper.R`
  reads them. They are also the inputs to the figures, so a figure can never
  drift out of step with the table it summarises.
- `../tables/tab_*.tex`: LaTeX table fragments, `\input` by the manuscript.
- `../tables/values_*.tex`: `\newcommand` definitions for every number quoted in
  the running text: the true index values in Section 5.1, the coverage averages
  in Section 5.3, the small-sample coverages in Section 5.4, and all of the
  estimates, standard errors, translation constants and *p*-values in Sections
  6.1-6.3. The manuscript `\input`s these six files in its preamble and refers
  to the numbers by name (`\bbAlphaHat`, `\htCPACIa`, and so on), so the prose
  updates together with the tables and cannot fall out of step with them.
- `../figures/`: figures in both EPS and PDF. The manuscript refers to them
  without an extension, so `pdflatex` picks the PDF and a `latex`/`dvips` route
  picks the EPS.

Re-running any script overwrites its fragments; recompiling the manuscript
(three times, for cross-references and float placement) then picks up the new
numbers everywhere they appear.

`results/` is shipped with this package. The three grid simulations take hours,
and a reader who wants to check a table against the raw output, or to rebuild
the figures and the text macros, should not have to run them first. Deleting the
directory and running `Rscript run_all.R` rebuilds it from nothing.

### Nothing in the manuscript is typed by hand

That claim is meant literally, and it is now true of the running text as well as
of the tables. Every number the prose quotes about a result comes from a macro
that one of these scripts writes: the coverages averaged over the design and the
counts of discarded replications from `09_values_text.R`, the peaks and crossing
values of Figure 3 and the shape at which the median equals the scale from
`01_figures.R`, the censored sample of Section 6.1 and the complete-data fit from
`10_realdata1_ballbearings.R`. The only decimals left in the source are design
constants the author chose, namely the quantile levels 0.25, 0.5 and 0.75, the nominal
level 0.95 and the tolerance 0.01, and equation and section numbers inside
citations. If a simulation is re-run and a number moves, every sentence that
quotes it moves with it.

## Reproducibility notes

- `MASTER_SEED` is defined in `00_gie_pffc.R` and every script derives its own
  seed from it. `05_sim_ci.R` and `06_sim_heavytail.R` assign one seed per
  configuration of the design grid, so their results do not depend on the order
  of traversal and are identical whether run serially or in parallel. The RNG
  kind is named explicitly in those two scripts: `furrr` installs an
  L'Ecuyer-CMRG seed in each worker, and `set.seed()` without a `kind` argument
  would re-seed *that* generator, so a parallel run would silently differ from a
  serial one.
- The bootstrap in both data analyses uses a fixed seed. Bootstrap interval
  endpoints therefore reproduce exactly.
- The censored sample of the first application is **constructed by the script**,
  not typed in. `10_realdata1_ballbearings.R` forms `n` groups of `k` units from
  the 23 endurance times, applies the censoring plan, and writes the resulting
  sample to `../tables/values_realdata1.tex` as `\bbSample`, which the manuscript
  prints. This matters because the construction has to obey the scheme it
  illustrates: with `k = 2` the design consumes `n * k` units, which cannot
  exceed 23, and the first observation must be the smallest lifetime among all
  units on test, since every observation is a group minimum. The sample used in
  the original submission satisfied neither condition: it began at 41.52 when
  the smallest observation is 17.88, and its plan `R = (3, 0^8)` called for 24
  units. `07_sim_smallm.R` now reads the design back from that macro file
  through `read_macros()` instead of repeating the constants, so the two cannot
  drift apart.
- `04` to `07` report the number of replications in which the fit converged and, for
  the asymptotic interval, the number in which the observed information was
  positive definite. These counts are the numerical check on nonsingularity
  referred to in Section 3 and Appendix B.
- The simulation settings in the scripts are the ones stated in the paper:
  2000 replications for point estimation, 500 for intervals, and `B = 250`
  bootstrap resamples (4000 / 2000 / 250 in `07_sim_smallm.R`, where `m` is
  small and the extra replications are cheap).

## Every closed form in the paper is used, not just stated

The paper displays a number of closed-form results. Each is implemented and sits
on the path that produces the reported numbers; none is present only as a
side-check. `00_gie_pffc.R` opens with the full correspondence, and the four
that matter most are:

- **Eq. (3), the closed-form moment for integer shape.** `gie_moment()`
  dispatches exactly as Section 2.1 describes: the closed form when `alpha` is a
  positive integer with `r < alpha`, the quadrature of Eq. (4) otherwise. The
  true index values that define the simulation targets therefore go through the
  closed form, while every estimate goes through the quadrature, since `alpha`
  estimates are never exactly integer. The closed form is capped at
  `alpha <= 20`: see `ALPHA_EXACT_MAX` for why.
- **Eq. (15), the fixed-point iteration for the MLE of lambda**, together with
  the explicit expression for `alpha_hat` given `lambda`. This is the *primary*
  estimation route in `fit_mle()`, not a footnote. A direct maximisation is used
  only where the iteration fails to converge or does not land on an interior
  maximum; the `route` field records which was used, and `04_sim_bias_mse.R`
  reports the proportion of fits the published iteration handled on its own.
- **The Hessian displayed in Section 3.** `loglik_hessian()` and
  `observed_information()` evaluate those three second derivatives directly, and
  every standard error in the paper comes from them. `optim`'s numerical Hessian
  is computed only as a cross-check; the two agree to five or six significant
  figures.
- **The derivatives in the Remark of Section 3.** `grad_CL_quantile()` is the
  analytical gradient used by the delta method for the quantile-based index.
  Numerical differentiation is used for one thing only, the gradient of the
  moment-based index, whose value depends on quadrature, which is exactly the
  fallback the Remark itself describes.

The `loglik_score()` function is the same idea one step earlier: the estimating
equations before they are rearranged. `fit_mle()` uses it to confirm that the
point returned by the fixed-point iteration really is a stationary point.

## Estimation details

The maximum likelihood estimates are obtained by profiling `alpha` out
analytically,

```
alpha_hat(lambda) = -m / (k * sum((R_i + 1) * log(1 - exp(-lambda/x_i))))
```

and maximising the resulting one-dimensional profile likelihood over a
log-spaced grid in `lambda`. That grid maximum is not the estimate; it is the
starting value for the fixed-point iteration of Eq. (15), which is the primary
route and the one that converges in the great majority of fits. When the
iteration fails to converge, `optim(..., method = "L-BFGS-B", hessian = TRUE)`
is used as a fallback from the same starting value. `fit_mle()` records which
route was taken in its `route` element, and the scripts report how often the
fixed-point route sufficed.

Standard errors do not come from `optim`. They are computed from the analytic
observed information of Eq. (17), `observed_information()`, on whichever
parameter value the estimation returned, combined with the analytical gradient
of the index through the delta method; the gradient of the moment-based index,
which involves numerical quadrature, is obtained with `numDeriv::grad`. The
numerical Hessian from `optim` is retained only when `numeric_hessian = TRUE`,
and only so that the two can be compared: the applications print the largest
relative difference between them, and `validation/v5_check_information_and_mle.R`
checks the analytic form against numerical differentiation directly.

## Changes relative to the scripts used for the original submission

The original scripts have been reorganised into the structure above. Four
substantive corrections were made along the way; they are listed here so that
any change in a reported number can be traced.

1. **Quadrature.** The moment integral was previously evaluated by mapping
   `x = t/(1-t)` and applying Gauss-Legendre. That map does not carry `lambda`,
   so accuracy degraded as `lambda` grew, and it leaves an endpoint singularity
   of order `alpha-1-r`, so the rule converged only algebraically when `alpha`
   was close to `r`: the error at `alpha = 2.1`, `r = 2` was about 40%. The
   substitutions documented in `gie_moment_gl()` remove both problems. This
   affects the left edge of Figure 3, the `alpha = 2.5` row of Table 1, and
   `C_L^M` in Section 6.1 in its fourth decimal.
2. **Simulation settings.** The confidence-interval script had been left at
   `ds <- 20`, `B <- 20` after a debugging run and did not reproduce the
   published tables; it is now at the 500 / 250 stated in the paper.
3. **Seeds.** The two data-analysis scripts called `set.seed(NULL)`, so their
   bootstrap intervals were not reproducible. All seeds are now fixed and
   derived from `MASTER_SEED`.
4. **Figures from data.** Figures 4 and 5 previously had their numbers typed in
   by hand from the tables, rounded to two decimals in the case of Figure 5.
   They now read `results/*.csv`.

Two further changes are numerical hygiene rather than corrections:
`log(1 - exp(z))` is computed as `log(-expm1(z))` throughout, since the
cancellation as `z -> 0` was severe enough to drop replications silently; and
the maximum likelihood starting value is now obtained by profiling `alpha` out
and maximising the one-dimensional profile likelihood on a grid, rather than
being supplied by hand.
