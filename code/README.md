# Reproducibility material

**A General Quantile-Based Lifetime Performance Index and its Application to the GIE Distribution under Progressive First-Failure Censoring**

This directory contains all `R` code needed to reproduce every table, figure and
numerical statement in the paper.

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

R ≥ 4.1 is assumed. No other software is required.

## How to run

```
cd code
Rscript run_all.R quick    # validation + figures 1-3 + both data analyses (minutes)
Rscript run_all.R          # everything, including the simulations (hours)
```

Individual scripts can also be run on their own; each sources `00_gie_pffc.R`
and sets its own seed, so results do not depend on execution order.

After a full run, recompile the manuscript from the parent directory:

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
| `00_gie_pffc.R` | — | — | — | shared definitions only |
| `12_table1_translation.R` | 1 | — | `values_translation.tex` | seconds |
| `01_figures.R` | — | 1, 2, 3 | — | seconds |
| `10_realdata1_ballbearings.R` | 9 | — | `values_realdata1.tex` | < 1 min |
| `11_realdata2_guineapig.R` | 10 | — | `values_realdata2.tex` | < 1 min |
| `04_sim_bias_mse.R` | 2, 3, 4 | — | `values_sim.tex` | 20–40 min |
| `05_sim_ci.R` | 5, 11 (App. C) | — | — | ~1–3 h, parallel |
| `06_sim_heavytail.R` | 6, 7 | — | `values_heavytail.tex` | ~1–2 h, parallel |
| `07_sim_smallm.R` | 8 | — | `values_smallm.tex` | 30–60 min |
| `08_figures_simulation.R` | — | 4, 5 | — | seconds (needs 04, 05) |

The first four rows are the `quick` set; the rest are the simulations.

### Validation scripts

These do not produce output for the paper; they check that the analytical
results in the paper agree with independent numerical calculations. Each stops
with an error if a check fails.

| File | Checks |
|---|---|
| `validation/v1_check_derivatives.R` | the analytical derivatives in the Remark of Section 3, against `numDeriv` |
| `validation/v2_check_invariance.R` | location–scale invariance of both indexes |
| `validation/v3_check_moment_formula.R` | the closed-form moment of Eq. (3), against quadrature and `integrate()` |
| `validation/v4_check_information.R` | the inequalities and bounds of Appendix B, including heavy-tailed cases with `alpha <= 2` |
| `validation/v5_check_information_and_mle.R` | the closed-form Hessian of Section 3 against `numDeriv`, the explicit MLE of `alpha` given `lambda`, and the fixed-point iteration of Eq. (15) against direct maximisation |

## Output

- `results/` — raw numbers as `.csv` / `.rds`. These are the inputs to the
  figures, so a figure can never drift out of step with the table it summarises.
- `../tables/tab_*.tex` — LaTeX table fragments, `\input` by the manuscript.
- `../tables/values_*.tex` — `\newcommand` definitions for every number quoted in
  the running text: the true index values in Section 5.1, the coverage averages
  in Section 5.3, the small-sample coverages in Section 5.4, and all of the
  estimates, standard errors, translation constants and *p*-values in Sections
  6.1–6.3. The manuscript `\input`s these six files in its preamble and refers
  to the numbers by name (`\bbAlphaHat`, `\htCPACIa`, and so on), so the prose
  updates together with the tables and cannot fall out of step with them.
- `../figures/` — figures in both EPS and PDF. The manuscript refers to them
  without an extension, so `pdflatex` picks the PDF and a `latex`/`dvips` route
  picks the EPS.

Re-running any script overwrites its fragments; recompiling the manuscript
(three times, for cross-references and float placement) then picks up the new
numbers everywhere they appear.

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
- `04`–`07` report the number of replications in which the fit converged and, for
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
  closed form, while every estimate — `alpha` estimates are never exactly
  integer — goes through the quadrature. The closed form is capped at
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
  Numerical differentiation is used for one thing only — the gradient of the
  moment-based index, whose value depends on quadrature — which is exactly the
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

maximising the resulting one-dimensional profile likelihood over a log-spaced
grid in `lambda`, and using that as the starting value for
`optim(..., method = "L-BFGS-B", hessian = TRUE)`. Standard errors come from the
observed information returned by `optim`, combined with the analytical gradient
of the index through the delta method; the gradient of the moment-based index,
which involves numerical quadrature, is obtained with `numDeriv::grad`.

## Changes relative to the scripts used for the original submission

The original scripts have been reorganised into the structure above. Four
substantive corrections were made along the way; they are listed here so that
any change in a reported number can be traced.

1. **Quadrature.** The moment integral was previously evaluated by mapping
   `x = t/(1-t)` and applying Gauss–Legendre. That map does not carry `lambda`,
   so accuracy degraded as `lambda` grew, and it leaves an endpoint singularity
   of order `alpha-1-r`, so the rule converged only algebraically when `alpha`
   was close to `r` — the error at `alpha = 2.1`, `r = 2` was about 40%. The
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

## The marked-up manuscript

`../paper_QREI_diff.tex` is the revised manuscript with everything new or
changed printed in red, for the referee. It is produced by `../make_markup.py`:

```
cd ..
python make_markup.py paper_QREI.tex      # the original submission
pdflatex paper_QREI_diff.tex              # three times
```

`latexdiff` is the usual tool for this, but two things about this manuscript
make a direct comparison useless without preprocessing: the revision refers to
its numbers through macros (`\bbAlphaHat` and the rest) and pulls its tables in
with `\input`, so a naive diff would compare macro names against numbers and
report the whole of every table as changed. `make_markup.py` resolves both
first, then marks at the block level, refining to word level inside ordinary
paragraphs. Deleted material is not shown; only additions and changes are
marked, which is what the marked-up copy is for.
