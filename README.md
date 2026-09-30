# A General Quantile-Based Lifetime Performance Index: reproduction code

`R` code reproducing every table, figure and reported number in

> Alkhuffash, O. and Kuş, C. *A General Quantile-Based Lifetime Performance
> Index and its Application to the GIE Distribution under Progressive
> First-Failure Censoring.* Submitted to the *South African Journal of Science*.

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
Rscript run_all.R
```

About an hour and a quarter on the machine the paper was produced on, almost all
of it in the four Monte Carlo studies. That runs the two data analyses, the
eight validation scripts, all eleven tables, Figures 1-5, and every number the
running text quotes.

## Requirements

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
`run_all.R` ends by printing `sessionInfo()`, so any run records its own
versions in its own output and the block above can be checked against it.

## One entry point, and why

`run_all.R` is the only entry point, and it is deliberate. No script in this
package reads a stored result: not a LaTeX file, not a csv, not an rds. Every
number in the paper is therefore the output of a computation carried out in the
run that reports it, and never of a file being reread. Two things follow.

The four Monte Carlo studies of Sections 5.1-5.4 are always run, so a
reproduction takes over an hour; the csv files in `results/` are a record of what
a run produced, not an input to the next one. Everything is seeded from
`MASTER_SEED`, and the two parallel studies assign one seed per configuration, so
the same numbers come back whether the run is serial or parallel.

The steps run in one session and in a fixed order, because three of them need
what earlier ones computed. `08_figures_simulation.R` draws Figures 4 and 5 from
the two baseline studies; `09_values_text.R` takes from all four studies every
number the running text quotes about a simulated result, and from the second
application the three interval lengths Section 6.2 compares. Each is handed the
objects the earlier steps left in the session, and each stops with an
explanation if it is run on its own.

The order is: the two applications, the eight validation scripts, `12`, `01`,
then `04` to `07`, then `08` and `09`. The applications come first because four
of the validation scripts check the design of the first one.

That design is the one thing that crosses between steps without going through
the session. `07_sim_smallm.R` simulates at it, and `v4`, `v5`, `v7` and `v8`
include it among the designs they check. None of them reads it back either: it
is rebuilt from the 23 published endurance times by `ball_bearing_design()` in
`00_gie_pffc.R`, which takes milliseconds and gives every caller the same sample
and the same fit. `run_all.R` compares the small-sample study with the
application once more at the end, and says so if they have drifted apart.

## What is here

```
00_gie_pffc.R         the model: density, quantiles, moments, likelihood,
                      MLE, delta method, bootstrap, PFFC generator
01,10,11,12_*.R       figures and the two data analyses
04,05,06,07_sim_*.R   the four Monte Carlo studies
08,09_*.R             figures and text macros, from the studies in memory
run_all.R             the entry point: everything, in order, in one session
validation/           eight scripts checking the analytical results in the
                      paper against independent numerical calculations
results/              a record of what the run produced, as csv
tables/               written by the scripts: LaTeX table fragments and
figures/              \newcommand macro files, and the figures
```

`tables/` and `figures/` start empty: the manuscript is not part of this
repository, and the fragments are regenerated rather than tracked. Copy them
into the manuscript's directory and compile it there with `pdflatex`, three
times: the table floats need a third pass to settle.

## Files

Every table and every figure in the manuscript is produced by one of these
scripts, and so is every number quoted in the running text. The submitted
manuscript carries them inline rather than `\input`-ing them, so that it
compiles on its own, but each can be compared with the fragment or the macro
that produced it.

| File | Tables | Figures | Text macros | Runtime |
|---|---|---|---|---|
| `00_gie_pffc.R` | - | - | - | shared definitions only |
| `12_table1_translation.R` | 1 | - | `values_translation.tex` | seconds |
| `01_figures.R` | - | 1, 2, 3 | `values_figures.tex` | seconds |
| `10_realdata1_ballbearings.R` | 9 | - | `values_realdata1.tex` | < 1 min |
| `11_realdata2_guineapig.R` | 10 | - | `values_realdata2.tex` | < 1 min |
| `04_sim_bias_mse.R` | 2, 3, 4 | - | `values_sim.tex` | 2 min |
| `05_sim_ci.R` | 5, 11 (App. C) | - | - | 52 min, parallel |
| `06_sim_heavytail.R` | 6, 7 | - | `values_heavytail.tex` | 9 min, parallel |
| `07_sim_smallm.R` | 8 | - | `values_smallm.tex` | 11 min (needs the design of 10) |
| `08_figures_simulation.R` | - | 4, 5 | - | seconds (needs 04 and 05 in the session) |
| `09_values_text.R` | - | - | `values_text.tex` | seconds (needs 04-07 and 11 in the session) |

`run_all.R` runs every row, in the order given above, in one session.

### Validation scripts

These eight scripts do not produce output for the paper; they check that the
analytical results in the paper agree with independent numerical calculations.
Each stops with an error if a check fails.

| File | Checks |
|---|---|
| `validation/v1_check_derivatives.R` | the analytical derivatives in the Remark of Section 3, against `numDeriv` |
| `validation/v2_check_invariance.R` | location-scale invariance of both indexes |
| `validation/v3_check_moment_formula.R` | the closed-form moment of Eq. (3), against quadrature and `integrate()` |
| `validation/v4_check_information.R` | the inequalities and bounds of Appendix B, including heavy-tailed cases with `alpha <= 2` |
| `validation/v5_check_information_and_mle.R` | the closed-form Hessian of Section 3 against `numDeriv`, the explicit MLE of `alpha` given `lambda`, and the fixed-point iteration of Eq. (16) against direct maximisation |
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
and not a width. It takes about half a minute and runs as part of `run_all.R`.

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
general code. Ten to twenty seconds; also part of `run_all.R`.

## Output

### Nothing here reads a stored result

`tables/`, `figures/` and `results/` hold output only. No script reads a LaTeX
file, a csv or an rds back, so nothing computed can depend on a file the
manuscript may since have been edited into, and no table can be a reprint of an
older run. What one step passes to the next it passes in memory, within the one
run; the only design that crosses otherwise is rebuilt from the published data
by `ball_bearing_design()`.

- `results/`: raw numbers as `.csv`, written by the four studies and by the two
  applications. **Tracked in this repository**, as a record of the run the paper
  reports, so that a reader can compare a table with the output it came from.
  Nothing reads them. Deleting the directory changes nothing, and
  `Rscript run_all.R` rewrites it.
- `tables/tab_*.tex`: the eleven LaTeX table fragments. These are the eleven
  tables of the manuscript and can be compared with it line by line.
- `tables/values_*.tex`: `\newcommand` definitions for every number quoted in
  the running text: the true index values in Section 5.1, the coverage averages
  in Sections 5.2 and 5.3, the small-sample coverages in Section 5.4, the counts
  behind every "in all N configurations" the prose asserts, and all of the
  estimates, standard errors, translation constants and *p*-values in Sections
  6.1-6.3. Eight such files are written, and every number in them is named
  (`\bbAlphaHat`, `\ctCPACIa`, and so on), so each sentence in the paper can be
  checked against the output it describes.
- `figures/`: figures in both EPS and PDF. The manuscript refers to them without
  an extension, so `pdflatex` picks the PDF and a `latex`/`dvips` route picks
  the EPS.

Re-running any script overwrites its fragments; recompiling the manuscript
(three times, for cross-references and float placement) then picks up the new
numbers everywhere they appear.

### Every number in the manuscript comes from here

The submitted manuscript carries its tables and its numbers inline rather than
`\input`-ing them, so that it compiles on its own. Every one of those numbers is
nevertheless produced here and can be checked against the file that produces it:
the coverages averaged over the design, the counts of discarded replications and
the counts behind every design-wide ordering from `09_values_text.R`; the peaks
and crossing values of Figure 3 and the shape at which the median equals the
scale from `01_figures.R`; the censored sample of Section 6.1 and the
complete-data fit from `10_realdata1_ballbearings.R`. The only decimals chosen by
the authors are design constants: the quantile levels 0.25, 0.5 and 0.75, the
nominal level 0.95, the one-sided level 0.05 and the standard
$C_{0}^{\xi} = 0.5$ of Section 6.3, and the two specification limits. If a
simulation is re-run and a number moves, `run_all.R` rewrites every table and
macro that quotes it, and the difference from the paper is then visible.

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
- `run_all.R` prints `sessionInfo()` when it finishes.
- The censored sample of the first application is **constructed, not typed in**.
  `ball_bearing_design()` in `00_gie_pffc.R` forms `n` groups of `k` units from
  the 23 endurance times, applies the censoring plan and fits the model; the
  analysis, `07_sim_smallm.R` and the four validation scripts that include this
  design all call it, so none of them can drift apart from the others. This
  matters because the construction has to obey the scheme it illustrates: with
  `k = 2` the design consumes `n * k` units, which cannot exceed 23, and the
  first observation must be the smallest lifetime among all units on test, since
  every observation is a group minimum. The sample used in the original
  submission satisfied neither condition: it began at 41.52 when the smallest
  observation is 17.88, and its plan `R = (3, 0^8)` called for 24 units.
- The censored sample of the **second** application is the one exception: it is
  given rather than constructed. It satisfies the scheme, but the draw that
  produced it is not recorded, and `build_pffc_sample()` does not reproduce it
  under this package's seed. It is the one construction here a reader cannot
  repeat, and `guinea_pig_design()` says so.
- `04` to `07` report the number of replications in which the fit converged and,
  for the asymptotic interval, the number in which the observed information was
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
  positive integer with `r < alpha`, the quadrature of Eq. (6) otherwise. The
  true index values that define the simulation targets therefore go through the
  closed form, while every estimate goes through the quadrature, since `alpha`
  estimates are never exactly integer. The closed form is capped at
  `alpha <= 20`: see `ALPHA_EXACT_MAX` for why.
- **Eq. (16), the fixed-point iteration for the MLE of lambda**, together with
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
starting value for the fixed-point iteration of Eq. (16), which is the primary
route and the one that converges in the great majority of fits. When the
iteration fails to converge, `optim(..., method = "L-BFGS-B", hessian = TRUE)`
is used as a fallback from the same starting value. `fit_mle()` records which
route was taken in its `route` element, and the scripts report how often the
fixed-point route sufficed.

Standard errors do not come from `optim`. They are computed from the analytic
observed information of Eq. (18), `observed_information()`, on whichever
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
   They are now drawn from the studies themselves, which
   `08_figures_simulation.R` takes from the session that produced them.

Two further changes are numerical hygiene rather than corrections:
`log(1 - exp(z))` is computed as `log(-expm1(z))` throughout, since the
cancellation as `z -> 0` was severe enough to drop replications silently; and
the maximum likelihood starting value is now obtained by profiling `alpha` out
and maximising the one-dimensional profile likelihood on a grid, rather than
being supplied by hand.

## Changes made for the revision

1. **One entry point, and nothing read back.** The package had two runners and
   a reuse mode: `run_paper.R` rebuilt the tables in two minutes by reading
   `results/*.csv` instead of simulating. A table produced that way is a
   reprint of an earlier run, not a reproduction, and several scripts took the
   design of the first application from a file as well. The reuse mode and both
   old runners are gone. `run_all.R` is the only entry point; it runs every
   step in one session and in order, the studies are always run, and the three
   steps that need an earlier result are handed it in memory.
2. **The one-sided coverage in Table 8.** `07_sim_smallm.R` recorded the
   coverage of the one-sided lower bound but did not print it. It is now a
   column of Table 8, because it is that bound, and not the two-sided interval
   beside it, that the decision rule of Section 6.3 uses.
3. **The translated threshold.** An earlier draft of Section 6.3 also tested
   against a standard inherited on the moment scale, carrying the uncertainty
   of the translated critical value through a scalar `G`. That calculation was
   dropped from the paper, and it has been removed from
   `10_realdata1_ballbearings.R` with the macros that served it. Section 2.3
   still reports `kappa` and `delta`, and `12_table1_translation.R` still
   verifies the affine map.
4. **Quadrature accuracy.** The header of `gie_moment_gl()` claimed a relative
   error of 1e-5 at 100 nodes and 1e-6 at 400. The measured values are 7.5e-5
   and 2.9e-6; the comment now states the two accuracies the paper states, and
   `validation/v3_check_moment_formula.R` tests exactly those.
5. **The design of the first application is rebuilt, not read.**
   `07_sim_smallm.R` and four of the validation scripts took it from the macro
   file `10_realdata1_ballbearings.R` writes into `tables/`. That made a
   computed result depend on a directory whose contents belong to the
   manuscript, and it failed outright on a fresh checkout, where `tables/` is
   empty. The construction now lives in `ball_bearing_design()` in
   `00_gie_pffc.R` and every caller rebuilds it from the 23 published endurance
   times, in milliseconds; `read_macros()` has been removed.

## License

MIT, see `LICENSE`.
