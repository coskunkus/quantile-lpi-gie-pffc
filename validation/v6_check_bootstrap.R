################################################################################
##  validation/v6_check_bootstrap.R
##
##  Why the percentile bootstrap undercovers while the ACI does not.
##
##  Table 5 shows the percentile interval (PB) covering well below nominal --
##  around 0.85-0.93 -- while the asymptotic interval (ACI) sits at nominal and
##  the bias-corrected normal bootstrap (NB) sits at or slightly above it.  That
##  ordering looks wrong at first sight, and the obvious suspicion is a coding
##  error in the bootstrap.  It is not one.  This script establishes that by
##  measuring, in a single configuration of the Table 5 design, the four
##  quantities that decide the matter:
##
##    (a) the Monte Carlo bias of Chat, and the bootstrap's estimate of it;
##    (b) which side each interval misses on;
##    (c) where each interval is centred relative to the true index;
##    (d) how often a bootstrap refit fails, which is what would make the
##        bootstrap distribution selectively truncated.
##
##  The explanation the numbers give is the textbook one.  Chat is biased; the
##  bootstrap measures that bias accurately; and the percentile interval is
##  built from the quantiles of a bootstrap distribution centred at Chat + bias,
##  which relative to the true C sits at roughly C + 2 * bias.  The bias
##  therefore enters the percentile interval twice.  NB subtracts the estimated
##  bias explicitly and so removes it; the ACI carries it once.  The percentile
##  interval is longer than the ACI, not shorter, so its lower coverage is a
##  displacement and not a width problem.
##
##  The interval functions exercised here -- `bootstrap_index`, `ci_percentile`,
##  `ci_normal_boot`, `ci_asymptotic` -- are the ones `05_sim_ci.R` calls to
##  build Table 5, not re-implementations of them, so a fault in any of them
##  would show up here.
##
##  Runtime: about five minutes on one core.
################################################################################

source("../00_gie_pffc.R")
need("numDeriv", "statmod")

set.seed(MASTER_SEED + 106L)

## One cell of the Table 5 design: the m = 25, Early cell, where the gap between
## the two coverages is widest.  Fewer replications than the paper uses -- this
## is a diagnostic, not a table -- but the effects it has to resolve are large
## compared with the Monte Carlo error at this size.
M       <- 25
K       <- 2
ALPHA   <- 5
LAMBDA  <- 2.0
L_LIMIT <- 2.0
SCHEME  <- "Early"
N_REP   <- 200
B_BOOT  <- 150
LEVEL   <- 0.95

R      <- make_scheme(M, SCHEME)
tp     <- c(ALPHA, LAMBDA)
C_true <- CL_quantile(tp, L_LIMIT)

cat(sprintf("Design: m=%d, k=%d, alpha=%g, lambda=%g, L=%g, %s scheme\n",
            M, K, ALPHA, LAMBDA, L_LIMIT, SCHEME))
cat(sprintf("True index C_L^xi = %.4f;  %d replications, B = %d\n\n",
            C_true, N_REP, B_BOOT))

Chat   <- rep(NA_real_, N_REP)
se_dm  <- rep(NA_real_, N_REP)
b_bias <- rep(NA_real_, N_REP)      # bootstrap estimate of the bias
b_sd   <- rep(NA_real_, N_REP)
lo     <- hi <- matrix(NA_real_, N_REP, 4,
                       dimnames = list(NULL, c("ACI", "PB", "NB", "basic")))
n_fail <- 0L                        # bootstrap refits that did not converge
n_try  <- 0L

for (i in seq_len(N_REP)) {
  if (i %% 25 == 0) { cat(sprintf("  %d/%d\r", i, N_REP)); utils::flush.console() }
  x  <- generate_pffc(M, K, R, ALPHA, LAMBDA)
  ft <- fit_mle(x, R, K)
  if (is.null(ft)) next
  Chat[i]  <- CL_quantile(ft$par, L_LIMIT)
  se_dm[i] <- delta_se(ft, L_LIMIT, "quantile")

  ## The bootstrap of Section 4, through the same function `05_sim_ci.R` uses:
  ## resample from the fitted parameters, keeping the group size and censoring
  ## plan, and refit.
  bq     <- bootstrap_index(ft$par, M, K, R, L_LIMIT, B_BOOT, index = "quantile")
  n_try  <- n_try + B_BOOT
  n_fail <- n_fail + sum(!is.finite(bq))
  if (sum(is.finite(bq)) < B_BOOT / 2 || !is.finite(Chat[i])) next

  b_bias[i] <- mean(bq, na.rm = TRUE) - Chat[i]
  b_sd[i]   <- sd(bq, na.rm = TRUE)

  ci <- ci_asymptotic(Chat[i], se_dm[i], LEVEL);      lo[i, "ACI"] <- ci[1]; hi[i, "ACI"] <- ci[2]
  ci <- ci_percentile(bq, LEVEL);                     lo[i, "PB"]  <- ci[1]; hi[i, "PB"]  <- ci[2]
  ci <- ci_normal_boot(Chat[i], bq, LEVEL);           lo[i, "NB"]  <- ci[1]; hi[i, "NB"]  <- ci[2]
  ## The basic (pivotal) interval is not used in the paper.  It is computed here
  ## only because it is the percentile interval with the bias reflected out, so
  ## the difference between the two isolates the quantity under discussion.
  qq <- ci_percentile(bq, LEVEL)
  lo[i, "basic"] <- 2 * Chat[i] - qq[2]; hi[i, "basic"] <- 2 * Chat[i] - qq[1]
}
cat(strrep(" ", 20), "\r", sep = "")

ok       <- is.finite(Chat)
mc_bias  <- mean(Chat[ok]) - C_true
mc_sd    <- sd(Chat[ok])
se_bias  <- mc_sd / sqrt(sum(ok))
boot_b   <- mean(b_bias, na.rm = TRUE)
fail_pct <- 100 * n_fail / max(1L, n_try)

## ---------------------------------------------------------------------------
## (a) does the bootstrap measure the bias correctly?
## ---------------------------------------------------------------------------
cat("(a) Bias of the estimator\n")
cat(sprintf("    Monte Carlo bias of Chat        %+.4f   (MC standard error %.4f)\n",
            mc_bias, se_bias))
cat(sprintf("    bootstrap estimate of the bias  %+.4f\n", boot_b))
cat(sprintf("    Monte Carlo sd of Chat          %.4f\n", mc_sd))
cat(sprintf("    mean delta-method se            %.4f\n", mean(se_dm[ok], na.rm = TRUE)))
cat(sprintf("    mean bootstrap sd               %.4f\n\n", mean(b_sd, na.rm = TRUE)))

## ---------------------------------------------------------------------------
## (b), (c) coverage, length, side of the miss, and centre
## ---------------------------------------------------------------------------
cat("(b) Coverage, and which side each interval misses on\n")
summ <- function(nm) {
  l <- lo[, nm]; h <- hi[, nm]
  m_ <- is.finite(l) & is.finite(h)
  c(CP    = mean(l[m_] < C_true & h[m_] > C_true),
    AL    = mean(h[m_] - l[m_]),
    ## "entirely below" means the whole interval lies below the true value, so
    ## the procedure has overstated the index; "entirely above" is the reverse.
    below = mean(h[m_] <= C_true),
    above = mean(l[m_] >= C_true),
    ctr   = mean(0.5 * (l[m_] + h[m_])) - C_true)
}
tab <- t(sapply(c("ACI", "PB", "NB", "basic"), summ))
cat("           CP      AL     entirely   entirely   centre\n")
cat("                          below      above      minus C\n")
for (nm in rownames(tab))
  cat(sprintf("    %-6s %.3f  %.4f    %.3f      %.3f     %+.4f\n",
              nm, tab[nm, "CP"], tab[nm, "AL"], tab[nm, "below"],
              tab[nm, "above"], tab[nm, "ctr"]))
cat("\n")

## ---------------------------------------------------------------------------
## (d) selective loss of bootstrap replicates
## ---------------------------------------------------------------------------
cat("(d) Bootstrap refits that failed to converge\n")
cat(sprintf("    %.2f%% of %d refits\n\n", fail_pct, n_try))

## ---------------------------------------------------------------------------
## The checks
## ---------------------------------------------------------------------------
fails <- character(0)
chk <- function(cond, msg) if (!isTRUE(cond)) fails <<- c(fails, msg)

## 1. The bootstrap's bias estimate must track the Monte Carlo bias.  If it did
##    not, the resampling itself would be wrong.
chk(abs(boot_b - mc_bias) < 0.25 * abs(mc_bias) + 3 * se_bias,
    sprintf("bootstrap bias estimate %+.4f does not track the Monte Carlo bias %+.4f",
            boot_b, mc_bias))

## 2. The bootstrap distribution must not be selectively truncated: refits have
##    to converge essentially always, or the percentile quantiles would be
##    measuring convergence rather than sampling variability.
chk(fail_pct < 2,
    sprintf("%.2f%% of bootstrap refits failed; the percentile quantiles are truncated",
            fail_pct))

## 3. PB must be *no shorter* than the ACI while covering less.  This is the
##    point of the whole script: the deficit is a displacement, not a width.
chk(tab["PB", "AL"] >= tab["ACI", "AL"],
    "the percentile interval is shorter than the ACI, so its undercoverage is a width problem")
chk(tab["PB", "CP"] < tab["ACI", "CP"],
    "the percentile interval did not undercover relative to the ACI in this run")

## 4. PB must miss overwhelmingly on one side, and that side must be the one the
##    sign of the bias predicts: a downward-biased Chat drags the interval down,
##    so the interval falls entirely below the true value.
pb_far  <- if (mc_bias < 0) tab["PB", "below"] else tab["PB", "above"]
pb_near <- if (mc_bias < 0) tab["PB", "above"] else tab["PB", "below"]
chk(pb_far - pb_near > 0.04,
    sprintf("the percentile interval's misses were not one-sided (%.3f vs %.3f)",
            pb_far, pb_near))

## 5. NB, which subtracts the estimated bias, must be centred much closer to the
##    truth than PB, which does not.
chk(abs(tab["NB", "ctr"]) < 0.5 * abs(tab["PB", "ctr"]),
    sprintf("the bias-corrected interval is not better centred than the percentile interval (%+.4f vs %+.4f)",
            tab["NB", "ctr"], tab["PB", "ctr"]))

## 6. The displacement of PB must be larger than one bias, which is what
##    distinguishes "the bias enters twice" from "the bias enters once".  The
##    factor exceeds two in practice, the excess being the skewness of the
##    bootstrap distribution, so 1.5 is a conservative threshold.
ratio <- tab["PB", "ctr"] / mc_bias
chk(ratio > 1.5,
    sprintf("the percentile interval is displaced by %.2f biases, not more than one", ratio))

cat(sprintf("Percentile interval displaced by %.2f x the bias; ACI by %.2f; NB by %.2f.\n\n",
            ratio, tab["ACI", "ctr"] / mc_bias, tab["NB", "ctr"] / mc_bias))

if (length(fails) == 0) {
  cat("PASS: the bootstrap is implemented correctly and the undercoverage of the\n",
      "percentile interval is a property of that interval, not a defect in the\n",
      "code.  Chat is biased; the bootstrap measures the bias accurately; the\n",
      "percentile interval inherits it twice over, the ACI once, and NB not at\n",
      "all, which is the ordering of coverages seen in Table 5.\n", sep = "")
} else {
  for (f in fails) cat("  - ", f, "\n", sep = "")
  stop("FAIL: the bootstrap diagnostic did not behave as the explanation requires.")
}
