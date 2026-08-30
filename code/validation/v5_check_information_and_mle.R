################################################################################
##  validation/v5_check_information_and_mle.R
##
##  Two claims are checked here, both of which concern closed forms that the
##  paper displays and the code now uses on the main path rather than in a side
##  calculation.
##
##  1.  The Hessian of Section 3 -- the three second derivatives displayed
##      immediately after Eq. (17) -- against numerical differentiation of the
##      log-likelihood.  These expressions supply the observed information used
##      for every standard error in the paper, so they had better be right.
##
##  2.  The estimation procedure of Section 3: the closed-form MLE of alpha
##      given lambda, and the fixed-point iteration of Eq. (15) for lambda.  We
##      check that the iteration, where it converges, lands on the same point as
##      a direct maximisation, and record how often it converges.
################################################################################

source("../00_gie_pffc.R")
need("numDeriv", "statmod")

set.seed(MASTER_SEED + 105L)

## ---------------------------------------------------------------------------
## 1.  The closed-form Hessian
## ---------------------------------------------------------------------------

cat("Checking the closed-form Hessian of Section 3 against numDeriv ...\n")

app_fit <- read_macros(file.path(TABLES_DIR, "values_realdata1.tex"))
designs <- list(
  ## the design of the first application, read from the macro file that
  ## 10_realdata1_ballbearings.R writes rather than copied here
  list(m = as.numeric(app_fit$bbM), k = as.numeric(app_fit$bbK),
       alpha = as.numeric(app_fit$bbAlphaHat),
       lambda = as.numeric(app_fit$bbLambdaHat), sch = "app",
       R1 = as.numeric(app_fit$bbRfirst)),
  list(m = 25,  k = 2, alpha = 5,    lambda = 2,     sch = "Early"),
  list(m = 25,  k = 5, alpha = 1,    lambda = 2,     sch = "Late"),
  list(m = 50,  k = 2, alpha = 2,    lambda = 2,     sch = "Middle"),
  list(m = 100, k = 5, alpha = 10,   lambda = 2,     sch = "Equal"))

worst_H <- 0
for (d in designs) {
  R <- if (identical(d$sch, "app")) c(d$R1, rep(0, d$m - 1))
         else make_scheme(d$m, d$sch)
  for (i in 1:40) {
    x   <- generate_pffc(d$m, d$k, R, d$alpha, d$lambda)
    par <- c(d$alpha, d$lambda)
    ana <- loglik_hessian(par, x, R, d$k)
    ## Pass the data by closure, not through `...`: numDeriv::hessian has a
    ## formal argument named `x`, so `hessian(neg_loglik, par, x = x, ...)`
    ## binds the data to that formal and shifts `par` into `method`.
    num <- -numDeriv::hessian(function(p) neg_loglik(p, x, R, d$k), par)
    worst_H <- max(worst_H, max(abs(ana - num) / pmax(1, abs(num))))
  }
}
cat(sprintf("  largest relative discrepancy over %d samples: %.3e\n",
            length(designs) * 40, worst_H))

## The same comparison at the fitted values on the two real datasets, which is
## where the reported standard errors actually come from.
cat("\nAt the fitted values of the two applications:\n")
## The ball bearing sample is not written out here.  It is constructed by
## 10_realdata1_ballbearings.R and recorded in the macro file that script
## writes; a copy typed into this check would go stale the moment the
## construction changed, which is exactly what happened once already.
apps <- list(
  list(name = "ball bearings",
       x = as.numeric(strsplit(app_fit$bbSample, ",\\s*")[[1]]),
       R = c(as.numeric(app_fit$bbRfirst), rep(0, as.numeric(app_fit$bbM) - 1)),
       k = as.numeric(app_fit$bbK)),
  list(name = "guinea pigs",
       x = c(12, 15, 22, 24, 32, 34, 38, 38, 44, 53, 54, 54, 55, 56, 57, 58,
             65, 67, 70, 73, 81, 98, 109, 110, 131, 258),
       R = c(10, rep(0, 25)), k = 2))
for (a in apps) {
  f <- fit_mle(a$x, a$R, a$k, numeric_hessian = TRUE)
  if (is.null(f) || is.null(f$hessian_numeric)) {
    cat(sprintf("  %-14s no numerical Hessian available; skipped\n", a$name))
    next
  }
  rel <- max(abs(f$information - f$hessian_numeric) / pmax(1, abs(f$information)))
  cat(sprintf("  %-14s analytic vs optim numerical Hessian: %.2e\n", a$name, rel))
  worst_H <- max(worst_H, rel)
}

## ---------------------------------------------------------------------------
## 2.  The estimation procedure
## ---------------------------------------------------------------------------

cat("\nChecking the closed-form MLE of alpha given lambda ...\n")
worst_a <- 0
for (d in designs[1:3]) {
  R <- if (identical(d$sch, "app")) c(d$R1, rep(0, d$m - 1))
         else make_scheme(d$m, d$sch)
  for (i in 1:20) {
    x <- generate_pffc(d$m, d$k, R, d$alpha, d$lambda)
    for (lam in c(0.5, 1, 2, 10) * d$lambda) {
      a_closed <- alpha_hat_given_lambda(lam, x, R, d$k)
      ## maximise the likelihood over alpha alone, with lambda held fixed
      a_numeric <- optimise(function(a) neg_loglik(c(a, lam), x, R, d$k),
                            interval = c(1e-6, 50 * max(1, a_closed)),
                            tol = 1e-12)$minimum
      worst_a <- max(worst_a, abs(a_closed - a_numeric) / max(1, a_closed))
    }
  }
}
cat(sprintf("  largest relative discrepancy: %.3e\n", worst_a))

cat("\nChecking the fixed-point iteration of Eq. (15) ...\n")
n_tot <- n_fp <- 0L
worst_p <- 0
for (d in designs) {
  R <- if (identical(d$sch, "app")) c(d$R1, rep(0, d$m - 1))
         else make_scheme(d$m, d$sch)
  for (i in 1:60) {
    x  <- generate_pffc(d$m, d$k, R, d$alpha, d$lambda)
    n_tot <- n_tot + 1L
    ps <- profile_start(x, R, d$k)
    if (is.null(ps)) next
    fp <- lambda_fixed_point(x, R, d$k, ps$lambda)
    if (is.null(fp)) next
    n_fp <- n_fp + 1L
    op <- optim(fp$par, neg_loglik, x = x, R = R, k = d$k,
                method = "L-BFGS-B", lower = c(1e-8, 1e-8))
    ## where the iteration converges, it must agree with direct maximisation
    worst_p <- max(worst_p, max(abs(fp$par - op$par) / pmax(1, abs(op$par))))
  }
}
cat(sprintf("  converged on %d of %d samples (%.1f%%)\n",
            n_fp, n_tot, 100 * n_fp / n_tot))
cat(sprintf("  where it converged, largest relative difference from optim: %.3e\n",
            worst_p))
cat("  (fit_mle() accepts the fixed point only when the closed-form score\n",
    "   vanishes, the observed information is positive definite and the\n",
    "   likelihood has not decreased; otherwise it falls back to optim, so the\n",
    "   samples on which the iteration fails are still handled.)\n", sep = "")

if (is.finite(worst_H) && worst_H > 0 && worst_H < 1e-5 &&
    worst_a < 1e-6 && worst_p < 1e-5) {
  cat("\nPASS: the closed-form Hessian of Section 3, the explicit MLE of alpha,\n",
      "and the fixed-point iteration of Eq. (15) all agree with independent\n",
      "numerical computation.\n", sep = "")
} else {
  stop("FAIL: a closed form used on the main path disagrees with numerical computation.")
}
