################################################################################
##  validation/v8_check_score_moments.R
##
##  Are the moments of the score finite when the lifetime has no moments?
##
##  This is the question the referee asked under M5, and Section 3 now answers
##  it by an identity rather than by an inequality.  The observations are not
##  independent draws from F: they are progressively first-failure-censored
##  order statistics.  Applying the probability integral transform to F, as
##  though they were independent, gives the wrong answer, and an earlier draft
##  of this paper did exactly that.  The correct route goes through the
##  representation the bootstrap of Section 4 already uses.  Writing
##
##      E_i = -log{1 - F_k(X_i)},      F_k(x) = 1 - (1 - e^{-lambda/x})^{alpha k},
##
##  the algorithm of Balakrishnan and Sandhu (1995) gives
##
##      E_i = sum_{j<=i} W_j / (n - sum_{l<j} R_l - j + 1),   W_j iid Exp(1),
##
##  and therefore the identity on which the whole argument rests,
##
##      log(1 - e^{-lambda/X_i}) = -E_i / (alpha k).                       (*)
##
##  Since each E_i is a fixed positive linear combination of independent
##  standard exponential variables, it has moments of every order.  The two
##  score components then read
##
##      dl/dalpha  = m/alpha - alpha^{-1} sum_i (R_i + 1) E_i,
##      dl/dlambda = m/lambda - lambda^{-1} sum_i Z_i
##                   + lambda^{-1} sum_i {alpha k (R_i+1) - 1} g(E_i),
##
##  with Z_i = lambda / X_i = -log{1 - e^{-E_i/(alpha k)}} and
##  g(E) = Z(E) {e^{E/(alpha k)} - 1}.  The first is linear in the E_i.  In the
##  second, g is bounded by 1 and Z has only a logarithmic singularity at
##  E = 0, so both have moments of every order.
##
##  Nothing in that chain refers to E(X) or E(X^2), which is the point: the
##  score can have finite moments of every order while the lifetime itself has
##  none.  This script checks each link, and then exhibits the contrast
##  directly by letting the sample size grow: the running average of X^2
##  wanders off when alpha <= 2, while the running average of the squared score
##  settles down.
##
##  Runtime: two to three minutes.
################################################################################

source("../00_gie_pffc.R")

set.seed(MASTER_SEED + 108L)

fails <- character(0)
chk <- function(cond, msg) if (!isTRUE(cond)) fails <<- c(fails, msg)

## ---------------------------------------------------------------------------
## 1.  g(E) is bounded by 1, with the two limits the paper claims
## ---------------------------------------------------------------------------
## Written in u = E/(alpha k), so the check is free of alpha, k and lambda.
## This is bound 1 of Appendix B in another guise: with t = e^{-u},
## g = -log(1-t)(1-t)/t, which increases to 1 and never reaches it.

u <- exp(seq(log(1e-12), log(60), length.out = 400000))
g <- -log1p(-exp(-u)) * expm1(u)

cat(sprintf("g(u) = -log(1-e^-u)(e^u - 1):  max = %.12f   (must not exceed 1)\n",
            max(g)))
cat(sprintf("  limit as u -> 0   : %.3e   (claimed 0)\n", g[1]))
cat(sprintf("  limit as u -> inf : %.12f   (claimed 1)\n", g[length(g)]))
chk(max(g) <= 1 + 1e-9, "g(E) exceeds its claimed bound of 1")
chk(g[1] < 1e-8, "g(E) does not vanish as E -> 0")
chk(abs(g[length(g)] - 1) < 1e-9, "g(E) does not tend to 1 as E -> infinity")

## ---------------------------------------------------------------------------
## 2.  Identity (*), and the exponential spacings behind it
## ---------------------------------------------------------------------------
## E_i is recovered from the sample by (*) and then differenced.  If the
## representation is right, the normalised spacings must be independent
## standard exponential variables.  The samples come from generate_pffc(),
## the function the simulations actually use, so this checks the production
## code and not a private copy of the algorithm.

spacing_denominators <- function(m, R) {
  n <- sum(R) + m
  n - c(0, cumsum(R)[-m]) - seq_len(m) + 1
}

check_identity <- function(m, k, R, alpha, lambda, N = 10000, label = "") {
  d <- spacing_denominators(m, R)
  D <- matrix(NA_real_, N, m)
  for (i in seq_len(N)) {
    x  <- generate_pffc(m, k, R, alpha, lambda)
    Ei <- -alpha * k * log1p(-exp(-lambda / x))     # identity (*), solved for E
    chk(all(diff(Ei) > 0), sprintf("%s: recovered E is not increasing", label))
    D[i, ] <- d * diff(c(0, Ei))                    # normalised spacings
  }
  p  <- vapply(seq_len(m), function(j) suppressWarnings(
                 ks.test(D[, j], "pexp", 1)$p.value), 0)
  cr <- max(abs(cor(D)[upper.tri(diag(m))]))
  cat(sprintf("%-22s m=%3d k=%d alpha=%-5g : min KS p = %.4f, max |cor| = %.4f, mean = %.4f\n",
              label, m, k, alpha, min(p), cr, mean(D)))
  chk(min(p) > 1e-4,
      sprintf("%s: normalised spacings are not standard exponential (KS p = %.5f)",
              label, min(p)))
  chk(cr < 4 / sqrt(N),
      sprintf("%s: normalised spacings are correlated (max |cor| = %.4f)", label, cr))
  invisible(NULL)
}

cat("\nIdentity log(1 - e^{-lambda/x}) = -E/(alpha k), and the spacings behind it:\n")
app <- read_macros(file.path(TABLES_DIR, "values_realdata1.tex"))
designs <- list(
  list(m = 25, k = 2, sch = "Early", alpha = 5,   lambda = 2,   lab = "Early, light tail"),
  list(m = 12, k = 3, sch = "Late",  alpha = 1,   lambda = 2,   lab = "Late, alpha = 1"),
  list(m = 12, k = 1, sch = "Equal", alpha = 0.5, lambda = 100, lab = "Equal, alpha = 0.5"),
  list(m = as.numeric(app$bbM), k = as.numeric(app$bbK), sch = "app",
       alpha = as.numeric(app$bbAlphaHat), lambda = as.numeric(app$bbLambdaHat),
       R = c(as.numeric(app$bbRfirst), rep(0, as.numeric(app$bbM) - 1)),
       lab = "first application"))
for (dd in designs) {
  R <- if (is.null(dd$R)) make_scheme(dd$m, dd$sch) else dd$R
  check_identity(dd$m, dd$k, R, dd$alpha, dd$lambda, label = dd$lab)
}

## ---------------------------------------------------------------------------
## 3.  The two closed forms for the score
## ---------------------------------------------------------------------------
## The expressions Section 3 differentiates are compared with loglik_score(),
## the function the estimation code uses.  They are algebraically the same
## thing; if they disagree, one of the two is wrong.

cat("\nClosed forms of the score against loglik_score():\n")
worst <- 0
for (dd in designs) {
  R <- if (is.null(dd$R)) make_scheme(dd$m, dd$sch) else dd$R
  a <- dd$alpha; lam <- dd$lambda; k <- dd$k; m <- dd$m
  for (i in 1:200) {
    x  <- generate_pffc(m, k, R, a, lam)
    Ei <- -a * k * log1p(-exp(-lam / x))
    Zi <- lam / x
    gi <- Zi * expm1(Ei / (a * k))
    mine <- c(m / a - sum((R + 1) * Ei) / a,
              m / lam - sum(Zi) / lam + sum((a * k * (R + 1) - 1) * gi) / lam)
    ref  <- loglik_score(c(a, lam), x, R, k)
    worst <- max(worst, max(abs(mine - ref) / pmax(1, abs(ref))))
  }
  cat(sprintf("  %-22s largest relative difference so far: %.3e\n", dd$lab, worst))
}
chk(worst < 1e-8,
    sprintf("the closed forms of Section 3 disagree with loglik_score() (%.3e)", worst))

## ---------------------------------------------------------------------------
## 4.  Score moments where the lifetime has none
## ---------------------------------------------------------------------------
## Var(X) is infinite for alpha <= 2 and E(X) as well for alpha <= 1.  The
## claim under test is that the score is untouched by this.  Both are estimated
## from the same simulated samples, at two sample sizes, so that a quantity
## with no population moment is seen to move while one with finite moments is
## seen to stay put.  A sample average of X^2 drawn from a distribution with
## infinite second moment is dominated by its largest draw, so it grows with N
## instead of converging.

score_moments <- function(m, k, R, alpha, lambda, N) {
  S <- matrix(NA_real_, N, 2); x2 <- numeric(N)
  for (i in seq_len(N)) {
    x <- generate_pffc(m, k, R, alpha, lambda)
    S[i, ] <- loglik_score(c(alpha, lambda), x, R, k)
    x2[i]  <- mean(x^2)
  }
  c(mean_a = mean(S[, 1]), mean_l = mean(S[, 2]),
    m2_a = mean(S[, 1]^2), m2_l = mean(S[, 2]^2),
    m4_a = mean(S[, 1]^4), m4_l = mean(S[, 2]^4),
    se_a = sd(S[, 1]) / sqrt(N), se_l = sd(S[, 2]) / sqrt(N),
    x2 = mean(x2))
}

cat("\nFourth moments of the two score components, at two Monte Carlo sizes.\n")
cat("A finite moment barely moves when N is quadrupled; a nonexistent one does\n")
cat("not settle.  The last column is the sample average of X^2, shown for\n")
cat("contrast: it has no population counterpart once alpha <= 2.\n\n")
cat(sprintf("%-7s %-6s %12s %12s %12s %12s\n",
            "alpha", "N", "E[S_a^4]", "E[S_l^4]", "mean score", "avg X^2"))

R25 <- make_scheme(25, "Equal")
for (alpha in c(5, 2, 1, 0.5)) {
  prev <- NULL
  for (N in c(5000, 20000)) {
    r <- score_moments(25, 2, R25, alpha, 2, N)
    cat(sprintf("%-7g %-6d %12.4g %12.4g %12.4f %12.4g\n",
                alpha, N, r["m4_a"], r["m4_l"],
                max(abs(r[c("mean_a", "mean_l")])), r["x2"]))
    ## The score has mean zero at the true parameter.  This is the strongest
    ## single check here: it holds only if the likelihood, its derivative and
    ## the generator all describe the same model.
    z <- abs(r[c("mean_a", "mean_l")]) / r[c("se_a", "se_l")]
    chk(all(z < 4), sprintf("alpha = %g, N = %d: the score does not have mean zero (|z| = %.2f)",
                            alpha, N, max(z)))
    chk(all(is.finite(r[c("m4_a", "m4_l")])),
        sprintf("alpha = %g: a fourth moment of the score is not finite", alpha))
    if (!is.null(prev)) {
      ratio <- max(r[c("m4_a", "m4_l")] / prev[c("m4_a", "m4_l")])
      chk(ratio < 3 && ratio > 1 / 3,
          sprintf("alpha = %g: the fourth moment of the score changed by a factor of %.2f when N was quadrupled, so it does not look finite",
                  alpha, ratio))
    }
    prev <- r
  }
}

## ---------------------------------------------------------------------------
cat("\n")
if (length(fails) == 0) {
  cat("PASS: the score-moment argument of Section 3 holds as written.  The\n",
      "identity log(1 - e^{-lambda/x}) = -E/(alpha k) is satisfied by the\n",
      "samples the simulations use, the normalised spacings behind it are\n",
      "standard exponential and independent, the closed forms of the two score\n",
      "components agree with the estimation code, and the fourth moments of\n",
      "both components are finite and stable at alpha = 0.5, 1 and 2, where the\n",
      "lifetime has no variance and, at alpha <= 1, no mean either.\n", sep = "")
} else {
  for (f in fails) cat("  - ", f, "\n", sep = "")
  stop("FAIL: the score-moment argument of Section 3 does not hold as written.")
}
