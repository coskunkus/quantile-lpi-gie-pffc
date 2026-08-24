################################################################################
##  validation/v3_check_moment_formula.R
##
##  Checks the closed-form moment expression of Eq. (3), which is proved in
##  Appendix A, against two independent numerical calculations: Gauss-Legendre
##  quadrature (Eq. 4) and R's adaptive integrate().
##
##  The formula is valid only for integer alpha with r < alpha, which is why the
##  paper also provides the quadrature route.
################################################################################

source("../00_gie_pffc.R")
need("statmod")

integrate_moment <- function(r, alpha, lambda)
  integrate(function(x) x^r * gie_pdf(x, alpha, lambda),
            lower = 0, upper = Inf, rel.tol = 1e-10)$value

## Two comparisons are made, and they test different things.
##   worst_ad : closed form vs adaptive integration.  This is the actual test of
##              Eq. (3), and both sides are accurate to near machine precision,
##              so the tolerance can be tight.
##   worst_gl : closed form vs the Gauss-Legendre rule of Eq. (4) at the number
##              of nodes used by CL_moment().  This measures the accuracy of the
##              quadrature, not the correctness of Eq. (3), and its tolerance is
##              set from the accuracy claimed in Section 2.1.
rows <- list()
worst_ad <- worst_gl <- 0
for (alpha in 2:8) {
  for (r in 1:(alpha - 1)) {
    for (lambda in c(0.5, 2, 10)) {
      ex <- gie_moment_exact(r, alpha, lambda)
      gl <- gie_moment_gl(r, alpha, lambda, 400)
      ad <- integrate_moment(r, alpha, lambda)
      worst_ad <- max(worst_ad, abs(ex - ad) / max(1, abs(ex)))
      worst_gl <- max(worst_gl, abs(ex - gl) / max(1, abs(ex)))
      if (lambda == 2 && r <= 2)
        rows[[length(rows) + 1]] <-
          data.frame(alpha, r, lambda, exact = ex, quadrature = gl, adaptive = ad)
    }
  }
}

cat("Comparison for lambda = 2 (Eq. 3 vs Eq. 4 vs integrate()):\n")
print(do.call(rbind, rows), row.names = FALSE, digits = 10)
cat("\nLargest relative discrepancy over all integer alpha in 2..8, r < alpha,\n")
cat("and lambda in {0.5, 2, 10}:\n")
cat(sprintf("  Eq. (3) vs adaptive integration : %.3e   (tests Eq. 3 itself)\n", worst_ad))
cat(sprintf("  Eq. (3) vs quadrature, N = 400  : %.3e   (tests the accuracy of Eq. 4)\n", worst_gl))

if (worst_ad < 1e-9 && worst_gl < 1e-5) {
  cat("\nPASS: the closed form of Eq. (3) reproduces the adaptive integral to near\n",
      "machine precision, and the Gauss-Legendre rule of Eq. (4) attains the\n",
      "accuracy claimed in Section 2.1.\n", sep = "")
} else {
  stop("FAIL: Eq. (3) does not match the numerical integrals.")
}

## The identity on which Appendix A turns, checked directly.  An (alpha-1)-st
## forward difference annihilates every polynomial of degree below alpha-1, so
## the sums below vanish for every r < alpha.  In the proof this is what makes
## the polynomial part of the partial-fraction expansion disappear and what
## leaves the residues summing to zero, which in turn is why the integral in
## Eq. (25) of the appendix converges.
cat("\nForward-difference cancellation the proof relies on (should be ~0):\n")
for (alpha in 3:8) for (r in 1:(alpha - 1)) {
  j <- 0:(alpha - 1)
  s <- sum((-1)^j * choose(alpha - 1, j) * (j + 1)^(r - 1))
  cat(sprintf("  alpha=%d r=%d : %+.3e\n", alpha, r, s))
}
