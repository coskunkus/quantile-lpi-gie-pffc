################################################################################
##  validation/v1_check_derivatives.R
##
##  Checks the analytical partial derivatives given in the Remark of Section 3
##  against numerical differentiation.  Both the component derivatives of the
##  quantile function and the assembled gradient of C_L^xi are tested, over a
##  grid of parameter values rather than at a single point.
################################################################################

source("../00_gie_pffc.R")
need("numDeriv")

cat("Checking d xi(p)/d alpha and d xi(p)/d lambda ...\n")
worst_xi <- 0
for (p in c(0.25, 0.5, 0.75, 0.05, 0.95)) {
  for (alpha in c(0.5, 1, 2, 3.5, 10, 50)) {
    for (lambda in c(0.1, 1, 100)) {
      num <- numDeriv::grad(function(th) gie_quantile(p, th[1], th[2]),
                            c(alpha, lambda))
      ana <- c(dxi_dalpha(p, alpha, lambda), dxi_dlambda(p, alpha))
      worst_xi <- max(worst_xi, max(abs(num - ana) / pmax(1, abs(num))))
    }
  }
}
cat(sprintf("  largest relative discrepancy: %.3e\n", worst_xi))

cat("Checking the gradient of C_L^xi ...\n")
worst_g <- 0
detail  <- NULL
for (alpha in c(0.5, 1, 2, 3.5, 10, 50)) {
  for (lambda in c(0.1, 1, 100)) {
    for (L in c(0, 0.5 * lambda, 2 * lambda)) {
      num <- numDeriv::grad(CL_quantile, c(alpha, lambda), L = L)
      ana <- grad_CL_quantile(c(alpha, lambda), L)
      rel <- max(abs(num - ana) / pmax(1, abs(num)))
      if (rel > worst_g) { worst_g <- rel; detail <- c(alpha, lambda, L) }
    }
  }
}
cat(sprintf("  largest relative discrepancy: %.3e  at alpha=%.2f lambda=%.2f L=%.2f\n",
            worst_g, detail[1], detail[2], detail[3]))

if (worst_xi < 1e-6 && worst_g < 1e-6) {
  cat("\nPASS: the analytical derivatives in the Remark of Section 3 agree with\n",
      "numerical differentiation throughout the grid tested.\n", sep = "")
} else {
  stop("FAIL: discrepancy between analytical and numerical derivatives.")
}
