################################################################################
##  validation/v4_check_information.R
##
##  Checks the two elementary inequalities on which Appendix B rests, and then
##  confirms numerically that the resulting bounds on the observed information
##  hold for simulated PFFC samples -- including heavy-tailed cases with
##  alpha <= 2, where Var(X) is infinite.
##
##  The point of the exercise is the one made in Appendix B: the bounds involve
##  only alpha, lambda, k and R, and not the observed lifetimes, so they cannot
##  be affected by the weight of the tail.
################################################################################

source("../00_gie_pffc.R")

## ---------------------------------------------------------------------------
## 1.  The two scalar inequalities
## ---------------------------------------------------------------------------

z <- exp(seq(log(1e-8), log(50), length.out = 200000))

g1 <- z / expm1(z)                         # must lie in (0, 1]
g2 <- z^2 * exp(z) / expm1(z)^2            # must lie in (0, 1]

cat(sprintf("max over z of z/(e^z - 1)          = %.12f   (bound 1)\n", max(g1)))
cat(sprintf("max over z of z^2 e^z/(e^z - 1)^2  = %.12f   (bound 1)\n", max(g2)))
stopifnot(max(g1) <= 1 + 1e-9, max(g2) <= 1 + 1e-9, min(g1) > 0, min(g2) > 0)
cat("Both inequalities of Appendix B hold.\n\n")

## ---------------------------------------------------------------------------
## 2.  The bounds on the observed information, on simulated samples
## ---------------------------------------------------------------------------

set.seed(MASTER_SEED + 104L)

check_one <- function(m, k, R, alpha, lambda) {
  x  <- generate_pffc(m, k, R, alpha, lambda)
  zi <- lambda / x
  ## -expm1(-zi) rather than 1 - exp(-zi): the cancellation as zi -> 0 is
  ## severe enough that the naive form breaches the bound by rounding error
  ## alone, which would make this check fail at random.
  t1 <- exp(-zi) / (x * -expm1(-zi))             # bounded by 1/lambda
  t2 <- exp(-zi) / (x^2 * expm1(-zi)^2)          # bounded by 1/lambda^2

  H_aa <- -m / alpha^2
  H_al <- sum(k * (R + 1) * t1)
  H_ll <- -m / lambda^2 - sum((alpha * k * (R + 1) - 1) * t2)

  b_aa <- m / alpha^2
  b_al <- (k / lambda) * sum(R + 1)
  b_ll <- m / lambda^2 + (1 / lambda^2) * sum(abs(alpha * k * (R + 1) - 1))

  tol <- 1e-8                                    # relative, not absolute
  c(t1_ok = all(t1 > 0 & t1 <= (1 / lambda) * (1 + tol)),
    t2_ok = all(t2 > 0 & t2 <= (1 / lambda^2) * (1 + tol)),
    aa_ok = abs(H_aa) <= b_aa * (1 + tol),
    al_ok = abs(H_al) <= b_al * (1 + tol),
    ll_ok = abs(H_ll) <= b_ll * (1 + tol),
    finite = all(is.finite(c(H_aa, H_al, H_ll))))
}

grid <- expand.grid(alpha = c(0.5, 1, 1.5, 2, 5, 10),
                    lambda = c(0.5, 2, 100),
                    m = c(9, 25, 100), k = c(1, 2, 5),
                    scheme = c("Early", "Middle", "Late", "Equal"),
                    stringsAsFactors = FALSE)

res <- t(vapply(seq_len(nrow(grid)), function(i) {
  R <- make_scheme(grid$m[i], grid$scheme[i])
  check_one(grid$m[i], grid$k[i], R, grid$alpha[i], grid$lambda[i])
}, numeric(6)))

cat("Configurations tested:", nrow(grid),
    "(including alpha <= 2, where Var(X) is infinite)\n")
cat("Failures per check:\n")
print(colSums(res == 0))

if (all(res == 1)) {
  cat("\nPASS: every element of the observed information is finite and obeys the\n",
      "bounds of Appendix B, for all shape values tested including alpha <= 2.\n", sep = "")
} else {
  stop("FAIL: a bound of Appendix B was violated.")
}

## ---------------------------------------------------------------------------
## 3.  Nonsingularity in the configurations actually used in the paper
## ---------------------------------------------------------------------------

cat("\nSmallest eigenvalue of the observed information over 200 simulated\n",
    "samples in each of the reported designs (should be strictly positive):\n", sep = "")
designs <- list(
  list(m = 25,  k = 2, alpha = 5,   lambda = 2,      sch = "Early"),
  list(m = 25,  k = 5, alpha = 1,   lambda = 2,      sch = "Late"),
  list(m = 100, k = 2, alpha = 2,   lambda = 2,      sch = "Equal"),
  list(m = 9,   k = 2, alpha = 3.42, lambda = 155.3, sch = "Early"))
for (d in designs) {
  R <- if (d$m == 9) c(3, rep(0, d$m - 1)) else make_scheme(d$m, d$sch)
  mn <- Inf
  for (i in 1:200) {
    x  <- generate_pffc(d$m, d$k, R, d$alpha, d$lambda)
    fx <- fit_mle(x, R, d$k, start = c(d$alpha, d$lambda))
    if (is.null(fx)) next
    mn <- min(mn, min(eigen(fx$information, symmetric = TRUE, only.values = TRUE)$values))
  }
  cat(sprintf("  m=%3d k=%d alpha=%.2f lambda=%.1f %-6s : %.4g\n",
              d$m, d$k, d$alpha, d$lambda, d$sch, mn))
}
