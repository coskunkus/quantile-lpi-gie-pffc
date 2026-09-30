################################################################################
##  validation/v2_check_invariance.R
##
##  Checks the location-scale invariance claimed at the end of Section 2.2.
##  For Y = aX + b with a > 0 and L' = aL + b, both indexes must be unchanged.
##
##  The check has to be done carefully to mean anything.  Writing
##      E[Y] = a E[X] + b,   sd(Y) = a sd(X),   xi_p(Y) = a xi_p(X) + b
##  and substituting reduces the transformed index to the original one as an
##  algebraic identity, so a "check" built that way cannot fail and measures
##  only rounding error.  Instead we compute the quantities for Y from the
##  distribution of Y itself, without using any of those three relations:
##    - the quantiles of Y by solving F_Y(y) = p, with F_Y(y) = F_X((y-b)/a),
##      using uniroot on the GIE distribution function;
##    - the mean and second moment of Y by numerically integrating against the
##      density f_Y(y) = f_X((y-b)/a)/a over the support (b, Inf).
##  Only then are the two index values compared.
################################################################################

source("../00_gie_pffc.R")
need("statmod")

set.seed(MASTER_SEED + 102L)

## Quantile of Y obtained by inverting F_Y numerically.
quantile_Y <- function(p, alpha, lambda, a, b) {
  f <- function(y) gie_cdf((y - b) / a, alpha, lambda) - p
  ## bracket: start from the transformed GIE quantile only as a search seed,
  ## then widen until the sign change is captured.
  seed <- a * gie_quantile(p, alpha, lambda) + b
  lo <- b + (seed - b) / 100
  hi <- b + (seed - b) * 100
  while (f(lo) > 0) lo <- b + (lo - b) / 10
  while (f(hi) < 0) hi <- b + (hi - b) * 10
  uniroot(f, c(lo, hi), tol = .Machine$double.eps^0.75)$root
}

## Raw moments of Y by direct numerical integration of y^r f_Y(y).
moment_Y <- function(r, alpha, lambda, a, b) {
  dens <- function(y) gie_pdf((y - b) / a, alpha, lambda) / a
  integrate(function(y) y^r * dens(y), lower = b, upper = Inf,
            rel.tol = 1e-10, subdivisions = 2000L)$value
}

worst_m <- worst_q <- 0
rows <- list()
for (i in 1:25) {
  alpha  <- runif(1, 3.0, 12)          # alpha > 2 so that C_L^M exists
  lambda <- runif(1, 0.5, 20)
  L      <- runif(1, 0, 2 * lambda)
  a      <- exp(runif(1, -1.5, 1.5))
  b      <- runif(1, -2, 20)

  ## --- original scale, computed from the GIE distribution ---
  mu  <- gie_moment_gl(1, alpha, lambda, 800)
  ex2 <- gie_moment_gl(2, alpha, lambda, 800)
  sdv <- sqrt(ex2 - mu^2)
  med <- gie_quantile(0.50, alpha, lambda)
  q1  <- gie_quantile(0.25, alpha, lambda)
  q3  <- gie_quantile(0.75, alpha, lambda)
  CM  <- (mu - L) / sdv
  CX  <- (med - L) / (q3 - q1)

  ## --- transformed scale, computed from the distribution of Y ---
  Lp   <- a * L + b
  muY  <- moment_Y(1, alpha, lambda, a, b)
  ex2Y <- moment_Y(2, alpha, lambda, a, b)
  sdY  <- sqrt(ex2Y - muY^2)
  medY <- quantile_Y(0.50, alpha, lambda, a, b)
  q1Y  <- quantile_Y(0.25, alpha, lambda, a, b)
  q3Y  <- quantile_Y(0.75, alpha, lambda, a, b)
  CMp  <- (muY - Lp) / sdY
  CXp  <- (medY - Lp) / (q3Y - q1Y)

  worst_m <- max(worst_m, abs(CM - CMp) / max(1, abs(CM)))
  worst_q <- max(worst_q, abs(CX - CXp) / max(1, abs(CX)))
  if (i <= 5) rows[[length(rows) + 1]] <-
    data.frame(alpha, lambda, L, a, b, C_M = CM, C_M_trans = CMp,
               C_xi = CX, C_xi_trans = CXp)
}

cat("First five of 25 random transformations.\n")
cat("The '_trans' columns are computed from the distribution of Y = aX + b,\n")
cat("not from the invariance relations being tested.\n\n")
print(do.call(rbind, rows), row.names = FALSE, digits = 7)

cat(sprintf("\nLargest relative discrepancy over 25 transformations:\n"))
cat(sprintf("  moment-based index   : %.3e\n", worst_m))
cat(sprintf("  quantile-based index : %.3e\n", worst_q))

if (worst_m < 1e-6 && worst_q < 1e-8) {
  cat("\nPASS: both indexes are location-scale invariant.\n")
} else {
  stop("FAIL: invariance violated.")
}
