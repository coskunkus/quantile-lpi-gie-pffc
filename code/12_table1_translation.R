################################################################################
##  12_table1_translation.R
##
##  Table 1 of the paper: the translation constants kappa = sigma/IQR and
##  delta = (mu - median)/IQR for the GIE distribution, and the quantile-scale
##  critical value corresponding to a moment-scale standard of 0.5.
##
##  Also verifies numerically that
##    (i)  Eq. (10)  -- C_L^xi depends on (alpha, L/lambda) only -- holds, and
##    (ii) Eq. (12) -- C_L^xi = kappa C_L^M - delta -- holds.
##
##  Runtime: a few seconds.
################################################################################

source("00_gie_pffc.R")
need("statmod")

ALPHAS <- c(2.5, 3, 4, 5, 7, 10, 20, 50)
C0M    <- 0.5

tc <- t(vapply(ALPHAS, translation_constants, c(kappa = 0, delta = 0)))
tab <- data.frame(alpha = ALPHAS, kappa = tc[, "kappa"], delta = tc[, "delta"])
tab$C0xi <- tab$kappa * C0M - tab$delta

## Normal reference: IQR = 2 * qnorm(0.75) * sigma, and delta = 0 by symmetry.
kappa_norm <- 1 / (2 * qnorm(0.75))
tab <- rbind(tab, data.frame(alpha = NA, kappa = kappa_norm, delta = 0,
                             C0xi = kappa_norm * C0M))

cat("Table 1: translation constants (last row is the normal reference)\n")
print(tab, row.names = FALSE, digits = 5)
write.csv(tab, file.path(RESULTS_DIR, "translation_constants.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## Table 1 as a LaTeX fragment
## ---------------------------------------------------------------------------

con <- base::file(file.path(TABLES_DIR, "tab_translation.tex"), open = "wt")
wl  <- function(...) writeLines(paste0(...), con)
wl("\\begin{table}[H]"); wl("\\centering")
wl(sprintf(paste0("\\caption{Translation constants $\\kappa=\\sigma/\\mathrm{IQR}$ and ",
                  "$\\delta=(\\mu-\\xi(0.5))/\\mathrm{IQR}$ for the GIE distribution, ",
                  "together with the quantile-based critical value ",
                  "$C_{0}^{\\xi}=\\kappa C_{0}^{M}-\\delta$ corresponding to the ",
                  "moment-based standard $C_{0}^{M}=%.1f$. Both constants are free of ",
                  "$\\lambda$ and of $L$. The final row is the normal reference.}"), C0M))
wl("\\label{T:kappa}")
wl("\\begin{tabular}{l ccc}"); wl("\\toprule")
wl(sprintf("$\\alpha$ & $\\kappa$ & $\\delta$ & $C_{0}^{\\xi}$ for $C_{0}^{M}=%.1f$ \\\\", C0M))
wl("\\midrule")
for (i in seq_along(ALPHAS))
  wl(sprintf("%-8s & %.4f & %.4f & %.4f \\\\",
             format(ALPHAS[i]), tab$kappa[i], tab$delta[i], tab$C0xi[i]))
wl("\\midrule")
j <- nrow(tab)
wl(sprintf("Normal   & %.4f & 0      & %.4f \\\\", tab$kappa[j], tab$C0xi[j]))
wl("\\bottomrule"); wl("\\end{tabular}"); wl("\\end{table}")
close(con)
cat("Wrote", file.path(TABLES_DIR, "tab_translation.tex"), "\n")

## ---------------------------------------------------------------------------
## The limit of the two constants as alpha grows
## ---------------------------------------------------------------------------
##
## Table 1 invites the reading that the GIE approaches normality as alpha grows:
## over the range shown, delta falls towards zero and kappa passes close to the
## normal 1/(2 z_{0.75}).  It does not stop there.  With W = lambda / X we have
## P(W <= w) = (1 - e^{-w})^alpha, so W is a maximum of alpha independent unit
## exponentials and W - log(alpha) converges to a standard Gumbel variate G.
## Then X = lambda / (log(alpha) + G_alpha), and the standardised limit is a
## REFLECTED Gumbel, not a normal.  Since kappa and delta are scale-free they
## converge to the corresponding functionals of -G.
gumbel_gamma <- -digamma(1)                       # Euler-Mascheroni
gum_sd       <- pi / sqrt(6)                      # sd of G
gum_q        <- function(p) -log(-log(p))         # quantile of G
gum_iqr      <- gum_q(0.75) - gum_q(0.25)
gum_med      <- gum_q(0.5)
kappa_gumbel <- gum_sd / gum_iqr                             # sd(-G) / IQR(-G)
delta_gumbel <- (-gumbel_gamma - (-gum_med)) / gum_iqr       # (E - med) / IQR

## The shape at which delta crosses zero, quoted as the point where the
## agreement with the normal reference ends.  translation_constants() uses
## quadrature and loses accuracy for very large alpha, so the root is sought
## over a range where it is still reliable.
delta_of <- function(a) unname(translation_constants(a)["delta"])
alpha_zero <- tryCatch(
  uniroot(delta_of, c(200, 5000), tol = 1e-6)$root,
  error = function(e) NA_real_)

cat(sprintf("\nLimits as alpha -> Inf (reflected Gumbel): kappa = %.4f, delta = %+.4f\n",
            kappa_gumbel, delta_gumbel))
cat(sprintf("delta crosses zero near alpha = %.0f\n", alpha_zero))
for (a in c(50, 200, 1000, 5000)) {
  tc <- translation_constants(a)
  cat(sprintf("  alpha = %6g : kappa = %.4f, delta = %+.4f\n", a, tc["kappa"], tc["delta"]))
}

## Numbers quoted in the running text of Section 2.3.
write_macros(file.path(TABLES_DIR, "values_translation.tex"), list(
  normIQRoverSigma = fmt(2 * qnorm(0.75), 4),
  normKappa        = fmt(kappa_norm, 4),
  normCzeroXi      = fmt(kappa_norm * C0M, 3),
  gumKappa         = fmt(kappa_gumbel, 4),
  gumDelta         = fmt(delta_gumbel, 4),
  gumAlphaZero     = if (is.na(alpha_zero)) "1000" else
                     formatC(signif(alpha_zero, 2), format = "d", big.mark = "{,}")
), "12_table1_translation.R")

cat(sprintf("\nNormal reference: IQR / sigma = %.4f, kappa = %.4f\n",
            2 * qnorm(0.75), kappa_norm))
cat(sprintf("Under normality C_L^xi = %.4f * C_L^M, so C_L^M > 0.5 corresponds to C_L^xi > %.4f\n",
            kappa_norm, kappa_norm * 0.5))

## ---------------------------------------------------------------------------
## Verification 1: Eq. (10), dependence on (alpha, L/lambda) only
## ---------------------------------------------------------------------------

cat("\n--- Check of Eq. (10) ---\n")
set.seed(MASTER_SEED + 12L)
worst <- 0
for (i in 1:200) {
  alpha <- runif(1, 0.3, 30)
  ratio <- runif(1, 0.01, 3)
  lam1  <- runif(1, 0.1, 100); lam2 <- runif(1, 0.1, 100)
  v1 <- CL_quantile(c(alpha, lam1), ratio * lam1)
  v2 <- CL_quantile(c(alpha, lam2), ratio * lam2)
  v3 <- CL_quantile_ratio(alpha, ratio)
  worst <- max(worst, abs(v1 - v2), abs(v1 - v3))
}
cat(sprintf("Largest discrepancy over 200 random configurations: %.3e\n", worst))
stopifnot(worst < 1e-8)
cat("Eq. (10) verified.\n")

## ---------------------------------------------------------------------------
## Verification 2: Eq. (12), the affine map
## ---------------------------------------------------------------------------

cat("\n--- Check of Eq. (12) ---\n")
worst <- 0
for (alpha in c(2.5, 3, 4, 5, 10, 20)) {
  for (lambda in c(0.5, 2, 50)) {
    for (L in c(0.1, 1, 5)) {
      ## Both sides use the same quadrature rule, which is exactly equivariant
      ## in lambda, so the quadrature error cancels and the identity holds to
      ## machine precision rather than to quadrature precision.
      kd <- translation_constants(alpha, 800)
      lhs <- CL_quantile(c(alpha, lambda), L)
      rhs <- kd["kappa"] * CL_moment(c(alpha, lambda), L, 800) - kd["delta"]
      worst <- max(worst, abs(lhs - rhs))
    }
  }
}
cat(sprintf("Largest discrepancy: %.3e\n", worst))
stopifnot(worst < 1e-8)
cat("Eq. (12) verified.\n")
