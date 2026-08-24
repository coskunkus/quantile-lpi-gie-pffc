################################################################################
##  10_realdata1_ballbearings.R
##
##  First application, Sections 6.1 and 6.3: deep groove ball bearing endurance
##  data (Lieblein & Zelen 1956; Lawless 2011).
##
##  Reproduces Table 9 and every number quoted in Sections 6.1 and 6.3,
##  including the one-sided decision rule of Eq. (22).
##
##  The bootstrap seed is fixed.  The earlier version of this script called
##  set.seed(NULL), which made the bootstrap intervals irreproducible.
##
##  Runtime: under a minute.
################################################################################

source("00_gie_pffc.R")
need("numDeriv", "statmod")

set.seed(MASTER_SEED + 10L)

B_BOOT <- 1000
LEVEL  <- 0.95

## ---------------------------------------------------------------------------
## Data
## ---------------------------------------------------------------------------

## Complete data: 23 endurance times, in millions of revolutions.
full <- c(17.88, 28.92, 33.00, 41.52, 42.12, 45.60, 48.40, 51.84, 51.96,
          54.12, 55.56, 67.80, 68.64, 68.64, 68.88, 84.12, 93.12, 98.64,
          105.12, 105.84, 127.92, 128.04, 173.40)

## PFFC sample generated from the complete data with m = 9, k = 2,
## R = (3, 0^8), as described in Section 6.1.
x <- c(41.52, 42.12, 48.40, 51.84, 51.96, 55.56, 68.64, 127.92, 173.40)
m <- length(x)
k <- 2
R <- c(3, rep(0, m - 1))
L <- 20

cat("Ball bearing data:  m =", m, " k =", k, " n =", m + sum(R), " L =", L, "\n\n")

## ---------------------------------------------------------------------------
## Point estimates and asymptotic standard errors
## ---------------------------------------------------------------------------

fit <- fit_mle(x, R, k, numeric_hessian = TRUE)
stopifnot(!is.null(fit))
par_hat <- fit$par
I_obs   <- fit$information          # analytic, Eq. (16)

cat(sprintf("Estimation route: %s", fit$route))
if (identical(fit$route, "fixedpoint"))
  cat(sprintf(" (Eq. (15) converged in %d iterations)", fit$iterations))
cat(sprintf("\nMaximised log-likelihood: %.8f\n", fit$loglik))

ev <- eigen(I_obs, symmetric = TRUE, only.values = TRUE)$values
cat("Observed information eigenvalues:", sprintf("%.4g", ev),
    if (min(ev) > 0) " (positive definite)\n" else " (NOT positive definite)\n")
if (is.null(fit$hessian_numeric)) {
  cat("optim did not return a numerical Hessian; the analytic information of\n",
      "  Eq. (17) is used regardless, and validation/v5_ checks it.\n", sep = "")
} else {
  cat(sprintf("Largest relative difference between the analytic information of Eq. (17)\n  and the numerical Hessian from optim: %.2e\n\n",
              max(abs(I_obs - fit$hessian_numeric) / pmax(1, abs(I_obs)))))
}

V  <- solve(I_obs)
se_alpha  <- sqrt(V[1, 1])
se_lambda <- sqrt(V[2, 2])

Cq <- CL_quantile(par_hat, L)
Cm <- CL_moment(par_hat, L)
se_q <- delta_se(fit, L, "quantile")
se_m <- delta_se(fit, L, "moment")

cat("Point estimates and asymptotic standard errors\n")
print(data.frame(
  Parameter = c("alpha", "lambda", "C_L^xi", "C_L^M"),
  Estimate  = c(par_hat[1], par_hat[2], Cq, Cm),
  SE        = c(se_alpha, se_lambda, se_q, se_m)),
  row.names = FALSE, digits = 5)

## Approximate 95% interval for alpha, quoted in Section 6.1 in support of the
## statement that the finiteness of the second moment is not established.
ci_alpha <- par_hat[1] + c(-1, 1) * qnorm(0.975) * se_alpha
cat(sprintf("\n95%% CI for alpha: (%.3f, %.3f)  -- contains values below 2: %s\n",
            ci_alpha[1], ci_alpha[2], ci_alpha[1] < 2))

## Translation constants at the fitted shape, quoted in Section 6.1.
## Same number of nodes on both sides, so that the identity of Eq. (12) is
## displayed as the exact relation it is rather than up to quadrature error.
kd <- translation_constants(par_hat[1], n_nodes = 400)
cat(sprintf("Translation constants at alpha-hat: kappa = %.4f, delta = %.4f\n",
            kd["kappa"], kd["delta"]))
cat(sprintf("Check of Eq. (12):  kappa*C_M - delta = %.4f   vs   C_xi = %.4f\n",
            kd["kappa"] * Cm - kd["delta"], Cq))
cat(sprintf("L / lambda-hat = %.4f\n", L / par_hat[2]))

## ---------------------------------------------------------------------------
## Bootstrap
## ---------------------------------------------------------------------------

cat("\nParametric bootstrap, B =", B_BOOT, "...\n")
bq <- bm <- rep(NA_real_, B_BOOT)
for (b in seq_len(B_BOOT)) {
  y  <- generate_pffc(m, k, R, par_hat[1], par_hat[2])
  fb <- fit_mle(y, R, k, start = par_hat)
  if (is.null(fb) || fb$convergence != 0) next
  bq[b] <- CL_quantile(fb$par, L)
  bm[b] <- CL_moment(fb$par, L)
}
cat("Successful bootstrap refits:", sum(is.finite(bq)), "of", B_BOOT, "\n")

## ---------------------------------------------------------------------------
## Table 9
## ---------------------------------------------------------------------------

tab <- rbind(
  ACI = c(ci_asymptotic(Cm, se_m, LEVEL), ci_asymptotic(Cq, se_q, LEVEL)),
  PB  = c(ci_percentile(bm, LEVEL),        ci_percentile(bq, LEVEL)),
  NB  = c(ci_normal_boot(Cm, bm, LEVEL),   ci_normal_boot(Cq, bq, LEVEL)))
colnames(tab) <- c("M_Lower", "M_Upper", "Q_Lower", "Q_Upper")
tab <- cbind(tab,
             M_Length = tab[, 2] - tab[, 1],
             Q_Length = tab[, 4] - tab[, 3])

cat("\n95% confidence intervals (Table 9)\n")
print(round(tab[, c("M_Lower", "M_Upper", "M_Length",
                    "Q_Lower", "Q_Upper", "Q_Length")], 4))

## ---------------------------------------------------------------------------
## Section 6.3: the one-sided decision rule, Eq. (22)
## ---------------------------------------------------------------------------

cat("\n--- Decision rule of Section 6.3 (quantile-based index only) ---\n")
lb <- lower_bound_one_sided(Cq, se_q, LEVEL)
cat(sprintf("One-sided 95%% lower bound: %.4f - %.4f * %.4f = %.4f\n",
            Cq, qnorm(LEVEL), se_q, lb))

for (C0M in c(0.5)) {
  C0xi <- kd["kappa"] * C0M - kd["delta"]
  cat(sprintf("\nPolicy stated on the moment scale as C_0^M = %.2f\n", C0M))
  cat(sprintf("  translated threshold C_0^xi = %.4f\n", C0xi))
  cat(sprintf("  reject H0 (process capable)? %s   (z = %.4f, p = %.4f)\n",
              lb > C0xi, (Cq - C0xi) / se_q, 1 - pnorm((Cq - C0xi) / se_q)))
}
for (C0xi in c(0.5)) {
  cat(sprintf("\nPolicy stated directly on the quantile scale as C_0^xi = %.2f\n", C0xi))
  cat(sprintf("  reject H0 (process capable)? %s   (z = %.4f, p = %.4f)\n",
              lb > C0xi, (Cq - C0xi) / se_q, 1 - pnorm((Cq - C0xi) / se_q)))
}

## One-sided bootstrap analogues, quoted at the end of Section 6.3.
lb_pb <- unname(quantile(bq[is.finite(bq)], 1 - LEVEL))
lb_nb <- (Cq - (mean(bq, na.rm = TRUE) - Cq)) - qnorm(LEVEL) * sd(bq, na.rm = TRUE)
cat(sprintf("\nOne-sided 95%% lower bounds from the bootstrap:\n"))
cat(sprintf("  percentile (5th pct of replicates): %.4f\n", lb_pb))
cat(sprintf("  bias-corrected normal:             %.4f\n", lb_nb))

saveRDS(list(par = par_hat, se = c(se_alpha, se_lambda, se_q, se_m),
             Cq = Cq, Cm = Cm, table = tab, boot_q = bq, boot_m = bm),
        file.path(RESULTS_DIR, "realdata1.rds"))
cat("\nSaved", file.path(RESULTS_DIR, "realdata1.rds"), "\n")

## ---------------------------------------------------------------------------
## Table 9 as a LaTeX fragment
## ---------------------------------------------------------------------------

con <- base::file(file.path(TABLES_DIR, "tab_realdata1.tex"), open = "wt")
wl  <- function(...) writeLines(paste0(...), con)
wl("\\begin{table}[H]"); wl("\\centering")
wl(sprintf(paste0("\\caption{Point estimates and 95\\%% confidence intervals for ",
                  "$C_{L}^{M}$ and $C_{L}^{\\xi}$ using the ball bearings data ",
                  "($L=%g$), with $B=%d$ bootstrap resamples.}"), L, B_BOOT))
wl("\\label{T3}")
wl("\\begin{tabular}{l ccc c ccc}"); wl("\\toprule")
wl("\\multirow{2}{*}{CI Method} & \\multicolumn{3}{c}{Moment-based ($C_{L}^{M}$)} & ",
   "\\phantom{abc} & \\multicolumn{3}{c}{Quantile-based ($C_{L}^{\\xi}$)} \\\\")
wl("\\cmidrule(r){2-4} \\cmidrule(l){6-8}")
wl(" & Lower & Upper & Length && Lower & Upper & Length \\\\")
wl("\\midrule")
wl(sprintf("Point Estimate & \\multicolumn{3}{c}{%.4f} && \\multicolumn{3}{c}{%.4f} \\\\",
           Cm, Cq))
wl("\\midrule")
for (i in seq_len(nrow(tab)))
  wl(sprintf("%-3s & %.4f & %.4f & %.4f && %.4f & %.4f & %.4f \\\\",
             rownames(tab)[i], tab[i, "M_Lower"], tab[i, "M_Upper"], tab[i, "M_Length"],
             tab[i, "Q_Lower"], tab[i, "Q_Upper"], tab[i, "Q_Length"]))
wl("\\bottomrule"); wl("\\end{tabular}"); wl("\\end{table}")
close(con)
cat("Wrote", file.path(TABLES_DIR, "tab_realdata1.tex"), "\n")

## ---------------------------------------------------------------------------
## Numbers quoted in the running text of Sections 6.1 and 6.3
## ---------------------------------------------------------------------------

lbound  <- lower_bound_one_sided(Cq, se_q, LEVEL)
CzeroM  <- 0.5
CzeroXi <- unname(kd["kappa"] * CzeroM - kd["delta"])
pTrans  <- 1 - pnorm((Cq - CzeroXi) / se_q)
pDirect <- 1 - pnorm((Cq - CzeroM)  / se_q)

write_macros(file.path(TABLES_DIR, "values_realdata1.tex"), list(
  bbM            = as.character(m),
  bbK            = as.character(k),
  bbL            = as.character(L),
  bbAlphaHat     = fmt(par_hat[1], 4),
  bbAlphaSE      = fmt(se_alpha, 4),
  bbLambdaHat    = fmt(par_hat[2], 2),
  bbLambdaSE     = fmt(se_lambda, 4),
  bbCM           = fmt(Cm, 4),
  bbCMSE         = fmt(se_m, 4),
  bbCX           = fmt(Cq, 4),
  bbCXSE         = fmt(se_q, 4),
  bbKappa        = fmt(kd["kappa"], 4),
  bbDelta        = fmt(kd["delta"], 4),
  bbRatio        = fmt(L / par_hat[2], 3),
  bbCzeroM       = fmt(CzeroM, 1),
  bbCzeroXi      = fmt(CzeroXi, 4),
  bbZlevel       = fmt(qnorm(LEVEL), 4),
  bbLowerBound   = fmt(lbound, 3),
  bbPBlower      = fmt(lb_pb, 3),
  bbNBlower      = fmt(lb_nb, 3),
  bbPtranslated  = fmt(pTrans, 3),
  bbPdirect      = fmt(pDirect, 3),
  bbBootB        = as.character(B_BOOT)
), "10_realdata1_ballbearings.R")
