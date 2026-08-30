################################################################################
##  11_realdata2_guineapig.R
##
##  Second application, Section 6.2: survival times of guinea pigs injected with
##  tubercle bacilli (Kinaci, Wu & Kus 2019).
##
##  Here alpha-hat < 2, so the moment-based index does not exist and only the
##  quantile-based index is analysed.  Reproduces Table 10.
##
##  Runtime: under a minute.
################################################################################

source("00_gie_pffc.R")
need("numDeriv", "statmod")

set.seed(MASTER_SEED + 11L)

B_BOOT <- 1000
LEVEL  <- 0.95

## Complete data: 72 survival times in days.
full <- c(12, 15, 22, 24, 24, 32, 32, 33, 34, 38, 38, 43, 44, 48, 52, 53, 54,
          54, 55, 56, 57, 58, 58, 59, 60, 60, 60, 60, 61, 62, 63, 65, 65, 67,
          68, 70, 70, 72, 73, 75, 76, 76, 81, 83, 84, 85, 87, 91, 95, 96, 98,
          99, 109, 110, 121, 127, 129, 131, 143, 146, 146, 175, 175, 211, 233,
          258, 258, 263, 297, 341, 341, 376)

## PFFC sample with m = 26, k = 2, R = (10, 0^25).
x <- c(12, 15, 22, 24, 32, 34, 38, 38, 44, 53, 54, 54, 55, 56, 57, 58,
       65, 67, 70, 73, 81, 98, 109, 110, 131, 258)
m <- length(x)
k <- 2
R <- c(10, rep(0, m - 1))
L <- 10

cat("Guinea pig data:  m =", m, " k =", k, " n =", m + sum(R), " L =", L, "\n\n")

fit <- fit_mle(x, R, k, numeric_hessian = TRUE)
stopifnot(!is.null(fit))
par_hat <- fit$par
I_obs   <- fit$information          # analytic, Eq. (17)

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

V <- solve(I_obs)
se_alpha <- sqrt(V[1, 1]); se_lambda <- sqrt(V[2, 2])
Cq   <- CL_quantile(par_hat, L)
se_q <- delta_se(fit, L, "quantile")

cat(sprintf("alpha-hat = %.4f (SE %.4f);  lambda-hat = %.4f (SE %.4f)\n",
            par_hat[1], se_alpha, par_hat[2], se_lambda))
ci_alpha <- par_hat[1] + c(-1, 1) * qnorm(0.975) * se_alpha
cat(sprintf("95%% CI for alpha: (%.3f, %.3f)\n", ci_alpha[1], ci_alpha[2]))
cat(sprintf("Moment-based index is undefined here (alpha-hat = %.4f <= 2).\n",
            par_hat[1]))
cat(sprintf("C_L^xi = %.4f (SE %.4f);  L / lambda-hat = %.4f\n\n",
            Cq, se_q, L / par_hat[2]))

cat("Parametric bootstrap, B =", B_BOOT, "...\n")
bq <- rep(NA_real_, B_BOOT)
for (b in seq_len(B_BOOT)) {
  y  <- generate_pffc(m, k, R, par_hat[1], par_hat[2])
  fb <- fit_mle(y, R, k, start = par_hat)
  if (is.null(fb) || fb$convergence != 0) next
  bq[b] <- CL_quantile(fb$par, L)
}
cat("Successful bootstrap refits:", sum(is.finite(bq)), "of", B_BOOT, "\n")

tab <- rbind(ACI = ci_asymptotic(Cq, se_q, LEVEL),
             PB  = ci_percentile(bq, LEVEL),
             NB  = ci_normal_boot(Cq, bq, LEVEL))
colnames(tab) <- c("Lower", "Upper")
tab <- cbind(tab, Length = tab[, 2] - tab[, 1])

cat("\n95% confidence intervals for C_L^xi (Table 10)\n")
print(round(tab, 4))

cat(sprintf("\nOne-sided 95%% lower bound (Eq. 22): %.4f\n",
            lower_bound_one_sided(Cq, se_q, LEVEL)))

saveRDS(list(par = par_hat, se = c(se_alpha, se_lambda, se_q),
             Cq = Cq, table = tab, boot_q = bq),
        file.path(RESULTS_DIR, "realdata2.rds"))
cat("Saved", file.path(RESULTS_DIR, "realdata2.rds"), "\n")

## ---------------------------------------------------------------------------
## Table 10 as a LaTeX fragment
## ---------------------------------------------------------------------------

con <- base::file(file.path(TABLES_DIR, "tab_realdata2.tex"), open = "wt")
wl  <- function(...) writeLines(paste0(...), con)
wl("\\begin{table}[H]"); wl("\\centering")
wl(sprintf(paste0("\\caption{Point estimate and 95\\%% confidence intervals for ",
                  "$C_{L}^{\\xi}$ using the guinea pig survival data ($L=%g$), with ",
                  "$B=%d$ bootstrap resamples. The moment-based index is undefined ",
                  "here because $\\widehat{\\alpha}<2$.}"), L, B_BOOT))
wl("\\label{T4}")
wl("\\begin{tabular}{l ccc}"); wl("\\toprule")
wl("\\multirow{2}{*}{CI Method} & \\multicolumn{3}{c}{Quantile-based ($C_{L}^{\\xi}$)} \\\\")
wl("\\cmidrule(r){2-4}")
wl(" & Lower & Upper & Length \\\\")
wl("\\midrule")
wl(sprintf("Point Estimate & \\multicolumn{3}{c}{%.4f} \\\\", Cq))
wl("\\midrule")
for (i in seq_len(nrow(tab)))
  wl(sprintf("%-3s & %.4f & %.4f & %.4f \\\\",
             rownames(tab)[i], tab[i, "Lower"], tab[i, "Upper"], tab[i, "Length"]))
wl("\\bottomrule"); wl("\\end{tabular}"); wl("\\end{table}")
close(con)
cat("Wrote", file.path(TABLES_DIR, "tab_realdata2.tex"), "\n")

## Numbers quoted in the running text of Section 6.2.
write_macros(file.path(TABLES_DIR, "values_realdata2.tex"), list(
  gpM         = as.character(m),
  gpK         = as.character(k),
  gpL         = as.character(L),
  gpAlphaHat  = fmt(par_hat[1], 4),
  gpAlphaSE   = fmt(se_alpha, 4),
  gpLambdaHat = fmt(par_hat[2], 4),
  gpLambdaSE  = fmt(se_lambda, 4),
  gpAlphaCIlo = fmt(ci_alpha[1], 2),
  gpAlphaCIhi = fmt(ci_alpha[2], 2),
  gpCX        = fmt(Cq, 4),
  gpCXSE      = fmt(se_q, 4),
  gpRatio     = fmt(L / par_hat[2], 3),
  gpBootB     = as.character(B_BOOT),
  ## the three interval lengths, so that Section 6.2 can quote how far apart
  ## they are without the figure being typed into the manuscript
  gpALACI     = fmt(tab["ACI", "Length"], 4),
  gpALPB      = fmt(tab["PB",  "Length"], 4),
  gpALNB      = fmt(tab["NB",  "Length"], 4)
), "11_realdata2_guineapig.R")
