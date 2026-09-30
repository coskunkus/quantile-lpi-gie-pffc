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

## The complete data and the censored sample both come from
## guinea_pig_design() in 00_gie_pffc.R.  It is defined there, and not here,
## because validation/v5_check_information_and_mle.R checks the same fit and
## rebuilds it the same way; no script in this package reads a stored result.
## The construction touches no random stream, so the bootstrap below draws
## exactly what it drew before.
app  <- guinea_pig_design()
full <- app$full
x    <- app$x
m <- app$m; k <- app$k; R <- app$R; L <- app$L

cat("Guinea pig data:  m =", m, " k =", k, " n =", m + sum(R), " L =", L, "\n\n")

fit <- fit_mle(x, R, k, numeric_hessian = TRUE)
stopifnot(!is.null(fit))
par_hat <- fit$par
I_obs   <- fit$information          # analytic, Eq. (18)

cat(sprintf("Estimation route: %s", fit$route))
if (identical(fit$route, "fixedpoint"))
  cat(sprintf(" (Eq. (16) converged in %d iterations)", fit$iterations))
cat(sprintf("\nMaximised log-likelihood: %.8f\n", fit$loglik))

ev <- eigen(I_obs, symmetric = TRUE, only.values = TRUE)$values
cat("Observed information eigenvalues:", sprintf("%.4g", ev),
    if (min(ev) > 0) " (positive definite)\n" else " (NOT positive definite)\n")
if (is.null(fit$hessian_numeric)) {
  cat("optim did not return a numerical Hessian; the analytic information of\n",
      "  Eq. (18) is used regardless, and validation/v5_ checks it.\n", sep = "")
} else {
  cat(sprintf("Largest relative difference between the analytic information of Eq. (18)\n  and the numerical Hessian from optim: %.2e\n\n",
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

cat(sprintf("\nOne-sided 95%% lower bound (Eq. (24)): %.4f\n",
            lower_bound_one_sided(Cq, se_q, LEVEL)))

saveRDS(list(par = par_hat, se = c(se_alpha, se_lambda, se_q),
             Cq = Cq, table = tab, boot_q = bq),
        file.path(RESULTS_DIR, "realdata2.rds"))
cat("Saved", file.path(RESULTS_DIR, "realdata2.rds"), "\n")

## The design, the fit and the three interval lengths, as one row of csv.
## 09_values_text.R takes the spread of those lengths from here, so that the
## sentence Section 6.2 quotes is not typed in.  As elsewhere in this package,
## what one script passes to another goes through results/ and never through
## tables/, which holds output only.
write.csv(data.frame(
  m = m, k = k, R1 = R[1], L = L, n_full = length(full),
  sample = paste(fmt(x, 0), collapse = ", "),
  alpha_hat = par_hat[1], alpha_se = se_alpha,
  lambda_hat = par_hat[2], lambda_se = se_lambda,
  C_xi = Cq, C_xi_se = se_q,
  AL_ACI = tab["ACI", "Length"], AL_PB = tab["PB", "Length"],
  AL_NB = tab["NB", "Length"],
  B_boot = B_BOOT, level = LEVEL),
  file.path(RESULTS_DIR, "realdata2.csv"), row.names = FALSE)
cat("Saved", file.path(RESULTS_DIR, "realdata2.csv"), "\n")

## The three interval lengths, left in the session for 09_values_text.R,
## which quotes how far apart they are in Section 6.2.
APP_GUINEAPIG <- list(AL = as.numeric(tab[, "Length"]))

## ---------------------------------------------------------------------------
## Table 10 as a LaTeX fragment
## ---------------------------------------------------------------------------

con <- base::file(file.path(TABLES_DIR, "tab_realdata2.tex"), open = "wt")
wl  <- function(...) writeLines(paste0(...), con)
wl("\\begin{table}[H]"); wl("\\centering")
wl(sprintf(paste0("\\caption{Point estimate and 95\\%% CIs for ",
                  "$C_{L}^{\\xi}$ using the guinea pig survival data ($L=%g$), with ",
                  "$B=%d$ bootstrap resamples. The moment-based index is not ",
                  "reported here, since $\\widehat{\\alpha}<2$.}"), L, B_BOOT))
wl("\\label{T4}"); wl("\\scriptsize")
wl("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}l ccc}"); wl("\\toprule")
wl("\\multirow{2}{*}{CI Method} & \\multicolumn{3}{c}{Quantile-based ($C_{L}^{\\xi}$)} \\\\")
wl("\\cmidrule(r){2-4}")
wl(" & Lower & Upper & Length \\\\")
wl("\\midrule")
wl(sprintf("Point Estimate & \\multicolumn{3}{c}{%.4f} \\\\", Cq))
wl("\\midrule")
for (i in seq_len(nrow(tab)))
  wl(sprintf("%-3s & %.4f & %.4f & %.4f \\\\",
             rownames(tab)[i], tab[i, "Lower"], tab[i, "Upper"], tab[i, "Length"]))
wl("\\bottomrule"); wl("\\end{tabular*}"); wl("\\end{table}")
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
