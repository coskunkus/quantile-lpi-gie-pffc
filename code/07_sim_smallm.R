################################################################################
##  07_sim_smallm.R
##
##  Small-sample study matched to the first application, Section 5.4.
##  Produces Table 8 and writes it as a LaTeX fragment into ../tables/.
##
##  The design copies the ball bearing application exactly: the group size,
##  the censoring plan and the parameters are read from the macro file that
##  10_realdata1_ballbearings.R writes, so this block cannot drift out of step
##  with the application it mirrors.  m is varied over {9,12,15,20,25} so that the transition from
##  the sample size actually used (m = 9) up to the smallest size in the main
##  design (m = 25) can be seen.
##
##  Runtime: roughly 30-60 minutes on a single core; m is small here.
################################################################################

source("00_gie_pffc.R")
need("numDeriv", "statmod", "progressr")

set.seed(MASTER_SEED + 7L)

N_REP_PT <- 4000
N_REP_CI <- 2000
B_BOOT   <- 250
## The design constants are read from the macro file written by
## 10_realdata1_ballbearings.R rather than copied here, so that this block
## cannot drift out of step with the application it is supposed to mirror.
## Run 10_realdata1_ballbearings.R first.
app <- read_macros(file.path(TABLES_DIR, "values_realdata1.tex"))
ALPHA    <- as.numeric(app$bbAlphaHat)
LAMBDA   <- as.numeric(app$bbLambdaHat)
L_LIMIT  <- as.numeric(app$bbL)
K        <- as.numeric(app$bbK)
R_FIRST  <- as.numeric(app$bbRfirst)
cat(sprintf("Design read from the application: alpha = %g, lambda = %g, L = %g, k = %g, R_1 = %g\n",
            ALPHA, LAMBDA, L_LIMIT, K, R_FIRST))
M_SET    <- c(9, 12, 15, 20, 25)
LEVEL    <- 0.95
N_TICK_PT <- 5
N_TICK_CI <- 20
TICK_PT   <- max(1L, N_REP_PT %/% N_TICK_PT)
TICK_CI   <- max(1L, N_REP_CI %/% N_TICK_CI)
N_TICK_PT <- N_REP_PT %/% TICK_PT   # counts actually delivered
N_TICK_CI <- N_REP_CI %/% TICK_CI

C_true <- CL_quantile(c(ALPHA, LAMBDA), L_LIMIT)
cat(sprintf("True index at the application design: C_L^xi = %.4f  (L/lambda = %.4f)\n",
            C_true, L_LIMIT / LAMBDA))

res <- with_progress_bar(length(M_SET) * (N_TICK_PT + N_TICK_CI), function(p) {
lapply(M_SET, function(m) {
  R <- c(R_FIRST, rep(0, m - 1))
  true <- c(ALPHA, LAMBDA)

  Ch <- rep(NA_real_, N_REP_PT)
  for (i in seq_len(N_REP_PT)) {
    tick_if(p, i, TICK_PT)
    x  <- generate_pffc(m, K, R, ALPHA, LAMBDA)
    ft <- fit_mle(x, R, K, start = true)
    if (is.null(ft)) next
    Ch[i] <- CL_quantile(ft$par, L_LIMIT)
  }

  cov <- len <- matrix(NA_real_, N_REP_CI, 3,
                       dimnames = list(NULL, c("ACI", "PB", "NB")))
  ## The decision rule of Section 6.3 is one-sided, and the coverage of a
  ## two-sided interval says nothing about the level of a one-sided test: a
  ## two-sided 95% interval whose misses are lopsided can have a lower endpoint
  ## that exceeds the true value far more, or far less, than 5% of the time.
  ## So the one-sided lower bound at level 1 - LEVEL is recorded separately.
  ## one[, j] = 1 when the bound lies below the true index, i.e. when the
  ## one-sided interval covers; its complement is the type-I error rate of the
  ## corresponding test at a threshold equal to the true value.
  one <- matrix(NA_real_, N_REP_CI, 3,
                dimnames = list(NULL, c("ACI", "PB", "NB")))
  for (i in seq_len(N_REP_CI)) {
    tick_if(p, i, TICK_CI)
    x  <- generate_pffc(m, K, R, ALPHA, LAMBDA)
    ft <- fit_mle(x, R, K, start = true)
    if (is.null(ft)) next
    Chat <- CL_quantile(ft$par, L_LIMIT)
    if (!is.finite(Chat)) next

    se <- delta_se(ft, L_LIMIT, "quantile")
    if (is.finite(se)) {
      ci <- ci_asymptotic(Chat, se, LEVEL)
      cov[i, "ACI"] <- as.numeric(ci[1] < C_true && ci[2] > C_true)
      len[i, "ACI"] <- diff(ci)
      one[i, "ACI"] <- as.numeric(lower_bound_one_sided(Chat, se, LEVEL) < C_true)
    }
    bt <- bootstrap_index(ft$par, m, K, R, L_LIMIT, B_BOOT, index = "quantile")
    if (sum(is.finite(bt)) >= B_BOOT / 2) {
      ci <- ci_percentile(bt, LEVEL)
      cov[i, "PB"] <- as.numeric(ci[1] < C_true && ci[2] > C_true)
      len[i, "PB"] <- diff(ci)
      ci <- ci_normal_boot(Chat, bt, LEVEL)
      cov[i, "NB"] <- as.numeric(ci[1] < C_true && ci[2] > C_true)
      len[i, "NB"] <- diff(ci)
      ## one-sided bootstrap bounds, built the same way as in Section 6.3
      bt_ok <- bt[is.finite(bt)]
      one[i, "PB"] <- as.numeric(unname(quantile(bt_ok, 1 - LEVEL)) < C_true)
      one[i, "NB"] <- as.numeric((Chat - (mean(bt_ok) - Chat)) -
                                 qnorm(LEVEL) * sd(bt_ok) < C_true)
    }
  }

  data.frame(
    m = m, C_true = C_true,
    Bias = mean(Ch, na.rm = TRUE) - C_true,
    MSE  = mean((Ch - C_true)^2, na.rm = TRUE),
    AL_ACI = mean(len[, "ACI"], na.rm = TRUE), CP_ACI = mean(cov[, "ACI"], na.rm = TRUE),
    AL_PB  = mean(len[, "PB"],  na.rm = TRUE), CP_PB  = mean(cov[, "PB"],  na.rm = TRUE),
    AL_NB  = mean(len[, "NB"],  na.rm = TRUE), CP_NB  = mean(cov[, "NB"],  na.rm = TRUE),
    ## one-sided lower-bound coverage; 1 - this is the type-I error rate of the
    ## one-sided test of Section 6.3 at a threshold equal to the true index
    CP1_ACI = mean(one[, "ACI"], na.rm = TRUE),
    CP1_PB  = mean(one[, "PB"],  na.rm = TRUE),
    CP1_NB  = mean(one[, "NB"],  na.rm = TRUE),
    ## As in 05_sim_ci.R: replications that failed the positive-definiteness
    ## check or the bootstrap-convergence check are excluded from the means, so
    ## these counts are needed to read the coverages correctly.
    n_pt = sum(is.finite(Ch)), n_ACI = sum(is.finite(cov[, "ACI"])),
    n_boot = sum(is.finite(cov[, "PB"])), n_rep_ci = N_REP_CI)
})
})
out <- do.call(rbind, res)

for (i in seq_len(nrow(out)))
  cat(sprintf("m=%2d | bias %+.4f mse %.4f | ACI %.4f/%.4f PB %.4f/%.4f NB %.4f/%.4f\n",
              out$m[i], out$Bias[i], out$MSE[i], out$AL_ACI[i], out$CP_ACI[i],
              out$AL_PB[i], out$CP_PB[i], out$AL_NB[i], out$CP_NB[i]))
cat("\nOne-sided 95% lower-bound coverage (1 - type-I error of the Section 6.3 rule):\n")
for (i in seq_len(nrow(out)))
  cat(sprintf("m=%2d | ACI %.4f  PB %.4f  NB %.4f\n",
              out$m[i], out$CP1_ACI[i], out$CP1_PB[i], out$CP1_NB[i]))
write.csv(out, file.path(RESULTS_DIR, "sim_smallm.csv"), row.names = FALSE)

con <- base::file(file.path(TABLES_DIR, "tab_smallm.tex"), open = "wt")
wl <- function(...) writeLines(paste0(...), con)
wl("\\begin{table}[H]"); wl("\\centering")
wl(sprintf(paste0("\\caption{Small-sample performance at the design of the first ",
                  "application: $k=%d$, $\\mathbf{R}=(%g,0^{m-1})$, $\\alpha=%.4f$, ",
                  "$\\lambda=%.2f$, $L=%g$, so that $L/\\lambda=%.3f$ and ",
                  "$C_L^{\\xi}=%.4f$. Bias and MSE from %d replications; AL and CP ",
                  "from %d replications with $B=%d$.}"),
           K, R_FIRST, ALPHA, LAMBDA, L_LIMIT, L_LIMIT / LAMBDA, C_true,
           N_REP_PT, N_REP_CI, B_BOOT))
wl("\\label{T:smallm}"); wl("\\small")
wl("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}l rr rr rr rr}"); wl("\\toprule")
wl(" & \\multicolumn{2}{c}{Point estimator} & \\multicolumn{2}{c}{ACI} & ",
   "\\multicolumn{2}{c}{PB} & \\multicolumn{2}{c}{NB} \\\\")
wl("\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9}")
wl("$m$ & \\multicolumn{1}{c}{Bias} & \\multicolumn{1}{c}{MSE}",
   paste(rep(" & \\multicolumn{1}{c}{AL} & \\multicolumn{1}{c}{CP}", 3), collapse = ""), " \\\\")
wl("\\midrule")
for (i in seq_len(nrow(out)))
  wl(sprintf("%d & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
             out$m[i], out$Bias[i], out$MSE[i],
             out$AL_ACI[i], out$CP_ACI[i], out$AL_PB[i], out$CP_PB[i],
             out$AL_NB[i], out$CP_NB[i]))
wl("\\bottomrule"); wl("\\end{tabular*}"); wl("\\end{table}")
close(con)

## Numbers quoted in the running text of Sections 5.4 and 6.1.
r9 <- out[out$m == min(M_SET), ]
write_macros(file.path(TABLES_DIR, "values_smallm.tex"), list(
  smNrepPt   = as.character(N_REP_PT),
  smNrepCi   = as.character(N_REP_CI),
  smBootB    = as.character(B_BOOT),
  smAlpha    = fmt(ALPHA, 4),
  smLambda   = fmt(LAMBDA, 2),
  smL        = as.character(L_LIMIT),
  smK        = as.character(K),
  smRfirst   = as.character(R_FIRST),
  smRatio    = fmt(L_LIMIT / LAMBDA, 3),
  smCtrue    = fmt(C_true, 4),
  smMmin     = as.character(min(M_SET)),
  smMmax     = as.character(max(M_SET)),
  smMset     = paste(M_SET, collapse = ","),
  smCPACI    = fmt(r9$CP_ACI, 3),
  smCPPB     = fmt(r9$CP_PB, 3),
  smCPNB     = fmt(r9$CP_NB, 3),
  smOneACI   = fmt(r9$CP1_ACI, 3),
  smOnePB    = fmt(r9$CP1_PB, 3),
  smOneNB    = fmt(r9$CP1_NB, 3)
), "07_sim_smallm.R")

cat("\nWrote", file.path(TABLES_DIR, "tab_smallm.tex"), "\n")
