################################################################################
##  04_sim_bias_mse.R
##
##  Baseline simulation study, Section 5.1.
##  Produces Tables 2 and 3 (Bias and MSE of the two indexes) and Table 4
##  (the scale-free comparison between them).
##
##  Design:  m in {25,50,100}, k in {2,5}, alpha in {5,10}, lambda = 2, L = 2,
##           four censoring plans, sum(R) = m so that n = 2m.
##  2000 Monte Carlo replications per configuration.
##
##  Runtime: roughly 20-40 minutes on a single core.  Set N_REP smaller for a
##  quick check; the numbers in the paper use N_REP = 2000.
##
##  REUSE MODE.  The tables and the macro file are written from the data frame
##  the simulation produces, and that data frame is also saved to
##  results/sim_bias_mse.csv.  If the option gie.reuse_results is set, the saved
##  file is read instead of re-running the study, and everything after it is
##  unchanged.  That is how run_paper.R rebuilds Tables 2 to 4 without spending
##  half an hour on a simulation whose output is already shipped.  Running the
##  script on its own simulates, as before.
################################################################################

source("00_gie_pffc.R")
need("numDeriv", "statmod", "progressr")

set.seed(MASTER_SEED + 4L)

N_REP     <- 2000
LAMBDA    <- 2.0
L_LIMIT   <- 2.0
M_SET     <- c(25, 50, 100)
K_SET     <- c(2, 5)
ALPHA_SET <- c(5, 10)
SCHEMES   <- c("Early", "Middle", "Late", "Equal")
N_TICK    <- 20                  # progress-bar ticks per configuration
TICK_EVERY <- max(1L, N_REP %/% N_TICK)
N_TICK    <- N_REP %/% TICK_EVERY   # the count actually delivered

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------

REUSE <- isTRUE(getOption("gie.reuse_results", FALSE))
CSV   <- file.path(RESULTS_DIR, "sim_bias_mse.csv")
if (REUSE && !file.exists(CSV))
  stop("04_sim_bias_mse.R: reuse was asked for but ", CSV, " does not exist.\n",
       "  Run this script without the option, or run run_simulations.R.")

n_cell <- length(M_SET) * length(K_SET) * length(ALPHA_SET) * length(SCHEMES)
if (!REUSE)
cat(sprintf("%d configurations, %d replications each.\n", n_cell, N_REP))
utils::flush.console()
t_start <- Sys.time()

run_cell <- function(m, k, alpha, sch, p = NULL) {
  R    <- make_scheme(m, sch)
  true <- c(alpha, LAMBDA)
  Cq_true <- CL_quantile(true, L_LIMIT)
  Cm_true <- CL_moment(true, L_LIMIT)

  Cq <- Cm <- rep(NA_real_, N_REP)
  n_fp <- 0L                       # replications solved by Eq. (15) alone
  for (i in seq_len(N_REP)) {
    tick_if(p, i, TICK_EVERY)
    x  <- generate_pffc(m, k, R, alpha, LAMBDA)
    ft <- fit_mle(x, R, k, start = true)
    if (is.null(ft)) next
    if (identical(ft$route, "fixedpoint")) n_fp <- n_fp + 1L
    Cq[i] <- CL_quantile(ft$par, L_LIMIT)
    Cm[i] <- CL_moment(ft$par, L_LIMIT)
  }

  data.frame(
    m = m, k = k, alpha = alpha, scheme = sch,
    Cq_true = Cq_true, Cm_true = Cm_true,
    n_ok_q  = sum(is.finite(Cq)), n_ok_m = sum(is.finite(Cm)),
    n_fixedpoint = n_fp,
    Bias_q  = mean(Cq, na.rm = TRUE) - Cq_true,
    MSE_q   = mean((Cq - Cq_true)^2, na.rm = TRUE),
    Bias_m  = mean(Cm, na.rm = TRUE) - Cm_true,
    MSE_m   = mean((Cm - Cm_true)^2, na.rm = TRUE),
    stringsAsFactors = FALSE)
}

cells <- expand.grid(scheme = SCHEMES, alpha = ALPHA_SET, k = K_SET, m = M_SET,
                     stringsAsFactors = FALSE)

if (REUSE) {

  cat("Reading the saved study from ", CSV, " rather than re-running it.\n",
      sep = "")
  out <- read.csv(CSV, stringsAsFactors = FALSE)
  ## The saved file must describe the design this script declares, or the
  ## tables written below would carry the wrong row labels.
  have <- out[, c("m", "k", "alpha", "scheme")]
  want <- cells[, c("m", "k", "alpha", "scheme")]
  key  <- function(d) sort(do.call(paste, c(d[c("m", "k", "alpha", "scheme")],
                                            sep = "|")))
  if (!identical(key(have), key(want)))
    stop("04_sim_bias_mse.R: ", basename(CSV), " does not match the design in ",
         "this script.\n  Delete it and run run_simulations.R.")
  out <- out[order(out$m, out$k, out$alpha, match(out$scheme, SCHEMES)), ]

} else {

  res <- with_progress_bar(n_cell * N_TICK, function(p) {
    lapply(seq_len(nrow(cells)), function(i)
      run_cell(cells$m[i], cells$k[i], cells$alpha[i], cells$scheme[i], p))
  })

  cat(sprintf("\nFinished in %.1f minutes.\n",
              as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

  out <- do.call(rbind, res)
  out <- out[order(out$m, out$k, out$alpha, match(out$scheme, SCHEMES)), ]

  for (i in seq_len(nrow(out)))
    cat(sprintf("m=%3d k=%d alpha=%2d %-6s | q: %+.4f %.4f | m: %+.4f %.4f\n",
                out$m[i], out$k[i], out$alpha[i], out$scheme[i],
                out$Bias_q[i], out$MSE_q[i], out$Bias_m[i], out$MSE_m[i]))

  ## Relative (scale-free) measures used in Table 4 and Figure 4.
  out$RB_q   <- out$Bias_q / abs(out$Cq_true)
  out$RMSE_q <- out$MSE_q  / out$Cq_true^2
  out$RB_m   <- out$Bias_m / abs(out$Cm_true)
  out$RMSE_m <- out$MSE_m  / out$Cm_true^2
  out$ratio  <- out$RMSE_m / out$RMSE_q

  write.csv(out, CSV, row.names = FALSE)
}

cat("\nTrue index values used:\n")
print(unique(out[, c("alpha", "Cm_true", "Cq_true")]))
cat("\nQuantile index has smaller relative MSE in",
    sum(out$ratio > 1), "of", nrow(out), "configurations.\n")
cat("Quantile index has smaller |relative bias| in",
    sum(abs(out$RB_q) < abs(out$RB_m)), "of", nrow(out), "configurations.\n")
fp_rate <- sum(out$n_fixedpoint) / (nrow(out) * N_REP)
cat(sprintf("The fixed-point iteration of Eq. (15) attained the maximum on its own in %.1f%% of the %d fits.\n",
            100 * fp_rate, nrow(out) * N_REP))

## ---------------------------------------------------------------------------
## LaTeX fragments
## ---------------------------------------------------------------------------

ord <- order(out$m, out$k, out$alpha,
             match(out$scheme, SCHEMES))
o   <- out[ord, ]

## ---- Tables 2 and 3: Bias and MSE, one table per index ---------------------
## Layout: rows are (m, k, alpha); the four censoring schemes run across.

write_bias_mse_table <- function(o, bias_col, mse_col, path, label, caption) {
  con <- base::file(path, open = "wt"); wl <- function(...) writeLines(paste0(...), con)
  wl("\\begin{table}[H]"); wl("\\centering")
  wl("\\caption{", caption, "}"); wl("\\label{", label, "}"); wl("\\small")
  wl("\\begin{tabular}{lll rr rr rr rr}"); wl("\\toprule")
  wl(" & & & \\multicolumn{2}{c}{Early} & \\multicolumn{2}{c}{Middle} & ",
     "\\multicolumn{2}{c}{Late} & \\multicolumn{2}{c}{Equal} \\\\")
  wl("\\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9} \\cmidrule(lr){10-11}")
  wl("$m$ & $k$ & $\\alpha$",
     paste(rep(" & \\multicolumn{1}{c}{Bias} & \\multicolumn{1}{c}{MSE}", 4), collapse = ""),
     " \\\\")
  wl("\\midrule")
  first_m <- TRUE
  for (m in M_SET) {
    if (!first_m) wl("\\midrule")
    first_m <- FALSE
    for (ki in seq_along(K_SET)) {
      k <- K_SET[ki]
      if (ki > 1) wl("\\addlinespace")
      for (ai in seq_along(ALPHA_SET)) {
        a <- ALPHA_SET[ai]
        cells <- vapply(SCHEMES, function(s) {
          r <- o[o$m == m & o$k == k & o$alpha == a & o$scheme == s, ]
          sprintf("%.3f & %.3f", r[[bias_col]], r[[mse_col]])
        }, character(1))
        wl(sprintf("%s & %s & %d & %s \\\\",
                   if (ki == 1 && ai == 1) m else "",
                   if (ai == 1) k else "", a, paste(cells, collapse = " & ")))
      }
    }
  }
  wl("\\bottomrule"); wl("\\end{tabular}"); wl("\\end{table}")
  close(con)
}

write_bias_mse_table(o, "Bias_m", "MSE_m",
                     file.path(TABLES_DIR, "tab_sim_moment.tex"),
                     "T:sim_moment",
                     "Bias and MSE of the MLE for the Moment-based Index.")
write_bias_mse_table(o, "Bias_q", "MSE_q",
                     file.path(TABLES_DIR, "tab_sim_quantile.tex"),
                     "T:sim_quantile",
                     "Bias and MSE of the MLE for the Quantile-based Index.")

## ---- Table 4: the scale-free comparison ------------------------------------

con <- base::file(file.path(TABLES_DIR, "tab_relative.tex"), open = "wt")
wl  <- function(...) writeLines(paste0(...), con)
wl("\\begin{table}[H]")
wl("\\centering")
wl("\\caption{Scale-free comparison of the two estimators: relative bias ",
   "$\\mathrm{RB}=\\mathrm{Bias}(\\widehat{C})/|C|$ and relative mean squared ",
   "error $\\mathrm{RelMSE}=\\mathrm{MSE}(\\widehat{C})/C^{2}$, computed from ",
   "Tables \\ref{T:sim_moment} and \\ref{T:sim_quantile}. Smaller values ",
   "indicate better relative performance. The final column is the ratio of the ",
   "two relative MSEs; values above one favour the quantile-based index.}")
wl("\\label{T:relative}")
wl("\\scriptsize")
wl("\\begin{tabular}{lll l rr rr r}")
wl("\\toprule")
wl(" & & & & \\multicolumn{2}{c}{Moment-based ($C_L^{M}$)} & ",
   "\\multicolumn{2}{c}{Quantile-based ($C_L^{\\xi}$)} & \\\\")
wl("\\cmidrule(lr){5-6} \\cmidrule(lr){7-8}")
wl("$m$ & $k$ & $\\alpha$ & Scheme & \\multicolumn{1}{c}{RB} & ",
   "\\multicolumn{1}{c}{RelMSE} & \\multicolumn{1}{c}{RB} & ",
   "\\multicolumn{1}{c}{RelMSE} & \\multicolumn{1}{c}{Ratio} \\\\")
wl("\\midrule")
blk <- 0
for (i in seq_len(nrow(o))) {
  first <- (i - 1) %% 4 == 0
  if (first && i > 1) {
    blk <- blk + 1
    wl(if (blk %% 2 == 0) "\\midrule" else "\\addlinespace")
  }
  c1 <- if (first && (i == 1 || o$m[i] != o$m[i - 1])) o$m[i] else ""
  c2 <- if (first && (i == 1 || o$k[i] != o$k[i - 1] || o$m[i] != o$m[i - 1])) o$k[i] else ""
  c3 <- if (first) o$alpha[i] else ""
  wl(sprintf("%s & %s & %s & %s & %.4f & %.4f & %.4f & %.4f & %.2f \\\\",
             c1, c2, c3, o$scheme[i],
             o$RB_m[i], o$RMSE_m[i], o$RB_q[i], o$RMSE_q[i], o$ratio[i]))
}
wl("\\bottomrule"); wl("\\end{tabular}"); wl("\\end{table}")
close(con)

## ---- Numbers quoted in the running text of Section 5.1 --------------------

write_macros(file.path(TABLES_DIR, "values_sim.tex"), list(
  simNrep        = as.character(N_REP),
  simLambda      = as.character(LAMBDA),
  simL           = as.character(L_LIMIT),
  simRatio       = as.character(L_LIMIT / LAMBDA),
  simCMfive      = fmt(out$Cm_true[out$alpha == 5][1], 4),
  simCMten       = fmt(out$Cm_true[out$alpha == 10][1], 4),
  simCXfive      = fmt(out$Cq_true[out$alpha == 5][1], 4),
  simCXten       = fmt(out$Cq_true[out$alpha == 10][1], 4),
  simNconfig     = as.character(nrow(out)),
  simNbetterMSE  = as.character(sum(out$ratio > 1)),
  simNbetterBias = as.character(sum(abs(out$RB_q) < abs(out$RB_m))),
  simRatioLo     = fmt(min(out$ratio), 1),
  simRatioHi     = fmt(max(out$ratio), 1),
  simFPrate      = fmt(100 * fp_rate, 1)
), "04_sim_bias_mse.R")

cat("\nWrote tab_sim_moment.tex, tab_sim_quantile.tex and tab_relative.tex to",
    TABLES_DIR, "\n")
cat(if (REUSE) "Read" else "Wrote", CSV, "\n")
