################################################################################
##  05_sim_ci.R
##
##  Confidence interval simulation for the baseline design, Section 5.2.
##  Produces Table 5 (quantile-based index, main text) and the table of
##  Appendix C (moment-based index).
##
##  NOTE ON A CHANGE FROM THE EARLIER CODE.  An earlier version of this script
##  was left with ds <- 20 and B <- 20 after a debugging run, which does not
##  reproduce the published tables.  The values below (500 replications, 250
##  bootstrap resamples) are the ones stated in the paper and used to produce
##  Table 5.
##
##  Runtime: this is the heaviest of the four scripts, because the moment-based
##  index has to be evaluated for every bootstrap replicate as well.  It runs in
##  parallel by default; set PARALLEL <- FALSE for a single core.  Progress is
##  A progress bar is shown while it runs.
################################################################################

source("00_gie_pffc.R")
need("numDeriv", "statmod", "progressr")

PARALLEL <- TRUE           # set FALSE to run on a single core
N_REP    <- 500
B_BOOT   <- 250
LAMBDA   <- 2.0
L_LIMIT  <- 2.0
M_SET     <- c(25, 50, 100)
K_SET     <- c(2, 5)
ALPHA_SET <- c(5, 10)
SCHEMES   <- c("Early", "Middle", "Late", "Equal")
LEVEL     <- 0.95
N_TICK    <- 20                  # progress-bar ticks per configuration
TICK_EVERY <- max(1L, N_REP %/% N_TICK)
N_TICK     <- N_REP %/% TICK_EVERY   # the count actually delivered

## ---------------------------------------------------------------------------
## One configuration
## ---------------------------------------------------------------------------

run_cell <- function(m, k, alpha, sch, seed, p = NULL) {
  ## The RNG kind must be named explicitly.  furrr_options(seed = TRUE) installs
  ## an L'Ecuyer-CMRG .Random.seed in the worker before each element is
  ## evaluated, and set.seed() with no `kind` re-seeds whatever generator
  ## .Random.seed currently names.  Without this argument the workers would draw
  ## from L'Ecuyer while a serial run drew from Mersenne-Twister, and the two
  ## would silently produce different tables.
  set.seed(seed, kind = "Mersenne-Twister")
  R  <- make_scheme(m, sch)
  tp <- c(alpha, LAMBDA)
  Cq_true <- CL_quantile(tp, L_LIMIT)
  Cm_true <- CL_moment(tp, L_LIMIT)

  nm  <- c("ACI", "PB", "NB")
  cq  <- lq <- matrix(NA_real_, N_REP, 3, dimnames = list(NULL, nm))
  cm  <- lm <- matrix(NA_real_, N_REP, 3, dimnames = list(NULL, nm))

  for (i in seq_len(N_REP)) {
    ## Ticked at the top of the body, before any `next`, so that the bar always
    ## receives exactly N_TICK ticks per configuration and reaches 100%.
    tick_if(p, i, TICK_EVERY)
    x  <- generate_pffc(m, k, R, alpha, LAMBDA)
    fx <- fit_mle(x, R, k, start = tp)
    if (is.null(fx)) next
    Chq <- CL_quantile(fx$par, L_LIMIT)
    Chm <- CL_moment(fx$par, L_LIMIT)

    ## --- asymptotic intervals ---
    seq_ <- delta_se(fx, L_LIMIT, "quantile")
    if (is.finite(Chq) && is.finite(seq_)) {
      ci <- ci_asymptotic(Chq, seq_, LEVEL)
      cq[i, "ACI"] <- as.numeric(ci[1] < Cq_true && ci[2] > Cq_true); lq[i, "ACI"] <- diff(ci)
    }
    sem <- delta_se(fx, L_LIMIT, "moment")
    if (is.finite(Chm) && is.finite(sem)) {
      ci <- ci_asymptotic(Chm, sem, LEVEL)
      cm[i, "ACI"] <- as.numeric(ci[1] < Cm_true && ci[2] > Cm_true); lm[i, "ACI"] <- diff(ci)
    }

    ## --- bootstrap: refit once per resample, evaluate both indexes ---
    bq <- bm <- rep(NA_real_, B_BOOT)
    for (b in seq_len(B_BOOT)) {
      y  <- generate_pffc(m, k, R, fx$par[1], fx$par[2])
      fb <- fit_mle(y, R, k, start = fx$par)
      if (is.null(fb)) next
      bq[b] <- CL_quantile(fb$par, L_LIMIT)
      bm[b] <- CL_moment(fb$par, L_LIMIT)
    }
    if (sum(is.finite(bq)) >= B_BOOT / 2 && is.finite(Chq)) {
      ci <- ci_percentile(bq, LEVEL)
      cq[i, "PB"] <- as.numeric(ci[1] < Cq_true && ci[2] > Cq_true); lq[i, "PB"] <- diff(ci)
      ci <- ci_normal_boot(Chq, bq, LEVEL)
      cq[i, "NB"] <- as.numeric(ci[1] < Cq_true && ci[2] > Cq_true); lq[i, "NB"] <- diff(ci)
    }
    if (sum(is.finite(bm)) >= B_BOOT / 2 && is.finite(Chm)) {
      ci <- ci_percentile(bm, LEVEL)
      cm[i, "PB"] <- as.numeric(ci[1] < Cm_true && ci[2] > Cm_true); lm[i, "PB"] <- diff(ci)
      ci <- ci_normal_boot(Chm, bm, LEVEL)
      cm[i, "NB"] <- as.numeric(ci[1] < Cm_true && ci[2] > Cm_true); lm[i, "NB"] <- diff(ci)
    }
  }

  data.frame(m = m, k = k, alpha = alpha, scheme = sch,
             AL_ACI_q = mean(lq[, "ACI"], na.rm = TRUE), CP_ACI_q = mean(cq[, "ACI"], na.rm = TRUE),
             AL_PB_q  = mean(lq[, "PB"],  na.rm = TRUE), CP_PB_q  = mean(cq[, "PB"],  na.rm = TRUE),
             AL_NB_q  = mean(lq[, "NB"],  na.rm = TRUE), CP_NB_q  = mean(cq[, "NB"],  na.rm = TRUE),
             AL_ACI_m = mean(lm[, "ACI"], na.rm = TRUE), CP_ACI_m = mean(cm[, "ACI"], na.rm = TRUE),
             AL_PB_m  = mean(lm[, "PB"],  na.rm = TRUE), CP_PB_m  = mean(cm[, "PB"],  na.rm = TRUE),
             AL_NB_m  = mean(lm[, "NB"],  na.rm = TRUE), CP_NB_m  = mean(cm[, "NB"],  na.rm = TRUE),
             ## Number of replications contributing to each summary.  Values
             ## below N_REP mean that some replications were discarded -- the
             ## observed information was not positive definite, or alpha-hat
             ## fell below 2 so that C_L^M was undefined, or fewer than half the
             ## bootstrap refits converged.  The reported coverage is then
             ## conditional on success, so these counts must be inspected before
             ## the corresponding CP is read as an unconditional coverage.
             n_ACI_q = sum(is.finite(cq[, "ACI"])), n_boot_q = sum(is.finite(cq[, "PB"])),
             n_ACI_m = sum(is.finite(cm[, "ACI"])), n_boot_m = sum(is.finite(cm[, "PB"])),
             n_rep = N_REP,
             stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## Drive
## ---------------------------------------------------------------------------

grid <- expand.grid(scheme = SCHEMES, alpha = ALPHA_SET, k = K_SET, m = M_SET,
                    stringsAsFactors = FALSE)
## One distinct seed per configuration.  Because each configuration owns its
## seed, results do not depend on the order in which the grid is traversed and
## are identical whether or not PARALLEL is used.
grid$seed <- MASTER_SEED + 5000L + seq_len(nrow(grid))

## ---------------------------------------------------------------------------
## Drive, with a progress bar
## ---------------------------------------------------------------------------

n_cell <- nrow(grid)
cat(sprintf("%d configurations, %d replications each with B = %d.\n",
            n_cell, N_REP, B_BOOT))

if (PARALLEL) {
  need("future", "furrr")
  n_work <- max(1, future::availableCores() - 1)
  future::plan(future::multisession, workers = n_work)
  cat(sprintf("Running on %d workers.\n", n_work))
} else {
  cat("Running on a single core.\n")
}
utils::flush.console()

t_start <- Sys.time()

## One future_pmap over the whole grid rather than batches: the progress bar
## now supplies the feedback that batching was there to provide, and dropping
## the batch boundaries removes the barrier at which every worker had to wait
## for the slowest configuration in its round.
parts <- with_progress_bar(n_cell * N_TICK, function(p) {
  if (PARALLEL) {
    furrr::future_pmap(
      list(grid$m, grid$k, grid$alpha, grid$scheme, grid$seed),
      function(m, k, alpha, sch, seed) run_cell(m, k, alpha, sch, seed, p),
      .options = furrr::furrr_options(seed = TRUE))
  } else {
    lapply(seq_len(n_cell), function(i)
      run_cell(grid$m[i], grid$k[i], grid$alpha[i], grid$scheme[i],
               grid$seed[i], p))
  }
})

if (PARALLEL) future::plan(future::sequential)

out <- do.call(rbind, parts)
cat(sprintf("\nFinished in %.1f minutes.\n",
            as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

## Per-configuration summary, printed after the bar so as not to disturb it.
for (i in seq_len(nrow(out)))
  cat(sprintf("m=%3d k=%d alpha=%2d %-6s | ACI %.4f/%.4f PB %.4f/%.4f NB %.4f/%.4f\n",
              out$m[i], out$k[i], out$alpha[i], out$scheme[i],
              out$AL_ACI_q[i], out$CP_ACI_q[i], out$AL_PB_q[i], out$CP_PB_q[i],
              out$AL_NB_q[i], out$CP_NB_q[i]))

out <- out[order(out$m, out$k, out$alpha, match(out$scheme, SCHEMES)), ]
write.csv(out, file.path(RESULTS_DIR, "sim_ci.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## LaTeX fragments
## ---------------------------------------------------------------------------

write_ci_table <- function(o, suffix, path, label, caption) {
  con <- base::file(path, open = "wt"); wl <- function(...) writeLines(paste0(...), con)
  wl("\\begin{table}[H]"); wl("\\centering")
  wl("\\caption{", caption, "}"); wl("\\label{", label, "}"); wl("\\scriptsize")
  wl("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lll l rr rr rr}"); wl("\\toprule")
  wl(" & & & & \\multicolumn{2}{c}{ACI} & \\multicolumn{2}{c}{PB} & \\multicolumn{2}{c}{NB} \\\\")
  wl("\\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}")
  wl("$m$ & $k$ & $\\alpha$ & Scheme",
     paste(rep(" & \\multicolumn{1}{c}{AL} & \\multicolumn{1}{c}{CP}", 3), collapse = ""), " \\\\")
  wl("\\midrule")
  blk <- 0
  for (i in seq_len(nrow(o))) {
    first <- (i - 1) %% 4 == 0
    if (first && i > 1) { blk <- blk + 1
      wl(if (blk %% 2 == 0) "\\midrule" else "\\addlinespace") }
    wl(sprintf("%s & %s & %s & %s & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
               if (first && (i == 1 || o$m[i] != o$m[i - 1])) o$m[i] else "",
               if (first && (i == 1 || o$k[i] != o$k[i - 1] || o$m[i] != o$m[i - 1])) o$k[i] else "",
               if (first) o$alpha[i] else "", o$scheme[i],
               o[[paste0("AL_ACI_", suffix)]][i], o[[paste0("CP_ACI_", suffix)]][i],
               o[[paste0("AL_PB_",  suffix)]][i], o[[paste0("CP_PB_",  suffix)]][i],
               o[[paste0("AL_NB_",  suffix)]][i], o[[paste0("CP_NB_",  suffix)]][i]))
  }
  wl("\\bottomrule"); wl("\\end{tabular*}"); wl("\\end{table}")
  close(con)
}

write_ci_table(out, "q", file.path(TABLES_DIR, "tab_ci_quantile.tex"),
               "T:ci_quantile",
               "ALs and CPs of 95\\% CIs for the Quantile-based Index.")
write_ci_table(out, "m", file.path(TABLES_DIR, "tab_ci_moment.tex"),
               "T:ci_moment",
               "ALs and CPs of 95\\% CIs for the Moment-based Index.")

cat("\nWrote sim_ci.csv, tab_ci_quantile.tex and tab_ci_moment.tex\n")
