################################################################################
##  06_sim_heavytail.R
##
##  Heavy-tailed simulation study, Section 5.3.
##  Produces Tables 6 (Bias/MSE) and 7 (AL/CP) and writes them as LaTeX
##  fragments into ../tables/, from where the manuscript \input's them.
##
##  Design:  alpha in {1,2} -- the moment-based index does not exist here --
##           lambda = 2 and L = 0.2, so that L/lambda = 0.1, close to the ratio
##           of 0.154 and 0.125 seen in the applications.  Everything else matches
##           the baseline design of Section 5.1.
##
##  Runtime: the interval part is the expensive one, roughly 2-5 hours on a
##  single core.  Reduce N_REP_CI or B for a quick check.
################################################################################

source("00_gie_pffc.R")
need("numDeriv", "statmod", "progressr")

## No global seed is needed: run_cell() seeds itself per configuration, so
## results do not depend on the order of traversal.

N_REP_PT  <- 2000      # replications for Bias and MSE
N_REP_CI  <- 500       # replications for AL and CP
B_BOOT    <- 250       # bootstrap resamples
LAMBDA    <- 2.0
L_RATIO   <- 0.1
L_LIMIT   <- LAMBDA * L_RATIO
M_SET     <- c(25, 50, 100)
K_SET     <- c(2, 5)
ALPHA_SET <- c(1, 2)
SCHEMES   <- c("Early", "Middle", "Late", "Equal")
LEVEL     <- 0.95
PARALLEL  <- TRUE       # set FALSE to run on a single core
N_TICK_PT <- 5          # progress ticks from the point-estimation loop
N_TICK_CI <- 20         # progress ticks from the interval loop
TICK_PT   <- max(1L, N_REP_PT %/% N_TICK_PT)
TICK_CI   <- max(1L, N_REP_CI %/% N_TICK_CI)
N_TICK_PT <- N_REP_PT %/% TICK_PT   # counts actually delivered
N_TICK_CI <- N_REP_CI %/% TICK_CI

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
  R      <- make_scheme(m, sch)
  true   <- c(alpha, LAMBDA)
  C_true <- CL_quantile(true, L_LIMIT)

  ## ---- point estimation ---------------------------------------------------
  Ch <- rep(NA_real_, N_REP_PT)
  for (i in seq_len(N_REP_PT)) {
    tick_if(p, i, TICK_PT)
    x  <- generate_pffc(m, k, R, alpha, LAMBDA)
    ft <- fit_mle(x, R, k, start = true)
    if (is.null(ft)) next
    Ch[i] <- CL_quantile(ft$par, L_LIMIT)
  }

  ## ---- intervals ----------------------------------------------------------
  cov <- len <- matrix(NA_real_, N_REP_CI, 3,
                       dimnames = list(NULL, c("ACI", "PB", "NB")))
  for (i in seq_len(N_REP_CI)) {
    ## Ticked before any `next`, so the tick count per configuration is exact.
    tick_if(p, i, TICK_CI)
    x  <- generate_pffc(m, k, R, alpha, LAMBDA)
    ft <- fit_mle(x, R, k, start = true)
    if (is.null(ft)) next
    Chat <- CL_quantile(ft$par, L_LIMIT)
    if (!is.finite(Chat)) next

    se <- delta_se(ft, L_LIMIT, "quantile")
    if (is.finite(se)) {
      ci <- ci_asymptotic(Chat, se, LEVEL)
      cov[i, "ACI"] <- as.numeric(ci[1] < C_true && ci[2] > C_true)
      len[i, "ACI"] <- diff(ci)
    }

    bt <- bootstrap_index(ft$par, m, k, R, L_LIMIT, B_BOOT, index = "quantile")
    if (sum(is.finite(bt)) >= B_BOOT / 2) {
      ci <- ci_percentile(bt, LEVEL)
      cov[i, "PB"] <- as.numeric(ci[1] < C_true && ci[2] > C_true)
      len[i, "PB"] <- diff(ci)
      ci <- ci_normal_boot(Chat, bt, LEVEL)
      cov[i, "NB"] <- as.numeric(ci[1] < C_true && ci[2] > C_true)
      len[i, "NB"] <- diff(ci)
    }
  }

  data.frame(
    m = m, k = k, alpha = alpha, scheme = sch, C_true = C_true,
    Bias = mean(Ch, na.rm = TRUE) - C_true,
    MSE  = mean((Ch - C_true)^2, na.rm = TRUE),
    AL_ACI = mean(len[, "ACI"], na.rm = TRUE), CP_ACI = mean(cov[, "ACI"], na.rm = TRUE),
    AL_PB  = mean(len[, "PB"],  na.rm = TRUE), CP_PB  = mean(cov[, "PB"],  na.rm = TRUE),
    AL_NB  = mean(len[, "NB"],  na.rm = TRUE), CP_NB  = mean(cov[, "NB"],  na.rm = TRUE),
    n_pt = sum(is.finite(Ch)),
    n_ACI = sum(is.finite(cov[, "ACI"])), n_boot = sum(is.finite(cov[, "PB"])),
    stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## Drive
## ---------------------------------------------------------------------------

grid <- expand.grid(scheme = SCHEMES, alpha = ALPHA_SET, k = K_SET, m = M_SET,
                    stringsAsFactors = FALSE)
## One distinct seed per configuration, so the results do not depend on the
## order of traversal and are the same serially or in parallel.
grid$seed <- MASTER_SEED + 6000L + seq_len(nrow(grid))
n_cell <- nrow(grid)

cat(sprintf("%d configurations, %d replications for Bias/MSE and %d with B = %d for the intervals.\n",
            n_cell, N_REP_PT, N_REP_CI, B_BOOT))

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

parts <- with_progress_bar(n_cell * (N_TICK_PT + N_TICK_CI), function(p) {
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

for (i in seq_len(nrow(out)))
  cat(sprintf("alpha=%d m=%3d k=%d %-6s | bias %+.4f mse %.4f | ACI %.4f/%.4f PB %.4f/%.4f NB %.4f/%.4f\n",
              out$alpha[i], out$m[i], out$k[i], out$scheme[i], out$Bias[i], out$MSE[i],
              out$AL_ACI[i], out$CP_ACI[i], out$AL_PB[i], out$CP_PB[i],
              out$AL_NB[i], out$CP_NB[i]))

out <- out[order(out$m, out$k, out$alpha, match(out$scheme, SCHEMES)), ]
write.csv(out, file.path(RESULTS_DIR, "sim_heavytail.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## LaTeX fragments
## ---------------------------------------------------------------------------

C1 <- out$C_true[out$alpha == 1][1]
C2 <- out$C_true[out$alpha == 2][1]
pick <- function(a, m, k, s, col)
  out[[col]][out$alpha == a & out$m == m & out$k == k & out$scheme == s]

## ---- Table 6: Bias and MSE -------------------------------------------------
con <- base::file(file.path(TABLES_DIR, "tab_heavytail_point.tex"), open = "wt")
wl <- function(...) writeLines(paste0(...), con)
wl("\\begin{table}[H]"); wl("\\centering")
wl(sprintf(paste0("\\caption{Bias and MSE of the MLE of the quantile-based index ",
                  "$C_L^{\\xi}$ in the heavy-tailed designs, $\\lambda=%g$ and $L=%g$ ",
                  "(so $L/\\lambda=%g$). The moment-based index does not exist for ",
                  "these shape values. True index values: $C_L^{\\xi}=%.4f$ for ",
                  "$\\alpha=1$ and $C_L^{\\xi}=%.4f$ for $\\alpha=2$. Based on %d ",
                  "Monte Carlo replications.}"),
           LAMBDA, L_LIMIT, L_RATIO, C1, C2, N_REP_PT))
wl("\\label{T:ht_point}"); wl("\\scriptsize")
wl("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lll rr rr rr rr}"); wl("\\toprule")
wl(" & & & \\multicolumn{2}{c}{Early} & \\multicolumn{2}{c}{Middle} & ",
   "\\multicolumn{2}{c}{Late} & \\multicolumn{2}{c}{Equal} \\\\")
wl("\\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9} \\cmidrule(lr){10-11}")
wl("$m$ & $k$ & $\\alpha$", paste(rep(" & \\multicolumn{1}{c}{Bias} & \\multicolumn{1}{c}{MSE}", 4), collapse = ""), " \\\\")
wl("\\midrule")
for (mi in seq_along(M_SET)) {
  m <- M_SET[mi]
  if (mi > 1) wl("\\midrule")
  for (ki in seq_along(K_SET)) {
    k <- K_SET[ki]
    if (ki > 1) wl("\\addlinespace")
    for (ai in seq_along(ALPHA_SET)) {
      a <- ALPHA_SET[ai]
      cells <- unlist(lapply(SCHEMES, function(s)
        sprintf("%.4f & %.4f", pick(a, m, k, s, "Bias"), pick(a, m, k, s, "MSE"))))
      wl(sprintf("%s & %s & %d & %s \\\\",
                 if (ki == 1 && ai == 1) m else "",
                 if (ai == 1) k else "", a, paste(cells, collapse = " & ")))
    }
  }
}
wl("\\bottomrule"); wl("\\end{tabular*}"); wl("\\end{table}")
close(con)

## ---- Table 7: AL and CP ----------------------------------------------------
con <- base::file(file.path(TABLES_DIR, "tab_heavytail_ci.tex"), open = "wt")
wl <- function(...) writeLines(paste0(...), con)
wl("\\begin{table}[H]"); wl("\\centering")
wl(sprintf(paste0("\\caption{Average lengths (AL) and coverage probabilities (CP) of ",
                  "the 95\\%% confidence intervals for $C_L^{\\xi}$ in the heavy-tailed ",
                  "designs, $\\lambda=%g$ and $L=%g$. Based on %d Monte Carlo ",
                  "replications with $B=%d$ bootstrap resamples.}"),
           LAMBDA, L_LIMIT, N_REP_CI, B_BOOT))
wl("\\label{T:ht_ci}"); wl("\\scriptsize")
wl("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lll l rr rr rr}"); wl("\\toprule")
wl(" & & & & \\multicolumn{2}{c}{ACI} & \\multicolumn{2}{c}{PB} & \\multicolumn{2}{c}{NB} \\\\")
wl("\\cmidrule(lr){5-6} \\cmidrule(lr){7-8} \\cmidrule(lr){9-10}")
wl("$m$ & $k$ & $\\alpha$ & Scheme", paste(rep(" & \\multicolumn{1}{c}{AL} & \\multicolumn{1}{c}{CP}", 3), collapse = ""), " \\\\")
wl("\\midrule")
for (mi in seq_along(M_SET)) {
  m <- M_SET[mi]; nblk <- 0
  if (mi > 1) wl("\\midrule")
  for (k in K_SET) for (a in ALPHA_SET) {
    if (nblk > 0) wl("\\addlinespace")
    for (si in seq_along(SCHEMES)) {
      s <- SCHEMES[si]
      wl(sprintf("%s & %s & %s & %s & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\",
                 if (nblk == 0 && si == 1) m else "",
                 if (a == ALPHA_SET[1] && si == 1) k else "",
                 if (si == 1) a else "", s,
                 pick(a, m, k, s, "AL_ACI"), pick(a, m, k, s, "CP_ACI"),
                 pick(a, m, k, s, "AL_PB"),  pick(a, m, k, s, "CP_PB"),
                 pick(a, m, k, s, "AL_NB"),  pick(a, m, k, s, "CP_NB")))
    }
    nblk <- nblk + 1
  }
}
wl("\\bottomrule"); wl("\\end{tabular*}"); wl("\\end{table}")
close(con)

## ---------------------------------------------------------------------------
## Numbers quoted in the running text of Section 5.3
## ---------------------------------------------------------------------------

cpm <- function(m, col) mean(out[[col]][out$m == m])
shortest <- sum(out$AL_ACI <= pmin(out$AL_PB, out$AL_NB))

write_macros(file.path(TABLES_DIR, "values_heavytail.tex"), list(
  htNrepPt      = as.character(N_REP_PT),
  htNrepCi      = as.character(N_REP_CI),
  htBootB       = as.character(B_BOOT),
  htLambda      = as.character(LAMBDA),
  htL           = as.character(L_LIMIT),
  htRatio       = as.character(L_RATIO),
  htCone        = fmt(C1, 4),
  htCtwo        = fmt(C2, 4),
  htCratioOne   = fmt(CL_quantile(c(2, LAMBDA), LAMBDA), 3),   # alpha = 2 at L/lambda = 1
  htNconfig     = as.character(nrow(out)),
  htNshortest   = as.character(shortest),
  htCPACIa      = fmt(cpm(25, "CP_ACI"), 3),
  htCPACIb      = fmt(cpm(50, "CP_ACI"), 3),
  htCPACIc      = fmt(cpm(100, "CP_ACI"), 3),
  htCPACImin    = fmt(min(out$CP_ACI), 3),
  htCPACImax    = fmt(max(out$CP_ACI), 3),
  htCPNBsmall   = fmt(cpm(25, "CP_NB"), 3),
  htCPPBsmall   = fmt(cpm(25, "CP_PB"), 3),
  htCPPBlarge   = fmt(cpm(100, "CP_PB"), 3)
), "06_sim_heavytail.R")

cat("\nWrote tab_heavytail_point.tex and tab_heavytail_ci.tex to", TABLES_DIR, "\n")
