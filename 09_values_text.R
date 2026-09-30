################################################################################
##  09_values_text.R
##
##  Every summary number quoted in the running text of Sections 5.1 to 5.4,
##  Appendix B and Appendix C, written out as \newcommand macros.
##
##  The tables have always been generated; these are the figures the prose
##  quotes about them -- the coverage averaged over the design, how many
##  replications were discarded, in how many configurations an ordering holds
##  without exception -- and until this script existed they were typed by hand.
##  Typed numbers drift: the simulations are re-run, a table changes, and a
##  sentence three sections away still describes the old one.  Everything the
##  running text asserts about a simulated result is computed here, so that any
##  such sentence can be checked against the raw output it describes.
##
##  It reads only the csv files the simulations leave in results/ and so takes
##  seconds; it does not re-run anything.  If a csv is missing it says which
##  simulation produces it and stops.
##
##  Runtime: seconds.
################################################################################

source("00_gie_pffc.R")

## The four studies are taken from the session, not from results/: nothing in
## this package reads a stored result back.  This step is therefore not a
## script to run on its own; run_all.R runs the studies and this step in one
## session, in order.
from_step <- function(obj, produced_by) {
  if (!exists(obj, envir = globalenv()))
    stop("09_values_text.R: ", obj, " is not in this session.\n",
         "  ", produced_by, " has to have run in it first. Use run_all.R.")
  get(obj, envir = globalenv())
}

ci <- from_step("SIM_CI",        "05_sim_ci.R")
bm <- from_step("SIM_BIAS_MSE",  "04_sim_bias_mse.R")
sm <- from_step("SIM_SMALLM",    "07_sim_smallm.R")
ht <- from_step("SIM_HEAVYTAIL", "06_sim_heavytail.R")

## ---------------------------------------------------------------------------
## Counting the design-wide orderings the prose asserts
## ---------------------------------------------------------------------------
##
## Sections 5.1 to 5.3 each claim that an ordering holds in every one of N
## configurations.  Three groupings are involved and each appears more than
## once, so they are written as functions rather than repeated: the censoring
## plan is compared within (m, k, alpha); the group size within
## (m, alpha, plan); and the sample size within (k, alpha, plan).  Each returns
## the number of groups in which the claim holds together with the number of
## groups, so a sentence quoting "all N" can be checked against both.
##
## `cols` is a named logical vector: the column to look at, and whether the
## comparison is on the absolute value, as it is for a bias.

group_keys <- function(d, by) unique(d[, by, drop = FALSE])

group_rows <- function(d, key, i, by) {
  sel <- rep(TRUE, nrow(d))
  for (v in by) sel <- sel & d[[v]] == key[i, v]
  d[sel, , drop = FALSE]
}

## Early censoring smallest and late censoring largest, within (m, k, alpha).
plan_ordering <- function(d, cols) {
  by  <- c("m", "k", "alpha"); key <- group_keys(d, by); n <- 0L
  for (i in seq_len(nrow(key))) {
    g <- group_rows(d, key, i, by); good <- TRUE
    for (cc in names(cols)) {
      v <- if (cols[[cc]]) abs(g[[cc]]) else g[[cc]]
      good <- good && g$scheme[which.min(v)] == "Early" &&
                      g$scheme[which.max(v)] == "Late"
    }
    if (good) n <- n + 1L
  }
  c(ok = n, total = nrow(key))
}

## Every column strictly larger at the larger k, within (m, alpha, plan).
k_ordering <- function(d, cols) {
  by  <- c("m", "alpha", "scheme"); key <- group_keys(d, by); n <- 0L
  for (i in seq_len(nrow(key))) {
    g <- group_rows(d, key, i, by); g <- g[order(g$k), ]; good <- TRUE
    for (cc in names(cols)) {
      v <- if (cols[[cc]]) abs(g[[cc]]) else g[[cc]]
      good <- good && all(diff(v) > 0)
    }
    if (good) n <- n + 1L
  }
  c(ok = n, total = nrow(key))
}

## Every column strictly falling with m, within (k, alpha, plan).
m_ordering <- function(d, cols) {
  by  <- c("k", "alpha", "scheme"); key <- group_keys(d, by); n <- 0L
  for (i in seq_len(nrow(key))) {
    g <- group_rows(d, key, i, by); g <- g[order(g$m), ]; good <- TRUE
    for (cc in names(cols)) {
      v <- if (cols[[cc]]) abs(g[[cc]]) else g[[cc]]
      good <- good && all(diff(v) < 0)
    }
    if (good) n <- n + 1L
  }
  c(ok = n, total = nrow(key))
}

## ---------------------------------------------------------------------------
## Section 5.1: the baseline point-estimation study
## ---------------------------------------------------------------------------

BOTH <- c(Bias_q = TRUE, MSE_q = FALSE, Bias_m = TRUE, MSE_m = FALSE)
pt_plan <- plan_ordering(bm, BOTH)
pt_k    <- k_ordering(bm, BOTH)
pt_m    <- m_ordering(bm, BOTH)
bias_neg <- all(bm$Bias_q < 0) && all(bm$Bias_m < 0)

## Replications discarded.  04_sim_bias_mse.R does not record its replication
## count in the csv, so it is recovered as the largest count observed: no cell
## can exceed the number run.  A replication is lost to the quantile-based
## index only when the fit fails to converge; the moment-based index loses
## those as well, and in addition every replication with alpha-hat <= 2.
n_rep_pt   <- max(bm$n_ok_q, bm$n_ok_m)
pt_total   <- n_rep_pt * nrow(bm)
lost_q     <- pt_total - sum(bm$n_ok_q)
lost_m     <- pt_total - sum(bm$n_ok_m)
lost_alpha <- lost_m - lost_q
worst      <- bm[which.min(bm$n_ok_m), ]
worst_lost <- n_rep_pt - worst$n_ok_m

## The scale-free comparison, and the factorisation of its ratio.  The second
## factor is a property of the two targets alone and so is fixed within a
## shape; the first compares the two estimators.  Section 5.1 quotes both, by
## shape, because the difference between the shapes lies in the targets.
alphas  <- sort(unique(bm$alpha))
fac_two <- vapply(alphas, function(a) {
  r <- bm[bm$alpha == a, ][1, ]; (r$Cq_true / r$Cm_true)^2
}, 0)
fac_one <- vapply(alphas, function(a)
  mean(bm$MSE_m[bm$alpha == a] / bm$MSE_q[bm$alpha == a]), 0)
rat_by_alpha <- vapply(alphas, function(a) mean(bm$ratio[bm$alpha == a]), 0)

## Setting the two estimators against one common target through Eq. (13): the
## quantity kappa*C^M_hat - delta estimates C^xi and has mean squared error
## kappa^2 * MSE(C^M_hat), so no normalisation is involved at all.
kappa_alpha <- vapply(alphas, function(a)
  unname(translation_constants(a)["kappa"]), 0)
common_ok <- sum(vapply(seq_along(alphas), function(j) {
  s <- bm$alpha == alphas[j]
  sum(kappa_alpha[j]^2 * bm$MSE_m[s] > bm$MSE_q[s])
}, 0L))

## ---------------------------------------------------------------------------
## Section 5.2 and Appendix B: the interval study on the same design
## ---------------------------------------------------------------------------

n_rep_ci  <- ci$n_rep[1]
ci_total  <- n_rep_ci * nrow(ci)
ci_lost_q <- ci_total - sum(ci$n_ACI_q)
ci_lost_m <- ci_total - sum(ci$n_ACI_m)
ci_min_q  <- min(ci$n_ACI_q)
ci_min_m  <- min(ci$n_ACI_m)

ms  <- sort(unique(ci$m))
avg <- function(col, m) mean(ci[[col]][ci$m == m])

## Lengths: the ordering of the three procedures, and their response to the
## design.  The three are compared per configuration; the plan, group size and
## sample size are compared as in Section 5.1, on all three lengths at once.
ALLAL_q <- c(AL_ACI_q = FALSE, AL_PB_q = FALSE, AL_NB_q = FALSE)
ALLAL_m <- c(AL_ACI_m = FALSE, AL_PB_m = FALSE, AL_NB_m = FALSE)
order_q <- sum(ci$AL_ACI_q < ci$AL_PB_q & ci$AL_PB_q < ci$AL_NB_q)
order_m <- sum(ci$AL_ACI_m < ci$AL_PB_m & ci$AL_PB_m < ci$AL_NB_m)
ci_m_q  <- m_ordering(ci, ALLAL_q)
ci_plan_q <- plan_ordering(ci, ALLAL_q)
ci_k_q    <- k_ordering(ci, ALLAL_q)
ci_plan_m <- plan_ordering(ci, ALLAL_m)

## The NB interval is centred correctly and still over-covers, because its
## width comes from the bootstrap standard error rather than the delta-method
## one.  Both intervals are z times their standard error wide, so the ratio of
## the lengths is the ratio of the standard errors.
se_excess <- vapply(ms, function(m)
  100 * (mean(ci$AL_NB_q[ci$m == m] / ci$AL_ACI_q[ci$m == m]) - 1), 0)

## ---------------------------------------------------------------------------
## Appendix C: how the moment-based index differs
## ---------------------------------------------------------------------------
nearer_nb_m <- sum(abs(ci$CP_NB_m - 0.95) < abs(ci$CP_ACI_m - 0.95))
nb_m_by_m   <- vapply(ms, function(m) mean(ci$CP_NB_m[ci$m == m]), 0)
nb_m_cons   <- ms[nb_m_by_m > 0.95]

## ---------------------------------------------------------------------------
## Section 5.3: the heavy-tailed designs
## ---------------------------------------------------------------------------

ht_plan_mse  <- plan_ordering(ht, c(MSE = FALSE))
ht_plan_bias <- plan_ordering(ht, c(Bias = TRUE))
ht_k_mse     <- k_ordering(ht, c(MSE = FALSE))
ht_k_bias    <- k_ordering(ht, c(Bias = TRUE))

ht_ms   <- sort(unique(ht$m))
ht_cpm  <- function(m, col) mean(ht[[col]][ht$m == m])
ht_ord  <- sum(ht$AL_ACI < ht$AL_PB & ht$AL_PB < ht$AL_NB)
HTAL    <- c(AL_ACI = FALSE, AL_PB = FALSE, AL_NB = FALSE)
ht_m_al    <- m_ordering(ht, HTAL)
ht_plan_al <- plan_ordering(ht, HTAL)
ht_k_al    <- k_ordering(ht, HTAL)

## How far apart the four plans are, averaged over the design points at each m.
ht_spread <- vapply(ht_ms, function(mm) {
  key <- unique(ht[ht$m == mm, c("k", "alpha")])
  mean(vapply(seq_len(nrow(key)), function(i) {
    g <- ht[ht$m == mm & ht$k == key$k[i] & ht$alpha == key$alpha[i], ]
    max(g$Bias) - min(g$Bias)
  }, 0))
}, 0)

## ---------------------------------------------------------------------------
## Section 6.2: how far apart the three interval lengths are in Table 10
## ---------------------------------------------------------------------------
gp_len <- from_step("APP_GUINEAPIG", "11_realdata2_guineapig.R")$AL
if (anyNA(gp_len) || length(gp_len) != 3L)
  stop("09_values_text.R: the second application did not leave its three ",
       "interval lengths in the session.")

## ---------------------------------------------------------------------------
## Where the prose brackets a set of values ("between x and y"), the endpoints
## are rounded outwards rather than to nearest, so that the bracket printed is
## always true of every value in the set.
## ---------------------------------------------------------------------------
fmt_lo <- function(x, d) fmt(floor(x * 10^d + 1e-9) / 10^d, d)
fmt_hi <- function(x, d) fmt(ceiling(x * 10^d - 1e-9) / 10^d, d)

vals <- list(
  ## ---- Section 5.1 --------------------------------------------------------
  ctNrepPt      = as.character(n_rep_pt),
  ctNconfig     = as.character(nrow(bm)),
  ctPtTotal     = formatC(pt_total, big.mark = "{,}", format = "d"),
  ctLostQ       = as.character(lost_q),
  ctLostAlpha   = as.character(lost_alpha),
  ctLostM       = as.character(lost_m),
  ctWorstLost   = as.character(worst_lost),
  ctWorstM      = as.character(worst$m),
  ctWorstK      = as.character(worst$k),
  ctWorstAlpha  = as.character(worst$alpha),
  ctWorstScheme = tolower(worst$scheme),
  ## the three design-wide orderings, each holding for both indexes
  ctPlanOK      = as.character(pt_plan["ok"]),
  ctPlanCells   = as.character(pt_plan["total"]),
  ctKOK         = as.character(pt_k["ok"]),
  ctKCells      = as.character(pt_k["total"]),
  ctMOK         = as.character(pt_m["ok"]),
  ctMCells      = as.character(pt_m["total"]),
  ctBiasNeg     = if (bias_neg) "all" else "not all",
  ## the scale-free comparison, by shape
  ctAlphaLo     = as.character(alphas[1]),
  ctAlphaHi     = as.character(alphas[2]),
  ctFacTwoLo    = fmt(fac_two[1], 2),
  ctFacTwoHi    = fmt(fac_two[2], 2),
  ctFacOneLo    = fmt(fac_one[1], 2),
  ctFacOneHi    = fmt(fac_one[2], 2),
  ctRatioLo     = fmt(rat_by_alpha[1], 1),
  ctRatioHi     = fmt(rat_by_alpha[2], 1),
  ctCommonOK    = as.character(common_ok),
  ## ---- Section 5.2 and Appendix B ----------------------------------------
  ctNrepCi      = as.character(n_rep_ci),
  ctCiTotal     = formatC(ci_total, big.mark = "{,}", format = "d"),
  ctCiLostQ     = as.character(ci_lost_q),
  ctCiLostQpct  = fmt(100 * ci_lost_q / ci_total, 2),
  ctCiLostM     = as.character(ci_lost_m),
  ctCiWorstCell = as.character(n_rep_ci - ci_min_q),
  ctCiWorstCellM= as.character(n_rep_ci - ci_min_m),
  ## lengths
  ctOrderACIPBNB= as.character(order_q),
  ctOrderMoment = as.character(order_m),
  ctCiMonoM     = as.character(ci_m_q["ok"]),
  ctCiMonoCells = as.character(ci_m_q["total"]),
  ctCiPlanOK    = as.character(ci_plan_q["ok"]),
  ctCiPlanCells = as.character(ci_plan_q["total"]),
  ctCiKOK       = as.character(ci_k_q["ok"]),
  ctCiKCells    = as.character(ci_k_q["total"]),
  ctCiPlanOKm   = as.character(ci_plan_m["ok"]),
  ## coverage averages quoted in Section 5.2, quantile-based index
  ctCPACIa      = fmt(avg("CP_ACI_q", ms[1]), 3),
  ctCPACIb      = fmt(avg("CP_ACI_q", ms[2]), 3),
  ctCPACIc      = fmt(avg("CP_ACI_q", ms[3]), 3),
  ctCPACIqMin   = fmt_lo(min(ci$CP_ACI_q), 3),
  ctCPACIqMax   = fmt_hi(max(ci$CP_ACI_q), 3),
  ctCPNBa       = fmt(avg("CP_NB_q",  ms[1]), 3),
  ctCPNBc       = fmt(avg("CP_NB_q",  ms[3]), 3),
  ctCPPBa       = fmt(avg("CP_PB_q",  ms[1]), 3),
  ctCPPBc       = fmt(avg("CP_PB_q",  ms[3]), 3),
  ctSEexcA      = fmt(se_excess[1], 0),
  ctSEexcB      = fmt(se_excess[2], 0),
  ctSEexcC      = fmt(se_excess[3], 0),
  ## ---- Appendix C ---------------------------------------------------------
  ctCPACIqAll   = fmt(mean(ci$CP_ACI_q), 3),
  ctCPACImAll   = fmt(mean(ci$CP_ACI_m), 3),
  ctCPACImMin   = fmt_lo(min(ci$CP_ACI_m), 3),
  ctNearerNB    = as.character(nearer_nb_m),
  ctNBconsM     = paste(nb_m_cons, collapse = " and "),
  ## ---- Section 5.3 --------------------------------------------------------
  ctHtCells     = as.character(ht_plan_mse["total"]),
  ctHtMseOK     = as.character(ht_plan_mse["ok"]),
  ctHtBiasOK    = as.character(ht_plan_bias["ok"]),
  ctHtKmse      = as.character(ht_k_mse["ok"]),
  ctHtKbias     = as.character(ht_k_bias["ok"]),
  ctHtKCells    = as.character(ht_k_mse["total"]),
  ctHtSpreadA   = fmt(ht_spread[1], 3),
  ctHtSpreadC   = fmt(ht_spread[3], 3),
  ctHtOrder     = as.character(ht_ord),
  ctHtMonoM     = as.character(ht_m_al["ok"]),
  ctHtPlanAL    = as.character(ht_plan_al["ok"]),
  ctHtKAL       = as.character(ht_k_al["ok"]),
  ctHtCPACIa    = fmt(ht_cpm(ht_ms[1], "CP_ACI"), 3),
  ctHtCPACIb    = fmt(ht_cpm(ht_ms[2], "CP_ACI"), 3),
  ctHtCPACIc    = fmt(ht_cpm(ht_ms[3], "CP_ACI"), 3),
  ctHtCPACImin  = fmt_lo(min(ht$CP_ACI), 3),
  ctHtCPACImax  = fmt_hi(max(ht$CP_ACI), 3),
  ctHtCPNBa     = fmt(ht_cpm(ht_ms[1], "CP_NB"), 3),
  ctHtCPNBc     = fmt(ht_cpm(ht_ms[3], "CP_NB"), 3),
  ctHtCPPBa     = fmt(ht_cpm(ht_ms[1], "CP_PB"), 3),
  ctHtCPPBc     = fmt(ht_cpm(ht_ms[3], "CP_PB"), 3),
  ## ---- Section 5.4 --------------------------------------------------------
  ctSmCPACImin  = fmt_lo(min(sm$CP_ACI), 4),
  ctSmCPACImax  = fmt_hi(max(sm$CP_ACI), 4),
  ctSmCPACInine = fmt(sm$CP_ACI[sm$m == min(sm$m)], 4),
  ctSmOneACImin = fmt_lo(min(sm$CP1_ACI), 3),
  ctSmOneACImax = fmt_hi(max(sm$CP1_ACI), 3),
  ctSmOneACInine= fmt(sm$CP1_ACI[sm$m == min(sm$m)], 3),
  ctSmOnePBnine = fmt(sm$CP1_PB[sm$m == min(sm$m)], 3),
  ctSmOneNBnine = fmt(sm$CP1_NB[sm$m == min(sm$m)], 3),
  ctSmErrACI    = fmt(100 * (1 - sm$CP1_ACI[sm$m == min(sm$m)]), 1),
  ctSmErrPB     = fmt(100 * (1 - sm$CP1_PB[sm$m == min(sm$m)]), 1),
  ctSmErrNB     = fmt(100 * (1 - sm$CP1_NB[sm$m == min(sm$m)]), 1),
  ## ---- Section 6.2 --------------------------------------------------------
  gpSpread      = fmt(100 * (max(gp_len) - min(gp_len)) / min(gp_len), 1)
)

write_macros(file.path(TABLES_DIR, "values_text.tex"), vals, "09_values_text.R")

cat("\nNumbers quoted in the running text, now generated:\n")
for (nm in names(vals)) cat(sprintf("  %-14s %s\n", nm, vals[[nm]]))

## ---------------------------------------------------------------------------
## The orderings the prose asserts, printed as counts out of the number of
## groups, so that a sentence saying "in all N" can be read off directly.
## ---------------------------------------------------------------------------
cat("\nDesign-wide orderings (count / groups):\n")
cat(sprintf("  5.1 plan, both indexes            %d / %d\n", pt_plan[1], pt_plan[2]))
cat(sprintf("  5.1 group size, both indexes      %d / %d\n", pt_k[1],    pt_k[2]))
cat(sprintf("  5.1 sample size, both indexes     %d / %d\n", pt_m[1],    pt_m[2]))
cat(sprintf("  5.1 common target, Eq. (13)       %d / %d\n", common_ok,  nrow(bm)))
cat(sprintf("  5.2 ACI < PB < NB, quantile       %d / %d\n", order_q,    nrow(ci)))
cat(sprintf("  5.2 ACI < PB < NB, moment         %d / %d\n", order_m,    nrow(ci)))
cat(sprintf("  5.2 lengths fall with m           %d / %d\n", ci_m_q[1],  ci_m_q[2]))
cat(sprintf("  5.2 lengths by plan               %d / %d\n", ci_plan_q[1], ci_plan_q[2]))
cat(sprintf("  5.2 lengths by group size         %d / %d\n", ci_k_q[1],  ci_k_q[2]))
cat(sprintf("  5.3 plan, MSE                     %d / %d\n", ht_plan_mse[1],  ht_plan_mse[2]))
cat(sprintf("  5.3 plan, bias                    %d / %d\n", ht_plan_bias[1], ht_plan_bias[2]))
cat(sprintf("  5.3 group size, MSE               %d / %d\n", ht_k_mse[1],  ht_k_mse[2]))
cat(sprintf("  5.3 group size, bias              %d / %d\n", ht_k_bias[1], ht_k_bias[2]))
cat(sprintf("  5.3 ACI < PB < NB                 %d / %d\n", ht_ord,     nrow(ht)))
cat(sprintf("  C   lengths by plan, moment       %d / %d\n", ci_plan_m[1], ci_plan_m[2]))
cat(sprintf("\n  every bias negative in Section 5.1: %s\n", bias_neg))
