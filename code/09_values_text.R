################################################################################
##  09_values_text.R
##
##  Every summary number quoted in the running text of Sections 5.1, 5.2, 5.4,
##  Appendix B and Appendix C, written out as \newcommand macros.
##
##  The tables have always been generated; these are the figures the prose
##  quotes about them -- the coverage averaged over the design, how many
##  replications were discarded, how often one interval is nearer to nominal
##  than another -- and until this script existed they were typed into the
##  manuscript by hand.  Typed numbers drift: the simulations are re-run, a
##  table changes, and a sentence three sections away still describes the old
##  one.  This script removes that possibility for the last of them.
##
##  It reads only the csv files the simulations leave in results/ and so takes
##  seconds; it does not re-run anything.  If a csv is missing it says which
##  simulation produces it and stops.
##
##  Runtime: seconds.
################################################################################

source("00_gie_pffc.R")

need_csv <- function(f, produced_by) {
  p <- file.path(RESULTS_DIR, f)
  if (!file.exists(p))
    stop("09_values_text.R: ", p, " not found. Run ", produced_by, " first.")
  read.csv(p, stringsAsFactors = FALSE)
}

ci <- need_csv("sim_ci.csv",       "05_sim_ci.R")
bm <- need_csv("sim_bias_mse.csv", "04_sim_bias_mse.R")
sm <- need_csv("sim_smallm.csv",   "07_sim_smallm.R")

## ---------------------------------------------------------------------------
## Section 5.1: replications discarded in the point-estimation study
## ---------------------------------------------------------------------------
## 04_sim_bias_mse.R does not record its replication count in the csv, so it is
## recovered as the largest count observed: no cell can exceed the number run.
n_rep_pt <- max(bm$n_ok_q, bm$n_ok_m)
pt_total <- n_rep_pt * nrow(bm)
lost_q   <- pt_total - sum(bm$n_ok_q)
lost_m   <- pt_total - sum(bm$n_ok_m)
worst    <- bm[which.min(bm$n_ok_m), ]
worst_lost <- n_rep_pt - worst$n_ok_m

## ---------------------------------------------------------------------------
## Section 5.2 and Appendix B: replications discarded in the interval study
## ---------------------------------------------------------------------------
n_rep_ci <- ci$n_rep[1]
ci_total <- n_rep_ci * nrow(ci)
ci_lost_q <- ci_total - sum(ci$n_ACI_q)
ci_lost_m <- ci_total - sum(ci$n_ACI_m)
ci_min_q  <- min(ci$n_ACI_q)
ci_min_m  <- min(ci$n_ACI_m)
ci_worst_cell_q <- n_rep_ci - ci_min_q

## ---------------------------------------------------------------------------
## Section 5.2: coverage averaged over the design, by m
## ---------------------------------------------------------------------------
avg <- function(col, m) mean(ci[[col]][ci$m == m])
ms  <- sort(unique(ci$m))

## ---------------------------------------------------------------------------
## Appendix C: how the moment-based index differs
## ---------------------------------------------------------------------------
nearer_nb_m <- sum(abs(ci$CP_NB_m - 0.95) < abs(ci$CP_ACI_m - 0.95))
nb_m_by_m   <- vapply(ms, function(m) mean(ci$CP_NB_m[ci$m == m]), 0)
nb_m_cons   <- ms[nb_m_by_m > 0.95]

## ---------------------------------------------------------------------------
## Section 5.4: the range of ACI coverage in the small-sample block
## ---------------------------------------------------------------------------

vals <- list(
  ## ---- Section 5.1 --------------------------------------------------------
  ctNrepPt      = as.character(n_rep_pt),
  ctNconfig     = as.character(nrow(bm)),
  ctPtTotal     = formatC(pt_total, big.mark = "{,}", format = "d"),
  ctLostQ       = as.character(lost_q),
  ctLostQpct    = fmt(100 * lost_q / pt_total, 2),
  ctLostM       = as.character(lost_m),
  ctLostMpct    = fmt(100 * lost_m / pt_total, 2),
  ctWorstLost   = as.character(worst_lost),
  ctWorstM      = as.character(worst$m),
  ctWorstK      = as.character(worst$k),
  ctWorstScheme = tolower(worst$scheme),
  ## ---- Section 5.2 and Appendix B ----------------------------------------
  ctNrepCi      = as.character(n_rep_ci),
  ctCiTotal     = formatC(ci_total, big.mark = "{,}", format = "d"),
  ctCiLostQ     = as.character(ci_lost_q),
  ctCiLostQpct  = fmt(100 * ci_lost_q / ci_total, 2),
  ctCiLostM     = as.character(ci_lost_m),
  ctCiMinQ      = as.character(ci_min_q),
  ctCiMinM      = as.character(ci_min_m),
  ctCiWorstCell = as.character(ci_worst_cell_q),
  ## coverage averages quoted in Section 5.2, quantile-based index
  ctCPACIa      = fmt(avg("CP_ACI_q", ms[1]), 3),
  ctCPACIb      = fmt(avg("CP_ACI_q", ms[2]), 3),
  ctCPACIc      = fmt(avg("CP_ACI_q", ms[3]), 3),
  ctCPNBa       = fmt(avg("CP_NB_q",  ms[1]), 3),
  ctCPNBc       = fmt(avg("CP_NB_q",  ms[3]), 3),
  ctCPPBa       = fmt(avg("CP_PB_q",  ms[1]), 3),
  ctCPPBc       = fmt(avg("CP_PB_q",  ms[3]), 3),
  ## ---- Appendix C ---------------------------------------------------------
  ctCPACIqAll   = fmt(mean(ci$CP_ACI_q), 3),
  ctCPACImAll   = fmt(mean(ci$CP_ACI_m), 3),
  ctCPACImMin   = fmt(min(ci$CP_ACI_m),  3),
  ctNearerNB    = as.character(nearer_nb_m),
  ctNBconsM     = paste(nb_m_cons, collapse = " and "),
  ctShortestACI = as.character(sum(ci$AL_ACI_q < ci$AL_PB_q & ci$AL_ACI_q < ci$AL_NB_q)),
  ctOrderACIPBNB= as.character(sum(ci$AL_ACI_q < ci$AL_PB_q & ci$AL_PB_q < ci$AL_NB_q)),
  ## ---- Section 5.4 --------------------------------------------------------
  ctSmCPACImin  = fmt(min(sm$CP_ACI), 3),
  ctSmCPACImax  = fmt(max(sm$CP_ACI), 3)
)

write_macros(file.path(TABLES_DIR, "values_text.tex"), vals, "09_values_text.R")

cat("\nNumbers quoted in the running text, now generated:\n")
for (nm in names(vals)) cat(sprintf("  %-14s %s\n", nm, vals[[nm]]))
