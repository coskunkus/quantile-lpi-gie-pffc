################################################################################
##  run_all.R
##
##  The reproduction of the paper, from nothing.
##
##      Rscript run_all.R
##
##  or, from the RStudio console, after Session > Set Working Directory > To
##  Source File Location:
##
##      source("run_all.R")
##
##  This is the only entry point, and it is deliberate.  No script in this
##  package reads a stored result: not a LaTeX file, not a csv, not an rds.
##  Every number in the paper is therefore the output of a computation carried
##  out in this run, and never of a file being reread.  Two consequences follow.
##
##  First, the four Monte Carlo studies are always run.  They take over an hour
##  between them, and there is no shortcut: the csv files in results/ are a
##  record of what a run produced, not an input to the next one.
##
##  Second, the steps must run in one session and in this order, because three
##  of them need what earlier ones computed.  08_figures_simulation.R draws
##  Figures 4 and 5 from the two baseline studies, 09_values_text.R takes from
##  all four studies every number the running text quotes about a simulated
##  result, and 09 also needs the interval lengths of the second application.
##  Each is handed the objects the earlier steps left in the session, and each
##  stops with an explanation if it is run on its own.
##
##  The design of the first application is the one thing that crosses between
##  steps without going through the session: 07_sim_smallm.R simulates at it,
##  and the validation scripts v4, v5, v7 and v8 include it among the designs
##  they check.  None of them reads it back either.  It is rebuilt from the 23
##  published endurance times by ball_bearing_design() in 00_gie_pffc.R, which
##  takes milliseconds and gives every caller the same sample and the same fit.
##
##  ORDER
##
##      10_realdata1_ballbearings.R   Sections 6.1 and 6.3, Table 9
##      11_realdata2_guineapig.R      Section 6.2, Table 10
##      the eight validation scripts  the analytical results of the paper
##                                    against independent numerical checks
##      12_table1_translation.R       Table 1
##      01_figures.R                  Figures 1-3
##      04_sim_bias_mse.R             Tables 2-4        (2 min)
##      05_sim_ci.R                   Tables 5 and 11   (52 min, parallel)
##      06_sim_heavytail.R            Tables 6 and 7    (9 min, parallel)
##      07_sim_smallm.R               Table 8           (11 min)
##      08_figures_simulation.R       Figures 4 and 5
##      09_values_text.R              every number the running text quotes
##
##  The applications come first because the validation scripts check the design
##  of the first one, and because 07 simulates at it.
##
##  Runtime: an hour and a quarter, almost all of it in 04 to 07.  The stated
##  times are the ones measured on the machine the paper was produced on; a
##  slower machine, or one with fewer cores for 05 and 06, will take longer.
##
##  Output goes to tables/ and figures/: the LaTeX table fragments and the
##  \newcommand macro files the manuscript reports, and the figures.  The
##  raw output of the four studies is left in results/ as a record of the run.
##  Copy the fragments into the manuscript's directory and compile it there
##  with pdflatex, three times: the table floats need a third pass to settle.
################################################################################

if (!file.exists("00_gie_pffc.R"))
  stop("run_all.R must be run from the root of this package, the directory ",
       "that holds 00_gie_pffc.R; the working directory is currently '",
       getwd(), "'.")

source("00_gie_pffc.R")

applications <- c("10_realdata1_ballbearings.R",
                  "11_realdata2_guineapig.R")

checks <- file.path("validation",
                    c("v1_check_derivatives.R",
                      "v2_check_invariance.R",
                      "v3_check_moment_formula.R",
                      "v4_check_information.R",
                      "v5_check_information_and_mle.R",
                      "v6_check_bootstrap.R",
                      "v7_check_pffc_generator.R",
                      "v8_check_score_moments.R"))

analytic <- c("12_table1_translation.R",
              "01_figures.R")

studies <- c("04_sim_bias_mse.R",
             "05_sim_ci.R",
             "06_sim_heavytail.R",
             "07_sim_smallm.R")

summaries <- c("08_figures_simulation.R",
               "09_values_text.R")

run <- function(f) {
  cat("\n", strrep("=", 78), "\n== ", f, "\n", strrep("=", 78), "\n", sep = "")
  t0 <- Sys.time()
  if (startsWith(f, "validation/")) {
    ## the validation scripts source ../00_gie_pffc.R, so run them from there
    owd <- getwd(); setwd("validation")
    tryCatch(source(basename(f), echo = FALSE), finally = setwd(owd))
  } else {
    ## sourced into the global environment on purpose: the studies leave their
    ## results there for 08 and 09, which is how this package passes a result
    ## from one step to the next without writing and rereading a file
    source(f, echo = FALSE)
  }
  cat(sprintf("\n-- %s finished in %.1f minutes\n", f,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cat("run_all.R: the two applications, the eight validation scripts, all eleven\n",
    "           tables, Figures 1-5, and every number quoted in the running\n",
    "           text.  The four Monte Carlo studies are run, not read back, so\n",
    "           this takes over an hour.\n", sep = "")

t_all <- Sys.time()
for (f in c(applications, checks, analytic, studies, summaries)) run(f)

## ---------------------------------------------------------------------------
## Does the small-sample study still describe the application it mirrors?
## ---------------------------------------------------------------------------
##
## 07_sim_smallm.R takes its shape, scale, limit, group size and censoring plan
## from ball_bearing_design(), so the two cannot drift apart within a run.  The
## comparison is made anyway, and on the objects this run produced, because it
## is the one place where a change to the application would silently change
## what Section 5.4 describes.
app <- ball_bearing_design()
want <- CL_quantile(c(app$alpha_hat, app$lambda_hat), app$L)
have <- SIM_SMALLM$C_true[1]

drift <- character(0)
if (!isTRUE(all.equal(want, have, tolerance = 5e-4)))
  drift <- c(drift, sprintf("true index: study %.4f, application %.4f", have, want))
if (!isTRUE(all.equal(SIM_SMALLM$k[1], app$k)))
  drift <- c(drift, sprintf("group size k: study %s, application %s",
                            SIM_SMALLM$k[1], app$k))
if (!isTRUE(all.equal(SIM_SMALLM$R1[1], app$R[1])))
  drift <- c(drift, sprintf("censoring plan R_1: study %s, application %s",
                            SIM_SMALLM$R1[1], app$R[1]))

cat("\n", strrep("=", 78), "\n", sep = "")
if (length(drift)) {
  cat("WARNING: the small-sample study of Section 5.4 does not match the first\n")
  cat("         application it claims to mirror.\n")
  for (d in drift) cat("  - ", d, "\n", sep = "")
  cat("  Table 8 and the numbers Section 5.4 quotes are therefore not\n")
  cat("  descriptions of the experiment reported in Section 6.1.\n")
} else {
  cat("Section 5.4 matches the first application: C_L^xi = ",
      sprintf("%.4f", have), ", k = ", app$k, ", R_1 = ", app$R[1], ".\n",
      sep = "")
}

cat(sprintf("\nDone in %.1f minutes.  All eleven table fragments, all macro files\n",
            as.numeric(difftime(Sys.time(), t_all, units = "mins"))))
cat("and all five figures were written from this run.\n")
cat("They are in ", TABLES_DIR, " and the figures in ", FIGURES_DIR, ".\n",
    "The raw output of the four studies is in ", RESULTS_DIR,
    " as a record of the run.\n\n", sep = "")
print(sessionInfo())
