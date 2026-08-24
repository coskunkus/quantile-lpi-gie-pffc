################################################################################
##  run_simulations.R
##
##  The four Monte Carlo studies, and nothing else.
##
##      cd code
##      Rscript run_simulations.R
##
##  or, from the RStudio console, after Session > Set Working Directory > To
##  Source File Location:
##
##      source("run_simulations.R")
##
##  What this does, in dependency order:
##
##      04_sim_bias_mse.R    Section 5.1, Tables 2, 3, 4      20-40 min
##      05_sim_ci.R          Section 5.2, Tables 5 and 11     1-3 h, parallel
##      06_sim_heavytail.R   Section 5.3, Tables 6 and 7      1-2 h, parallel
##      07_sim_smallm.R      Section 5.4, Table 8             30-60 min
##
##  Total: several hours.  05 and 06 run in parallel by default; set
##  PARALLEL <- FALSE at the top of either to use a single core.  All four
##  report progress with an estimate of the time remaining.
##
##  Each writes its raw output to results/ as a csv and then builds its own
##  LaTeX tables from it.  Re-running overwrites both.
##
##  ORDER MATTERS FOR 07.  The small-sample block of Section 5.4 mirrors the
##  first application exactly, and reads its shape, scale, limit and group size
##  from ../tables/values_realdata1.tex, which 10_realdata1_ballbearings.R
##  writes.  That file ships with the package, so this script runs as it
##  stands; but if the application has been changed, run
##
##      Rscript 10_realdata1_ballbearings.R
##
##  before this, or Section 5.4 will mirror the old design.  run_paper.R does
##  that as a matter of course and checks the two agree afterwards.
##
##  When this finishes, run run_paper.R: Figures 4 and 5 and the summary
##  numbers Sections 5.1, 5.2, 5.4 and Appendices B and C quote are built from
##  the csv files this script produces, and are not rebuilt here.
################################################################################

if (!file.exists("00_gie_pffc.R"))
  stop("run_simulations.R must be run from the code/ directory; the working ",
       "directory is currently '", getwd(), "'.")

source("00_gie_pffc.R")

sims <- c("04_sim_bias_mse.R",
          "05_sim_ci.R",
          "06_sim_heavytail.R",
          "07_sim_smallm.R")

if (!file.exists(file.path(TABLES_DIR, "values_realdata1.tex")))
  stop("run_simulations.R: ../tables/values_realdata1.tex not found. ",
       "07_sim_smallm.R reads the design of the first application from it. ",
       "Run 10_realdata1_ballbearings.R first.")

run <- function(f) {
  cat("\n", strrep("=", 78), "\n== ", f, "\n", strrep("=", 78), "\n", sep = "")
  t0 <- Sys.time()
  source(f, echo = FALSE)
  cat(sprintf("\n-- %s finished in %.1f minutes\n", f,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cat("run_simulations.R: the four Monte Carlo studies of Sections 5.1-5.4.\n",
    "                   This takes several hours.  Progress is reported as it\n",
    "                   goes, and each study writes its results to results/ as\n",
    "                   soon as it finishes, so an interrupted run keeps what\n",
    "                   it had completed.\n", sep = "")

t_all <- Sys.time()
for (f in sims) run(f)

cat("\n", strrep("=", 78), "\n", sep = "")
cat(sprintf("All four finished in %.1f hours.\n",
            as.numeric(difftime(Sys.time(), t_all, units = "hours"))))
cat("\nNow run\n    Rscript run_paper.R\n",
    "to rebuild Figures 4 and 5 and the numbers the running text quotes from\n",
    "these results, and then recompile the manuscript.\n\n", sep = "")
print(sessionInfo())
