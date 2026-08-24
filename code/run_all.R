################################################################################
##  run_all.R
##
##  Runs everything, in dependency order, from this directory.
##
##  From a terminal, with the working directory set to this folder:
##      Rscript run_all.R quick        # validation + Table 1 + Figures 1-3 +
##                                     # both data analyses          (minutes)
##      Rscript run_all.R              # the above plus the simulations (hours)
##
##  From the RStudio console, set the working directory to this folder first
##  (Session > Set Working Directory > To Source File Location) and then:
##      QUICK <- TRUE;  source("run_all.R")     # quick pass
##      QUICK <- FALSE; source("run_all.R")     # full run
##  The QUICK object is checked before the command-line argument, because
##  source() passes no arguments and the full run would otherwise start by
##  accident.
##
##  The long simulations are 04, 05, 06 and 07; 05 and 06 dominate and both run
##  in parallel by default (set PARALLEL <- FALSE inside either for a single
##  core).  All four report progress as they go, with an estimate of the time
##  remaining.
##
##  Output locations:
##      results/    csv and rds files with the raw numbers
##      ../tables/  LaTeX table fragments, \input by the manuscript
##      ../figures/ figures, in both EPS and PDF
##
##  After a full run, recompile the manuscript:
##      cd .. && pdflatex paper_QREI_revised.tex   (three times)
################################################################################

args  <- commandArgs(trailingOnly = TRUE)
quick <- if (exists("QUICK", inherits = TRUE)) {
  isTRUE(QUICK)
} else {
  length(args) > 0 && args[1] == "quick"
}

if (!file.exists("00_gie_pffc.R"))
  stop("run_all.R must be run from the code/ directory; the working directory ",
       "is currently '", getwd(), "'.")

cat(if (quick) "QUICK pass: validation, Table 1, Figures 1-3, both data analyses.\n"
    else "FULL run: everything, including the simulations. This takes hours.\n")

fast <- c("12_table1_translation.R",
          "01_figures.R",
          "10_realdata1_ballbearings.R",
          "11_realdata2_guineapig.R")
slow <- c("04_sim_bias_mse.R",
          "05_sim_ci.R",
          "06_sim_heavytail.R",
          "07_sim_smallm.R")
after <- c("08_figures_simulation.R")
checks <- file.path("validation",
                    c("v1_check_derivatives.R", "v2_check_invariance.R",
                      "v3_check_moment_formula.R", "v4_check_information.R",
                      "v5_check_information_and_mle.R"))

run <- function(f) {
  cat("\n", strrep("=", 78), "\n== ", f, "\n", strrep("=", 78), "\n", sep = "")
  t0 <- Sys.time()
  if (startsWith(f, "validation/")) {
    ## the validation scripts source ../00_gie_pffc.R, so run them from there
    owd <- getwd()
    setwd("validation")
    tryCatch(source(basename(f), echo = FALSE), finally = setwd(owd))
  } else {
    source(f, echo = FALSE)
  }
  cat(sprintf("\n-- %s finished in %.1f minutes\n", f,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

todo <- c(checks, fast, if (!quick) c(slow, after))
for (f in todo) run(f)

cat("\n", strrep("=", 78), "\n", sep = "")
cat("Done.  Session information:\n\n")
print(sessionInfo())
if (quick)
  cat("\nNOTE: run without the 'quick' argument to regenerate the simulation",
      "tables and Figures 4-5.\n")
