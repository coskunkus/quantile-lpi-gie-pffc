################################################################################
##  run_paper.R
##
##  Everything in the paper except the simulations.
##
##      cd code
##      Rscript run_paper.R
##
##  or, from the RStudio console, after Session > Set Working Directory > To
##  Source File Location:
##
##      source("run_paper.R")
##
##  What this does, in dependency order:
##
##      the seven validation scripts   analytical results against independent
##                                     numerical calculations
##      12_table1_translation.R        Table 1
##      01_figures.R                   Figures 1-3 and the constants Sections
##                                     2.2 and 2.3 quote about them
##      10_realdata1_ballbearings.R    Section 6.1 and 6.3, Table 9; builds the
##                                     censored sample from the 23 endurance
##                                     times and writes it into the manuscript
##      11_realdata2_guineapig.R       Section 6.2, Table 10
##      04_sim_bias_mse.R              Tables 2-4 and the numbers Section 5.1
##         (in reuse mode)             quotes, rebuilt from results/ without
##                                     re-running the study
##      08_figures_simulation.R        Figures 4 and 5, from results/
##      09_values_text.R               the numbers Sections 5.1, 5.2, 5.4 and
##                                     Appendices B and C quote, from results/
##
##  Runtime: roughly ten minutes, almost all of it in the validation scripts.
##
##  WHAT THIS DOES NOT DO.  It does not run the Monte Carlo studies of
##  Sections 5.1 to 5.4.  Those take hours; run_simulations.R does them.  The
##  csv files they leave behind are shipped with this package in results/, so
##  the scripts above that read them work without re-running anything, and any
##  table can be checked against the raw output it came from.  Tables 2 to 4
##  are rewritten here from the saved csv, so their layout and their captions
##  are always current; Tables 5 to 8 and Table 11 are left as they are, since
##  the scripts that write them do not have a reuse mode.
##
##  It leaves its output in ../tables/ and ../figures/: the LaTeX table
##  fragments and \newcommand macro files the manuscript \inputs, and the
##  figures. Copy those into the manuscript's directory and recompile it there
##  (pdflatex, three times).
################################################################################

if (!file.exists("00_gie_pffc.R"))
  stop("run_paper.R must be run from the code/ directory; the working ",
       "directory is currently '", getwd(), "'.")

source("00_gie_pffc.R")

checks <- file.path("validation",
                    c("v1_check_derivatives.R",
                      "v2_check_invariance.R",
                      "v3_check_moment_formula.R",
                      "v4_check_information.R",
                      "v5_check_information_and_mle.R",
                      "v6_check_bootstrap.R",
                      "v7_check_pffc_generator.R"))

main <- c("12_table1_translation.R",
          "01_figures.R",
          "10_realdata1_ballbearings.R",
          "11_realdata2_guineapig.R")

from_results <- c("04_sim_bias_mse.R",
                  "08_figures_simulation.R",
                  "09_values_text.R")

## 04_sim_bias_mse.R simulates when run on its own.  This option tells it to
## read results/sim_bias_mse.csv instead and go straight to the tables.
old_opt <- options(gie.reuse_results = TRUE)

## ---------------------------------------------------------------------------
## The simulation output this script consumes rather than produces
## ---------------------------------------------------------------------------
needed <- c("sim_bias_mse.csv", "sim_ci.csv", "sim_heavytail.csv",
            "sim_smallm.csv")
missing <- needed[!file.exists(file.path(RESULTS_DIR, needed))]
if (length(missing))
  stop("run_paper.R: results/ is missing ", paste(missing, collapse = ", "),
       ".\n  These are shipped with the package. If you have deleted them, ",
       "run\n  Rscript run_simulations.R first -- it takes several hours.")

run <- function(f) {
  cat("\n", strrep("=", 78), "\n== ", f, "\n", strrep("=", 78), "\n", sep = "")
  t0 <- Sys.time()
  if (startsWith(f, "validation/")) {
    ## the validation scripts source ../00_gie_pffc.R, so run them from there
    owd <- getwd(); setwd("validation")
    tryCatch(source(basename(f), echo = FALSE), finally = setwd(owd))
  } else {
    source(f, echo = FALSE)
  }
  cat(sprintf("\n-- %s finished in %.1f minutes\n", f,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

cat("run_paper.R: validation, Tables 1-4, Figures 1-5, both applications, and\n",
    "             the numbers quoted in the running text.  The Monte Carlo\n",
    "             studies are not re-run: Tables 2-4 are rebuilt from the saved\n",
    "             csv, and Tables 5-8 and 11 are left as they are.\n",
    sep = "")

## The option is restored even if a script stops, so that running
## 04_sim_bias_mse.R by hand afterwards in the same session simulates, as its
## header promises, rather than silently reading the saved csv.
tryCatch(for (f in c(checks, main, from_results)) run(f),
         finally = options(old_opt))

## ---------------------------------------------------------------------------
## Is the shipped Section 5.4 block still the right one?
## ---------------------------------------------------------------------------
##
## 07_sim_smallm.R mirrors the first application: it takes its shape, scale,
## limit and group size from the macro file that 10_realdata1_ballbearings.R
## writes.  If that application changes -- a different censoring plan, a
## different lower limit -- the shipped Table 8 no longer describes it, and no
## error would otherwise be raised, because the table is simply not rewritten
## here.  So the true index the shipped block was built at is compared with the
## one the current application implies.
app <- read_macros(file.path(TABLES_DIR, "values_realdata1.tex"))
sm  <- read.csv(file.path(RESULTS_DIR, "sim_smallm.csv"))
## On a fresh checkout tables/ is empty until the scripts have written into it,
## and values_smallm.tex is written by 07_sim_smallm.R, which this script does
## not run.  Its absence is normal, not an error; the comparison then rests on
## the true index alone, and says so.
sm_file <- file.path(TABLES_DIR, "values_smallm.tex")
smv <- if (file.exists(sm_file)) read_macros(sm_file) else list()
want <- CL_quantile(c(as.numeric(app$bbAlphaHat), as.numeric(app$bbLambdaHat)),
                    as.numeric(app$bbL))
have <- sm$C_true[1]

## The index pins down alpha, lambda and L, but not the censoring plan or the
## group size: two designs can share a true index and still not be the same
## experiment.  Those are compared separately.  Older runs of 07_sim_smallm.R
## did not record them, in which case the comparison is skipped and said so.
drift <- character(0)
if (!isTRUE(all.equal(want, have, tolerance = 5e-4)))
  drift <- c(drift, sprintf("true index: block %.4f, application %.4f", have, want))
for (nm in list(c("smRfirst", "bbRfirst", "censoring plan R_1"),
                c("smK", "bbK", "group size k"))) {
  if (!file.exists(sm_file)) {
    ## nothing to compare against; reported once, below
  } else if (is.null(smv[[nm[1]]])) {
    drift <- c(drift, paste0(nm[3], ": not recorded by the run that produced ",
                             "Table 8, so it cannot be checked"))
  } else if (!identical(as.numeric(smv[[nm[1]]]), as.numeric(app[[nm[2]]]))) {
    drift <- c(drift, sprintf("%s: block %s, application %s",
                              nm[3], smv[[nm[1]]], app[[nm[2]]]))
  }
}

cat("\n", strrep("=", 78), "\n", sep = "")
if (length(drift)) {
  cat("WARNING: the small-sample block of Section 5.4 does not match the first\n")
  cat("         application it claims to mirror.\n")
  for (d in drift) cat("  - ", d, "\n", sep = "")
  cat("  Table 8 and the numbers Section 5.4 quotes are therefore not\n")
  cat("  descriptions of the experiment reported in Section 6.1.\n")
  cat("  Fix: Rscript 07_sim_smallm.R   (30-60 minutes), then re-run this.\n")
} else if (!file.exists(sm_file)) {
  cat("Section 5.4 matches the first application on the true index (C_L^xi = ",
      sprintf("%.4f", have), ").\n", sep = "")
  cat("  The group size and censoring plan could not be compared: ",
      basename(sm_file), " does not\n  exist yet. It is written by ",
      "07_sim_smallm.R, which this script does not run.\n", sep = "")
} else {
  cat("Section 5.4 matches the first application: C_L^xi = ",
      sprintf("%.4f", have), ", k = ", app$bbK, ", R_1 = ", app$bbRfirst, ".\n",
      sep = "")
}

cat("\nDone.  Tables 5-8 and 11 were not regenerated; they do not depend on\n",
    "anything this script computes, except Table 8, which is checked above.\n",
    "The LaTeX fragments and macro files are in ", TABLES_DIR, "\n",
    "and the figures in ", FIGURES_DIR, ".\n\n", sep = "")
print(sessionInfo())
