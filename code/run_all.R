################################################################################
##  run_all.R
##
##  Everything, from nothing: the four Monte Carlo studies followed by the rest
##  of the paper.  Several hours.
##
##      cd code
##      Rscript run_all.R
##
##  This is a convenience wrapper around the two scripts that do the work, and
##  it is almost never the one you want:
##
##      Rscript run_simulations.R   the Monte Carlo studies of Sections 5.1-5.4
##                                  Tables 2-8 and 11.  Hours.
##
##      Rscript run_paper.R         everything else -- validation, Table 1,
##                                  Figures 1-5, both applications, and every
##                                  number the running text quotes.  Minutes.
##
##  The simulations do not depend on anything else in the paper, and their raw
##  output ships with this package in results/.  So unless you mean to rebuild
##  that output from scratch, run_paper.R on its own reproduces the manuscript.
##
##  One dependency runs the other way: the small-sample block of Section 5.4
##  mirrors the first application and reads its design from the macro file
##  10_realdata1_ballbearings.R writes.  This script therefore builds the
##  applications first, then simulates, then rebuilds what reads the results.
################################################################################

if (!file.exists("00_gie_pffc.R"))
  stop("run_all.R must be run from the code/ directory; the working directory ",
       "is currently '", getwd(), "'.")

if (length(commandArgs(trailingOnly = TRUE)))
  message("run_all.R takes no arguments; it runs everything. ",
          "For part of it, use run_simulations.R or run_paper.R.")

cat(strrep("=", 78), "\n",
    "run_all.R: the applications, then the simulations, then the rest.\n",
    "           Several hours.  Ctrl-C is safe: each study writes its results\n",
    "           as soon as it finishes.\n", strrep("=", 78), "\n", sep = "")

## The applications first, so that Section 5.4 mirrors the current one.
source("10_realdata1_ballbearings.R", echo = FALSE)
source("11_realdata2_guineapig.R",    echo = FALSE)

source("run_simulations.R", echo = FALSE)
source("run_paper.R",       echo = FALSE)
