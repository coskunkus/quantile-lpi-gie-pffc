################################################################################
##  08_figures_simulation.R
##
##  Figures 4 and 5.
##
##  Figure 4  fig4_relmse   relative MSE of the two indexes, with relative bias
##                          printed above each point (k = 2, alpha = 5)
##  Figure 5  fig5_ci       average length of the three intervals for the
##                          quantile-based index, with CP printed above each
##                          point (k = 2, alpha = 5)
##
##  IMPORTANT.  The earlier versions of these two scripts had the numbers typed
##  in by hand, copied from the tables, and in the case of Figure 5 rounded to
##  two decimals.  Both scripts now read the simulation output from
##  results/*.csv, so the figures cannot drift out of step with the tables.
##  Run 04_sim_bias_mse.R and 05_sim_ci.R first.
##
##  Runtime: a few seconds.
################################################################################

source("00_gie_pffc.R")
need("ggplot2", "dplyr")

f_pt <- file.path(RESULTS_DIR, "sim_bias_mse.csv")
f_ci <- file.path(RESULTS_DIR, "sim_ci.csv")
if (!file.exists(f_pt)) stop("Run 04_sim_bias_mse.R first: ", f_pt, " not found.")
if (!file.exists(f_ci)) stop("Run 05_sim_ci.R first: ", f_ci, " not found.")

SCHEME_ORDER <- c("Early", "Middle", "Equal", "Late")

save_both <- function(plot, stem, width, height) {
  for (dev in c("eps", "pdf"))
    ggplot2::ggsave(file.path(FIGURES_DIR, paste0(stem, ".", dev)),
                    plot = plot, device = dev, width = width, height = height,
                    units = "in")
  cat("Wrote ", stem, ".eps / .pdf\n", sep = "")
}

## ---------------------------------------------------------------------------
## Figure 4
## ---------------------------------------------------------------------------

pt <- read.csv(f_pt)
d4 <- subset(pt, k == 2 & alpha == 5)
d4 <- rbind(
  data.frame(m = d4$m, scheme = d4$scheme, index = "Moment",
             RMSE = d4$RMSE_m, RB = d4$RB_m),
  data.frame(m = d4$m, scheme = d4$scheme, index = "Quantile",
             RMSE = d4$RMSE_q, RB = d4$RB_q))
d4$scheme  <- factor(d4$scheme, levels = SCHEME_ORDER)
d4$index   <- factor(d4$index, levels = c("Moment", "Quantile"))
d4$m_facet <- factor(paste("m =", d4$m), levels = paste("m =", c(25, 50, 100)))

p4 <- ggplot2::ggplot(d4, ggplot2::aes(scheme, RMSE, shape = index,
                                       linetype = index, group = index)) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 3, fill = "white") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", RB)),
                     size = 2.8, vjust = -1.0, fontface = "italic") +
  ggplot2::facet_wrap(~ m_facet, scales = "free_y") +
  ggplot2::scale_shape_manual(values = c(Moment = 16, Quantile = 17)) +
  ggplot2::scale_linetype_manual(values = c(Moment = "solid", Quantile = "dashed")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.12, 0.20))) +
  ggplot2::labs(x = "Censoring scheme",
                y = expression(paste("Relative MSE, ", MSE(hat(C)) / C^2)),
                shape = "Index type", linetype = "Index type") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 strip.background = ggplot2::element_rect(fill = "grey90", colour = "grey20"),
                 strip.text = ggplot2::element_text(face = "bold"))
save_both(p4, "fig4_relmse", 8, 5)

## ---------------------------------------------------------------------------
## Figure 5
## ---------------------------------------------------------------------------

ci <- read.csv(f_ci)
d5r <- subset(ci, k == 2 & alpha == 5)
d5 <- rbind(
  data.frame(m = d5r$m, scheme = d5r$scheme, method = "ACI", AL = d5r$AL_ACI_q, CP = d5r$CP_ACI_q),
  data.frame(m = d5r$m, scheme = d5r$scheme, method = "PB",  AL = d5r$AL_PB_q,  CP = d5r$CP_PB_q),
  data.frame(m = d5r$m, scheme = d5r$scheme, method = "NB",  AL = d5r$AL_NB_q,  CP = d5r$CP_NB_q))
d5$scheme  <- factor(d5$scheme, levels = SCHEME_ORDER)
d5$method  <- factor(d5$method, levels = c("ACI", "PB", "NB"))
d5$m_facet <- factor(paste("m =", d5$m), levels = paste("m =", c(25, 50, 100)))

p5 <- ggplot2::ggplot(d5, ggplot2::aes(scheme, AL, shape = method,
                                       linetype = method, group = method)) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::geom_point(size = 3, fill = "white") +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", CP)), size = 2.6,
                     position = ggplot2::position_dodge(width = 0.55), vjust = -1.0) +
  ggplot2::facet_wrap(~ m_facet, scales = "free_y") +
  ggplot2::scale_shape_manual(values = c(ACI = 16, PB = 17, NB = 15)) +
  ggplot2::scale_linetype_manual(values = c(ACI = "solid", PB = "dashed", NB = "dotted")) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.22))) +
  ggplot2::labs(x = "Censoring scheme", y = "Average length (AL)",
                shape = "CI method", linetype = "CI method") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "top",
                 strip.background = ggplot2::element_rect(fill = "grey90", colour = "grey20"),
                 strip.text = ggplot2::element_text(face = "bold"))
save_both(p5, "fig5_ci", 8, 5)
