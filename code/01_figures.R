################################################################################
##  01_figures.R
##
##  Figures 1, 2 and 3 of the paper.  Each figure is written to ../figures/ in
##  both EPS (for the journal) and PDF (so that the manuscript compiles with
##  pdflatex without a ghostscript conversion step).  The manuscript refers to
##  the files without an extension, so whichever format the compiler prefers is
##  picked up automatically.
##
##  Figure 1  fig1_pdf          GIE densities, lambda = 1, several alpha
##  Figure 2  fig2_momconv      Gauss-Legendre convergence, alpha = 4, lambda = 5
##  Figure 3  fig3_comparison   C_L^M and C_L^xi against alpha and against lambda
##
##  Runtime: a few seconds.
################################################################################

source("00_gie_pffc.R")
need("statmod", "ggplot2", "tidyr", "patchwork")

save_both <- function(plot, stem, width, height) {
  for (dev in c("eps", "pdf"))
    ggplot2::ggsave(file.path(FIGURES_DIR, paste0(stem, ".", dev)),
                    plot = plot, device = dev,
                    width = width, height = height, units = "in")
  cat("Wrote ", stem, ".eps / .pdf\n", sep = "")
}

## ---------------------------------------------------------------------------
## Figure 1: densities
## ---------------------------------------------------------------------------

alpha_vals <- c(0.5, 0.9, 2, 5, 10)
xs <- seq(0.01, 8, length.out = 800)
d1 <- do.call(rbind, lapply(alpha_vals, function(a)
  data.frame(x = xs, y = gie_pdf(xs, a, 1), alpha = factor(a, levels = alpha_vals))))

p1 <- ggplot2::ggplot(d1, ggplot2::aes(x, y, linetype = alpha)) +
  ggplot2::geom_line(linewidth = 0.9, colour = "black") +
  ggplot2::scale_linetype_manual(
    values = c("solid", "dashed", "dotted", "dotdash", "longdash"),
    labels = lapply(alpha_vals, function(a) bquote(alpha == .(a)))) +
  ggplot2::coord_cartesian(xlim = c(0, 5)) +
  ggplot2::labs(x = "x", y = "PDF", linetype = NULL) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(legend.position = "inside",
                 legend.position.inside = c(0.88, 0.75),
                 legend.background = ggplot2::element_blank())
save_both(p1, "fig1_pdf", 9, 4.2)

## ---------------------------------------------------------------------------
## Figure 2: quadrature convergence
## ---------------------------------------------------------------------------

a0 <- 4; l0 <- 5
ref <- c(gie_moment_gl(1, a0, l0, 800), gie_moment_gl(2, a0, l0, 800))
cat(sprintf("Reference values (N = 800):  E[X] = %.6f   E[X^2] = %.6f\n", ref[1], ref[2]))

Ns <- 2:10
d2 <- rbind(
  data.frame(N = Ns, value = vapply(Ns, function(n) gie_moment_gl(1, a0, l0, n), 0),
             moment = "Mean (E[X])"),
  data.frame(N = Ns, value = vapply(Ns, function(n) gie_moment_gl(2, a0, l0, n), 0),
             moment = "Second moment (E[X^2])"))
refs <- data.frame(moment = c("Mean (E[X])", "Second moment (E[X^2])"), y = ref)

p2 <- ggplot2::ggplot(d2, ggplot2::aes(N, value)) +
  ggplot2::geom_line(ggplot2::aes(linetype = moment), linewidth = 0.8, colour = "black") +
  ggplot2::geom_point(ggplot2::aes(shape = moment), size = 2.2, fill = "white") +
  ggplot2::geom_hline(data = refs, ggplot2::aes(yintercept = y),
                      linetype = "dashed", colour = "gray40", linewidth = 0.7) +
  ggplot2::facet_wrap(~ moment, scales = "free_y") +
  ggplot2::scale_shape_manual(values = c(16, 17)) +
  ggplot2::scale_linetype_manual(values = c("solid", "longdash")) +
  ggplot2::labs(x = "Number of nodes (N)", y = "Moment value") +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(legend.position = "none",
                 strip.background = ggplot2::element_rect(fill = "gray90", colour = "black"),
                 strip.text = ggplot2::element_text(face = "bold"))
save_both(p2, "fig2_momconv", 9.5, 4.0)

## ---------------------------------------------------------------------------
## Figure 3: behaviour of the two indexes
## ---------------------------------------------------------------------------
## Left panel:  lambda fixed at 5, alpha varying.
## Right panel: alpha fixed at 4, lambda varying.
## Both panels: L = 1.  These three constants are stated in the caption.

L_FIG <- 1
LAMBDA_FIXED <- 5
ALPHA_FIXED  <- 4

al <- seq(2.1, 100, by = 0.5)
dA <- data.frame(
  alpha  = al,
  C_L_M  = vapply(al, function(a) CL_moment(c(a, LAMBDA_FIXED), L_FIG, 400), 0),
  C_L_xi = vapply(al, function(a) CL_quantile(c(a, LAMBDA_FIXED), L_FIG), 0))
lm_ <- seq(1, 20, by = 0.5)
dL <- data.frame(
  lambda = lm_,
  C_L_M  = vapply(lm_, function(l) CL_moment(c(ALPHA_FIXED, l), L_FIG, 400), 0),
  C_L_xi = vapply(lm_, function(l) CL_quantile(c(ALPHA_FIXED, l), L_FIG), 0))

styl <- list(
  ggplot2::scale_colour_manual(values = c(C_L_M = "black", C_L_xi = "gray40"),
                               labels = c(C_L_M = "Moment-based", C_L_xi = "Quantile-based")),
  ggplot2::scale_linetype_manual(values = c(C_L_M = "solid", C_L_xi = "dashed"),
                                 labels = c(C_L_M = "Moment-based", C_L_xi = "Quantile-based")),
  ggplot2::theme_bw(base_size = 13),
  ggplot2::theme(legend.title = ggplot2::element_blank(),
                 plot.title = ggplot2::element_text(hjust = 0.5)))

pA <- ggplot2::ggplot(tidyr::pivot_longer(dA, -alpha, names_to = "Index"),
                      ggplot2::aes(alpha, value, colour = Index, linetype = Index)) +
  ggplot2::geom_line(linewidth = 1) + styl +
  ggplot2::labs(title = "Effect of shape parameter",
                x = expression(paste("Shape parameter (", alpha, ")")),
                y = "Performance index value")
pB <- ggplot2::ggplot(tidyr::pivot_longer(dL, -lambda, names_to = "Index"),
                      ggplot2::aes(lambda, value, colour = Index, linetype = Index)) +
  ggplot2::geom_line(linewidth = 1) + styl +
  ggplot2::labs(title = "Effect of scale parameter",
                x = expression(paste("Scale parameter (", lambda, ")")), y = NULL)

p3 <- (pA + pB) + patchwork::plot_layout(guides = "collect") &
  ggplot2::theme(legend.position = "bottom")
save_both(p3, "fig3_comparison", 10, 5.2)

cat("\nFigure 3 constants (quoted in the caption): lambda =", LAMBDA_FIXED,
    "in the left panel, alpha =", ALPHA_FIXED, "in the right panel, L =", L_FIG, "\n")
