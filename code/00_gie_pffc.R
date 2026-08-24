################################################################################
##  00_gie_pffc.R
##
##  Shared functions for
##    "A General Quantile-Based Lifetime Performance Index and its Application
##     to the GIE Distribution under Progressive First-Failure Censoring"
##
##  Every closed-form expression displayed in the paper is implemented here and
##  is on the path that produces the reported numbers, not merely checked in a
##  side script.  The correspondence is:
##
##    Eq. (1)  density                      gie_pdf, gie_log_pdf
##    Eq. (2)  distribution function        gie_cdf, gie_survival, gie_log_survival
##    Eq. (3)  closed-form moment, integer  gie_moment_exact   (used by gie_moment)
##    Eq. (4)  Gauss-Legendre moment        gie_moment_gl      (used by gie_moment)
##    Eq. (6)  moment-based index           CL_moment
##    Eq. (7)  quantile function            gie_quantile
##    Eq. (8)  quantile-based index         CL_quantile
##    Eq. (9)  h_p(alpha)                   gie_h
##    Eq. (10) index through (alpha, L/lam) CL_quantile_ratio  (the implementation
##                                          CL_quantile actually calls)
##    Eq. (11) translation constants        translation_constants
##    Eq. (12) affine map between indexes   checked in 12_ and reported in 10_
##    Eq. (14) log-likelihood               neg_loglik
##    Eq. (15) fixed-point iteration for    lambda_fixed_point  (the primary
##             the MLE of lambda            estimation route in fit_mle)
##             explicit MLE of alpha        alpha_hat_given_lambda
##             score / estimating equations loglik_score
##    Eq. (17) Fisher information; the      loglik_hessian, observed_information
##             Hessian displayed below it
##    Eq. (18) delta-method variance        delta_se
##    Remark   derivatives of xi and C      dxi_dalpha, dxi_dlambda, grad_CL_quantile
##    Eq. (19) asymptotic interval          ci_asymptotic
##    Eq. (21) bias-corrected normal boot   ci_normal_boot
##    Eq. (22) percentile bootstrap         ci_percentile
##    Eq. (23) one-sided decision rule      lower_bound_one_sided
##    App. B   bounds on the information    checked in validation/v4_
##
##  Sourcing this file has no side effects other than defining objects and
##  creating the output directories.
##
##  Notation follows the paper:
##    alpha  shape parameter          (alpha > 0)
##    lambda scale parameter          (lambda > 0)
##    m      number of observed failures
##    k      group size
##    R      censoring plan, a vector of length m of non-negative integers
##    L      lower specification limit
##    n      number of groups, n = m + sum(R)
################################################################################

## ---------------------------------------------------------------------------
## 1.  The GIE distribution
## ---------------------------------------------------------------------------
##
## Throughout, 1 - exp(-lambda/x) is computed as -expm1(-lambda/x).  The
## cancellation as lambda/x -> 0 is severe enough that the naive form loses all
## precision, which in the likelihood shows up as replications being silently
## discarded.

#' Density of GIE(alpha, lambda).  Eq. (1).
gie_pdf <- function(x, alpha, lambda) {
  out <- rep(0, length(x))
  ok  <- x > 0
  if (any(ok)) out[ok] <- exp(gie_log_pdf(x[ok], alpha, lambda))
  out
}

#' log density.  This is the form used inside the likelihood.
gie_log_pdf <- function(x, alpha, lambda) {
  z <- -lambda / x
  log(alpha) + log(lambda) - 2 * log(x) + z + (alpha - 1) * log(-expm1(z))
}

#' Distribution function of GIE(alpha, lambda).  Eq. (2).
gie_cdf <- function(x, alpha, lambda) {
  out <- rep(0, length(x))
  ok  <- x > 0
  if (any(ok)) out[ok] <- 1 - (-expm1(-lambda / x[ok]))^alpha
  out
}

#' Survival function, S(x) = (1 - exp(-lambda/x))^alpha.
gie_survival <- function(x, alpha, lambda) 1 - gie_cdf(x, alpha, lambda)

#' log survival function.  Used inside the likelihood.
gie_log_survival <- function(x, alpha, lambda) alpha * log(-expm1(-lambda / x))

#' h_p(alpha) of Eq. (9).  The quantile function factorises as
#' xi(p) = lambda * h_p(alpha), which is why the scale parameter enters the
#' index only through the ratio L/lambda.
gie_h <- function(p, alpha) -1 / log1p(-(1 - p)^(1 / alpha))

#' Quantile function xi(p).  Eq. (7), written through Eq. (9).
gie_quantile <- function(p, alpha, lambda) lambda * gie_h(p, alpha)

## ---------------------------------------------------------------------------
## 2.  Moments
## ---------------------------------------------------------------------------

#' Largest integer shape for which the closed form of Eq. (3) may be used.
#'
#' The sum in Eq. (3) is an (alpha-1)-st alternating finite difference: the
#' individual terms grow like the binomial coefficient C(alpha-1, j) while the
#' result stays O(1), so in double precision it loses roughly one digit per unit
#' of alpha.  Measured relative error against adaptive integration at lambda = 2:
#'
#'      alpha      10       15       20       25       30       50
#'      r = 2    6e-10    4e-09    3e-08    2e-04    2e-03    5e+03
#'
#' Beyond about alpha = 22 the expression is unusable and can even return a
#' negative "second moment".  We therefore use it only up to alpha = 20 and fall
#' back to the quadrature of Eq. (4)
#' above that.  At alpha = 20 the measured relative error is about 1e-9 for
#' the mean and 1e-7 for the second moment; the quadrature is accurate to about
#' 1e-6 uniformly in alpha, so
#' nothing is lost.  This is a property of the arithmetic, not of the formula:
#' Eq. (3) is exact, as validation/v3_ confirms for alpha in 2..8.
ALPHA_EXACT_MAX <- 20

#' r-th raw moment in closed form, integer alpha only.  Eq. (3), proved in
#' Appendix A.
gie_moment_exact <- function(r, alpha, lambda) {
  if (r >= alpha || r <= 0 || alpha != floor(alpha))
    stop("gie_moment_exact(): need integer alpha and 0 < r < alpha.")
  j <- 0:(alpha - 1)
  s <- sum((-1)^j * (j + 1)^(r - 1) * ifelse(j == 0, 0, log(j + 1)) /
             (gamma(j + 1) * gamma(alpha - j)))
  ((-lambda)^r / factorial(r - 1)) * gamma(alpha + 1) * s
}

#' r-th raw moment by Gauss-Legendre quadrature.  Eq. (4).
#'
#' The substitution used matters.  Mapping x = t/(1-t) directly places the nodes
#' at fixed absolute values of x while the distribution lives on the scale
#' lambda, so accuracy degrades as lambda grows; and it leaves an endpoint
#' singularity of order alpha-1-r, so the rule converges only algebraically when
#' alpha is close to r.  We therefore substitute u = lambda/x, giving
#'
#'   E[X^r] = alpha * lambda^r * int_0^Inf u^(-r) e^(-u) (1-e^(-u))^(alpha-1) du,
#'
#' which is manifestly proportional to lambda^r; then w = 1 - e^(-u) maps the
#' range to (0,1), and finally w = s^(1/(alpha-r)) absorbs the singularity:
#'
#'   E[X^r] = (alpha * lambda^r / (alpha-r)) * int_0^1 phi(s^(1/(alpha-r))) ds,
#'            phi(w) = (w / -log(1-w))^r,   phi(0) = 1.
#'
#' Accurate to about 1e-5 with 100 nodes and 1e-6 with 400, uniformly in alpha.
#' Cache for the Gauss-Legendre rule.
#'
#'
#' statmod::gauss.quad(n) computes the nodes and weights from the eigenvalues of
#' an n x n tridiagonal matrix, which for n = 400 costs milliseconds.  The rule
#' depends only on n, not on any parameter, but CL_moment() is called once per
#' bootstrap replicate -- about 130,000 times per configuration in
#' 05_sim_ci.R, and so a quarter of a million rule constructions -- so
#' recomputing it each time dominated the entire simulation, costing far more
#' than the estimation it was serving.  Caching by n reduces that to one
#' evaluation per rule per R process.
.gq_cache <- new.env(parent = emptyenv())

gauss_legendre <- function(n) {
  key <- as.character(n)
  if (!exists(key, envir = .gq_cache, inherits = FALSE))
    assign(key, statmod::gauss.quad(n, kind = "legendre"), envir = .gq_cache)
  get(key, envir = .gq_cache, inherits = FALSE)
}

gie_moment_gl <- function(r, alpha, lambda, n_nodes = 400) {
  if (alpha <= r || r < 0)
    stop("gie_moment_gl(): need alpha > r >= 0 for the moment to exist.")
  gq  <- gauss_legendre(n_nodes)
  s   <- (gq$nodes + 1) / 2
  wt  <- gq$weights / 2
  p   <- alpha - r
  w   <- s^(1 / p)
  phi <- ifelse(w < 1e-300, 1, (w / (-log1p(-w)))^r)
  alpha * lambda^r / p * sum(wt * phi)
}

#' r-th raw moment, dispatching exactly as Section 2.1 of the paper describes:
#' the closed form of Eq. (3) when alpha is a positive integer with r < alpha,
#' and the quadrature of Eq. (4) otherwise.
#'
#' The test is exact equality with round(alpha), not a tolerance: a shape
#' estimate produced by optimisation is never an exact integer, so estimates go
#' through the quadrature while the true parameter values used to define the
#' targets of the simulation go through the closed form.  That is the intended
#' behaviour and it is what the paper claims.
gie_moment <- function(r, alpha, lambda, n_nodes = 400) {
  if (alpha > 0 && alpha == round(alpha) && alpha <= ALPHA_EXACT_MAX &&
      r == round(r) && r > 0 && r < alpha)
    gie_moment_exact(r, alpha, lambda)
  else
    gie_moment_gl(r, alpha, lambda, n_nodes)
}

## ---------------------------------------------------------------------------
## 3.  The two performance indexes
## ---------------------------------------------------------------------------

#' Quantile-based index written through Eq. (10), as a function of the shape and
#' the ratio L/lambda only.  This is the implementation; CL_quantile is a
#' wrapper on it, so Eq. (10) is on the path of every reported number.
CL_quantile_ratio <- function(alpha, L_over_lambda) {
  (gie_h(0.50, alpha) - L_over_lambda) /
    (gie_h(0.75, alpha) - gie_h(0.25, alpha))
}

#' Quantile-based index C_L^xi.  Eq. (8).
#' `par` is c(alpha, lambda); this signature is what numDeriv::grad expects.
CL_quantile <- function(par, L) {
  alpha <- par[1]; lambda <- par[2]
  if (anyNA(par) || alpha <= 0 || lambda <= 0) return(NA_real_)
  v <- CL_quantile_ratio(alpha, L / lambda)
  if (!is.finite(v)) return(NA_real_)
  v
}

#' Moment-based index C_L^M.  Eq. (6).  NA when alpha <= 2, where the second
#' moment does not exist and the index is undefined.
CL_moment <- function(par, L, n_nodes = 400) {
  alpha <- par[1]; lambda <- par[2]
  if (anyNA(par) || alpha <= 2 || lambda <= 0) return(NA_real_)
  mu  <- tryCatch(gie_moment(1, alpha, lambda, n_nodes), error = function(e) NA_real_)
  ex2 <- tryCatch(gie_moment(2, alpha, lambda, n_nodes), error = function(e) NA_real_)
  if (is.na(mu) || is.na(ex2)) return(NA_real_)
  v <- ex2 - mu^2
  if (!is.finite(v) || v <= 0) return(NA_real_)
  (mu - L) / sqrt(v)
}

#' Translation constants kappa and delta of Eq. (11).  Defined only for
#' alpha > 2; both are free of lambda, so they are evaluated at lambda = 1.
translation_constants <- function(alpha, n_nodes = 800) {
  if (alpha <= 2) return(c(kappa = NA_real_, delta = NA_real_))
  mu  <- gie_moment(1, alpha, 1, n_nodes)
  ex2 <- gie_moment(2, alpha, 1, n_nodes)
  sdv <- sqrt(ex2 - mu^2)
  iqr <- gie_quantile(0.75, alpha, 1) - gie_quantile(0.25, alpha, 1)
  med <- gie_quantile(0.50, alpha, 1)
  c(kappa = sdv / iqr, delta = (mu - med) / iqr)
}

## ---------------------------------------------------------------------------
## 4.  Derivatives (the Remark in Section 3, and Eq. (16))
## ---------------------------------------------------------------------------

#' d xi(p) / d lambda = h_p(alpha).
dxi_dlambda <- function(p, alpha) gie_h(p, alpha)

#' d xi(p) / d alpha.
dxi_dalpha <- function(p, alpha, lambda) {
  q  <- (1 - p)^(1 / alpha)
  lg <- log1p(-q)
  lambda * log(1 - p) * q / (alpha^2 * (1 - q) * lg^2)
}

#' Gradient of C_L^xi with respect to (alpha, lambda), by the quotient rule.
#' This analytical gradient, not a numerical one, is what the delta method uses
#' for the quantile-based index throughout the paper.
grad_CL_quantile <- function(par, L) {
  alpha <- par[1]; lambda <- par[2]
  S <- gie_quantile(0.50, alpha, lambda) - L
  D <- gie_quantile(0.75, alpha, lambda) - gie_quantile(0.25, alpha, lambda)
  dS_a <- dxi_dalpha(0.50, alpha, lambda)
  dD_a <- dxi_dalpha(0.75, alpha, lambda) - dxi_dalpha(0.25, alpha, lambda)
  dS_l <- dxi_dlambda(0.50, alpha)
  dD_l <- dxi_dlambda(0.75, alpha) - dxi_dlambda(0.25, alpha)
  c((dS_a * D - S * dD_a) / D^2,
    (dS_l * D - S * dD_l) / D^2)
}

#' Score (gradient of the log-likelihood).  These are the estimating equations
#' of Section 3, before they are rearranged into the explicit expression for
#' alpha_hat and the fixed-point iteration of Eq. (15):
#'
#'   dl/dalpha  = m/alpha + k sum (R_i+1) log(1 - exp(-lambda/x_i))
#'   dl/dlambda = m/lambda - sum 1/x_i
#'                + sum (alpha k(R_i+1) - 1) x_i^-1 e^(-lam/x_i)
#'                      (1 - e^(-lam/x_i))^-1
#'
#' Used by fit_mle() to confirm that the point returned by the fixed-point
#' iteration really is a stationary point of the likelihood.
loglik_score <- function(par, x, R, k) {
  alpha <- par[1]; lambda <- par[2]
  m  <- length(x)
  om <- -expm1(-lambda / x)
  lw <- log(om)
  d  <- (1 / x) * exp(-lambda / x) / om
  c(m / alpha + k * sum((R + 1) * lw),
    m / lambda - sum(1 / x) + sum((alpha * k * (R + 1) - 1) * d))
}

#' Hessian of the log-likelihood, from the closed forms displayed in Section 3
#' immediately after Eq. (17):
#'
#'   d2 l / d alpha^2        = -m / alpha^2
#'   d2 l / d alpha d lambda = sum k(R_i+1) e^(-lam/x_i) / (x_i (1 - e^(-lam/x_i)))
#'   d2 l / d lambda^2       = -m/lam^2
#'                             - sum (alpha k(R_i+1) - 1) e^(-lam/x_i)
#'                                   / (x_i^2 (1 - e^(-lam/x_i))^2)
loglik_hessian <- function(par, x, R, k) {
  alpha <- par[1]; lambda <- par[2]
  m <- length(x)
  z  <- -lambda / x
  ez <- exp(z)
  om <- -expm1(z)                        # 1 - exp(-lambda/x), stable
  t1 <- ez / (x * om)
  t2 <- ez / (x^2 * om^2)
  H <- matrix(NA_real_, 2, 2)
  H[1, 1] <- -m / alpha^2
  H[1, 2] <- H[2, 1] <- sum(k * (R + 1) * t1)
  H[2, 2] <- -m / lambda^2 - sum((alpha * k * (R + 1) - 1) * t2)
  H
}

#' Observed information, Eq. (17) evaluated at the estimates: minus the Hessian
#' above.  This replaces the numerical Hessian returned by optim(); the two
#' agree to five or six significant figures, which validation/v5_ checks.
observed_information <- function(par, x, R, k) -loglik_hessian(par, x, R, k)

## ---------------------------------------------------------------------------
## 5.  PFFC data generation
## ---------------------------------------------------------------------------

#' Build one of the four censoring plans used in the paper.  The total number of
#' censored groups is `val`, which defaults to m, so that n = 2m for every plan.
make_scheme <- function(m, type, val = m) {
  R <- rep(0, m)
  switch(type,
         "Early"  = { R[1] <- val },
         "Middle" = { R[max(1L, floor(m / 2))] <- val },
         "Late"   = { R[m] <- val },
         "Equal"  = { R[] <- val / m },
         stop("make_scheme(): unknown scheme type '", type, "'"))
  R
}

#' Generate one PFFC sample of size m from GIE(alpha, lambda).
#'
#' Algorithm of Balakrishnan & Aggarwala (2000): generate a progressively
#' Type-II censored uniform sample, then transform through the quantile function
#' of F_k(x) = 1 - (1 - F(x))^k, which by Wu & Kus (2009) is the correct marginal
#' for a first-failure-censored observation.  Since
#' F_k(x) = 1 - (1 - e^(-lam/x))^(alpha k), its quantile function is the GIE
#' quantile function of Eq. (6) with alpha replaced by alpha*k.
generate_pffc <- function(m, k, R, alpha, lambda) {
  n     <- sum(R) + m
  cumR  <- c(0, cumsum(R)[-m])
  denom <- n - cumR - seq_len(m) + 1
  E     <- cumsum(rexp(m) / denom)
  u     <- pmin(-expm1(-E), 1 - 1e-14)
  gie_quantile(u, alpha * k, lambda)
}

## ---------------------------------------------------------------------------
## 6.  Likelihood and maximum likelihood estimation
## ---------------------------------------------------------------------------

#' Negative log-likelihood.  Built directly from the density and survival
#' function, in the form displayed in Section 3:
#'   l = m log k + sum log f(x_i) + sum (k(R_i+1) - 1) log S(x_i).
neg_loglik <- function(par, x, R, k) {
  alpha <- par[1]; lambda <- par[2]
  if (!is.finite(alpha) || !is.finite(lambda) || alpha <= 0 || lambda <= 0)
    return(Inf)
  if (anyNA(x) || any(x <= 0)) return(Inf)
  lf <- gie_log_pdf(x, alpha, lambda)
  ls <- gie_log_survival(x, alpha, lambda)
  if (any(!is.finite(lf)) || any(!is.finite(ls))) return(Inf)
  ll <- length(x) * log(k) + sum(lf) + sum((k * (R + 1) - 1) * ls)
  if (!is.finite(ll)) return(Inf)
  -ll
}

#' Closed-form MLE of alpha given lambda, from Section 3:
#'   alpha_hat = -m / (k sum (R_i+1) log(1 - exp(-lambda/x_i))).
alpha_hat_given_lambda <- function(lambda, x, R, k) {
  lw <- log(-expm1(-lambda / x))
  -length(x) / (k * sum((R + 1) * lw))
}

#' One step of the fixed-point iteration of Eq. (15):
#'   lambda^(h+1) = m { sum 1/x_i
#'                      - sum (alpha k(R_i+1) - 1) x_i^-1 e^(-lam/x_i)
#'                            (1 - e^(-lam/x_i))^-1 }^-1
lambda_fixed_point_step <- function(lambda, alpha, x, R, k) {
  d <- (1 / x) * exp(-lambda / x) / (-expm1(-lambda / x))
  D <- sum(1 / x) - sum((alpha * k * (R + 1) - 1) * d)
  if (!is.finite(D) || D <= 0) return(NA_real_)
  length(x) / D
}

#' The estimation procedure of Section 3: alternate the fixed-point step for
#' lambda with the closed-form expression for alpha until both stabilise.
#'
#' Returns NULL if the iteration leaves the parameter space or fails to converge
#' within `maxit`.  The cap is deliberately modest: the iteration is only
#' linearly convergent, so a sample that has not settled within a few hundred
#' steps is better handed to the direct maximisation than ground through
#' thousands more.  Both outcomes occur for a small fraction of samples -- see the
#' `route` field of fit_mle() and the counts the simulation scripts report --
#' which is why fit_mle() falls back to a direct maximisation in those cases.
lambda_fixed_point <- function(x, R, k, lambda0, tol = 1e-11, maxit = 500) {
  lambda <- lambda0
  alpha  <- alpha_hat_given_lambda(lambda, x, R, k)
  if (!is.finite(alpha) || alpha <= 0) return(NULL)
  for (it in seq_len(maxit)) {
    lam_new <- lambda_fixed_point_step(lambda, alpha, x, R, k)
    if (!is.finite(lam_new) || lam_new <= 0) return(NULL)
    a_new <- alpha_hat_given_lambda(lam_new, x, R, k)
    if (!is.finite(a_new) || a_new <= 0) return(NULL)
    done <- abs(lam_new - lambda) <= tol * max(1, lambda) &&
            abs(a_new  - alpha)  <= tol * max(1, alpha)
    lambda <- lam_new; alpha <- a_new
    if (done) return(list(par = c(alpha, lambda), iterations = it))
  }
  NULL
}

#' Starting value: profile alpha out analytically and maximise the resulting
#' one-dimensional profile likelihood over a log-spaced grid in lambda.  Cheap,
#' and it removes the dependence of the answer on a hand-chosen starting point.
profile_start <- function(x, R, k, n_grid = 200) {
  grid <- exp(seq(log(1e-3 * stats::median(x)),
                  log(1e3 * stats::median(x)), length.out = n_grid))
  prof <- vapply(grid, function(l) {
    a <- alpha_hat_given_lambda(l, x, R, k)
    if (!is.finite(a) || a <= 0) return(Inf)
    neg_loglik(c(a, l), x, R, k)
  }, numeric(1))
  if (all(!is.finite(prof))) return(NULL)
  i <- which.min(prof)
  list(lambda = grid[i], loglik = -prof[i])
}

#' Maximum likelihood estimation, following Section 3.
#'
#' (i)   A profile-likelihood grid supplies a starting value: alpha is profiled
#'       out with the closed-form expression, leaving a one-dimensional search
#'       in lambda.
#' (ii)  The fixed-point iteration of Eq. (15), alternated with the explicit
#'       expression for alpha_hat, is run from there.  This is the procedure the
#'       paper describes and it is the primary route.
#' (iii) The result is accepted only if it is a genuine interior maximum: the
#'       closed-form score must vanish, the observed information of Eq. (17)
#'       must be positive definite, and the likelihood must be no worse than at
#'       the starting value.  If the iteration fails to converge or any of these
#'       checks fails, a direct maximisation with optim()/L-BFGS-B is used
#'       instead.
#'
#' The `route` field records which of the two produced the answer, so the
#' scripts can report how often the published iteration sufficed on its own.
#' Because optim() is only invoked on the residual cases, this is faster than
#' maximising numerically for every sample.
#'
#' Returns a list with par, loglik, convergence (0 = success), information (the
#' analytic observed information of Eq. (17)), route and iterations; or NULL if
#' both routes fail.  With numeric_hessian = TRUE the numerical Hessian from
#' optim is also returned, for the cross-check in validation/v5_.
fit_mle <- function(x, R, k, start = NULL, numeric_hessian = FALSE,
                    score_tol = 1e-6) {
  if (is.null(start)) {
    ps <- profile_start(x, R, k)
    if (is.null(ps)) return(NULL)
    lam0 <- ps$lambda; ref_ll <- ps$loglik
  } else {
    lam0 <- start[2]
    a0 <- if (length(start) >= 1 && is.finite(start[1]) && start[1] > 0)
      start[1] else alpha_hat_given_lambda(lam0, x, R, k)
    ref_ll <- -neg_loglik(c(a0, lam0), x, R, k)
  }
  if (!is.finite(lam0) || lam0 <= 0) return(NULL)

  ## ---- (ii) the paper's fixed-point iteration ----
  fp <- lambda_fixed_point(x, R, k, lam0)
  if (!is.null(fp)) {
    ll <- -neg_loglik(fp$par, x, R, k)
    I  <- observed_information(fp$par, x, R, k)
    sc <- loglik_score(fp$par, x, R, k)
    m <- length(x)
    ok <- is.finite(ll) && (!is.finite(ref_ll) || ll >= ref_ll - 1e-8) &&
          is_posdef(I) &&
          abs(sc[1]) <= score_tol * (m / fp$par[1]) &&
          abs(sc[2]) <= score_tol * (m / fp$par[2])
    if (ok)
      return(c(list(par = fp$par, loglik = ll, value = -ll, convergence = 0L,
                    information = I, route = "fixedpoint",
                    iterations = fp$iterations, score = sc),
               if (numeric_hessian) list(hessian_numeric =
                 numeric_hessian_from_optim(fp$par, x, R, k)) else NULL))
  }

  ## ---- (iii) fallback: direct maximisation ----
  st <- if (is.null(fp)) c(alpha_hat_given_lambda(lam0, x, R, k), lam0) else fp$par
  if (!is.finite(st[1]) || st[1] <= 0) st <- c(1, lam0)
  op <- tryCatch(
    stats::optim(par = st, fn = neg_loglik, x = x, R = R, k = k,
                 method = "L-BFGS-B", lower = c(1e-8, 1e-8),
                 hessian = numeric_hessian),
    error = function(e) NULL)
  if (is.null(op) || op$convergence != 0 || !is.finite(op$value)) return(NULL)
  par <- op$par
  out <- list(par = par, loglik = -op$value, value = op$value, convergence = 0L,
              information = observed_information(par, x, R, k),
              route = "optim", iterations = NA_integer_,
              score = loglik_score(par, x, R, k))
  if (numeric_hessian) out$hessian_numeric <- op$hessian
  out
}

#' The numerical Hessian optim would have produced at `par`, used only when
#' fit_mle() is asked for it but took the fixed-point route and so never called
#' optim.  One extra optim() call from the optimum; a few evaluations.
numeric_hessian_from_optim <- function(par, x, R, k) {
  op <- tryCatch(
    stats::optim(par = par, fn = neg_loglik, x = x, R = R, k = k,
                 method = "L-BFGS-B", lower = c(1e-8, 1e-8), hessian = TRUE,
                 control = list(maxit = 1)),
    error = function(e) NULL)
  if (is.null(op)) NULL else op$hessian
}

## ---------------------------------------------------------------------------
## 7.  Standard errors and confidence intervals
## ---------------------------------------------------------------------------

#' Is a 2x2 symmetric matrix positive definite?  This is the numerical check on
#' the nonsingularity condition discussed in Section 3 and Appendix B.
is_posdef <- function(M) {
  ev <- tryCatch(eigen(M, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NA_real_)
  !anyNA(ev) && all(is.finite(ev)) && min(ev) > 0
}

#' Delta-method standard error of an index, Eq. (18).
#'
#' `fit` is the object returned by fit_mle(); its `information` component is the
#' analytic observed information of Eq. (16).  The gradient of the quantile-based
#' index is the analytical one from the Remark in Section 3; for the
#' moment-based index, whose value depends on numerical quadrature, the gradient
#' is obtained with numDeriv::grad, exactly as that Remark says.
#'
#' Returns NA unless the observed information is positive definite.
delta_se <- function(fit, L, index = c("quantile", "moment")) {
  index <- match.arg(index)
  if (is.null(fit)) return(NA_real_)
  I <- fit$information
  if (!is_posdef(I)) return(NA_real_)
  V <- tryCatch(solve(I), error = function(e) NULL)
  if (is.null(V)) return(NA_real_)
  g <- if (index == "quantile") {
    grad_CL_quantile(fit$par, L)
  } else {
    tryCatch(numDeriv::grad(CL_moment, fit$par, L = L), error = function(e) NULL)
  }
  if (is.null(g) || anyNA(g) || any(!is.finite(g))) return(NA_real_)
  v <- drop(t(g) %*% V %*% g)
  if (!is.finite(v) || v <= 0) return(NA_real_)
  sqrt(v)
}

#' Asymptotic (ACI) confidence interval, Eq. (19).
ci_asymptotic <- function(Chat, se, level = 0.95) {
  z <- qnorm(1 - (1 - level) / 2)
  c(Chat - z * se, Chat + z * se)
}

#' Parametric bootstrap replicates of an index.  Section 4, Steps 2-3.
bootstrap_index <- function(par_hat, m, k, R, L, B,
                            index = c("quantile", "moment")) {
  index <- match.arg(index)
  f <- if (index == "quantile") CL_quantile else CL_moment
  out <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    y  <- generate_pffc(m, k, R, par_hat[1], par_hat[2])
    fb <- fit_mle(y, R, k, start = par_hat)
    if (is.null(fb)) next
    out[b] <- f(fb$par, L)
  }
  out
}

#' Percentile bootstrap interval, Eq. (22).  R's default quantile type (7).
ci_percentile <- function(boot, level = 0.95) {
  a <- (1 - level) / 2
  boot <- boot[is.finite(boot)]
  if (length(boot) < 2) return(c(NA_real_, NA_real_))
  unname(quantile(boot, c(a, 1 - a)))
}

#' Bias-corrected normal-approximation bootstrap interval, Eq. (21).
ci_normal_boot <- function(Chat, boot, level = 0.95) {
  boot <- boot[is.finite(boot)]
  if (length(boot) < 2) return(c(NA_real_, NA_real_))
  z    <- qnorm(1 - (1 - level) / 2)
  bias <- mean(boot) - Chat
  se   <- sd(boot)
  centre <- Chat - bias
  c(centre - z * se, centre + z * se)
}

#' One-sided lower confidence bound used by the decision rule of Eq. (23).
lower_bound_one_sided <- function(Chat, se, level = 0.95) {
  Chat - qnorm(level) * se
}

## ---------------------------------------------------------------------------
## 8.  Housekeeping
## ---------------------------------------------------------------------------

#' Load a package, installing it first if necessary.
need <- function(...) {
  for (p in c(...)) {
    if (!requireNamespace(p, quietly = TRUE))
      install.packages(p, repos = "https://cloud.r-project.org")
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
  invisible(NULL)
}

#' Resolve output directories relative to the code/ directory, so that scripts
#' work whether they are run from code/ or from code/validation/.
CODE_DIR    <- if (file.exists("00_gie_pffc.R")) "." else ".."
RESULTS_DIR <- file.path(CODE_DIR, "results")
TABLES_DIR  <- file.path(CODE_DIR, "..", "tables")
FIGURES_DIR <- file.path(CODE_DIR, "..", "figures")
for (d in c(RESULTS_DIR, TABLES_DIR, FIGURES_DIR))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

#' Master seed.  Every script derives its own seed from this one, so scripts can
#' be run individually and still reproduce the reported numbers.
MASTER_SEED <- 20260822L

## ---------------------------------------------------------------------------
## 9.  Writing numbers back into the manuscript
## ---------------------------------------------------------------------------

#' Write a set of LaTeX macros to a file.
#'
#' Every number that appears in the running text of the manuscript -- not only
#' the ones in tables -- is defined here as a macro and \input by the paper, so
#' that re-running a script updates the prose as well as the tables.
#'
#' `values` is a named list; names must consist of letters only, since that is
#' what \newcommand accepts.  Values are written verbatim, so format them with
#' fmt() or sprintf() before passing them in.
write_macros <- function(file, values, script) {
  bad <- names(values)[grepl("[^A-Za-z]", names(values))]
  if (length(bad))
    stop("write_macros(): macro names must be letters only: ",
         paste(bad, collapse = ", "))
  con <- base::file(file, open = "wt")
  on.exit(close(con))
  writeLines(paste0("%% Generated by ", script, " -- do not edit by hand."), con)
  for (nm in names(values))
    writeLines(sprintf("\\newcommand{\\%s}{%s}", nm, values[[nm]]), con)
  cat("Wrote", file, "\n")
}

#' Read back a macro file written by write_macros().
#'
#' Lets one script use a constant another script computed, instead of the
#' constant being copied by hand into both.  Returns a named list of strings.
read_macros <- function(file) {
  if (!file.exists(file))
    stop("read_macros(): ", file, " not found. Run the script that writes it first.")
  ln <- readLines(file, warn = FALSE)
  m  <- regmatches(ln, regexec("^\\\\newcommand\\{\\\\([A-Za-z]+)\\}\\{(.*)\\}\\s*$", ln))
  m  <- Filter(function(z) length(z) == 3L, m)
  setNames(lapply(m, `[[`, 3L), vapply(m, `[[`, "", 2L))
}

## ---------------------------------------------------------------------------
## 10.  Progress reporting
## ---------------------------------------------------------------------------
##
## The simulations run for hours, and in parallel the workers' own output is not
## shown, so a progress bar is the only signal that anything is happening.  We
## use `progressr`, which is the mechanism `furrr` relays progress through, and
## degrade gracefully to no bar at all if it is not installed.

#' Choose the nicest available progress handler and make sure progress is
#' actually reported.
#'
#' Two traps here.  First, progressr reports progress only in interactive mode
#' unless the option `progressr.enable` is set, so under `Rscript run_all.R`
#' every tick would be signalled and silently discarded -- the bar would simply
#' never appear.  Second, the bar from the `progress` package renders only to a
#' stream it recognises, so it can print nothing when output is redirected to a
#' log file; base R's txtprogressbar always writes.  We therefore use `progress`
#' only when running interactively.
setup_progress <- function() {
  if (!requireNamespace("progressr", quietly = TRUE)) return(FALSE)
  if (is.null(getOption("progressr.enable"))) options(progressr.enable = TRUE)
  if (interactive() && requireNamespace("progress", quietly = TRUE)) {
    progressr::handlers("progress")
  } else {
    progressr::handlers("txtprogressbar")
  }
  TRUE
}

#' Run `fn(p)` with a progress bar of `steps` ticks, where `p` is the tick
#' function.  If progressr is unavailable, `fn(NULL)` is run instead and the
#' tick calls inside it become no-ops.  Returns whatever `fn` returns.
with_progress_bar <- function(steps, fn) {
  if (setup_progress()) {
    progressr::with_progress({
      p <- progressr::progressor(steps = steps)
      fn(p)
    })
  } else {
    fn(NULL)
  }
}

#' Tick helper: call `p` once every `every` iterations.
#'
#' The `every > 0` guard matters.  Callers compute `every` as
#' `N_REP %/% N_TICK`, and the script headers invite reducing `N_REP` for a
#' quick check; if it drops below `N_TICK` then `every` is 0, `i %% 0` is NaN,
#' and the `if` would fail with "missing value where TRUE/FALSE needed" --
#' inside a parallel worker, taking the whole run with it.
tick_if <- function(p, i, every) {
  if (!is.null(p) && every > 0 && i %% every == 0) p()
  invisible(NULL)
}

#' Convenience formatter: fixed number of decimals, as a character string.
fmt <- function(x, d = 4) sprintf(paste0("%.", d, "f"), x)
