###############################################################################
## endpoint_augmented_extremes.R
## Complete companion implementation for:
##   "Endpoint-augmented symbolic likelihoods for extreme-value inference"
##   A. H. Jamaluddin, M. A. Mohd Safari, F. Z. Che Rose, M. J. Mohd Nasir,
##   H. S. Zulkafli
##
## DEPENDENCIES: Base R only (base + parallel + tools + utils). 
##               No compiled dependencies, no tidyverse, no rstan.
##
## HOW TO RUN:
##   1. Standard run: Rscript endpoint_augmented_extremes.R
##   2. Fast smoke test (reduced iterations/configs): 
##      ENDPOINT_SMOKE=1 Rscript endpoint_augmented_extremes.R
##   3. Fallback to synthetic data (if NASA POWER API is unreachable):
##      ENDPOINT_ALLOW_SURROGATE=1 Rscript endpoint_augmented_extremes.R
##
## MAP: MANUSCRIPT -> CODE
##   Def. def:symbol, Eq. eq:counts, Def. def:osbins ... build_symbol()
##   Eq. eq:gev / eq:gpd / eq:endpoints ................ dgev pgev qgev dgpd pgpd qgpd
##   Prop. prop:like, Eq. eq:like ...................... ll_symbolic()
##   Eq. eq:control .................................... ll_control()
##   Sec. sec:priors ................................... log_prior()
##   Alg. alg:symbolic_inference ....................... rwm()
##   Sec. sec:laplace .................................. laplace_full()
##   Eq. eq:rl, Sec. sec:functionals ................... gev_return_level(), summarise_draws()
##   Sec. sec:simulation ............................... simulate_one(), run_simulation(), summarise_simulation()
##   Sec. sec:application .............................. get_power(), thi(), monthly_maxima(), fit_gev_stationary(), moving_window_gev()
##   Sec. sec:software ................................. write_all_tables(), make_figures()
##
## OUTPUTS (written to 'outputs/' and 'figures/' directories):
##   outputs/tab_sim_gev.tex, tab_sim_gpd.tex, tab_endpoint.tex
##   outputs/tab_app_desc.tex, tab_app_gev.tex, tab_app_pot.tex, tab_app_win.tex
##   outputs/tab_sim_misspec.tex (shape-support sensitivity)
##   figures/fig_sim_main.pdf, fig_endpoint.pdf, fig_app_summary.pdf
##   figures/fig_pot.pdf, fig_real_main.pdf, fig_app_trace.pdf
##   outputs/sim_summary.csv, app_summary.csv, windows_*.csv, results.rds
##
## PROVENANCE: The application fetches NASA POWER daily T2M_MAX and RH2M for 
## five sites (1981-01-01 to 2025-12-31). It will ABORT with an error if the 
## API is unavailable, unless ENDPOINT_ALLOW_SURROGATE=1 is set.
###############################################################################

# =============================================================================
# 0. Working directory configuration
#    Defaults to the current working directory. Users can override this by 
#    setting the ENDPOINT_WORKING_DIR environment variable.
# =============================================================================
target_dir <- Sys.getenv("ENDPOINT_WORKING_DIR", unset = getwd())
if (dir.exists(target_dir)) {
  setwd(target_dir)
}

options(stringsAsFactors = FALSE)
options(dplyr.summarise.inform = FALSE)
options(future.globals.maxSize = 1024^3) # 1 GB

## =============================================================================
## 1. Global configuration
## =============================================================================
CONFIG <- list(
  ## ---- reproducibility ------------------------------------------------------
  seed                = 20261014L,
  n_cores             = as.integer(Sys.getenv("ENDPOINT_CORES", unset = max(1L, parallel::detectCores() - 1L))),
  out_dir             = Sys.getenv("ENDPOINT_OUT_DIR", unset = "outputs"),
  fig_dir             = Sys.getenv("ENDPOINT_FIG_DIR", unset = "figures"),
  cache_dir           = file.path(Sys.getenv("ENDPOINT_OUT_DIR", unset = "outputs"), "cache"),
  
  ## ---- model support and priors (Sec. sec:priors) ---------------------------
  xi_lower            = -0.5,     # finite-variance support, lower edge
  xi_upper            =  0.5,     # finite-variance support, upper edge
  prior_m_sd          =  2.0,     # m   ~ N(0, 2^2)
  prior_logscale_sd   =  1.0,     # log omega ~ N(0, 1)
  xi_lower_wide       = -1.0,     # sensitivity support (Remark rem:bound)
  
  ## ---- sampler (Sec. sec:sampler) ------------------------------------------
  n_chains            = 2L,
  n_iter              = 2400L,    # retained iterations per chain (simulation)
  n_burn              =  700L,    # burn-in per chain (simulation)
  n_iter_app          = 6000L,    # retained iterations per chain (application)
  n_burn_app          = 1500L,
  cov_window          = 150L,     # empirical covariance window (states)
  target_acc          = 0.234,    # Roberts-Rosenthal target for adaptive RWM
  cov_inflation       = 2.38^2,   # multiplied by Sigma/d, d = dimension
  rr_cap              = 0.02,
  
  ## ---- simulation grid (Sec. sec:sim-design) --------------------------------
  fams                = c("gev", "gpd"),
  regimes             = c("bounded", "exponential", "heavy"),
  m_list              = c(500L, 5000L, 50000L, 500000L),
  B_list              = c(5L, 10L, 20L, 50L),
  reps                = as.integer(Sys.getenv("ENDPOINT_REPS", unset = 20L)),
  truth_gev           = list(mu = 25, sigma = 5, xi = c(bounded = -0.2, exponential = 0, heavy = 0.2)),
  truth_gpd           = list(sigma = 5, xi = c(bounded = -0.25, exponential = 0, heavy = 0.25)),
  
  ## ---- application (Sec. sec:app-data) -------------------------------------
  app_B               = 10L,      # twelve stored numbers per site
  app_window_chains   = 1L,       # chains per MOVING-WINDOW fit (screening device)
  app_start           = as.Date("1981-01-01"),
  app_end             = as.Date("2025-12-31"),
  app_community       = "AG",
  app_threshold_p     = 0.95,     # site-specific 95th percentile of daily THI
  window_W            = 120L,     # moving-window length in monthly blocks
  window_step         = 6L,       # advance by six blocks
  app_trace_site      = "Melbourne"
)

CONFIG$smoke <- nzchar(Sys.getenv("ENDPOINT_SMOKE"))
CONFIG$allow_surrogate <- identical(Sys.getenv("ENDPOINT_ALLOW_SURROGATE"), "1")

if (CONFIG$smoke) {
  message("*** SMOKE TEST MODE: reduced settings. NOT for reported numbers. ***")
  CONFIG$reps       <- 2L
  CONFIG$m_list     <- c(500L, 50000L)
  CONFIG$B_list     <- c(5L, 20L)
  CONFIG$n_iter     <- 300L
  CONFIG$n_burn     <- 150L
  CONFIG$n_iter_app <- 300L
  CONFIG$n_burn_app <- 150L
  CONFIG$window_step <- 60L
}

## Application sites (Sec. sec:app-data)
LOCATIONS <- data.frame(
  location  = c("Singapore", "Yakutsk", "Death Valley", "Melbourne", "Chicago"),
  latitude  = c(1.35, 62.03, 36.51, -37.81, 41.88),
  longitude = c(103.82, 129.73, -117.30, 144.96, -87.63),
  climate   = c("Equatorial", "Subarctic continental", "Hot arid", "Temperate oceanic", "Humid continental"),
  stringsAsFactors = FALSE
)

## =============================================================================
## 2. Small utilities
## =============================================================================
prog <- function(...) cat(sprintf(...), "\n", file = stderr())

dir_setup <- function() {
  for (d in c(CONFIG$out_dir, CONFIG$fig_dir, CONFIG$cache_dir)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

`%||%` <- function(a, b) if (is.null(a)) b else a

safe_finite <- function(x) x[is.finite(x)]
safe_mean <- function(x, na.rm = TRUE) {
  x <- x[is.finite(x)]; if (length(x) == 0L) NA_real_ else mean(x, na.rm = na.rm)
}
safe_median <- function(x) {
  x <- x[is.finite(x)]; if (length(x) == 0L) NA_real_ else stats::median(x)
}
safe_sd <- function(x) {
  x <- x[is.finite(x)]; if (length(x) < 2L) NA_real_ else stats::sd(x)
}
safe_quant <- function(x, p) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(rep(NA_real_, length(p)))
  as.numeric(stats::quantile(x, p, names = FALSE, type = 7))
}
safe_ratio <- function(a, b) if (!is.finite(a) || !is.finite(b) || b == 0) NA_real_ else a / b

## LaTeX number formatting
fmt <- function(x, d = 3L, na = "--") {
  ifelse(is.finite(x), formatC(x, format = "f", digits = d), na)
}
fmt_int <- function(x, na = "--") {
  ifelse(is.finite(x), formatC(round(x), format = "d"), na)
}
tex_num <- function(x, d = 3L) ifelse(is.finite(x), sprintf(paste0("%.", d, "f"), x), "--")

## =============================================================================
## 3. Extreme-value families (Eq. eq:gev, eq:gpd, eq:endpoints)
## =============================================================================
EPS_XI <- 1e-8

pgev <- function(x, mu = 0, sigma = 1, xi = 0) {
  z <- (x - mu) / sigma
  if (abs(xi) < EPS_XI) return(exp(-exp(-z)))
  t <- 1 + xi * z
  out <- exp(-t^(-1 / xi))
  out[t <= 0] <- if (xi > 0) 0 else 1
  pmin(pmax(out, 0), 1)
}

dgev <- function(x, mu = 0, sigma = 1, xi = 0, log = FALSE) {
  if (!is.finite(sigma) || sigma <= 0) {
    out <- rep(-Inf, length(x)); return(if (log) out else exp(out))
  }
  z <- (x - mu) / sigma
  if (abs(xi) < EPS_XI) {
    lp <- -log(sigma) - z - exp(-z)
  } else {
    t <- 1 + xi * z
    lp <- ifelse(t <= 0, -Inf, -log(sigma) - (1 + 1 / xi) * log(pmax(t, 1e-300)) - t^(-1 / xi))
  }
  if (log) lp else exp(lp)
}

qgev <- function(p, mu = 0, sigma = 1, xi = 0) {
  if (abs(xi) < EPS_XI) return(mu - sigma * log(-log(p)))
  mu + sigma / xi * ((-log(p))^(-xi) - 1)
}

pgpd <- function(x, sigma = 1, xi = 0) {
  x <- as.numeric(x); out <- numeric(length(x)); ok <- x >= 0
  if (!any(ok)) return(out)
  if (abs(xi) < EPS_XI) {
    out[ok] <- 1 - exp(-x[ok] / sigma)
  } else {
    t <- 1 + xi * x[ok] / sigma
    out[ok] <- ifelse(t <= 0, 1, 1 - t^(-1 / xi))
  }
  pmin(pmax(out, 0), 1)
}

dgpd <- function(x, sigma = 1, xi = 0, log = FALSE) {
  if (!is.finite(sigma) || sigma <= 0) {
    out <- rep(-Inf, length(x)); return(if (log) out else exp(out))
  }
  if (abs(xi) < EPS_XI) {
    lp <- ifelse(x < 0, -Inf, -log(sigma) - x / sigma)
  } else {
    t <- 1 + xi * x / sigma
    lp <- ifelse(x < 0 | t <= 0, -Inf,
                 -log(sigma) - (1 + 1 / xi) * log(pmax(t, 1e-300)))
  }
  if (log) lp else exp(lp)
}

qgpd <- function(p, sigma = 1, xi = 0) {
  if (abs(xi) < EPS_XI) return(-sigma * log(1 - p))
  sigma / xi * ((1 - p)^(-xi) - 1)
}

## Fitted upper endpoint (Eq. eq:endpoints); +Inf when the tail is not bounded
gev_endpoint <- function(mu, sigma, xi) ifelse(xi < 0, mu - sigma / xi, Inf)
gpd_endpoint <- function(sigma, xi)     ifelse(xi < 0, -sigma / xi,    Inf)

## T-year return level of the GEV (Eq. eq:rl); Gumbel limit at xi = 0
gev_return_level <- function(mu, sigma, xi, T) {
  ifelse(abs(xi) < EPS_XI,
         mu - sigma * log(-log(1 - 1 / T)),
         mu - sigma / xi * (1 - (-log(1 - 1 / T))^(-xi)))
}

## Sampling: exact inverse-cdf draws
rgev <- function(n, mu = 0, sigma = 1, xi = 0) qgev(stats::runif(n), mu, sigma, xi)
rgpd <- function(n, sigma = 1, xi = 0)     qgpd(stats::runif(n), sigma, xi)

## =============================================================================
## 4. Data-generating truth of the simulation grid (Sec. sec:sim-design)
## =============================================================================
truth_of <- function(fam, regime) {
  if (fam == "gev") {
    list(mu = CONFIG$truth_gev$mu, sigma = CONFIG$truth_gev$sigma,
         xi = unname(CONFIG$truth_gev$xi[regime]))
  } else {
    list(mu = 0, sigma = CONFIG$truth_gpd$sigma,
         xi = unname(CONFIG$truth_gpd$xi[regime]))
  }
}

simulate_series <- function(fam, regime, m) {
  tr <- truth_of(fam, regime)
  if (fam == "gev") rgev(m, tr$mu, tr$sigma, tr$xi) else rgpd(m, tr$sigma, tr$xi)
}
## =============================================================================
## 5. The endpoint-augmented symbolic observation
##    Definition def:symbol (equal-width rule, Eq. eq:counts) and
##    Definition def:osbins (order-statistic rule)
## =============================================================================
## build_symbol() constructs the symbol; the two extrema are REMOVED from the
## counts (they are recorded exactly, so counting them would double-count them).
##   equal : t_b = t_0 + b*Delta,  Delta = (t_B - t_0)/B      (used everywhere)
##   os    : edges at the interior order statistics              (asymptotics only)
build_symbol <- function(y, B, rule = c("equal", "os")) {
  rule <- match.arg(rule)
  y <- as.numeric(y); n <- length(y)
  if (n < 3L) stop("build_symbol() needs at least three observations")
  if (B < 1L) stop("build_symbol() needs B >= 1")
  ys  <- sort(y)
  t0  <- ys[1L]; tB <- ys[n]
  if (rule == "os") {                       # Definition def:osbins (B = n - 2)
    B <- n - 2L
    edges  <- ys                          # edges[1] = y_(1), ..., edges[n] = y_(n)
    counts <- rep(1L, B)
  } else {                                  # Definition def:symbol
    edges <- seq(t0, tB, length.out = B + 1L)
    if (any(diff(edges) <= 0)) {           # degenerate sample (all values equal)
      edges <- seq(t0 - 1e-6, tB + 1e-6, length.out = B + 1L)
    }
    interior <- ys[2L:(n - 1L)]
    idx <- as.integer(cut(interior, breaks = edges, labels = FALSE,
                          right = TRUE, include.lowest = TRUE))
    idx[is.na(idx)] <- B
    idx <- pmax(1L, pmin(B, idx))
    counts <- tabulate(idx, nbins = B)
  }
  list(B = B, rule = rule, n = n, t0 = t0, tB = tB,
       edges = edges, counts = counts,
       interior = if (rule == "os") ys[2L:(n - 1L)] else ys[2L:(n - 1L)],
       stored_numbers = B + 2L, width = tB - t0)
}

## Standardisation (Sec. sec:priors): GEV uses the sample median and IQR/1.349,
## the GPD uses the sample mean excess.  The symbol and the raw data are mapped
## with the SAME constants so that the two likelihoods are comparable.
standardise <- function(y, fam) {
  if (fam == "gev") {
    c0 <- stats::median(y)
    s0 <- stats::IQR(y) / 1.349
    if (!is.finite(s0) || s0 <= 0) s0 <- stats::sd(y)
    if (!is.finite(s0) || s0 <= 0) s0 <- 1
  } else {
    c0 <- 0
    s0 <- mean(y)
    if (!is.finite(s0) || s0 <= 0) s0 <- 1
  }
  list(c0 = c0, s0 = s0)
}

## Standardised copies of the symbol and of the data
standardise_symbol <- function(sym, std) {
  list(B = sym$B, n = sym$n,
       t0   = (sym$t0 - std$c0) / std$s0,
       tB   = (sym$tB - std$c0) / std$s0,
       edges = (sym$edges - std$c0) / std$s0,
       counts = sym$counts, width = sym$width / std$s0)
}

## Parameter vector layout, standardised scale:
##   GEV : par = (m, log omega, xi),  m = (mu - c0)/s0, omega = sigma/s0
##   GPD : par = (log omega, xi),     omega = sigma/s0
n_par <- function(fam) if (fam == "gev") 3L else 2L

unpack_par <- function(par, fam, std) {
  if (fam == "gev") {
    list(mu = par[1L] * std$s0 + std$c0, sigma = exp(par[2L]) * std$s0,
         xi = par[3L], logscale = par[2L], loc = par[1L])
  } else {
    list(mu = std$c0, sigma = exp(par[1L]) * std$s0,
         xi = par[2L], logscale = par[1L], loc = 0)
  }
}

init_par <- function(y, fam, std) {
  if (fam == "gev") {
    c(z = 0, logomega = 0, xi = 0)
  } else {
    c(logomega = 0, xi = 0)
  }
}

## =============================================================================
## 6. Likelihoods
## =============================================================================
## Eq. eq:like :  ell(theta; S) = log f(t_0) + log f(t_B) + sum_b c_b log p_b
## No renormalisation of the interior mass; a non-positive bin probability
## returns -Inf (Sec. sec:priors(c): bin probabilities are never floored).
ll_symbolic <- function(par, fam, syms) {
  if (any(!is.finite(par))) return(-Inf)
  if (fam == "gev") {
    m <- par[1L]; om <- exp(par[2L]); xi <- par[3L]
    d_end <- dgev(c(syms$t0, syms$tB), m, om, xi, log = TRUE)
    Fedges <- pgev(syms$edges, m, om, xi)
  } else {
    om <- exp(par[1L]); xi <- par[2L]
    d_end <- dgpd(c(syms$t0, syms$tB), om, xi, log = TRUE)
    Fedges <- pgpd(syms$edges, om, xi)
  }
  p <- diff(Fedges)
  if (any(!is.finite(p)) || any(p <= 0)) return(-Inf)
  ll <- sum(d_end) + sum(syms$counts * log(p))
  if (!is.finite(ll)) return(-Inf)
  ll
}

log_combinatorial <- function(sym) {
  lfactorial(sym$n) - sum(lfactorial(sym$counts))
}

## Eq. eq:control : endpoint-free grouped likelihood with renormalised interior
ll_control <- function(par, fam, syms) {
  if (any(!is.finite(par))) return(-Inf)
  if (fam == "gev") {
    Fedges <- pgev(syms$edges, par[1L], exp(par[2L]), par[3L])
  } else {
    Fedges <- pgpd(syms$edges, exp(par[1L]), par[2L])
  }
  denom <- Fedges[length(Fedges)] - Fedges[1L]
  if (!is.finite(denom) || denom <= 0) return(-Inf)
  p <- diff(Fedges) / denom
  if (any(!is.finite(p)) || any(p <= 0)) return(-Inf)
  ll <- sum(syms$counts * log(p))
  if (!is.finite(ll)) return(-Inf)
  ll
}

## Full-data log-likelihood on the standardised scale
ll_full <- function(par, fam, zs) {
  if (any(!is.finite(par))) return(-Inf)
  if (fam == "gev") {
    sum(dgev(zs, par[1L], exp(par[2L]), par[3L], log = TRUE))
  } else {
    sum(dgpd(zs, exp(par[1L]), par[2L], log = TRUE))
  }
}

## =============================================================================
## 7. Priors (Sec. sec:priors(b))
## =============================================================================
##   m ~ N(0, 2^2),  log omega ~ N(0, 1),  xi ~ U(xi_lower, xi_upper)
log_prior <- function(par, fam, xi_lower = CONFIG$xi_lower,
                      xi_upper = CONFIG$xi_upper) {
  if (any(!is.finite(par))) return(-Inf)
  if (fam == "gev") {
    xi <- par[3L]
    l <- stats::dnorm(par[1L], 0, CONFIG$prior_m_sd, log = TRUE) +
         stats::dnorm(par[2L], 0, CONFIG$prior_logscale_sd, log = TRUE)
  } else {
    xi <- par[2L]
    l <- stats::dnorm(par[1L], 0, CONFIG$prior_logscale_sd, log = TRUE)
  }
  if (xi <= xi_lower || xi >= xi_upper) return(-Inf)
  l - log(xi_upper - xi_lower)
}

## Penalised full-data criterion (same weak priors, so the two routes target
## the same quantity; Sec. sec:laplace)
ll_full_penalised <- function(par, fam, zs, xi_lower, xi_upper) {
  l <- ll_full(par, fam, zs)
  if (!is.finite(l)) return(-Inf)
  l + log_prior(par, fam, xi_lower, xi_upper)
}

## =============================================================================
## 8. Posterior functionals (Sec. sec:functionals)
## =============================================================================
## summarise_draws() converts the standardised draws back to the original scale
## draw-by-draw, so non-linear functionals (endpoint, return levels) get exact
## posterior summaries without a delta-method step.
summarise_draws <- function(draws, fam, std, T_years = c(10, 100),
                            xi_lower = CONFIG$xi_lower, xi_upper = CONFIG$xi_upper) {
  if (is.null(draws) || nrow(draws) == 0L) return(NULL)
  d <- ncol(draws)
  mu    <- if (fam == "gev") draws[, 1L] * std$s0 + std$c0 else rep(std$c0, nrow(draws))
  sigma <- exp(draws[, d - 1L]) * std$s0
  xi    <- draws[, d]
  ep    <- if (fam == "gev") gev_endpoint(mu, sigma, xi) else gpd_endpoint(sigma, xi)
  rl    <- lapply(T_years, function(TT) gev_return_level(mu, sigma, xi, TT))
  names(rl) <- paste0("z", T_years)
  heavy <- mean(xi >= 0)
  at_lo <- mean(xi <= xi_lower + 1e-3)
  at_hi <- mean(xi >= xi_upper - 1e-3)
  out <- list(
    xi = c(mean = safe_mean(xi), sd = safe_sd(xi),
           q05 = safe_quant(xi, 0.05), q50 = safe_median(xi), q95 = safe_quant(xi, 0.95)),
    endpoint = c(mean = safe_mean(ep[is.finite(ep)]), median = safe_median(ep[is.finite(ep)]),
                 q05 = safe_quant(ep[is.finite(ep)], 0.05), q95 = safe_quant(ep[is.finite(ep)], 0.95)),
    pr_heavy = heavy, pr_at_lower = at_lo, pr_at_upper = at_hi,
    mu = c(mean = safe_mean(mu), sd = safe_sd(mu)),
    sigma = c(mean = safe_mean(sigma), sd = safe_sd(sigma)),
    n_draws = nrow(draws)
  )
  for (nm in names(rl)) out[[nm]] <- c(mean = safe_mean(rl[[nm]]), median = safe_median(rl[[nm]]),
                                       q05 = safe_quant(rl[[nm]], 0.05),
                                       q95 = safe_quant(rl[[nm]], 0.95))
  out
}
## =============================================================================
## 9. Adaptive random-walk Metropolis sampler (Alg. alg:symbolic_inference)
##    After an initial non-adaptive phase the proposal covariance is replaced by
##    the empirical covariance of the most recent CONFIG$cov_window states,
##    scaled by 2.38^2/d (d = dimension); the scalar multiplier is adapted
##    towards CONFIG$target_acc during burn-in only, so the post-burn-in kernel
##    is fixed and the chain targets the stated posterior.
## =============================================================================
rwm <- function(logpost, theta0, n_iter = CONFIG$n_iter, n_burn = CONFIG$n_burn,
                seed = NULL, cov_window = CONFIG$cov_window,
                target_acc = CONFIG$target_acc, init_scale = 1.0,
                adapt_from = 100L, freeze_frac = 1.0) {
  if (!is.null(seed)) set.seed(seed)
  d <- length(theta0)
  total <- n_burn + n_iter
  freeze_at <- max(adapt_from, as.integer(n_burn * freeze_frac))

  theta <- as.numeric(theta0)
  lp <- logpost(theta)
  if (!is.finite(lp)) {                       # jitter an inadmissible start
    for (k in 1:50) {
      cand <- theta + stats::rnorm(d, 0, 0.25)
      lpc <- logpost(cand)
      if (is.finite(lpc)) { theta <- cand; lp <- lpc; break }
    }
    if (!is.finite(lp)) stop("rwm(): initial state has non-finite log posterior")
  }

  chain <- matrix(NA_real_, total, d)
  colnames(chain) <- paste0("th", seq_len(d))
  Sigma <- diag(d) * (CONFIG$cov_inflation / d)
  scale_cur <- init_scale
  hist <- matrix(NA_real_, cov_window, d); hidx <- 0L; hfill <- 0L
  n_acc_burn <- 0L; n_acc_samp <- 0L

  for (it in seq_len(total)) {
    prop <- theta + scale_cur * as.numeric(chol(Sigma + diag(1e-10, d)) %*% stats::rnorm(d))
    lpp <- logpost(prop)
    log_alpha <- lpp - lp
    accept <- is.finite(log_alpha) && log(stats::runif(1L)) < log_alpha
    if (accept) {
      theta <- prop; lp <- lpp
      if (it <= n_burn) n_acc_burn <- n_acc_burn + 1L else n_acc_samp <- n_acc_samp + 1L
    }
    if (it <= freeze_at) {
      hidx <- if (hidx >= cov_window) 1L else hidx + 1L
      hist[hidx, ] <- theta
      hfill <- min(hfill + 1L, cov_window)
      if (hfill >= max(25L, 3L * d) && it %% 10L == 0L) {
        S <- stats::cov(hist[seq_len(hfill), , drop = FALSE])
        if (all(is.finite(S)) && all(diag(S) > 0)) Sigma <- S * (CONFIG$cov_inflation / d)
        if (it > adapt_from) {
          rate <- (n_acc_burn / it)
          scale_cur <- scale_cur * exp(min(CONFIG$rr_cap, 5 / it) * (rate - target_acc))
          scale_cur <- min(max(scale_cur, 1e-6), 1e4)
        }
      }
    }
    chain[it, ] <- theta
  }

  list(draws = chain[(n_burn + 1L):total, , drop = FALSE],
       acceptance = (n_acc_burn + n_acc_samp) / total,
       acceptance_burn = n_acc_burn / max(1L, n_burn),
       acceptance_sampling = n_acc_samp / max(1L, n_iter),
       final_scale = scale_cur)
}

## Convenience wrappers used throughout -----------------------------------------
run_symbolic_fit <- function(y, fam, B, std, sym = NULL, xi_lower = CONFIG$xi_lower,
                             xi_upper = CONFIG$xi_upper, n_iter = CONFIG$n_iter,
                             n_burn = CONFIG$n_burn, n_chains = CONFIG$n_chains,
                             seed = NULL, measure = TRUE) {
  if (is.null(sym)) sym <- build_symbol(y, B)
  syms <- standardise_symbol(sym, std)
  lp <- function(par) { v <- ll_symbolic(par, fam, syms); if (!is.finite(v)) return(-Inf); v + log_prior(par, fam, xi_lower, xi_upper) }
  t0 <- proc.time()[["elapsed"]]
  draws <- list(); acc <- numeric(0); sc <- numeric(0)
  for (c_id in seq_len(n_chains)) {
    s <- if (is.null(seed)) NULL else seed + 1013L * c_id
    r <- rwm(lp, init_par(y, fam, std), n_iter = n_iter, n_burn = n_burn, seed = s)
    draws[[c_id]] <- r$draws; acc <- c(acc, r$acceptance_sampling); sc <- c(sc, r$final_scale)
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  list(fit = "S", fam = fam, B = B, draws = do.call(rbind, draws), chains = draws,
       acceptance = mean(acc), final_scale = mean(sc), elapsed = elapsed,
       sym = sym, syms = syms, std = std, y = y)
}

run_control_fit <- function(y, fam, B, std, sym = NULL, xi_lower = CONFIG$xi_lower,
                            xi_upper = CONFIG$xi_upper, n_iter = CONFIG$n_iter,
                            n_burn = CONFIG$n_burn, n_chains = CONFIG$n_chains,
                            seed = NULL) {
  if (is.null(sym)) sym <- build_symbol(y, B)
  syms <- standardise_symbol(sym, std)
  lp <- function(par) { v <- ll_control(par, fam, syms); if (!is.finite(v)) return(-Inf); v + log_prior(par, fam, xi_lower, xi_upper) }
  t0 <- proc.time()[["elapsed"]]
  draws <- list(); acc <- numeric(0)
  for (c_id in seq_len(n_chains)) {
    s <- if (is.null(seed)) NULL else seed + 2017L * c_id
    r <- rwm(lp, init_par(y, fam, std), n_iter = n_iter, n_burn = n_burn, seed = s)
    draws[[c_id]] <- r$draws; acc <- c(acc, r$acceptance_sampling)
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  list(fit = "EF", fam = fam, B = B, draws = do.call(rbind, draws), chains = draws,
       acceptance = mean(acc), elapsed = elapsed, sym = sym, syms = syms, std = std, y = y)
}

## =============================================================================
## 10. Full-data Laplace baseline (Sec. sec:laplace)
## =============================================================================
num_hessian <- function(fn, par, h = 1e-4) {
  d <- length(par); H <- matrix(0, d, d)
  for (i in seq_len(d)) {
    for (j in seq_len(d)) {
      hi <- h * max(1, abs(par[i])); hj <- h * max(1, abs(par[j]))
      p_pp <- par; p_pp[i] <- p_pp[i] + hi; p_pp[j] <- p_pp[j] + hj
      p_pm <- par; p_pm[i] <- p_pm[i] + hi; p_pm[j] <- p_pm[j] - hj
      p_mp <- par; p_mp[i] <- p_mp[i] - hi; p_mp[j] <- p_mp[j] + hj
      p_mm <- par; p_mm[i] <- p_mm[i] - hi; p_mm[j] <- p_mm[j] - hj
      H[i, j] <- (fn(p_pp) - fn(p_pm) - fn(p_mp) + fn(p_mm)) / (4 * hi * hj)
    }
  }
  H
}

laplace_full <- function(y, fam, std, xi_lower = CONFIG$xi_lower,
                         xi_upper = CONFIG$xi_upper) {
  zs <- (y - std$c0) / std$s0
  d <- n_par(fam)
  obj <- function(par) {
    v <- ll_full_penalised(par, fam, zs, xi_lower, xi_upper)
    if (!is.finite(v)) return(1e10) else -v
  }
  par0 <- init_par(y, fam, std)
  lower <- c(-50, -10, xi_lower + 1e-6); upper <- c(50, 10, xi_upper - 1e-6)
  if (fam == "gpd") { lower <- c(-10, xi_lower + 1e-6); upper <- c(10, xi_upper - 1e-6) }
  t0 <- proc.time()[["elapsed"]]
  op <- try(stats::optim(par0, obj, method = "L-BFGS-B", lower = lower, upper = upper,
                         control = list(maxit = 500, factr = 1e6)), silent = TRUE)
  elapsed <- proc.time()[["elapsed"]] - t0
  if (inherits(op, "try-error") || op$convergence != 0 || !is.finite(op$value)) {
    return(list(fit = "F", fam = fam, par = rep(NA_real_, d), se = rep(NA_real_, d),
                xi = NA_real_, se_xi = NA_real_, mu = NA_real_, sigma = NA_real_,
                endpoint = NA_real_, z10 = NA_real_, z100 = NA_real_,
                at_bound = FALSE, converged = FALSE, elapsed = elapsed, hessian_ok = FALSE))
  }
  par_hat <- op$par
  f <- function(p) { v <- ll_full_penalised(p, fam, zs, xi_lower, xi_upper); if (!is.finite(v)) -1e10 else v }
  H <- num_hessian(f, par_hat)
  se <- rep(NA_real_, d); hess_ok <- FALSE
  ev <- try(eigen(-H, symmetric = TRUE, only.values = TRUE)$values, silent = TRUE)
  if (!inherits(ev, "try-error") && all(is.finite(ev)) && all(ev > 0)) {
    Vinv <- try(solve(-H), silent = TRUE)
    if (!inherits(Vinv, "try-error")) { se <- sqrt(pmax(diag(Vinv), 0)); hess_ok <- TRUE }
  }
  u <- unpack_par(par_hat, fam, std)
  xi_hat <- u$xi
  ep <- if (fam == "gev") gev_endpoint(u$mu, u$sigma, xi_hat) else gpd_endpoint(u$sigma, xi_hat)
  at_bound <- abs(xi_hat - xi_lower) < 1e-3 || abs(xi_hat - xi_upper) < 1e-3
  is_gev <- identical(fam, "gev")
  list(fit = "F", fam = fam, par = par_hat, se = se, xi = xi_hat, se_xi = se[d],
       mu = u$mu, sigma = u$sigma, endpoint = ep,
       z10  = if (is_gev) gev_return_level(u$mu, u$sigma, xi_hat, 10)  else NA_real_,
       z100 = if (is_gev) gev_return_level(u$mu, u$sigma, xi_hat, 100) else NA_real_,
       at_bound = at_bound, converged = TRUE, hessian_ok = hess_ok, elapsed = elapsed)
}

## =============================================================================
## 11. Convergence diagnostics (Sec. sec:sampler)
## =============================================================================
rank_normalise <- function(x) {
  n <- length(x)
  r <- rank(x, ties.method = "average")
  stats::qnorm((r - 0.375) / (n + 0.25))
}

## Rank-normalised split R-hat (Vehtari et al. 2021)
split_rhat <- function(chains) {
  chains <- lapply(chains, function(v) v[is.finite(v)])
  chains <- chains[vapply(chains, length, integer(1L)) >= 4L]
  if (length(chains) < 2L) return(NA_real_)
  z <- rank_normalise(unlist(chains))
  idx <- 1L; parts <- list()
  for (v in chains) {
    h <- floor(length(v) / 2L)
    if (h < 2L) next
    parts[[length(parts) + 1L]] <- z[idx:(idx + h - 1L)]
    parts[[length(parts) + 1L]] <- z[(idx + h):(idx + 2L * h - 1L)]
    idx <- idx + length(v)
  }
  if (length(parts) < 2L) return(NA_real_)
  m <- length(parts); n <- min(vapply(parts, length, integer(1L)))
  parts <- lapply(parts, function(v) v[seq_len(n)])
  W <- mean(vapply(parts, stats::var, numeric(1L)))
  B <- n * stats::var(vapply(parts, mean, numeric(1L)))
  if (!is.finite(W) || W <= 0) return(NA_real_)
  sqrt(((n - 1) / n * W + B / n) / W)
}

## Autocorrelation-based effective sample size, Geyer initial positive sequence
ess_chain <- function(x) {
  x <- x[is.finite(x)]; n <- length(x)
  if (n < 20L) return(NA_real_)
  if (stats::var(x) <= 0) return(NA_real_)
  ac <- as.numeric(stats::acf(x, lag.max = min(n - 1L, 1000L), plot = FALSE)$acf)[-1L]
  if (length(ac) < 2L) return(n)
  pair <- ac[seq(1L, length(ac) - 1L, by = 2L)] + ac[seq(2L, length(ac), by = 2L)]
  first_neg <- which(pair <= 0)[1L]
  if (!is.na(first_neg)) pair <- pair[seq_len(max(1L, first_neg - 1L))]
  tau <- -1 + 2 * sum(pair)
  if (!is.finite(tau) || tau <= 0) tau <- 1
  min(n, max(1, n / tau))
}

ess_bulk_tail <- function(chains) {
  v <- unlist(chains); v <- v[is.finite(v)]
  if (length(v) < 40L) return(c(bulk = NA_real_, tail = NA_real_))
  z <- rank_normalise(v)
  parts <- list(); idx <- 1L
  for (ch in chains) {
    ch <- ch[is.finite(ch)]; if (length(ch) < 20L) next
    q <- stats::quantile(ch, c(0.05, 0.95), names = FALSE, type = 7)
    parts[[length(parts) + 1L]] <- z[idx:(idx + length(ch) - 1L)]
    idx <- idx + length(ch)
    attr(parts[[length(parts)]], "hi") <- ch > q[2]
    attr(parts[[length(parts)]], "lo") <- ch < q[1]
  }
  bulk <- mean(vapply(parts, ess_chain, numeric(1L)), na.rm = TRUE)
  hi <- unlist(lapply(parts, function(p) as.numeric(attr(p, "hi"))))
  lo <- unlist(lapply(parts, function(p) as.numeric(attr(p, "lo"))))
  tail_ess <- min(ess_chain(hi), ess_chain(lo))
  c(bulk = bulk, tail = tail_ess)
}

## Diagnostics for a fitted object
fit_diagnostics <- function(fit) {
  d <- ncol(fit$chains[[1L]])
  xi_col <- d
  list(
    rhat_xi = split_rhat(lapply(fit$chains, function(m) m[, xi_col])),
    ess_xi = mean(vapply(fit$chains, function(m) ess_chain(m[, xi_col]), numeric(1L)), na.rm = TRUE),
    ess_xi_bulk_tail = ess_bulk_tail(lapply(fit$chains, function(m) m[, xi_col])),
    rhat_all = vapply(seq_len(d), function(j)
      split_rhat(lapply(fit$chains, function(m) m[, j])), numeric(1L)),
    acceptance = fit$acceptance
  )
}
## =============================================================================
## 12. Simulation study (Sec. sec:simulation)
## =============================================================================
sim_grid <- function() {
  g <- expand.grid(fam = CONFIG$fams, regime = CONFIG$regimes,
                   m = CONFIG$m_list, B = CONFIG$B_list,
                   stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  g[order(match(g$fam, CONFIG$fams), match(g$regime, CONFIG$regimes), g$m, g$B), ]
}

## Extract the endpoint draws of a fitted object, on the original scale
endpoint_draws <- function(fit, fam, std) {
  dr <- fit$draws; d <- ncol(dr)
  if (fam == "gev") {
    gev_endpoint(dr[, 1L] * std$s0 + std$c0, exp(dr[, d - 1L]) * std$s0, dr[, d])
  } else {
    gpd_endpoint(exp(dr[, 1L]) * std$s0, dr[, 2L])
  }
}

## One replication: F (full-data Laplace), S (symbolic), EF (endpoint-free)
simulate_one <- function(fam, regime, m, B, rep_id, seed = NULL,
                         n_iter = CONFIG$n_iter, n_burn = CONFIG$n_burn,
                         n_chains = CONFIG$n_chains, store_draws = FALSE,
                         xi_lower = CONFIG$xi_lower, xi_upper = CONFIG$xi_upper) {
  if (!is.null(seed)) set.seed(seed)
  tr  <- truth_of(fam, regime)
  y   <- simulate_series(fam, regime, m)
  std <- standardise(y, fam)
  sym <- build_symbol(y, B)

  f_fit <- laplace_full(y, fam, std, xi_lower, xi_upper)
  s_fit <- run_symbolic_fit(y, fam, B, std, sym = sym, xi_lower = xi_lower,
                            xi_upper = xi_upper, n_iter = n_iter, n_burn = n_burn,
                            n_chains = n_chains, seed = (seed %||% 1L) + 11L)
  e_fit <- run_control_fit(y, fam, B, std, sym = sym, xi_lower = xi_lower,
                           xi_upper = xi_upper, n_iter = n_iter, n_burn = n_burn,
                           n_chains = n_chains, seed = (seed %||% 1L) + 29L)

  s_sum <- summarise_draws(s_fit$draws, fam, std, xi_lower = xi_lower, xi_upper = xi_upper)
  e_sum <- summarise_draws(e_fit$draws, fam, std, xi_lower = xi_lower, xi_upper = xi_upper)
  s_diag <- fit_diagnostics(s_fit)

  true_ep <- if (fam == "gev") gev_endpoint(tr$mu, tr$sigma, tr$xi) else gpd_endpoint(tr$sigma, tr$xi)
  rng <- diff(range(y)); fail_cut <- min(y) + 2 * rng

  ## Endpoint error moments are taken over the draws that return a usable finite
  ## endpoint (|y^F| at or below the failure cut).  Draws with xi >= 0 give
  ## y^F = +Inf and draws beyond the failure cut are reported separately, in
  ## ep_fail and ep_pr_inf; mixing them into the mean square would let a single
  ## unbounded draw (of order 1/xi) dominate the statistic.
  ep_stats <- function(ep) {
    fin <- ep[is.finite(ep) & ep <= fail_cut]
    list(median_err = if (length(fin)) safe_median(abs(fin - true_ep)) else NA_real_,
         mse = if (length(fin)) safe_mean((fin - true_ep)^2) else NA_real_,
         rmse = if (length(fin)) sqrt(safe_mean((fin - true_ep)^2)) else NA_real_,
         frac_usable = mean(is.finite(ep) & ep <= fail_cut),
         fail_frac = as.numeric(mean(!is.finite(ep) | ep > fail_cut)),
         pr_inf = as.numeric(mean(!is.finite(ep))))
  }
  s_ep <- ep_stats(endpoint_draws(s_fit, fam, std))
  e_ep <- ep_stats(endpoint_draws(e_fit, fam, std))
  f_ok <- is.finite(f_fit$endpoint) && f_fit$endpoint <= fail_cut
  f_ep <- list(median_err = if (f_ok) abs(f_fit$endpoint - true_ep) else NA_real_,
               mse = if (f_ok) (f_fit$endpoint - true_ep)^2 else NA_real_,
               rmse = if (f_ok) abs(f_fit$endpoint - true_ep) else NA_real_,
               frac_usable = as.numeric(f_ok),
               fail_frac = as.numeric(!f_ok),
               pr_inf = as.numeric(!is.finite(f_fit$endpoint)))

  f_se_ok <- length(f_fit$se_xi) == 1L && is.finite(f_fit$se_xi) && is.finite(f_fit$xi)
  f_ci <- if (f_se_ok) f_fit$xi + c(-1, 1) * stats::qnorm(0.95) * f_fit$se_xi else c(NA_real_, NA_real_)
  s_ci <- c(s_sum$xi[["q05"]], s_sum$xi[["q95"]])
  e_ci <- c(e_sum$xi[["q05"]], e_sum$xi[["q95"]])

  row <- function(method, est, se_used, ci_lo, ci_hi, ep_tab, runtime,
                  rhat = NA_real_, ess = NA_real_, acc = NA_real_) {
    ci_ok <- is.finite(ci_lo) && is.finite(ci_hi)
    data.frame(
      fam = fam, regime = regime, m = m, B = B, rep_id = rep_id, method = method,
      xi_true = tr$xi, xi_est = est, bias = est - tr$xi, sq_err = (est - tr$xi)^2,
      sd_reported = se_used, ci_lo = ci_lo, ci_hi = ci_hi, ci_ok = as.integer(ci_ok),
      cover90 = if (ci_ok) as.integer(ci_lo <= tr$xi && ci_hi >= tr$xi) else NA_integer_,
      width90 = if (ci_ok) ci_hi - ci_lo else NA_real_,
      ep_median_err = ep_tab$median_err, ep_mse = ep_tab$mse, ep_rmse = ep_tab$rmse,
      ep_frac_usable = ep_tab$frac_usable, ep_fail = ep_tab$fail_frac, ep_pr_inf = ep_tab$pr_inf,
      runtime = runtime, rhat_xi = rhat, ess_xi = ess, acceptance = acc,
      stringsAsFactors = FALSE)
  }
  rows <- rbind(
    row("F",  f_fit$xi, f_fit$se_xi, f_ci[1], f_ci[2], f_ep, f_fit$elapsed),
    row("S",  s_sum$xi[["mean"]], s_sum$xi[["sd"]], s_ci[1], s_ci[2], s_ep,
        s_fit$elapsed, s_diag$rhat_xi, s_diag$ess_xi, s_fit$acceptance),
    row("EF", e_sum$xi[["mean"]], e_sum$xi[["sd"]], e_ci[1], e_ci[2], e_ep,
        e_fit$elapsed, NA_real_, NA_real_, e_fit$acceptance)
  )
  out <- list(rows = rows, f_fit = f_fit, s_sum = s_sum, e_sum = e_sum,
              s_diag = s_diag,
              meta = list(fam = fam, regime = regime, m = m, B = B))
  if (store_draws) out$draws <- list(f = f_fit, s = s_fit, ef = e_fit, y = y, sym = sym)
  out
}

## Run the 96-configuration grid (Sec. sec:sim-design)
run_simulation <- function(reps = CONFIG$reps, seed = CONFIG$seed,
                           n_iter = CONFIG$n_iter, n_burn = CONFIG$n_burn,
                           n_chains = CONFIG$n_chains, grid = sim_grid(),
                           n_cores = CONFIG$n_cores, store_example = TRUE,
                           keep_draws_rep = 1L) {
  prog("[SIM] %d configurations x %d replications = %d fits x 3 methods",
       nrow(grid), reps, nrow(grid) * reps)
  jobs <- list(); k <- 0L
  for (i in seq_len(nrow(grid))) for (r in seq_len(reps)) {
    k <- k + 1L
    jobs[[k]] <- list(fam = grid$fam[i], regime = grid$regime[i], m = grid$m[i],
                      B = grid$B[i], rep_id = r, seed = seed + 7919L * i + 31L * r,
                      store = (store_example && r == keep_draws_rep))
  }
  run_one <- function(j) {
    simulate_one(j$fam, j$regime, j$m, j$B, j$rep_id, seed = j$seed,
                 n_iter = n_iter, n_burn = n_burn, n_chains = n_chains,
                 store_draws = j$store)
  }
  t00 <- proc.time()[["elapsed"]]
  if (n_cores > 1L && .Platform$OS.type == "unix") {
    res <- parallel::mclapply(jobs, run_one, mc.cores = min(n_cores, length(jobs)))
  } else {
    res <- lapply(jobs, run_one)
  }
  elapsed <- proc.time()[["elapsed"]] - t00
  detail <- do.call(rbind, lapply(res, function(z) z$rows))
  example <- NULL
  for (z in res) if (!is.null(z$draws)) { example <- z; break }
  prog("[SIM] grid finished in %.1f s", elapsed)
  list(detail = detail, example = example, elapsed = elapsed, grid = grid, reps = reps)
}

## Configuration-level summaries with Monte-Carlo standard errors
summarise_simulation <- function(detail) {
  key <- interaction(detail$fam, detail$regime, detail$m, detail$B, detail$method, drop = TRUE)
  parts <- split(detail, key)
  rows <- lapply(parts, function(d) {
    R <- nrow(d)
    rmse <- sqrt(safe_mean(d$sq_err))
    se_rmse <- if (R > 1L) stats::sd(d$sq_err, na.rm = TRUE) / (2 * rmse * sqrt(R)) else NA_real_
    cov_ok <- d$cover90[!is.na(d$cover90)]
    data.frame(
      fam = d$fam[1L], regime = d$regime[1L], m = d$m[1L], B = d$B[1L], method = d$method[1L],
      reps = R, xi_true = d$xi_true[1L],
      bias = safe_mean(d$bias), se_bias = if (R > 1L) stats::sd(d$bias, na.rm = TRUE) / sqrt(R) else NA_real_,
      rmse = rmse, se_rmse = se_rmse,
      mean_sd = safe_mean(d$sd_reported),
      cover90 = if (length(cov_ok)) mean(cov_ok) else NA_real_,
      frac_ci_available = safe_mean(d$ci_ok),
      width90 = safe_mean(d$width90),
      ep_median_err = safe_mean(d$ep_median_err),
      ep_mse = safe_mean(d$ep_mse),
      ep_mse_se = if (R > 1L) stats::sd(d$ep_mse, na.rm = TRUE) / sqrt(R) else NA_real_,
      ep_rmse = sqrt(safe_mean(d$ep_mse)),
      ep_frac_usable = safe_mean(d$ep_frac_usable),
      ep_fail = safe_mean(d$ep_fail), ep_pr_inf = safe_mean(d$ep_pr_inf),
      runtime = safe_median(d$runtime),
      rhat_xi = safe_mean(d$rhat_xi), ess_xi = safe_mean(d$ess_xi),
      acceptance = safe_mean(d$acceptance),
      stringsAsFactors = FALSE)
  })
  s <- do.call(rbind, rows); rownames(s) <- NULL
  s <- s[order(match(s$fam, CONFIG$fams), match(s$regime, CONFIG$regimes),
               s$m, s$B, match(s$method, c("F", "S", "EF"))), ]

  ## Paired symbolic / full-data RMSE ratio with a delta-method Monte-Carlo s.e.
  s$rmse_ratio_vs_full <- NA_real_; s$rmse_ratio_se <- NA_real_
  for (i in seq_len(nrow(s))) {
    if (s$method[i] != "S") next
    sel <- function(mm) detail[detail$fam == s$fam[i] & detail$regime == s$regime[i] &
                                 detail$m == s$m[i] & detail$B == s$B[i] &
                                 detail$method == mm, ]
    dS <- sel("S"); dF <- sel("F")
    if (nrow(dS) == nrow(dF) && nrow(dS) > 1L) {
      oS <- dS[order(dS$rep_id), ]; oF <- dF[order(dF$rep_id), ]
      a <- oS$sq_err; b <- oF$sq_err
      A <- mean(a); Bm <- mean(b)
      if (is.finite(A) && is.finite(Bm) && Bm > 0) {
        s$rmse_ratio_vs_full[i] <- sqrt(A / Bm)
        g <- c(1 / (2 * sqrt(A * Bm)), -sqrt(A) / (2 * Bm^1.5))
        V <- stats::cov(cbind(a, b))
        s$rmse_ratio_se[i] <- sqrt(max(as.numeric(t(g) %*% V %*% g) / nrow(dS), 0))
      }
    }
  }
  s
}
## =============================================================================
## 13. LaTeX tables (booktabs; filenames and labels must match the .tex)
## =============================================================================
TEX_HEAD <- function() c("% ============================================================",
                         "% GENERATED AUTOMATICALLY by endpoint_augmented_extremes.R",
                         "% DO NOT EDIT BY HAND - re-run the script instead",
                         "% ============================================================")

write_lines_tex <- function(lines, file) {
  writeLines(lines, file.path(CONFIG$out_dir, file))
  prog("[TEX] wrote %s (%d lines)", file, length(lines))
}

METHOD_TEX <- c(F = "Full data", S = "Symbolic", EF = "Endpoint-free")

sub_set <- function(s, fam, regimes) s[s$fam == fam & s$regime %in% regimes, ]

## ---- Tables tab:sim-gev and tab:sim-gpd --------------------------------------
write_sim_table_family <- function(s, fam, regimes, file, label, caption, note) {
  d <- sub_set(s, fam, regimes)
  d <- d[order(d$m, d$B, match(d$method, c("F", "S", "EF"))), ]
  body <- character(0)
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    ratio <- if (r$method == "S" && is.finite(r$rmse_ratio_vs_full))
      sprintf("%.2f (%.2f)", r$rmse_ratio_vs_full, r$rmse_ratio_se) else "--"
    body <- c(body, sprintf(
      "%d & %d & %s & %s & %s & %s & %s & %s & %s \\\\",
      r$m, r$B, METHOD_TEX[[r$method]],
      tex_num(r$bias), tex_num(r$rmse),
      if (is.finite(r$se_rmse)) sprintf("(%.3f)", r$se_rmse) else "(--)",
      if (is.finite(r$cover90)) sprintf("%.3f", r$cover90) else "--",
      tex_num(r$width90), ratio))
  }
  write_lines_tex(c(TEX_HEAD(),
    "\\begin{table}[t]", "\\centering", "\\small",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    "\\begin{tabular}{rrlcccccc}", "\\toprule",
    "$m$ & $B$ & Method & Bias & RMSE & s.e.(RMSE) & $90\\%$ cov. & Width & RMSE ratio \\\\",
    "\\midrule", body, "\\bottomrule", "\\end{tabular}",
    sprintf("\\\\[2pt]\\footnotesize %s", note),
    "\\end{table}"), file)
}

## ---- Table tab:endpoint ------------------------------------------------------
write_endpoint_table <- function(s, file = "tab_endpoint.tex",
                                 label = "tab:endpoint") {
  d <- s[s$regime == "bounded", ]
  d <- d[order(match(d$fam, CONFIG$fams), d$m, d$B, match(d$method, c("F", "S", "EF"))), ]
  body <- character(0)
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    body <- c(body, sprintf("%s & %d & %d & %s & %s & %s & %s & %s & %s \\\\",
      ifelse(r$fam == "gev", "GEV", "GPD"), r$m, r$B, METHOD_TEX[[r$method]],
      tex_num(r$ep_median_err), tex_num(r$ep_mse, 2),
      if (is.finite(r$ep_mse_se)) sprintf("(%.2f)", r$ep_mse_se) else "(--)",
      sprintf("%.3f", r$ep_fail), sprintf("%.3f", r$ep_pr_inf),
      if (is.finite(r$rmse)) tex_num(r$rmse) else "--"))
  }
  write_lines_tex(c(TEX_HEAD(),
    "\\begin{table}[t]", "\\centering", "\\small",
    "\\caption{Finite-endpoint recovery in the bounded regimes of the $96$-configuration grid ($\\xi=-0.2$ for the GEV with $\\mu=25$, $\\sigma=5$, endpoint $y^{F}=50.0$; $\\xi=-0.25$ for the GPD with $\\sigma=5$, endpoint $y^{F}=20.0$). Errors are measured against the true endpoint of the generating model. The median error is the median over draws of $|\\hat y^{F}-y^{F}|$ (for the full-data estimator, which returns a single value, it is the absolute error itself). MSE is the mean squared endpoint error with its Monte-Carlo standard error over replications in parentheses. Failure fraction is the proportion of draws whose fitted endpoint exceeds twice the observed range; $\\Pr(\\xi\\ge0)$ is the posterior probability that the tail is not bounded, for which the endpoint functional is $+\\infty$ and no endpoint error is defined. The RMSE of $\\xi$ is repeated for reference.}",
    sprintf("\\label{%s}", label),
    "\\begin{tabular}{llrlcccccc}", "\\toprule",
    "Family & $m$ & $B$ & Method & Median $|\\hat y^{F}-y^{F}|$ & MSE & s.e.(MSE) & Fail. frac. & $\\Pr(\\xi\\ge0)$ & RMSE$(\\xi)$ \\\\",
    "\\midrule", body, "\\bottomrule", "\\end{tabular}",
    "\\end{table}"), file)
}

## ---- Table tab:app-desc ------------------------------------------------------
write_app_desc_table <- function(desc, file = "tab_app_desc.tex",
                                 label = "tab:app-desc") {
  body <- sprintf("%s & %.2f & %.2f & %s & %s & %s & %s & %s & %s & %s \\\\",
    desc$location, desc$latitude, desc$longitude, desc$climate,
    fmt_int(desc$n_days), tex_num(desc$mean_thi), tex_num(desc$median_thi),
    tex_num(desc$max_thi), tex_num(desc$q95_thi), fmt_int(desc$n_excess))
  write_lines_tex(c(TEX_HEAD(),
    "\\begin{table}[t]", "\\centering", "\\small",
    "\\caption{Application sites and their daily temperature--humidity index (THI) series, NASA POWER daily retrievals for 1 January 1981 to 31 December 2025. $q_{0.95}$ is the site-specific threshold used in the threshold-excess analysis and $n_{\\mathrm{exc}}$ the number of exceedances over it. Every site contributes $540$ monthly blocks to the GEV analysis and $n_{\\mathrm{exc}}$ daily excesses to the GPD analysis.}",
    sprintf("\\label{%s}", label),
    "\\begin{tabular}{lrrlcccccc}", "\\toprule",
    "Site & Lat. & Lon. & Climate & $N$ & Mean & Median & Max. & $q_{0.95}$ & $n_{\\mathrm{exc}}$ \\\\",
    "\\midrule", body, "\\bottomrule", "\\end{tabular}",
    "\\end{table}"), file)
}

## ---- Table tab:app-gev -------------------------------------------------------
write_app_gev_table <- function(app, file = "tab_app_gev.tex", label = "tab:app-gev") {
  rows <- lapply(seq_len(nrow(app)), function(i) {
    a <- app[i, ]
    sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
      a$location,
      sprintf("%.3f (%.3f)", a$f_xi, a$f_se_xi),
      sprintf("%.3f [%.3f, %.3f]", a$s_xi_mean, a$s_xi_lo, a$s_xi_hi),
      sprintf("%.3f [%.3f, %.3f]", a$ef_xi_mean, a$ef_xi_lo, a$ef_xi_hi),
      tex_num(a$f_endpoint, 2), tex_num(a$s_endpoint_mean, 2), tex_num(a$ef_endpoint_mean, 2),
      tex_num(a$f_z100, 2), tex_num(a$s_z100_mean, 2), tex_num(a$ef_z100_mean, 2),
      sprintf("%.3f", a$s_pr_at_lower))
  })
  write_lines_tex(c(TEX_HEAD(),
    "\\begin{table}[t]", "\\centering", "\\scriptsize",
    "\\caption{Monthly-maximum GEV fits at five sites ($540$ monthly blocks per site). F is the full-data penalised maximum likelihood estimate with its asymptotic standard error; S is the endpoint-augmented symbolic posterior at $B=10$ (twelve stored numbers); EF is the endpoint-free grouped likelihood of \\eqref{eq:control} on the same symbol. Entries for S and EF are posterior means with $90\\%$ credible intervals. $z_{100}$ is the monthly-maximum $100$-year level of \\eqref{eq:rl}; $y^{F}$ is the fitted upper endpoint of \\eqref{eq:endpoints}, reported as unavailable (--) when $\\xi\\ge0$ so that no finite endpoint exists. $\\Pr(\\xi\\le-1/2^{+})$ is the fraction of symbolic draws within $10^{-3}$ of the lower shape bound and measures how far the posterior is compressed against the truncation point; a large value must not be read as a sharp posterior.}",
    sprintf("\\label{%s}", label),
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{lcccccccccc}", "\\toprule",
    "Site & $\\hat\\xi_{\\mathrm{F}}$ (s.e.) & $\\bar\\xi_{\\mathrm{S}}$ $[q_{.05},q_{.95}]$ & $\\bar\\xi_{\\mathrm{EF}}$ $[q_{.05},q_{.95}]$ & $\\hat y^{F}_{\\mathrm{F}}$ & $\\bar y^{F}_{\\mathrm{S}}$ & $\\bar y^{F}_{\\mathrm{EF}}$ & $z^{\\mathrm{F}}_{100}$ & $\\bar z^{\\mathrm{S}}_{100}$ & $\\bar z^{\\mathrm{EF}}_{100}$ & $\\Pr(\\xi\\le-1/2^{+})$ \\\\",
    "\\midrule", unlist(rows), "\\bottomrule", "\\end{tabular}}",
    "\\end{table}"), file)
}

## ---- Table tab:app-pot -------------------------------------------------------
write_app_pot_table <- function(pot, file = "tab_app_pot.tex", label = "tab:app-pot") {
  rows <- sprintf("%s & %d & %.1f & %.3f (%.3f) & %.3f [%.3f, %.3f] & %.3f & %s & %.3f \\\\",
    pot$location, pot$n_excess, pot$threshold_above, pot$f_xi, pot$f_se_xi,
    pot$s_xi_mean, pot$s_xi_lo, pot$s_xi_hi,
    abs(pot$s_xi_mean - pot$f_xi), tex_num(pot$s_endpoint_mean, 2), pot$s_pr_infinite)
  write_lines_tex(c(TEX_HEAD(),
    "\\begin{table}[t]", "\\centering", "\\small",
    "\\caption{Threshold-excess (GPD) analysis of the daily THI series. The threshold is the site-specific empirical $95$th percentile of the full daily series; $n_{\\mathrm{exc}}$ is the number of exceedances and $y^{F}_{\\mathrm{S}}$ the posterior mean of the fitted endpoint above the threshold, reported as -- when the symbolic posterior places mass on $\\xi\\ge0$. F is the full-data maximum likelihood estimate with its asymptotic standard error and S the endpoint-augmented symbolic posterior at $B=10$, i.e. twelve stored numbers for the $822$ excesses at each site. The fourth column is the absolute difference between the symbolic posterior mean and the full-data estimate; the fifth column is the symbolic $90\\%$ credible interval.}",
    sprintf("\\label{%s}", label),
    "\\begin{tabular}{lrcccccc}", "\\toprule",
    "Site & $n_{\\mathrm{exc}}$ & Thr. above & $\\hat\\xi_{\\mathrm{F}}$ (s.e.) & $\\bar\\xi_{\\mathrm{S}}$ $[q_{.05},q_{.95}]$ & $|\\bar\\xi_{\\mathrm{S}}-\\hat\\xi_{\\mathrm{F}}|$ & $\\bar y^{F}_{\\mathrm{S}}$ & $\\Pr(\\xi\\ge0)_{\\mathrm{S}}$ \\\\",
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
    "\\end{table}"), file)
}

## ---- Table tab:app-win -------------------------------------------------------
write_app_win_table <- function(win, file = "tab_app_win.tex", label = "tab:app-win") {
  rows <- sprintf("%s & %d & %.3f & %.2f & %s & %s & %s & %.2f \\\\",
    win$location, win$n_windows, win$rho1, win$n_eff,
    tex_num(win$slope_xi), tex_num(win$se_slope_xi),
    tex_num(win$slope_z100), win$frac_S_in_F)
  write_lines_tex(c(TEX_HEAD(),
    "\\begin{table}[t]", "\\centering", "\\small",
    "\\caption{Overlapping moving-window analysis of the monthly-maximum GEV. Each window holds $W=120$ monthly blocks and is advanced by six blocks, giving $71$ windows per site. The trend slope is the least-squares slope of the window-level posterior mean of $\\xi$ (and of $z_{100}$) on the window start date, in units of the response per decade; its standard error uses the effective sample size $n_{\\mathrm{eff}}$ of \\eqref{eq:neff} with the lag-one autocorrelation $\\hat\\rho_1$ of the window-level sequence, and is a screening statistic rather than a test. The last column is the fraction of windows in which the symbolic posterior mean of $\\xi$ lies inside the asymptotic full-data interval, restricted to those windows in which that interval is available.}",
    sprintf("\\label{%s}", label),
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{lrrrrrrr}", "\\toprule",
    "Site & $N_w$ & $\\hat\\rho_1$ & $n_{\\mathrm{eff}}$ & Slope $\\xi$ /decade & s.e. & Slope $z_{100}$ /decade & Frac. S in F \\\\",
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}}",
    "\\end{table}"), file)
}

## ---- Table tab_sim_misspec (Remark rem:bound sensitivity check) --------------
write_misspec_table <- function(ms, file = "tab_sim_misspec.tex",
                                label = "tab:sim-misspec") {
  rows <- sprintf("%s & %s & %.3f & %s & %.3f & %s & %s \\\\",
    ms$location, tex_num(ms$xi_narrow_F), ms$pr_at_bound_narrow_F,
    tex_num(ms$xi_wide_F), ms$pr_at_bound_wide_F,
    tex_num(ms$xi_narrow_S), tex_num(ms$xi_wide_S))
  write_lines_tex(c(TEX_HEAD(),
    "\\begin{table}[t]", "\\centering", "\\small",
    "\\caption{Sensitivity of the monthly-maximum GEV shape estimate to the width of the shape support, the check claimed in Remark~\\ref{rem:bound}. Narrow is the support $(-1/2,1/2)$ used throughout the paper; wide is $(-1,1/2)$. F is the full-data penalised maximum likelihood estimate and S the symbolic posterior mean at $B=10$. The columns headed $\\Pr(\\mathrm{bound})$ give the fraction of symbolic draws within $10^{-3}$ of the lower edge of the corresponding support. Widening the support moves the estimates to the new lower bound rather than freeing them, which is the evidence that the data prefer a very short tail.}",
    sprintf("\\label{%s}", label),
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{lcccccccc}", "\\toprule",
    "Site & $\\hat\\xi_{\\mathrm{F}}^{\\mathrm{narrow}}$ & $\\Pr(\\mathrm{bound})$ & $\\hat\\xi_{\\mathrm{F}}^{\\mathrm{wide}}$ & $\\Pr(\\mathrm{bound})$ & $\\bar\\xi_{\\mathrm{S}}^{\\mathrm{narrow}}$ & $\\bar\\xi_{\\mathrm{S}}^{\\mathrm{wide}}$ \\\\",
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}}",
    "\\end{table}"), file)
}

## =============================================================================
## 14. Figures (base graphics; Okabe-Ito colour-blind safe palette)
## =============================================================================
OI <- c(black = "#000000", orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7")
FAM_COL <- c(gev = OI[["blue"]], gpd = OI[["vermillion"]])
METHOD_COL <- c(F = OI[["black"]], S = OI[["blue"]], EF = OI[["vermillion"]])
REGIME_COL <- c(bounded = OI[["blue"]], exponential = OI[["green"]], heavy = OI[["vermillion"]])

panel_tag <- function(txt, cex = 1.0) {
  mtext(txt, side = 3, line = 0.4, adj = 0, font = 2, cex = cex)
}

open_pdf <- function(file, width = 9, height = 7) {
  grDevices::pdf(file.path(CONFIG$fig_dir, file), width = width, height = height)
  graphics::par(mfrow = c(2, 2), mar = c(4.6, 4.6, 2.6, 1.2), mgp = c(2.7, 0.7, 0),
                las = 1, cex.axis = 0.9, cex.lab = 1.0, cex.main = 1.0,
                col.axis = "grey20", col.lab = "grey10", tcl = -0.3, family = "sans")
}

## ---- Figure fig_sim_main -----------------------------------------------------
make_fig_sim <- function(s) {
  open_pdf("fig_sim_main.pdf", 9.5, 7.5)
  ss <- s[s$method == "S", ]
  ms <- sort(unique(ss$m))

  ## (a) ratio of symbolic to full-data RMSE of xi against m at B = 10 and B = 50
  graphics::plot(NA, xlim = range(ms), ylim = c(0, 3.4), log = "x",
                 xlab = "Sample size m (log scale)", ylab = expression(RMSE[S]/RMSE[F] ~ "for" ~ xi),
                 main = "")
  graphics::abline(h = 1, lty = 3, col = "grey45")
  for (b in c(10L, 50L)) {
    for (rg in CONFIG$regimes) {
      d <- ss[ss$B == b & ss$regime == rg, ]
      if (nrow(d) == 0L) next
      d <- d[order(d$m), ]
      graphics::lines(d$m, d$rmse_ratio_vs_full, type = "o", pch = ifelse(b == 10L, 16, 17),
                      lty = match(rg, CONFIG$regimes), col = REGIME_COL[[rg]], lwd = 1.7, cex = 1.1)
    }
  }
  graphics::legend("topright", legend = c(paste(CONFIG$regimes, "tail"), "B = 10", "B = 50"),
                   lty = c(1, 2, 3, NA, NA), pch = c(NA, NA, NA, 16, 17),
                   col = c(REGIME_COL[CONFIG$regimes], "grey25", "grey25"),
                   bty = "n", cex = 0.8, lwd = 1.7)
  panel_tag("(a)")

  ## (b) coverage of the 90% symbolic interval
  graphics::plot(NA, xlim = range(ms), ylim = c(0.6, 1.0), log = "x",
                 xlab = "Sample size m (log scale)", ylab = "Coverage of 90% interval", main = "")
  graphics::abline(h = 0.90, lty = 3, col = "grey45")
  for (rg in CONFIG$regimes) {
    d <- ss[ss$regime == rg, ]
    d <- stats::aggregate(cover90 ~ m, data = d, FUN = function(v) mean(v, na.rm = TRUE))
    d <- d[order(d$m), ]
    graphics::lines(d$m, d$cover90, type = "o", pch = 16, lwd = 1.7, col = REGIME_COL[[rg]])
  }
  graphics::legend("bottomright", legend = CONFIG$regimes, col = REGIME_COL[CONFIG$regimes],
                   lty = 1, pch = 16, bty = "n", cex = 0.85)
  panel_tag("(b)")

  ## (c) cost of one symbolic fit against m (log-log), cost independent of m
  graphics::plot(NA, xlim = range(ms), ylim = c(0.1, 20), log = "xy",
                 xlab = "Sample size m (log scale)", ylab = "Cost of one fit (s)", main = "")
  for (rg in CONFIG$regimes) {
    d <- ss[ss$regime == rg, ]
    d <- stats::aggregate(runtime ~ m, data = d, FUN = safe_median)
    d <- d[order(d$m), ]
    graphics::lines(d$m, d$runtime, type = "o", pch = 16, lwd = 1.7, col = REGIME_COL[[rg]])
  }
  dF <- stats::aggregate(runtime ~ m, data = s[s$method == "F", ], FUN = safe_median)
  dF <- dF[order(dF$m), ]
  graphics::lines(dF$m, dF$runtime, type = "o", pch = 4, lwd = 1.7, col = "grey20")
  graphics::legend("topleft", legend = c(CONFIG$regimes, "Full-data Laplace"),
                   col = c(REGIME_COL[CONFIG$regimes], "grey20"), lty = 1,
                   pch = c(16, 16, 16, 4), bty = "n", cex = 0.85)
  panel_tag("(c)")

  ## (d) cost against resolution B
  graphics::plot(NA, xlim = range(ss$B), ylim = c(0.1, 20), log = "y",
                 xlab = "Resolution B", ylab = "Cost of one fit (s)", main = "")
  for (fam in CONFIG$fams) {
    d <- ss[ss$fam == fam, ]
    d <- stats::aggregate(runtime ~ B, data = d, FUN = safe_median)
    graphics::lines(d$B, d$runtime, type = "o", pch = 16, lwd = 1.7, col = FAM_COL[[fam]])
  }
  graphics::abline(h = safe_median(s$runtime[s$method == "F"]), lty = 2, col = "grey20")
  graphics::legend("topleft", legend = c("Symbolic, GEV", "Symbolic, GPD",
                                         "Full-data Laplace (median)"),
                   col = c(FAM_COL[["gev"]], FAM_COL[["gpd"]], "grey20"),
                   lty = c(1, 1, 2), pch = c(16, 16, NA), bty = "n", cex = 0.85)
  panel_tag("(d)")
  grDevices::dev.off()
  prog("[FIG] fig_sim_main.pdf")
}

## ---- Figure fig_endpoint -----------------------------------------------------
make_fig_endpoint <- function(s, m_sel = 200L) {
  d <- s[s$m == m_sel & s$regime == "bounded", ]
  open_pdf("fig_endpoint.pdf", 9.5, 7.5)
  Bs <- sort(unique(d$B))

  ## (a) RMSE of the endpoint estimate against resolution
  ymax <- max(c(d$ep_rmse[d$method %in% c("S", "EF")]), na.rm = TRUE)
  graphics::plot(NA, xlim = range(Bs), ylim = c(0, 1.15 * ymax),
                 xlab = "Resolution B", ylab = expression(RMSE(~ hat(y)^F)), main = "")
  for (fam in CONFIG$fams) {
    for (mm in c("S", "EF")) {
      dd <- d[d$fam == fam & d$method == mm, ]
      dd <- dd[order(dd$B), ]
      graphics::lines(dd$B, dd$ep_rmse, type = "o", lwd = 1.7, cex = 1.05,
                      pch = ifelse(mm == "S", 16, 17), lty = ifelse(fam == "gev", 1, 2),
                      col = METHOD_COL[[mm]])
    }
  }
  graphics::legend("topright", legend = c("Symbolic (S), GEV", "Symbolic (S), GPD",
                                          "Endpoint-free (EF), GEV", "Endpoint-free (EF), GPD"),
                   col = METHOD_COL[c("S", "S", "EF", "EF")], lty = c(1, 2, 1, 2),
                   pch = c(16, 16, 17, 17), bty = "n", cex = 0.8, lwd = 1.6)
  panel_tag("(a)")

  ## (b) coverage of the 90% interval for xi
  graphics::plot(NA, xlim = range(Bs), ylim = c(0.5, 1.0),
                 xlab = "Resolution B", ylab = "Coverage of 90% interval for " ~ xi, main = "")
  graphics::abline(h = 0.9, lty = 3, col = "grey45")
  for (fam in CONFIG$fams) for (mm in c("S", "EF")) {
    dd <- d[d$fam == fam & d$method == mm, ]; dd <- dd[order(dd$B), ]
    graphics::lines(dd$B, dd$cover90, type = "o", lwd = 1.7, pch = ifelse(mm == "S", 16, 17),
                    lty = ifelse(fam == "gev", 1, 2), col = METHOD_COL[[mm]])
  }
  panel_tag("(b)")

  ## (c) RMSE of xi against resolution
  ymax2 <- max(d$rmse[d$method %in% c("S", "EF")], na.rm = TRUE)
  graphics::plot(NA, xlim = range(Bs), ylim = c(0, 1.15 * ymax2),
                 xlab = "Resolution B", ylab = expression(RMSE(xi)), main = "")
  for (fam in CONFIG$fams) for (mm in c("S", "EF")) {
    dd <- d[d$fam == fam & d$method == mm, ]; dd <- dd[order(dd$B), ]
    graphics::lines(dd$B, dd$rmse, type = "o", lwd = 1.7, pch = ifelse(mm == "S", 16, 17),
                    lty = ifelse(fam == "gev", 1, 2), col = METHOD_COL[[mm]])
  }
  panel_tag("(c)")

  ## (d) endpoint accuracy on the log scale: the median-error comparison
  med <- s[s$regime == "bounded" & s$method %in% c("S", "EF"), ]
  med <- med[is.finite(med$ep_median_err) & med$ep_median_err > 0, ]
  if (nrow(med) > 0L) {
    graphics::boxplot(ep_median_err ~ method, data = med, log = "y",
                      names = c("Endpoint-free (EF)", "Symbolic (S)"),
                      col = c(METHOD_COL[["EF"]], METHOD_COL[["S"]]),
                      border = "grey20", ylab = "Median endpoint error (log scale)",
                      xlab = "", main = "")
    panel_tag("(d)")
  } else {
    graphics::plot.new(); panel_tag("(d)")
  }
  grDevices::dev.off()
  prog("[FIG] fig_endpoint.pdf")
}
## =============================================================================
## 15. NASA POWER retrieval, THI, monthly maxima and threshold excesses
##     (Sec. sec:app-data).  ABORTS if the API is unavailable, unless the user
##     explicitly sets ENDPOINT_ALLOW_SURROGATE=1.
## =============================================================================
POWER_URL <- function(lat, lon, start, end, community) {
  sprintf(paste0("https://power.larc.nasa.gov/api/temporal/daily/point?parameters=T2M_MAX,RH2M",
                 "&community=%s&longitude=%.4f&latitude=%.4f&start=%s&end=%s&format=CSV&time-standard=UTC"),
          community, lon, lat, format(start, "%Y%m%d"), format(end, "%Y%m%d"))
}

## Daily temperature-humidity index, Eq. eq:thi
thi <- function(tmax, rh) tmax - (0.55 - 0.0055 * rh) * (tmax - 14.5)

get_power <- function(location, lat, lon, start, end,
                      cache_dir = CONFIG$cache_dir,
                      community = CONFIG$app_community, retries = 4L) {
  key <- tolower(gsub("[^a-zA-Z0-9]+", "_", location))
  cache_file <- file.path(cache_dir, paste0(key, "_power_daily.csv"))
  url <- POWER_URL(lat, lon, start, end, community)
  if (file.exists(cache_file)) {
    df <- utils::read.csv(cache_file, stringsAsFactors = FALSE)
    return(list(data = df, url = url, source = "cache", retrieved = NA_character_))
  }
  tmp <- tempfile(fileext = ".csv")
  ok <- FALSE
  for (k in seq_len(retries)) {
    ok <- tryCatch({
      utils::download.file(url, tmp, quiet = TRUE, mode = "wb"); TRUE
    }, error = function(e) { prog("[POWER] attempt %d failed (%s): %s", k, location, conditionMessage(e)); FALSE })
    if (ok) break
    Sys.sleep(2 * k)
  }
  if (!ok || !file.exists(tmp) || file.info(tmp)$size < 500) {
    stop(sprintf("NASA POWER retrieval failed for %s after %d attempts.\n  URL: %s\n  Set ENDPOINT_ALLOW_SURROGATE=1 to fall back to a SYNTHETIC series (results are then NOT NASA POWER).",
                 location, retries, url))
  }
  txt <- readLines(tmp, warn = FALSE)
  i <- grep("^YEAR,", txt)[1L]
  if (is.na(i)) stop("NASA POWER CSV for ", location, " has no YEAR header row")
  df <- utils::read.csv(text = paste(txt[i:length(txt)], collapse = "\n"),
                        stringsAsFactors = FALSE)
  df$date <- as.Date(sprintf("%d-%03d", df$YEAR, df$DOY), format = "%Y-%j")
  df <- df[, c("date", "T2M_MAX", "RH2M")]
  df <- df[is.finite(df$T2M_MAX) & is.finite(df$RH2M) &
             df$T2M_MAX > -900 & df$RH2M > -900, ]
  df <- df[order(df$date), ]
  utils::write.csv(df, cache_file, row.names = FALSE)
  list(data = df, url = url, source = "api",
       retrieved = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
}

## Deterministic synthetic stand-in, used ONLY with ENDPOINT_ALLOW_SURROGATE=1.
## It reproduces the seasonal cycle, the site climatology and a bounded upper
## tail, so that the pipeline can be exercised without network access.
surrogate_power <- function(location, lat, lon, start, end, seed = 4242L) {
  set.seed(seed + abs(round(lat * 100)))
  d <- seq(start, end, by = "day")
  doy <- as.integer(format(d, "%j"))
  warm <- cos(2 * pi * (doy - 200) / 365.25)          # north-hemisphere warm season
  if (lat < 0) warm <- -warm
  amp <- 8 + 0.45 * abs(lat); base <- 30 - 0.32 * abs(lat)
  ar <- as.numeric(stats::arima.sim(list(ar = 0.7), n = length(d), sd = 1.6))
  tmax <- base + amp * (0.55 + 0.45 * warm) + ar
  rh <- pmin(97, pmax(18, 62 - 0.28 * abs(lat) + 22 * warm + stats::rnorm(length(d), 0, 7)))
  data.frame(date = d, T2M_MAX = round(tmax, 2), RH2M = round(rh, 2))
}

## Monthly maxima of a daily series: the block is the calendar month
monthly_maxima <- function(df, value = "THI") {
  v <- df[[value]]
  key <- format(df$date, "%Y-%m")
  agg <- stats::aggregate(v, by = list(block = key), FUN = max)
  names(agg)[2L] <- "value"
  agg$date <- as.Date(paste0(agg$block, "-01"), format = "%Y-%m-%d")
  agg <- agg[order(agg$date), ]
  agg$month <- as.integer(format(agg$date, "%m"))
  rownames(agg) <- NULL
  agg
}

## Overlapping moving windows over a block series (Sec. sec:windows)
make_windows <- function(blocks, W = CONFIG$window_W, step = CONFIG$window_step) {
  n <- length(blocks)
  if (n < W) return(list())
  starts <- seq(1L, n - W + 1L, by = step)
  lapply(starts, function(s) blocks[s:(s + W - 1L)])
}

## =============================================================================
## 16. Stationary GEV fit (F, S, EF) for the application
## =============================================================================
fit_gev_stationary <- function(x, B = CONFIG$app_B, n_iter = CONFIG$n_iter_app,
                               n_burn = CONFIG$n_burn_app, n_chains = CONFIG$n_chains,
                               seed = NULL, xi_lower = CONFIG$xi_lower,
                               xi_upper = CONFIG$xi_upper) {
  std <- standardise(x, "gev"); sym <- build_symbol(x, B)
  F <- laplace_full(x, "gev", std, xi_lower, xi_upper)
  S <- run_symbolic_fit(x, "gev", B, std, sym = sym, xi_lower = xi_lower,
                        xi_upper = xi_upper, n_iter = n_iter, n_burn = n_burn,
                        n_chains = n_chains, seed = (seed %||% 1L) + 7L)
  E <- run_control_fit(x, "gev", B, std, sym = sym, xi_lower = xi_lower,
                       xi_upper = xi_upper, n_iter = n_iter, n_burn = n_burn,
                       n_chains = n_chains, seed = (seed %||% 1L) + 13L)
  list(F = F, S = S, EF = E, std = std, sym = sym,
       s = summarise_draws(S$draws, "gev", std, xi_lower = xi_lower, xi_upper = xi_upper),
       e = summarise_draws(E$draws, "gev", std, xi_lower = xi_lower, xi_upper = xi_upper),
       diag = fit_diagnostics(S))
}

## Stationary GPD fit on threshold excesses
fit_gpd_excesses <- function(exc, B = CONFIG$app_B, n_iter = CONFIG$n_iter_app,
                             n_burn = CONFIG$n_burn_app, n_chains = CONFIG$n_chains,
                             seed = NULL, xi_lower = CONFIG$xi_lower,
                             xi_upper = CONFIG$xi_upper) {
  std <- standardise(exc, "gpd"); sym <- build_symbol(exc, B)
  F <- laplace_full(exc, "gpd", std, xi_lower, xi_upper)
  S <- run_symbolic_fit(exc, "gpd", B, std, sym = sym, xi_lower = xi_lower,
                        xi_upper = xi_upper, n_iter = n_iter, n_burn = n_burn,
                        n_chains = n_chains, seed = (seed %||% 1L) + 17L)
  E <- run_control_fit(exc, "gpd", B, std, sym = sym, xi_lower = xi_lower,
                       xi_upper = xi_upper, n_iter = n_iter, n_burn = n_burn,
                       n_chains = n_chains, seed = (seed %||% 1L) + 23L)
  dr <- S$draws; dd <- ncol(dr)
  gmean <- function(v) safe_mean(v[is.finite(v)])
  list(F = F, S = S, EF = E, std = std, sym = sym,
       s = summarise_draws(S$draws, "gpd", std, xi_lower = xi_lower, xi_upper = xi_upper),
       e = summarise_draws(E$draws, "gpd", std, xi_lower = xi_lower, xi_upper = xi_upper),
       diag = fit_diagnostics(S),
       s_endpoint_mean = gmean(gpd_endpoint(exp(dr[, 1L]) * std$s0, dr[, 2L])))
}

## Effective sample size and trend slope from overlapping window summaries
neff_trend <- function(values, dates) {
  ok <- is.finite(values)
  if (sum(ok) < 4L) return(list(rho1 = NA_real_, n_eff = NA_real_, slope = NA_real_, se = NA_real_))
  v <- values[ok]; t <- as.numeric(dates[ok]) / 365.25
  rho1 <- suppressWarnings(stats::cor(v[-1L], v[-length(v)]))
  Nw <- length(v)
  n_eff <- if (is.finite(rho1) && abs(rho1) < 1) Nw * (1 - rho1) / (1 + rho1) else NA_real_
  fit <- stats::lm(v ~ t)
  slope_decade <- stats::coef(fit)[2L] * 10
  resid_sd <- stats::sd(stats::resid(fit))
  se <- if (is.finite(n_eff) && n_eff > 2) resid_sd / sqrt(n_eff - 2) * 10 else NA_real_
  list(rho1 = rho1, n_eff = n_eff, slope = slope_decade, se = se)
}
## =============================================================================
## 17. The application itself (Sec. sec:application)
## =============================================================================
run_application <- function(loc_table = LOCATIONS, B = CONFIG$app_B,
                            n_iter = CONFIG$n_iter_app, n_burn = CONFIG$n_burn_app,
                            n_chains = CONFIG$n_chains, n_cores = CONFIG$n_cores,
                            seed = CONFIG$seed, wide_support = FALSE) {
  provenance <- list()
  desc_rows <- list(); gev_rows <- list(); pot_rows <- list()
  win_rows <- list(); trace_store <- list(); window_traj <- list()

  for (i in seq_len(nrow(loc_table))) {
    loc <- loc_table[i, ]
    prog("[APP] %s ...", loc$location)
    if (CONFIG$allow_surrogate) {
      df_raw <- tryCatch(get_power(loc$location, loc$latitude, loc$longitude,
                                   CONFIG$app_start, CONFIG$app_end),
                         error = function(e) {
                           prog("[APP] API unavailable (%s): synthetic series used", conditionMessage(e))
                           list(data = surrogate_power(loc$location, loc$latitude, loc$longitude,
                                                       CONFIG$app_start, CONFIG$app_end),
                                url = NA_character_, source = "SYNTHETIC", retrieved = NA_character_)
                         })
    } else {
      df_raw <- get_power(loc$location, loc$latitude, loc$longitude,
                          CONFIG$app_start, CONFIG$app_end)
    }
    d <- df_raw$data
    d$THI <- thi(d$T2M_MAX, d$RH2M)
    d <- d[is.finite(d$THI), ]
    provenance[[loc$location]] <- c(url = df_raw$url %||% "NA", source = df_raw$source,
                                    retrieved = df_raw$retrieved %||% NA_character_,
                                    n_days = as.character(nrow(d)))

    ## ---- Table tab:app-desc -------------------------------------------------
    thr <- unname(stats::quantile(d$THI, CONFIG$app_threshold_p, type = 7))
    exc <- d$THI[d$THI > thr] - thr
    mm  <- monthly_maxima(d, "THI")
    desc_rows[[loc$location]] <- data.frame(
      location = loc$location, latitude = loc$latitude, longitude = loc$longitude,
      climate = loc$climate, n_days = nrow(d),
      mean_thi = mean(d$THI), median_thi = stats::median(d$THI),
      max_thi = max(d$THI), q95_thi = thr, n_excess = length(exc),
      n_blocks = nrow(mm), stringsAsFactors = FALSE)

    ## ---- GEV on monthly maxima (F, S, EF) -----------------------------------
    xb <- mm$value
    g <- fit_gev_stationary(xb, B = B, n_iter = n_iter, n_burn = n_burn,
                            n_chains = n_chains, seed = seed + 101L * i,
                            xi_lower = if (wide_support) CONFIG$xi_lower_wide else CONFIG$xi_lower)
    gev_rows[[loc$location]] <- data.frame(
      location = loc$location, n_blocks = length(xb),
      f_xi = g$F$xi, f_se_xi = g$F$se_xi, f_endpoint = g$F$endpoint,
      f_z100 = g$F$z100, f_ci_lo = g$F$xi - qnorm(.95) * g$F$se_xi,
      f_ci_hi = g$F$xi + qnorm(.95) * g$F$se_xi,
      s_xi_mean = g$s$xi[["mean"]], s_xi_lo = g$s$xi[["q05"]], s_xi_hi = g$s$xi[["q95"]],
      s_endpoint_mean = g$s$endpoint[["mean"]], s_z100_mean = g$s$z100[["mean"]],
      s_z10_mean = g$s$z10[["mean"]],
      s_pr_at_lower = g$s$pr_at_lower, s_pr_at_upper = g$s$pr_at_upper,
      s_pr_heavy = g$s$pr_heavy, rhat_xi = g$diag$rhat_xi, ess_xi = g$diag$ess_xi,
      acceptance = g$S$acceptance,

      ef_xi_mean = g$e$xi[["mean"]], ef_xi_lo = g$e$xi[["q05"]], ef_xi_hi = g$e$xi[["q95"]],
      ef_endpoint_mean = g$e$endpoint[["mean"]], ef_z100_mean = g$e$z100[["mean"]],
      ef_runtime = g$EF$elapsed, s_runtime = g$S$elapsed, f_runtime = g$F$elapsed,
      stored_numbers = g$sym$stored_numbers,
      stringsAsFactors = FALSE)

    if (identical(loc$location, CONFIG$app_trace_site)) {
      trace_store[[loc$location]] <- list(fit = g$S, sum = g$s, diag = g$diag,
                                          xi_lower = CONFIG$xi_lower, xi_upper = CONFIG$xi_upper)
    }

    ## ---- GPD on threshold excesses -----------------------------------------
    p <- fit_gpd_excesses(exc, B = B, n_iter = n_iter, n_burn = n_burn,
                          n_chains = n_chains, seed = seed + 211L * i,
                          xi_lower = if (wide_support) CONFIG$xi_lower_wide else CONFIG$xi_lower)
    pot_rows[[loc$location]] <- data.frame(
      location = loc$location, n_excess = length(exc), threshold_above = 0,
      threshold = thr, f_xi = p$F$xi, f_se_xi = p$F$se_xi,
      s_xi_mean = p$s$xi[["mean"]], s_xi_lo = p$s$xi[["q05"]], s_xi_hi = p$s$xi[["q95"]],
      s_endpoint_mean = p$s_endpoint_mean, s_pr_infinite = p$s$pr_heavy,
      ef_xi_mean = p$e$xi[["mean"]], rhat_xi = p$diag$rhat_xi, ess_xi = p$diag$ess_xi,
      s_runtime = p$S$elapsed, f_runtime = p$F$elapsed, stored_numbers = B + 2L,
      stringsAsFactors = FALSE)

    ## ---- Moving windows (F and S per window) --------------------------------
    wins <- make_windows(xb)
    if (length(wins) > 0L) {
      one_win <- function(wv) {
        st <- standardise(wv, "gev"); sy <- build_symbol(wv, B)
        Fw <- laplace_full(wv, "gev", st)
        Sw <- run_symbolic_fit(wv, "gev", B, st, sym = sy, n_iter = n_iter,
                               n_burn = n_burn, n_chains = n_chains, seed = seed + 5L)
        sw <- summarise_draws(Sw$draws, "gev", st)
        list(xi_F = Fw$xi, xi_F_lo = Fw$xi - qnorm(.95) * Fw$se_xi,
             xi_F_hi = Fw$xi + qnorm(.95) * Fw$se_xi, ci_ok = is.finite(Fw$se_xi),
             xi_S = sw$xi[["median"]], xi_S_mean = sw$xi[["mean"]],
             xi_S_lo = sw$xi[["q05"]], xi_S_hi = sw$xi[["q95"]],
             z100_F = Fw$z100, z100_S = sw$z100[["median"]],
             z100_S_lo = sw$z100[["q05"]], z100_S_hi = sw$z100[["q95"]])
      }
      if (n_cores > 1L && .Platform$OS.type == "unix") {
        wres <- parallel::mclapply(wins, one_win, mc.cores = min(n_cores, length(wins)))
      } else {
        wres <- lapply(wins, one_win)
      }
      traj <- do.call(rbind, lapply(wres, function(z) as.data.frame(z, stringsAsFactors = FALSE)))
      traj$location <- loc$location
      traj$window_id <- seq_len(nrow(traj))
      traj$window_start <- mm$date[seq_along(wins)]
      window_traj[[loc$location]] <- traj
      utils::write.csv(traj, file.path(CONFIG$out_dir,
        sprintf("windows_%s.csv", tolower(gsub("[^a-zA-Z0-9]+", "_", loc$location)))),
        row.names = FALSE)

      tr_xi <- neff_trend(traj$xi_S_mean, traj$window_start)
      tr_z  <- neff_trend(traj$z100_S, traj$window_start)
      ok <- traj$ci_ok & is.finite(traj$xi_F_lo) & is.finite(traj$xi_F_hi)
      frac_in <- if (any(ok)) mean(traj$xi_S_mean[ok] >= traj$xi_F_lo[ok] &
                                     traj$xi_S_mean[ok] <= traj$xi_F_hi[ok]) else NA_real_
      win_rows[[loc$location]] <- data.frame(
        location = loc$location, n_windows = nrow(traj), rho1 = tr_xi$rho1,
        n_eff = tr_xi$n_eff, slope_xi = tr_xi$slope, se_slope_xi = tr_xi$se,
        slope_z100 = tr_z$slope, frac_S_in_F = frac_in, stringsAsFactors = FALSE)
    }
  }

  list(desc = do.call(rbind, desc_rows), gev = do.call(rbind, gev_rows),
       pot = do.call(rbind, pot_rows), win = do.call(rbind, win_rows),
       traj = window_traj, trace = trace_store, provenance = provenance)
}

## =============================================================================
## 18. Application figures
## =============================================================================
make_fig_app_summary <- function(app) {
  g <- app$gev; p <- app$pot
  sites <- g$location
  open_pdf("fig_app_summary.pdf", 9.5, 7.5)

  ## (a) shape parameter by site, monthly-maximum GEV vs daily-excess GPD
  ylim <- range(c(g$s_xi_mean, p$s_xi_mean, g$f_xi, p$f_xi), na.rm = TRUE)
  ylim <- ylim + c(-0.06, 0.06)
  graphics::plot(NA, xlim = c(0.6, length(sites) + 0.4), ylim = ylim, xaxt = "n",
                 xlab = "", ylab = expression(xi), main = "")
  graphics::axis(1, at = seq_along(sites), labels = sites, cex.axis = 0.85)
  graphics::abline(h = 0, lty = 3, col = "grey50")
  graphics::abline(h = c(CONFIG$xi_lower, CONFIG$xi_upper), lty = 2, col = "grey65")
  off <- 0.16
  graphics::points(seq_along(sites) - off, g$s_xi_mean, pch = 16, col = FAM_COL[["gev"]], cex = 1.2)
  graphics::segments(seq_along(sites) - off, g$s_xi_lo, seq_along(sites) - off, g$s_xi_hi,
                     col = FAM_COL[["gev"]], lwd = 2)
  graphics::points(seq_along(sites) + off, p$s_xi_mean, pch = 17, col = FAM_COL[["gpd"]], cex = 1.2)
  graphics::segments(seq_along(sites) + off, p$s_xi_lo, seq_along(sites) + off, p$s_xi_hi,
                     col = FAM_COL[["gpd"]], lwd = 2)
  graphics::legend("bottomright", legend = c("Monthly-maximum GEV (S)", "Daily-excess GPD (S)"),
                   col = FAM_COL[c("gev", "gpd")], pch = c(16, 17), bty = "n", cex = 0.85)
  panel_tag("(a)")

  ## (b) shape posterior: full-data estimate vs symbolic posterior, 90% intervals
  graphics::plot(NA, xlim = c(0.5, length(sites) + 0.5), ylim = ylim, xaxt = "n",
                 xlab = "", ylab = expression(xi), main = "")
  graphics::axis(1, at = seq_along(sites), labels = sites, cex.axis = 0.85)
  graphics::abline(h = 0, lty = 3, col = "grey50")
  graphics::abline(h = c(CONFIG$xi_lower, CONFIG$xi_upper), lty = 2, col = "grey65")
  for (j in seq_along(sites)) {
    graphics::segments(j - 0.18, g$f_ci_lo[j], j - 0.18, g$f_ci_hi[j], col = METHOD_COL[["F"]], lwd = 2)
    graphics::points(j - 0.18, g$f_xi[j], pch = 3, col = METHOD_COL[["F"]], cex = 1.2)
    graphics::segments(j + 0.12, g$s_xi_lo[j], j + 0.12, g$s_xi_hi[j], col = METHOD_COL[["S"]], lwd = 2)
    graphics::points(j + 0.12, g$s_xi_mean[j], pch = 16, col = METHOD_COL[["S"]], cex = 1.2)
  }
  graphics::legend("bottomright", legend = c("Full data (MLE, asymptotic 90%)",
                                             "Symbolic (mean, 90% CI)"),
                   col = METHOD_COL[c("F", "S")], pch = c(3, 16), bty = "n", cex = 0.85)
  panel_tag("(b)")

  ## (c) monthly-maximum 100-year level by site
  graphics::plot(NA, xlim = c(0.5, length(sites) + 0.5),
                 ylim = range(c(g$f_z100, g$s_z100_mean, app$desc$max_thi), na.rm = TRUE, finite = TRUE) * c(0.98, 1.0),
                 xaxt = "n", xlab = "", ylab = expression(z[100] ~ "(monthly-maximum THI)"), main = "")
  graphics::axis(1, at = seq_along(sites), labels = sites, cex.axis = 0.85)
  graphics::abline(h = 1, lty = 3)
  graphics::points(seq_along(sites) - 0.15, g$f_z100, pch = 3, col = METHOD_COL[["F"]], cex = 1.2)
  graphics::points(seq_along(sites) + 0.15, g$s_z100_mean, pch = 16, col = METHOD_COL[["S"]], cex = 1.2)
  graphics::points(seq_along(sites), app$desc$max_thi, pch = 4, col = "grey35", cex = 1.1)
  graphics::legend("topleft", legend = c("Full data", "Symbolic", "Observed maximum"),
                   col = c(METHOD_COL[c("F", "S")], "grey35"), pch = c(3, 16, 4), bty = "n", cex = 0.85)
  panel_tag("(c)")

  ## (d) symbolic posterior width vs stored numbers, all sites
  graphics::plot(NA, xlim = range(app$gev$stored_numbers), ylim = range(c(g$f_ci_hi - g$f_ci_lo,
                 g$s_xi_hi - g$s_xi_lo), na.rm = TRUE),
                 xlab = "Stored numbers per symbol", ylab = "Width of 90% interval for " ~ xi, main = "")
  graphics::points(app$gev$stored_numbers, g$f_ci_hi - g$f_ci_lo, pch = 3, col = METHOD_COL[["F"]], cex = 1.2)
  graphics::points(app$gev$stored_numbers, g$s_xi_hi - g$s_xi_lo, pch = 16, col = METHOD_COL[["S"]], cex = 1.2)
  graphics::text(app$gev$stored_numbers, g$s_xi_hi - g$s_xi_lo, labels = sites, pos = 3, cex = 0.7)
  panel_tag("(d)")
  grDevices::dev.off()
  prog("[FIG] fig_app_summary.pdf")
}

make_fig_pot <- function(app) {
  p <- app$pot
  open_pdf("fig_pot.pdf", 9.5, 7.5)

  ## (a) posterior shape parameter, full data vs symbolic, with 90% intervals
  ylim <- range(c(p$f_xi - qnorm(.95) * p$f_se_xi, p$f_xi + qnorm(.95) * p$f_se_xi,
                  p$s_xi_lo, p$s_xi_hi), na.rm = TRUE, finite = TRUE) + c(-0.03, 0.03)
  graphics::plot(NA, xlim = c(0.5, nrow(p) + 0.5), ylim = ylim, xaxt = "n",
                 xlab = "", ylab = expression(xi), main = "")
  graphics::axis(1, at = seq_len(nrow(p)), labels = p$location, cex.axis = 0.85)
  graphics::abline(h = 0, lty = 3, col = "grey50")
  for (j in seq_len(nrow(p))) {
    graphics::segments(j - 0.18, p$f_xi[j] - qnorm(.95) * p$f_se_xi[j],
                       j - 0.18, p$f_xi[j] + qnorm(.95) * p$f_se_xi[j], col = METHOD_COL[["F"]], lwd = 2)
    graphics::points(j - 0.18, p$f_xi[j], pch = 3, col = METHOD_COL[["F"]], cex = 1.2)
    graphics::segments(j + 0.12, p$s_xi_lo[j], j + 0.12, p$s_xi_hi[j], col = METHOD_COL[["S"]], lwd = 2)
    graphics::points(j + 0.12, p$s_xi_mean[j], pch = 16, col = METHOD_COL[["S"]], cex = 1.2)
  }
  graphics::legend("topright", legend = c("Full data (MLE, asymptotic 90%)", "Symbolic (mean, 90% CI)"),
                   col = METHOD_COL[c("F", "S")], pch = c(3, 16), bty = "n", cex = 0.85)
  panel_tag("(a)")

  ## (b) number of threshold excesses retained at each site
  graphics::barplot(p$n_excess, names.arg = p$location, las = 2, cex.names = 0.8,
                    col = OI[["sky"]], border = "grey25", ylim = c(0, max(p$n_excess) * 1.15),
                    ylab = "Threshold excesses", main = "")
  graphics::abline(h = mean(p$n_excess), lty = 2, col = "grey30")
  panel_tag("(b)")

  ## (c) absolute gap between the symbolic posterior mean and the full-data MLE
  gap <- abs(p$s_xi_mean - p$f_xi)
  graphics::barplot(gap, names.arg = p$location, las = 2, cex.names = 0.8,
                    col = OI[["green"]], border = "grey25",
                    ylim = c(0, max(0.02, max(gap) * 1.2)),
                    ylab = expression("|" * bar(xi)[S] - hat(xi)[F] * "|"), main = "")
  graphics::abline(h = 0.012, lty = 2, col = OI[["vermillion"]])
  graphics::legend("topleft", legend = "0.012 (text claim)", lty = 2,
                   col = OI[["vermillion"]], bty = "n", cex = 0.8)
  panel_tag("(c)")

  ## (d) posterior compression against the shape bounds
  graphics::barplot(rbind(p$s_pr_infinite, app$gev$s_pr_at_lower), beside = TRUE,
                    names.arg = p$location, las = 2, cex.names = 0.8,
                    col = c(OI[["vermillion"]], OI[["purple"]]), border = "grey25",
                    ylab = "Posterior probability", main = "")
  graphics::legend("topleft", legend = c(expression("GPD: " * Pr(xi >= 0)),
                                         expression("GEV: " * Pr(xi <= -0.5 + 10^-3))),
                   fill = c(OI[["vermillion"]], OI[["purple"]]), bty = "n", cex = 0.8)
  panel_tag("(d)")
  grDevices::dev.off()
  prog("[FIG] fig_pot.pdf")
}

make_fig_windows <- function(app, max_sites = 5L) {
  traj <- app$traj
  if (length(traj) == 0L) { prog("[FIG] no window trajectories; skipping fig_real_main.pdf"); return(invisible(NULL)) }
  open_pdf("fig_real_main.pdf", 9.5, 7.5)
  cols <- c(OI[["blue"]], OI[["vermillion"]], OI[["green"]], OI[["purple"]], OI[["orange"]])

  ## (a) posterior median of xi at each window start, with the 90% band
  graphics::plot(NA, xlim = range(unlist(lapply(traj, function(z) range(z$window_start)))),
                 ylim = range(unlist(lapply(traj, function(z) c(z$xi_S_lo, z$xi_S_hi))), na.rm = TRUE),
                 xlab = "Window start", ylab = expression("posterior median of " * xi), main = "")
  graphics::abline(h = 0, lty = 3, col = "grey50")
  for (j in seq_along(traj)) {
    z <- traj[[j]]
    graphics::polygon(c(z$window_start, rev(z$window_start)),
                      c(z$xi_S_lo, rev(z$xi_S_hi)),
                      col = grDevices::adjustcolor(cols[j], 0.12), border = NA)
    graphics::lines(z$window_start, z$xi_S, col = cols[j], lwd = 1.8)
  }
  graphics::legend("topright", legend = names(traj), col = cols[seq_along(traj)],
                   lty = 1, lwd = 1.8, bty = "n", cex = 0.75, ncol = 2)
  panel_tag("(a)")

  ## (b) posterior median of z100 at each window start, with the 90% band
  graphics::plot(NA, xlim = range(unlist(lapply(traj, function(z) range(z$window_start)))),
                 ylim = range(unlist(lapply(traj, function(z) c(z$z100_S_lo, z$z100_S_hi))), na.rm = TRUE),
                 xlab = "Window start", ylab = expression("posterior median of " * z[100]), main = "")
  for (j in seq_along(traj)) {
    z <- traj[[j]]
    graphics::polygon(c(z$window_start, rev(z$window_start)),
                      c(z$z100_S_lo, rev(z$z100_S_hi)),
                      col = grDevices::adjustcolor(cols[j], 0.12), border = NA)
    graphics::lines(z$window_start, z$z100_S, col = cols[j], lwd = 1.8)
  }
  panel_tag("(b)")

  ## (c) window-by-window agreement: symbolic vs full-data xi
  allsz <- do.call(rbind, lapply(traj, function(z) z))
  graphics::plot(allsz$xi_F, allsz$xi_S_mean, pch = 16, cex = 0.7,
                 col = grDevices::adjustcolor(cols[as.integer(factor(allsz$location))], 0.75),
                 xlab = expression(hat(xi)[F] ~ "(full data, window)"),
                 ylab = expression(bar(xi)[S] ~ "(symbolic, window)"), main = "")
  lim <- range(c(allsz$xi_F, allsz$xi_S_mean), na.rm = TRUE)
  graphics::abline(0, 1, lty = 2, col = "grey35"); graphics::abline(v = CONFIG$xi_lower, lty = 3, col = "grey65")
  panel_tag("(c)")

  ## (d) the overlap pathology: lag-one autocorrelation of the window sequence
  rr <- if (nrow(app$win) > 0L) app$win$rho1 else NA_real_
  graphics::barplot(rr, names.arg = app$win$location, las = 2, cex.names = 0.8,
                    col = OI[["yellow"]], border = "grey25", ylim = c(0, 1.05),
                    ylab = expression(hat(rho)[1] ~ "(window-level " * xi * ")"), main = "")
  graphics::abline(h = 0.95, lty = 2, col = OI[["vermillion"]])
  panel_tag("(d)")
  grDevices::dev.off()
  prog("[FIG] fig_real_main.pdf")
}

make_fig_trace <- function(app) {
  if (length(app$trace) == 0L) { prog("[FIG] no trace fit stored; skipping fig_app_trace.pdf"); return(invisible(NULL)) }
  nm <- names(app$trace)[1L]
  tr <- app$trace[[nm]]; fit <- tr$fit
  ch <- fit$chains[[1L]]; d <- ncol(ch)
  xi_col <- d; mu_col <- 1L
  open_pdf("fig_app_trace.pdf", 9.5, 7.5)

  graphics::plot(seq_len(nrow(ch)), ch[, xi_col], type = "l", col = OI[["blue"]], lwd = 0.7,
                 xlab = "Iteration (post burn-in)", ylab = expression(xi),
                 ylim = c(CONFIG$xi_lower, CONFIG$xi_upper), main = "")
  graphics::abline(h = c(CONFIG$xi_lower, CONFIG$xi_upper), lty = 2, col = OI[["vermillion"]])
  panel_tag("(a)")

  graphics::plot(seq_len(nrow(ch)), ch[, mu_col], type = "l", col = OI[["green"]], lwd = 0.7,
                 xlab = "Iteration (post burn-in)", ylab = expression(m), main = "")
  panel_tag("(b)")

  graphics::hist(ch[, xi_col], breaks = 40, col = grDevices::adjustcolor(OI[["blue"]], 0.55),
                 border = "white", probability = TRUE,
                 xlim = c(CONFIG$xi_lower, CONFIG$xi_upper), xlab = expression(xi),
                 main = "")
  graphics::abline(v = c(CONFIG$xi_lower, CONFIG$xi_upper), lty = 2, col = OI[["vermillion"]], lwd = 2)
  graphics::legend("topright", legend = "shape support bounds", lty = 2, col = OI[["vermillion"]],
                   bty = "n", cex = 0.8)
  panel_tag("(c)")

  graphics::hist(ch[, mu_col], breaks = 40, col = grDevices::adjustcolor(OI[["green"]], 0.55),
                 border = "white", probability = TRUE, xlab = expression(m), main = "")
  panel_tag("(d)")
  grDevices::dev.off()
  prog("[FIG] fig_app_trace.pdf")
}
## =============================================================================
## 19. Reproducibility record and provenance (Sec. sec:software)
## =============================================================================
write_reproducibility_record <- function(script_path = "endpoint_augmented_extremes.R") {
  md5 <- tryCatch(as.character(tools::md5sum(script_path)), error = function(e) NA_character_)
  git <- tryCatch({
    h <- system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE)
    if (length(h) == 1L && grepl("^[0-9a-f]{7,40}$", h[1L])) h[1L] else "not-a-git-repository"
  }, error = function(e) "git-unavailable")
  rec <- c(
    "REPRODUCIBILITY RECORD",
    sprintf("timestamp          : %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("R version          : %s", R.version.string),
    sprintf("platform           : %s", R.version$platform),
    sprintf("script MD5         : %s", ifelse(is.na(md5), "script-not-found", md5)),
    sprintf("git commit         : %s", git),
    sprintf("master seed        : %d", CONFIG$seed),
    sprintf("RNG kind           : %s", paste(RNGkind(), collapse = ", ")),
    sprintf("cores              : %d", CONFIG$n_cores),
    sprintf("shape support      : (%.2f, %.2f)", CONFIG$xi_lower, CONFIG$xi_upper),
    sprintf("sampler            : adaptive RWM, window %d, target accept %.3f",
            CONFIG$cov_window, CONFIG$target_acc),
    sprintf("chains x iterations: %d chains of %d retained (application %d)",
            CONFIG$n_chains, CONFIG$n_iter, CONFIG$n_iter_app),
    "",
    "Package versions:",
    paste0("  base R only; no contributed packages are required"),
    "",
    "Session info:",
    utils::capture.output(utils::sessionInfo()))
  writeLines(rec, file.path(CONFIG$out_dir, "reproducibility_record.txt"))
  prog("[OUT] reproducibility_record.txt")
}

write_provenance <- function(sim_out, app) {
  conn <- file(file.path(CONFIG$out_dir, "provenance.txt"), "w")
  writeLines(c(
    "PROVENANCE OF EVERY TABLE AND FIGURE",
    sprintf("generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "Simulation tables (tab_sim_gev, tab_sim_gpd, tab_endpoint) and",
    "fig_sim_main / fig_endpoint come from the SIMULATION study only:",
    sprintf("  grid        : %d configurations", nrow(sim_out$grid)),
    sprintf("  replications: %d per configuration", sim_out$reps),
    sprintf("  total fits  : %d datasets x 3 methods", nrow(sim_out$grid) * sim_out$reps),
    "",
    "Application tables (tab_app_desc, tab_app_gev, tab_app_pot, tab_app_win)",
    "and figures fig_app_summary / fig_pot / fig_real_main / fig_app_trace are",
    "computed from the retrieval recorded below.  No surrogate or simulated",
    "series contributes to them unless the source field says SYNTHETIC.",
    ""), conn)
  for (nm in names(app$provenance)) {
    p <- app$provenance[[nm]]
    writeLines(sprintf("%-14s source=%-9s days=%-6s retrieved=%s\n               url=%s",
                       nm, p[["source"]], p[["n_days"]],
                       p[["retrieved"]] %||% "NA", p[["url"]] %||% "NA"), conn)
  }
  writeLines(c("", "If every source above is 'api', every application table is",
               "reproducible from the NASA POWER daily point API and the request",
               "URLs printed above.  The observation count alone does not certify",
               "provenance; the source field does."), conn)
  close(conn)
  prog("[OUT] provenance.txt")
}

## =============================================================================
## 20. Master driver
## =============================================================================
run_all_analysis <- function(mode = c("main", "pilot"), do_sim = TRUE, do_app = TRUE,
                             reps = CONFIG$reps, seed = CONFIG$seed,
                             n_iter = CONFIG$n_iter, n_burn = CONFIG$n_burn,
                             n_iter_app = CONFIG$n_iter_app, n_burn_app = CONFIG$n_burn_app,
                             n_cores = CONFIG$n_cores, script_path = "endpoint_augmented_extremes.R") {
  mode <- match.arg(mode)
  dir_setup()
  t_start <- proc.time()[["elapsed"]]
  prog("================================================================")
  prog(" %s run started %s", toupper(mode), format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  prog(" outputs -> %s | figures -> %s", normalizePath(CONFIG$out_dir),
       normalizePath(CONFIG$fig_dir))
  prog("================================================================")

  writeLines(capture.output(str(CONFIG)), file.path(CONFIG$out_dir, "run_config.txt"))
  write_reproducibility_record(script_path)

  sim_out <- NULL
  if (do_sim) {
    sim_out <- run_simulation(reps = reps, seed = seed, n_iter = n_iter,
                              n_burn = n_burn, n_chains = CONFIG$n_chains,
                              n_cores = n_cores)
    s <- summarise_simulation(sim_out$detail)
    utils::write.csv(sim_out$detail, file.path(CONFIG$out_dir, "sim_detail.csv"), row.names = FALSE)
    utils::write.csv(s, file.path(CONFIG$out_dir, "sim_summary.csv"), row.names = FALSE)

    write_sim_table_family(s, "gev", c("bounded", "exponential"),
      "tab_sim_gev.tex", "tab:sim-gev",
      paste0("GEV results over the bounded ($\\xi=-0.2$) and exponential ($\\xi=0$) regimes. ",
             "$m$ is the number of blocks, $B$ the interior resolution. F is the full-data penalised ",
             "maximum likelihood estimate with its asymptotic standard error, S the endpoint-augmented ",
             "symbolic posterior of \\eqref{eq:like}, EF the endpoint-free grouped posterior of ",
             "\\eqref{eq:control}. Coverage is of the nominal $90\\%$ interval and is undefined (--) ",
             "where the full-data asymptotic standard error does not exist because the estimate is pinned ",
             "at the shape bound. Numbers in parentheses are Monte-Carlo standard errors over ",
             "replications. The RMSE ratio column (S only) is the paired ratio of symbolic to full-data ",
             "RMSE of $\\xi$ with its delta-method Monte-Carlo standard error."),
      sprintf("Bias and RMSE are of the shape parameter $\\xi$; true value %s. %d replications per configuration.",
              "$\\xi=-0.2$ (bounded) or $\\xi=0$ (exponential)", sim_out$reps))

    write_sim_table_family(s, "gpd", c("bounded", "exponential"),
      "tab_sim_gpd.tex", "tab:sim-gpd",
      paste0("GPD results over the bounded ($\\xi=-0.25$) and exponential ($\\xi=0$) regimes; ",
             "$m$ is the number of threshold excesses. Columns are as in Table~\\ref{tab:sim-gev}."),
      sprintf("True shape $\\xi=-0.25$ (bounded) or $\\xi=0$ (exponential); scale $\\sigma=5$ in both. %d replications per configuration.",
              sim_out$reps))

    write_endpoint_table(s)
    make_fig_sim(s)
    make_fig_endpoint(s)
  }

  app <- NULL
  if (do_app) {
    app <- run_application(B = CONFIG$app_B, n_iter = n_iter_app, n_burn = n_burn_app,
                           n_cores = n_cores, seed = seed)
    if (!is.null(app$desc) && nrow(app$desc) > 0L) write_app_desc_table(app$desc)
    if (!is.null(app$gev) && nrow(app$gev) > 0L) write_app_gev_table(app$gev)
    if (!is.null(app$pot) && nrow(app$pot) > 0L) write_app_pot_table(app$pot)
    if (!is.null(app$win) && nrow(app$win) > 0L) write_app_win_table(app$win)
    make_fig_app_summary(app)
    make_fig_pot(app)
    make_fig_windows(app)
    make_fig_trace(app)

    ## application summary csv, as promised in Sec. sec:software
    utils::write.csv(app$gev, file.path(CONFIG$out_dir, "app_summary.csv"), row.names = FALSE)
    utils::write.csv(app$pot, file.path(CONFIG$out_dir, "app_pot_summary.csv"), row.names = FALSE)
    utils::write.csv(app$desc, file.path(CONFIG$out_dir, "app_sites.csv"), row.names = FALSE)

    ## ---- Remark rem:bound sensitivity check: shape support (-1, 1/2) --------
    prog("[APP] shape-support sensitivity check (Remark rem:bound) ...")
    wide <- run_application(B = CONFIG$app_B, n_iter = n_iter_app, n_burn = n_burn_app,
                            n_cores = n_cores, seed = seed, wide_support = TRUE)
    ms <- data.frame(
      location = app$gev$location,
      xi_narrow_F = app$gev$f_xi,
      pr_at_bound_narrow_F = as.numeric(abs(app$gev$f_xi - CONFIG$xi_lower) < 1e-3),
      xi_wide_F = wide$gev$f_xi,
      pr_at_bound_wide_F = as.numeric(abs(wide$gev$f_xi - CONFIG$xi_lower_wide) < 1e-3),
      xi_narrow_S = app$gev$s_xi_mean,
      xi_wide_S = wide$gev$s_xi_mean,
      stringsAsFactors = FALSE)
    utils::write.csv(ms, file.path(CONFIG$out_dir, "sim_misspec.csv"), row.names = FALSE)
    write_misspec_table(ms)

    if (!is.null(sim_out)) write_provenance(sim_out, app)
    saveRDS(list(config = CONFIG, seed = seed, sim = sim_out, app = app,
                 wide = wide, misspec = ms, elapsed = proc.time()[["elapsed"]] - t_start),
            file.path(CONFIG$out_dir, "results.rds"))
    prog("[OUT] results.rds")
  }

  prog("================================================================")
  prog(" finished in %.1f s", proc.time()[["elapsed"]] - t_start)
  prog("================================================================")
  invisible(list(sim = sim_out, app = app))
}

## =============================================================================
## 21. Overrides: hardened sampler, corrected window indexing, richer
##     sensitivity table, and the master driver.
##     (These definitions supersede the earlier ones; the file is sourced after
##     them, so R resolves every call to the versions below.)
## =============================================================================

## ---------------------------------------------------------------------------
## 21a. Adaptive random-walk Metropolis.
##   Phase 1 (first 80 iterations or until a usable covariance exists): diagonal
##   proposal with fixed relative scales init_sd, which cannot be wildly
##   mis-scaled.  Phase 2: proposal covariance = (2.38^2/d) * Sigma_hat, the
##   empirical covariance of the most recent CONFIG$cov_window states, refreshed
##   every 10 iterations during burn-in.  A scalar multiplier is adapted by
##   Robbins-Monro towards CONFIG$target_acc on the windowed acceptance rate and
##   frozen at the end of burn-in, so the retained chain targets the stated
##   posterior.  A non-finite log-posterior is rejected (never clamped).
## ---------------------------------------------------------------------------
rwm <- function(logpost, theta0, n_iter, n_burn, seed = NULL, init_sd = NULL,
                cov_window = CONFIG$cov_window, target_acc = CONFIG$target_acc,
                init_scale = 1.0, adapt_from = 40L, window_acc = 50L) {
  if (!is.null(seed)) set.seed(seed)
  d <- length(theta0)
  if (is.null(init_sd)) init_sd <- rep(0.2, d)
  total <- n_burn + n_iter

  theta <- as.numeric(theta0); lp <- logpost(theta)
  if (!is.finite(lp)) {
    for (k in 1:100) {
      cand <- theta + stats::rnorm(d, 0, 0.25 * init_sd)
      lpc <- logpost(cand)
      if (is.finite(lpc)) { theta <- cand; lp <- lpc; break }
    }
    if (!is.finite(lp)) stop("rwm(): initial state has non-finite log-posterior")
  }

  chain <- matrix(NA_real_, total, d)
  colnames(chain) <- paste0("th", seq_len(d))
  Sigma <- diag(CONFIG$cov_inflation / d, d)
  has_cov <- FALSE
  sc <- init_scale
  hist <- matrix(NA_real_, cov_window, d); hidx <- 0L; hfill <- 0L
  accbuf <- rep(NA_integer_, window_acc); abidx <- 0L
  n_acc_burn <- 0L; n_acc_samp <- 0L

  for (it in seq_len(total)) {
    if (!has_cov || it <= 80L) {
      prop <- theta + sc * init_sd * stats::rnorm(d)
    } else {
      prop <- theta + sc * as.numeric(chol(Sigma + diag(1e-10, d)) %*% stats::rnorm(d))
    }
    lpp <- logpost(prop)
    log_alpha <- lpp - lp
    accept <- is.finite(log_alpha) && log(stats::runif(1L)) < log_alpha
    if (accept) {
      theta <- prop; lp <- lpp
      if (it <= n_burn) n_acc_burn <- n_acc_burn + 1L else n_acc_samp <- n_acc_samp + 1L
    }

    abidx <- if (abidx >= window_acc) 1L else abidx + 1L
    accbuf[abidx] <- as.integer(accept)

    if (it <= n_burn) {
      hidx <- if (hidx >= cov_window) 1L else hidx + 1L
      hist[hidx, ] <- theta
      hfill <- min(hfill + 1L, cov_window)
      if (it %% 10L == 0L && hfill >= max(30L, 5L * d)) {
        S <- stats::cov(hist[seq_len(hfill), , drop = FALSE])
        if (all(is.finite(S))) {
          ev <- eigen(S, symmetric = TRUE)
          ev$values <- pmax(ev$values, 1e-10)
          S <- ev$vectors %*% diag(ev$values, d) %*% t(ev$vectors)
          Sigma <- S * (CONFIG$cov_inflation / d)
          if (!has_cov) { has_cov <- TRUE; sc <- 1.0 }
        }
      }
      if (has_cov && it >= adapt_from && it %% 10L == 0L && sum(!is.na(accbuf)) >= 20L) {
        rate <- mean(accbuf, na.rm = TRUE)
        gamma_n <- min(0.5, max(0.05, 20 / it))
        sc <- min(max(sc * exp(gamma_n * (rate - target_acc)), 1e-4), 1e3)
      }
    }
    chain[it, ] <- theta
  }

  list(draws = chain[(n_burn + 1L):total, , drop = FALSE],
       acceptance = (n_acc_burn + n_acc_samp) / total,
       acceptance_burn = n_acc_burn / max(1L, n_burn),
       acceptance_sampling = n_acc_samp / max(1L, n_iter),
       final_scale = sc)
}

init_sd_of <- function(fam) if (fam == "gev") c(0.35, 0.25, 0.18) else c(0.25, 0.18)

## ---------------------------------------------------------------------------
## 21b. Fit wrappers: dispersed starts across chains (better R-hat) and the
##      hardened sampler.
## ---------------------------------------------------------------------------
run_symbolic_fit <- function(y, fam, B, std, sym = NULL, xi_lower = CONFIG$xi_lower,
                             xi_upper = CONFIG$xi_upper, n_iter = CONFIG$n_iter,
                             n_burn = CONFIG$n_burn, n_chains = CONFIG$n_chains,
                             seed = NULL, measure = TRUE) {
  if (is.null(sym)) sym <- build_symbol(y, B)
  syms <- standardise_symbol(sym, std)
  lp <- function(par) {
    v <- ll_symbolic(par, fam, syms)
    if (!is.finite(v)) return(-Inf)
    v + log_prior(par, fam, xi_lower, xi_upper)
  }
  t0 <- proc.time()[["elapsed"]]
  draws <- list(); acc <- numeric(0); sc <- numeric(0); sd0 <- init_sd_of(fam)
  for (c_id in seq_len(n_chains)) {
    th0 <- init_par(y, fam, std) + (c_id - 1L) * 0.6 * sd0[seq_len(length(init_par(y, fam, std)))]
    s <- if (is.null(seed)) NULL else seed + 1013L * c_id
    r <- rwm(lp, th0, n_iter = n_iter, n_burn = n_burn, seed = s, init_sd = sd0)
    draws[[c_id]] <- r$draws; acc <- c(acc, r$acceptance_sampling); sc <- c(sc, r$final_scale)
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  list(fit = "S", fam = fam, B = B, draws = do.call(rbind, draws), chains = draws,
       acceptance = mean(acc), final_scale = mean(sc), elapsed = elapsed,
       sym = sym, syms = syms, std = std, y = y)
}

run_control_fit <- function(y, fam, B, std, sym = NULL, xi_lower = CONFIG$xi_lower,
                            xi_upper = CONFIG$xi_upper, n_iter = CONFIG$n_iter,
                            n_burn = CONFIG$n_burn, n_chains = CONFIG$n_chains,
                            seed = NULL) {
  if (is.null(sym)) sym <- build_symbol(y, B)
  syms <- standardise_symbol(sym, std)
  lp <- function(par) {
    v <- ll_control(par, fam, syms)
    if (!is.finite(v)) return(-Inf)
    v + log_prior(par, fam, xi_lower, xi_upper)
  }
  t0 <- proc.time()[["elapsed"]]
  draws <- list(); acc <- numeric(0); sd0 <- init_sd_of(fam)
  for (c_id in seq_len(n_chains)) {
    th0 <- init_par(y, fam, std) + (c_id - 1L) * 0.6 * sd0[seq_len(length(init_par(y, fam, std)))]
    s <- if (is.null(seed)) NULL else seed + 2017L * c_id
    r <- rwm(lp, th0, n_iter = n_iter, n_burn = n_burn, seed = s, init_sd = sd0)
    draws[[c_id]] <- r$draws; acc <- c(acc, r$acceptance_sampling)
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  list(fit = "EF", fam = fam, B = B, draws = do.call(rbind, draws), chains = draws,
       acceptance = mean(acc), elapsed = elapsed, sym = sym, syms = syms, std = std, y = y)
}

## ---------------------------------------------------------------------------
## 21c. Windows carry their block indices, so the x axis of fig_real_main is the
##      true window start date (the first version used the first 71 months).
## ---------------------------------------------------------------------------
make_windows <- function(blocks, W = CONFIG$window_W, step = CONFIG$window_step) {
  n <- length(blocks)
  if (n < W) return(list())
  starts <- seq(1L, n - W + 1L, by = step)
  list(starts = starts, mats = lapply(starts, function(s) blocks[s:(s + W - 1L)]))
}

## ---------------------------------------------------------------------------
## 21d. The application, with corrected window indexing and full diagnostics.
## ---------------------------------------------------------------------------
run_application <- function(loc_table = LOCATIONS, B = CONFIG$app_B,
                            n_iter = CONFIG$n_iter_app, n_burn = CONFIG$n_burn_app,
                            n_chains = CONFIG$n_chains, n_cores = CONFIG$n_cores,
                            seed = CONFIG$seed, wide_support = FALSE,
                            do_pot = TRUE, do_windows = TRUE) {
  provenance <- list(); desc_rows <- list(); gev_rows <- list(); pot_rows <- list()
  win_rows <- list(); trace_store <- list(); window_traj <- list()
  xl <- if (wide_support) CONFIG$xi_lower_wide else CONFIG$xi_lower

  for (i in seq_len(nrow(loc_table))) {
    loc <- loc_table[i, ]
    prog("[APP] %s ...", loc$location)
    if (CONFIG$allow_surrogate) {
      df_raw <- tryCatch(get_power(loc$location, loc$latitude, loc$longitude,
                                   CONFIG$app_start, CONFIG$app_end),
                         error = function(e) {
                           prog("[APP] API unavailable (%s): SYNTHETIC series used", conditionMessage(e))
                           list(data = surrogate_power(loc$location, loc$latitude, loc$longitude,
                                                       CONFIG$app_start, CONFIG$app_end),
                                url = NA_character_, source = "SYNTHETIC", retrieved = NA_character_)
                         })
    } else {
      df_raw <- get_power(loc$location, loc$latitude, loc$longitude,
                          CONFIG$app_start, CONFIG$app_end)
    }
    d <- df_raw$data
    d$THI <- thi(d$T2M_MAX, d$RH2M)
    d <- d[is.finite(d$THI), ]
    provenance[[loc$location]] <- c(url = df_raw$url %||% "NA", source = df_raw$source,
                                    retrieved = df_raw$retrieved %||% NA_character_,
                                    n_days = as.character(nrow(d)))

    thr <- unname(stats::quantile(d$THI, CONFIG$app_threshold_p, type = 7))
    exc <- d$THI[d$THI > thr] - thr
    mm  <- monthly_maxima(d, "THI")
    desc_rows[[loc$location]] <- data.frame(
      location = loc$location, latitude = loc$latitude, longitude = loc$longitude,
      climate = loc$climate, n_days = nrow(d),
      mean_thi = mean(d$THI), median_thi = stats::median(d$THI),
      max_thi = max(d$THI), q95_thi = thr, n_excess = length(exc),
      n_blocks = nrow(mm), stringsAsFactors = FALSE)

    ## ---- monthly-maximum GEV: F, S, EF --------------------------------------
    xb <- mm$value
    g <- fit_gev_stationary(xb, B = B, n_iter = n_iter, n_burn = n_burn,
                            n_chains = n_chains, seed = seed + 101L * i, xi_lower = xl)
    gev_rows[[loc$location]] <- data.frame(
      location = loc$location, n_blocks = length(xb),
      f_xi = g$F$xi, f_se_xi = g$F$se_xi, f_endpoint = g$F$endpoint, f_z100 = g$F$z100,
      f_ci_lo = if (is.finite(g$F$se_xi)) g$F$xi - stats::qnorm(0.95) * g$F$se_xi else NA_real_,
      f_ci_hi = if (is.finite(g$F$se_xi)) g$F$xi + stats::qnorm(0.95) * g$F$se_xi else NA_real_,
      s_xi_mean = g$s$xi[["mean"]], s_xi_lo = g$s$xi[["q05"]], s_xi_hi = g$s$xi[["q95"]],
      s_endpoint_mean = g$s$endpoint[["mean"]], s_z100_mean = g$s$z100[["mean"]],
      s_z10_mean = g$s$z10[["mean"]],
      s_pr_at_lower = g$s$pr_at_lower, s_pr_at_upper = g$s$pr_at_upper,
      s_pr_heavy = g$s$pr_heavy, rhat_xi = g$diag$rhat_xi, ess_xi = g$diag$ess_xi,
      acceptance = g$S$acceptance,
      ef_xi_mean = g$e$xi[["mean"]], ef_xi_lo = g$e$xi[["q05"]], ef_xi_hi = g$e$xi[["q95"]],
      ef_endpoint_mean = g$e$endpoint[["mean"]], ef_z100_mean = g$e$z100[["mean"]],
      ef_runtime = g$EF$elapsed, s_runtime = g$S$elapsed, f_runtime = g$F$elapsed,
      stored_numbers = g$sym$stored_numbers, wide_support = wide_support,
      stringsAsFactors = FALSE)
    if (identical(loc$location, CONFIG$app_trace_site)) {
      trace_store[[loc$location]] <- list(fit = g$S, sum = g$s, diag = g$diag,
                                          xi_lower = CONFIG$xi_lower, xi_upper = CONFIG$xi_upper)
    }

    if (do_pot) {
    ## ---- daily threshold excesses: F, S, EF ---------------------------------
    p <- fit_gpd_excesses(exc, B = B, n_iter = n_iter, n_burn = n_burn,
                          n_chains = n_chains, seed = seed + 211L * i, xi_lower = xl)
    pot_rows[[loc$location]] <- data.frame(
      location = loc$location, n_excess = length(exc), threshold = thr,
      f_xi = p$F$xi, f_se_xi = p$F$se_xi,
      s_xi_mean = p$s$xi[["mean"]], s_xi_lo = p$s$xi[["q05"]], s_xi_hi = p$s$xi[["q95"]],
      s_endpoint_mean = p$s_endpoint_mean, s_pr_infinite = p$s$pr_heavy,
      ef_xi_mean = p$e$xi[["mean"]], rhat_xi = p$diag$rhat_xi, ess_xi = p$diag$ess_xi,
      s_runtime = p$S$elapsed, f_runtime = p$F$elapsed, stored_numbers = B + 2L,
      wide_support = wide_support, stringsAsFactors = FALSE)
    }

    if (do_windows) {
    ## ---- overlapping moving windows (F and S in every window) ---------------
    wins <- make_windows(xb)
    if (length(wins$mats) > 0L) {
      one_win <- function(wv) {
        st <- standardise(wv, "gev"); sy <- build_symbol(wv, B)
        Fw <- laplace_full(wv, "gev", st)
        Sw <- run_symbolic_fit(wv, "gev", B, st, sym = sy,
                              n_iter = CONFIG$n_iter, n_burn = CONFIG$n_burn,
                              n_chains = CONFIG$app_window_chains, seed = seed + 5L)
        sw <- summarise_draws(Sw$draws, "gev", st)
        data.frame(
          xi_F = Fw$xi,
          xi_F_lo = if (is.finite(Fw$se_xi)) Fw$xi - stats::qnorm(0.95) * Fw$se_xi else NA_real_,
          xi_F_hi = if (is.finite(Fw$se_xi)) Fw$xi + stats::qnorm(0.95) * Fw$se_xi else NA_real_,
          ci_ok = is.finite(Fw$se_xi),
          xi_S = sw$xi[["q50"]], xi_S_mean = sw$xi[["mean"]],
          xi_S_lo = sw$xi[["q05"]], xi_S_hi = sw$xi[["q95"]],
          z100_F = Fw$z100, z100_S = sw$z100[["median"]],
          z100_S_lo = sw$z100[["q05"]], z100_S_hi = sw$z100[["q95"]],
          stringsAsFactors = FALSE)
      }
      if (n_cores > 1L && .Platform$OS.type == "unix") {
        wres <- parallel::mclapply(wins$mats, one_win, mc.cores = min(n_cores, length(wins$mats)))
      } else {
        wres <- lapply(wins$mats, one_win)
      }
      traj <- do.call(rbind, wres)
      traj$location <- loc$location
      traj$window_id <- seq_len(nrow(traj))
      traj$window_start <- mm$date[wins$starts]
      window_traj[[loc$location]] <- traj
      utils::write.csv(traj, file.path(CONFIG$out_dir,
        sprintf("windows_%s.csv", tolower(gsub("[^a-zA-Z0-9]+", "_", loc$location)))),
        row.names = FALSE)

      tr_xi <- neff_trend(traj$xi_S_mean, traj$window_start)
      tr_z  <- neff_trend(traj$z100_S, traj$window_start)
      ok <- traj$ci_ok & is.finite(traj$xi_F_lo) & is.finite(traj$xi_F_hi)
      frac_in <- if (any(ok)) mean(traj$xi_S_mean[ok] >= traj$xi_F_lo[ok] &
                                     traj$xi_S_mean[ok] <= traj$xi_F_hi[ok]) else NA_real_
      win_rows[[loc$location]] <- data.frame(
        location = loc$location, n_windows = nrow(traj), rho1 = tr_xi$rho1,
        n_eff = tr_xi$n_eff, slope_xi = tr_xi$slope, se_slope_xi = tr_xi$se,
        slope_z100 = tr_z$slope, frac_S_in_F = frac_in, stringsAsFactors = FALSE)
    }
  }
    }

  list(desc = do.call(rbind, desc_rows), gev = do.call(rbind, gev_rows),
       pot = do.call(rbind, pot_rows), win = do.call(rbind, win_rows),
       traj = window_traj, trace = trace_store, provenance = provenance)
}

## ---------------------------------------------------------------------------
## 21e. Shape-support sensitivity table (Remark rem:bound), now reporting the
##      symbolic fraction of draws sitting on the edge of each support.
## ---------------------------------------------------------------------------
write_misspec_table <- function(app, wide, file = "tab_sim_misspec.tex",
                                label = "tab:sim-misspec") {
  ms <- merge(app$gev[, c("location", "f_xi", "s_xi_mean", "s_pr_at_lower")],
              wide$gev[, c("location", "f_xi", "s_xi_mean", "s_pr_at_lower")],
              by = "location", suffixes = c("_narrow", "_wide"))
  ms <- ms[match(LOCATIONS$location, ms$location), ]
  rows <- sprintf("%s & %s & %s & %s & %s & %s & %s \\\\",
    ms$location, tex_num(ms$f_xi_narrow), tex_num(ms$f_xi_wide),
    tex_num(ms$s_xi_mean_narrow), sprintf("%.3f", ms$s_pr_at_lower_narrow),
    tex_num(ms$s_xi_mean_wide), sprintf("%.3f", ms$s_pr_at_lower_wide))
  write_lines_tex(c(TEX_HEAD(),
    "\\begin{table}[t]", "\\centering", "\\small",
    "\\caption{Sensitivity of the monthly-maximum GEV shape estimate to the width of the shape support, the check claimed in Remark~\\ref{rem:bound}. \\emph{Narrow} is the support $(-1/2,1/2)$ used throughout the paper and \\emph{wide} is $(-1,1/2)$. $\\hat\\xi_{\\mathrm{F}}$ is the full-data penalised maximum likelihood estimate and $\\bar\\xi_{\\mathrm{S}}$ the endpoint-augmented symbolic posterior mean at $B=10$. The columns headed $\\Pr(\\mathrm{bound})$ give the fraction of symbolic draws within $10^{-3}$ of the lower edge of the corresponding support. Widening the support moves the estimates onto the new bound rather than freeing them, which is the evidence that the data prefer a tail shorter than the finite-variance boundary.}",
    sprintf("\\label{%s}", label),
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{lcccccc}", "\\toprule",
    "Site & $\\hat\\xi_{\\mathrm{F}}^{\\mathrm{narrow}}$ & $\\hat\\xi_{\\mathrm{F}}^{\\mathrm{wide}}$ & $\\bar\\xi_{\\mathrm{S}}^{\\mathrm{narrow}}$ & $\\Pr(\\mathrm{bound})$ & $\\bar\\xi_{\\mathrm{S}}^{\\mathrm{wide}}$ & $\\Pr(\\mathrm{bound})$ \\\\",
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}}",
    "\\end{table}"), file)
  utils::write.csv(ms, file.path(CONFIG$out_dir, "sim_misspec.csv"), row.names = FALSE)
  invisible(ms)
}

## ---------------------------------------------------------------------------
## 21f. Master driver (supersedes the earlier definition).
## ---------------------------------------------------------------------------
run_all_analysis <- function(mode = c("main", "pilot"), do_sim = TRUE, do_app = TRUE,
                             reps = CONFIG$reps, seed = CONFIG$seed,
                             n_iter = CONFIG$n_iter, n_burn = CONFIG$n_burn,
                             n_iter_app = CONFIG$n_iter_app, n_burn_app = CONFIG$n_burn_app,
                             n_cores = CONFIG$n_cores,
                             script_path = "endpoint_augmented_extremes.R") {
  mode <- match.arg(mode)
  dir_setup()
  t_start <- proc.time()[["elapsed"]]
  prog("================================================================")
  prog(" %s run started %s", toupper(mode), format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  prog(" outputs -> %s | figures -> %s", normalizePath(CONFIG$out_dir),
       normalizePath(CONFIG$fig_dir))
  prog("================================================================")
  writeLines(capture.output(str(CONFIG)), file.path(CONFIG$out_dir, "run_config.txt"))
  write_reproducibility_record(script_path)

  sim_out <- NULL
  if (do_sim) {
    sim_out <- run_simulation(reps = reps, seed = seed, n_iter = n_iter,
                              n_burn = n_burn, n_chains = CONFIG$n_chains,
                              n_cores = n_cores)
    s <- summarise_simulation(sim_out$detail)
    utils::write.csv(sim_out$detail, file.path(CONFIG$out_dir, "sim_detail.csv"), row.names = FALSE)
    utils::write.csv(s, file.path(CONFIG$out_dir, "sim_summary.csv"), row.names = FALSE)
    write_sim_table_family(s, "gev", c("bounded", "exponential"), "tab_sim_gev.tex", "tab:sim-gev",
      paste0("GEV results over the bounded ($\\xi=-0.2$) and exponential ($\\xi=0$) regimes. ",
             "$m$ is the number of blocks and $B$ the interior resolution. F is the full-data ",
             "penalised maximum likelihood estimate with its asymptotic standard error, S the ",
             "endpoint-augmented symbolic posterior of \\eqref{eq:like} and EF the endpoint-free ",
             "grouped posterior of \\eqref{eq:control}. Coverage is of the nominal $90\\%$ interval ",
             "and is reported as undefined (--) where the full-data asymptotic standard error does ",
             "not exist because the estimate is pinned at the shape bound. Numbers in parentheses ",
             "are Monte-Carlo standard errors over replications. The RMSE ratio column (S only) is ",
             "the paired ratio of symbolic to full-data RMSE of $\\xi$, with its delta-method ",
             "Monte-Carlo standard error."),
      sprintf("Bias and RMSE are of the shape parameter $\\xi$; true value $\\xi=-0.2$ (bounded) or $\\xi=0$ (exponential). %d replications per configuration.", sim_out$reps))
    write_sim_table_family(s, "gpd", c("bounded", "exponential"), "tab_sim_gpd.tex", "tab:sim-gpd",
      paste0("GPD results over the bounded ($\\xi=-0.25$) and exponential ($\\xi=0$) regimes; ",
             "$m$ is the number of threshold excesses. Columns are as in Table~\\ref{tab:sim-gev}."),
      sprintf("True shape $\\xi=-0.25$ (bounded) or $\\xi=0$ (exponential); scale $\\sigma=5$ in both. %d replications per configuration.", sim_out$reps))
    write_endpoint_table(s)
    make_fig_sim(s)
    make_fig_endpoint(s)
  }

  app <- NULL
  if (do_app) {
    app <- run_application(B = CONFIG$app_B, n_iter = n_iter_app, n_burn = n_burn_app,
                           n_cores = n_cores, seed = seed)
    if (!is.null(app$desc) && nrow(app$desc) > 0L) write_app_desc_table(app$desc)
    if (!is.null(app$gev) && nrow(app$gev) > 0L) write_app_gev_table(app$gev)
    if (!is.null(app$pot) && nrow(app$pot) > 0L) write_app_pot_table(app$pot)
    if (!is.null(app$win) && nrow(app$win) > 0L) write_app_win_table(app$win)
    make_fig_app_summary(app); make_fig_pot(app); make_fig_windows(app); make_fig_trace(app)
    utils::write.csv(app$gev, file.path(CONFIG$out_dir, "app_summary.csv"), row.names = FALSE)
    utils::write.csv(app$pot, file.path(CONFIG$out_dir, "app_pot_summary.csv"), row.names = FALSE)
    utils::write.csv(app$desc, file.path(CONFIG$out_dir, "app_sites.csv"), row.names = FALSE)

    prog("[APP] shape-support sensitivity check (Remark rem:bound) ...")
    wide <- run_application(B = CONFIG$app_B, n_iter = n_iter_app, n_burn = n_burn_app,
                            n_cores = n_cores, seed = seed, wide_support = TRUE)
    ms <- write_misspec_table(app, wide)

    if (!is.null(sim_out)) write_provenance(sim_out, app)
    saveRDS(list(config = CONFIG, seed = seed, sim = sim_out, app = app, wide = wide,
                 misspec = ms, elapsed = proc.time()[["elapsed"]] - t_start),
            file.path(CONFIG$out_dir, "results.rds"))
    prog("[OUT] results.rds")
  }

  prog("================================================================")
  prog(" finished in %.1f s", proc.time()[["elapsed"]] - t_start)
  prog("================================================================")
  invisible(list(sim = sim_out, app = app))
}

## ===========================================================================
## 22. Final overrides
##   (a) get_power(): cached CSV is re-read with an explicit Date column, and
##       the API text is stripped of CR so that the last field parses.
##   (b) write_app_pot_table(): prints the threshold itself.
##   (c) run_simulation(): reports the first failing job instead of dying
##       inside do.call(rbind, ...).
## ===========================================================================
get_power <- function(location, lat, lon, start, end, cache_dir = CONFIG$cache_dir,
                      community = CONFIG$app_community, retries = 4L) {
  key <- tolower(gsub("[^a-zA-Z0-9]+", "_", location))
  cache_file <- file.path(cache_dir, paste0(key, "_power_daily.csv"))
  url <- POWER_URL(lat, lon, start, end, community)
  if (file.exists(cache_file)) {
    df <- utils::read.csv(cache_file, stringsAsFactors = FALSE)
    df$date <- as.Date(df$date)
    df <- df[order(df$date), ]
    return(list(data = df, url = url, source = "cache", retrieved = NA_character_))
  }
  tmp <- tempfile(fileext = ".csv"); ok <- FALSE
  for (k in seq_len(retries)) {
    ok <- tryCatch({ utils::download.file(url, tmp, quiet = TRUE, mode = "wb"); TRUE },
      error = function(e) { prog("[POWER] attempt %d failed (%s): %s", k, location,
                                 conditionMessage(e)); FALSE })
    if (ok) break
    Sys.sleep(2 * k)
  }
  if (!ok || !file.exists(tmp) || file.info(tmp)$size < 500) {
    stop(sprintf(paste0("NASA POWER retrieval failed for %s after %d attempts.\n  URL: %s\n",
                        "  Set ENDPOINT_ALLOW_SURROGATE=1 to fall back to a SYNTHETIC series ",
                        "(every affected table is then labelled SYNTHETIC, not NASA POWER)."),
                 location, retries, url))
  }
  txt <- gsub("\r", "", readLines(tmp, warn = FALSE))
  i <- grep("^YEAR,", txt)[1L]
  if (is.na(i)) stop("NASA POWER CSV for ", location, " has no YEAR header row")
  df <- utils::read.csv(text = paste(txt[i:length(txt)], collapse = "\n"),
                        stringsAsFactors = FALSE)
  df$date <- as.Date(sprintf("%d-%03d", df$YEAR, df$DOY), format = "%Y-%j")
  keep <- c("date", "T2M_MAX", "RH2M")
  if (!all(keep %in% names(df))) stop("NASA POWER CSV for ", location,
                                     " is missing T2M_MAX and/or RH2M")
  df <- df[, keep]
  df <- df[is.finite(df$T2M_MAX) & is.finite(df$RH2M) &
             df$T2M_MAX > -900 & df$RH2M > -900, ]
  df <- df[order(df$date), ]
  utils::write.csv(df, cache_file, row.names = FALSE)
  list(data = df, url = url, source = "api",
       retrieved = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
}

write_app_pot_table <- function(pot, file = "tab_app_pot.tex", label = "tab:app-pot") {
  rows <- sprintf("%s & %d & %.2f & %.3f (%.3f) & %.3f [%.3f, %.3f] & %.3f & %s & %.3f \\\\",
    pot$location, pot$n_excess, pot$threshold, pot$f_xi, pot$f_se_xi,
    pot$s_xi_mean, pot$s_xi_lo, pot$s_xi_hi, abs(pot$s_xi_mean - pot$f_xi),
    tex_num(pot$s_endpoint_mean, 2), pot$s_pr_infinite)
  write_lines_tex(c(TEX_HEAD(),
    "\\begin{table}[t]", "\\centering", "\\small",
    "\\caption{Threshold-excess (GPD) analysis of the daily THI series. The threshold is the site-specific empirical $95$th percentile of the full daily series; $n_{\\mathrm{exc}}$ is the number of exceedances over it and $n_{\\mathrm{exc}}/N$ the fraction of the record retained as a value in the analysis. F is the full-data maximum likelihood estimate with its asymptotic standard error and S the endpoint-augmented symbolic posterior at $B=10$, that is twelve stored numbers per site. The column $|\\bar\\xi_{\\mathrm{S}}-\\hat\\xi_{\\mathrm{F}}|$ is the absolute difference between the symbolic posterior mean and the full-data estimate; the next column is the symbolic $90\\%$ credible interval for $\\xi$ and $\\bar y^{F}_{\\mathrm{S}}$ the posterior mean of the fitted endpoint above the threshold, reported as unavailable (--) when the symbolic posterior places mass on $\\xi\\ge0$, for which no finite endpoint exists.}",
    sprintf("\\label{%s}", label),
    "\\begin{tabular}{lrcccccc}", "\\toprule",
    "Site & $n_{\\mathrm{exc}}$ & Threshold & $\\hat\\xi_{\\mathrm{F}}$ (s.e.) & $\\bar\\xi_{\\mathrm{S}}$ $[q_{.05},q_{.95}]$ & $|\\bar\\xi_{\\mathrm{S}}-\\hat\\xi_{\\mathrm{F}}|$ & $\\bar y^{F}_{\\mathrm{S}}$ & $\\Pr(\\xi\\ge0)_{\\mathrm{S}}$ \\\\",
    "\\midrule", rows, "\\bottomrule", "\\end{tabular}",
    "\\end{table}"), file)
}

run_simulation <- function(reps = CONFIG$reps, seed = CONFIG$seed,
                           n_iter = CONFIG$n_iter, n_burn = CONFIG$n_burn,
                           n_chains = CONFIG$n_chains, grid = sim_grid(),
                           n_cores = CONFIG$n_cores, store_example = TRUE,
                           keep_draws_rep = 1L) {
  prog("[SIM] %d configurations x %d replications = %d datasets x 3 methods",
       nrow(grid), reps, nrow(grid) * reps)
  jobs <- list(); k <- 0L
  for (i in seq_len(nrow(grid))) for (r in seq_len(reps)) {
    k <- k + 1L
    jobs[[k]] <- list(fam = grid$fam[i], regime = grid$regime[i], m = grid$m[i],
                      B = grid$B[i], rep_id = r, seed = seed + 7919L * i + 31L * r,
                      store = (store_example && r == keep_draws_rep),
                      label = sprintf("%s/%s/m=%d/B=%d/rep=%d", grid$fam[i], grid$regime[i],
                                      grid$m[i], grid$B[i], r))
  }
  run_one <- function(j) {
    tryCatch(simulate_one(j$fam, j$regime, j$m, j$B, j$rep_id, seed = j$seed,
                          n_iter = n_iter, n_burn = n_burn, n_chains = n_chains,
                          store_draws = j$store),
             error = function(e) list(rows = NULL, error = conditionMessage(e), label = j$label))
  }
  t00 <- proc.time()[["elapsed"]]
  if (n_cores > 1L && .Platform$OS.type == "unix") {
    res <- parallel::mclapply(jobs, run_one, mc.cores = min(n_cores, length(jobs)),
                              mc.preschedule = FALSE, mc.set.seed = TRUE)
  } else {
    res <- lapply(jobs, run_one)
  }
  elapsed <- proc.time()[["elapsed"]] - t00
  bad <- !vapply(res, function(z) is.list(z) && !is.null(z$rows), logical(1L))
  if (any(bad)) {
    res <- lapply(res, function(z) if (is.list(z)) z else
      list(rows = NULL, label = "unknown",
           error = paste(utils::capture.output(print(z)), collapse = " | ")))
    prog("[SIM] %d of %d jobs failed. First failure: %s -> %s", sum(bad), length(res),
         res[[which(bad)[1L]]]$label, res[[which(bad)[1L]]]$error)
    res <- res[!bad]
    if (length(res) < 1L) stop("run_simulation(): every job failed")
  }
  detail <- do.call(rbind, lapply(res, function(z) z$rows))
  example <- NULL
  for (z in res) if (!is.null(z$draws)) { example <- z; break }
  prog("[SIM] grid finished in %.1f s (%d datasets)", elapsed, nrow(detail) / 3L)
  list(detail = detail, example = example, elapsed = elapsed, grid = grid, reps = reps)
}
## =============================================================================
## 22b. Final overrides
##   (a) write_endpoint_table(): reports the usable fraction alongside the
##       failure fraction (they sum to one).
##   (b) make_fig_endpoint() / make_fig_sim(): pick the largest available sample
##       size when the manuscript's m = 200 is absent (smoke runs) and refuse to
##       draw a panel whose limits are not finite.
## =============================================================================
write_endpoint_table <- function(s, file = "tab_endpoint.tex", label = "tab:endpoint") {
  d <- s[s$regime == "bounded", ]
  d <- d[order(match(d$fam, CONFIG$fams), d$m, d$B, match(d$method, c("F", "S", "EF"))), ]
  body <- character(0)
  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    body <- c(body, sprintf("%s & %d & %d & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
      ifelse(r$fam == "gev", "GEV", "GPD"), r$m, r$B, METHOD_TEX[[r$method]],
      tex_num(r$ep_median_err), tex_num(r$ep_mse, 2),
      if (is.finite(r$ep_mse_se)) sprintf("(%.2f)", r$ep_mse_se) else "(--)",
      sprintf("%.3f", r$ep_frac_usable),
      sprintf("%.3f", r$ep_fail), sprintf("%.3f", r$ep_pr_inf),
      if (is.finite(r$rmse)) tex_num(r$rmse) else "--"))
  }
  write_lines_tex(c(TEX_HEAD(),
    "\\begin{table}[t]", "\\centering", "\\small",
    "\\caption{Finite-endpoint recovery in the bounded regimes of the $96$-configuration grid ($\\xi=-0.2$ for the GEV with $\\mu=25$, $\\sigma=5$, endpoint $y^{F}=50.0$; $\\xi=-0.25$ for the GPD with $\\sigma=5$, endpoint $y^{F}=20.0$). Errors are measured against the true endpoint of the generating model. The median error is the median over draws of $|\\hat y^{F}-y^{F}|$. MSE is the mean squared endpoint error with its Monte-Carlo standard error over replications in parentheses; it is taken over the draws that return a usable finite endpoint, because a draw with $\\xi\\ge0$ returns $y^{F}=+\\infty$. \\emph{Usable} and \\emph{Fail.} are the complementary fractions of draws with a finite endpoint at or below the failure cut and of draws above it or at infinity, and $\\Pr(\\xi\\ge0)$ isolates the unbounded draws; the three columns together show how often the endpoint functional is informative at all. The RMSE of $\\xi$ is repeated for reference.}",
    sprintf("\\label{%s}", label),
    "\\begin{tabular}{llrlccccccc}", "\\toprule",
    "\\multicolumn{4}{c}{} & \\multicolumn{3}{c}{Endpoint error} & \\multicolumn{4}{c}{Diagnostics} \\\\",
    "\\cmidrule(lr){5-7}\\cmidrule(lr){8-11}",
    "Family & $m$ & $B$ & Method & Median $|\\hat y^{F}-y^{F}|$ & MSE & s.e.(MSE) & Usable & Fail. frac. & $\\Pr(\\xi\\ge0)$ & RMSE$(\\xi)$ \\\\",
    "\\midrule", body, "\\bottomrule", "\\end{tabular}",
    "\\end{table}"), file)
}

## m_sel: the manuscript panel uses m = 200; fall back to the largest m present.
pick_m <- function(vals, preferred = 200L) if (preferred %in% vals) preferred else max(vals)

make_fig_endpoint <- function(s, m_sel = NULL) {
  d0 <- s[s$regime == "bounded", ]
  if (nrow(d0) == 0L) { prog("[FIG] fig_endpoint.pdf skipped (no bounded cells)"); return(invisible(NULL)) }
  if (is.null(m_sel)) m_sel <- pick_m(unique(d0$m))
  d <- d0[d0$m == m_sel, ]
  Bs <- sort(unique(d$B))
  open_pdf("fig_endpoint.pdf", 9.5, 7.5)

  rng2 <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) c(0, 1) else range(x) }

  ## (a) RMSE of the endpoint estimate against resolution
  ymax <- rng2(d$ep_rmse[d$method %in% c("S", "EF")])[2]
  graphics::plot(NA, xlim = range(Bs), ylim = c(0, 1.15 * ymax),
                 xlab = "Resolution B", ylab = expression(RMSE(~ hat(y)^F)), main = "")
  for (fam in CONFIG$fams) for (mm in c("S", "EF")) {
    dd <- d[d$fam == fam & d$method == mm, ]; dd <- dd[order(dd$B), ]
    if (nrow(dd) == 0L) next
    graphics::lines(dd$B, dd$ep_rmse, type = "o", lwd = 1.7, cex = 1.05,
                    pch = ifelse(mm == "S", 16, 17), lty = ifelse(fam == "gev", 1, 2),
                    col = METHOD_COL[[mm]])
  }
  graphics::legend("topright", legend = c("Symbolic (S), GEV", "Symbolic (S), GPD",
                                          "Endpoint-free (EF), GEV", "Endpoint-free (EF), GPD"),
                   col = METHOD_COL[c("S", "S", "EF", "EF")], lty = c(1, 2, 1, 2),
                   pch = c(16, 16, 17, 17), bty = "n", cex = 0.8, lwd = 1.6)
  panel_tag("(a)")

  ## (b) coverage of the 90% interval for xi
  graphics::plot(NA, xlim = range(Bs), ylim = c(0.5, 1.0),
                 xlab = "Resolution B", ylab = "Coverage of 90% interval for xi", main = "")
  graphics::abline(h = 0.9, lty = 3, col = "grey45")
  for (fam in CONFIG$fams) for (mm in c("S", "EF")) {
    dd <- d[d$fam == fam & d$method == mm, ]; dd <- dd[order(dd$B), ]
    if (nrow(dd) == 0L) next
    graphics::lines(dd$B, dd$cover90, type = "o", lwd = 1.7, pch = ifelse(mm == "S", 16, 17),
                    lty = ifelse(fam == "gev", 1, 2), col = METHOD_COL[[mm]])
  }
  panel_tag("(b)")

  ## (c) RMSE of xi against resolution
  ymax2 <- rng2(d$rmse[d$method %in% c("S", "EF")])[2]
  graphics::plot(NA, xlim = range(Bs), ylim = c(0, 1.15 * ymax2),
                 xlab = "Resolution B", ylab = expression(RMSE(xi)), main = "")
  for (fam in CONFIG$fams) for (mm in c("S", "EF")) {
    dd <- d[d$fam == fam & d$method == mm, ]; dd <- dd[order(dd$B), ]
    if (nrow(dd) == 0L) next
    graphics::lines(dd$B, dd$rmse, type = "o", lwd = 1.7, pch = ifelse(mm == "S", 16, 17),
                    lty = ifelse(fam == "gev", 1, 2), col = METHOD_COL[[mm]])
  }
  panel_tag(sprintf("(c) m = %d", m_sel))

  ## (d) the median endpoint error, pooled over the bounded configurations
  med <- d0[d0$method %in% c("S", "EF"), ]
  med <- med[is.finite(med$ep_median_err) & med$ep_median_err > 0, ]
  if (nrow(med) > 0L) {
    graphics::boxplot(ep_median_err ~ method, data = med, log = "y",
                      names = c("Endpoint-free (EF)", "Symbolic (S)"),
                      col = c(METHOD_COL[["EF"]], METHOD_COL[["S"]]),
                      border = "grey20", ylab = "Median endpoint error (log scale)",
                      xlab = "", main = "")
  } else {
    graphics::plot.new()
  }
  panel_tag("(d)")
  grDevices::dev.off()
  prog("[FIG] fig_endpoint.pdf (m = %d)", m_sel)
}

make_fig_sim <- function(s) {
  open_pdf("fig_sim_main.pdf", 9.5, 7.5)
  ss <- s[s$method == "S", ]
  if (nrow(ss) == 0L) { grDevices::dev.off(); return(invisible(NULL)) }
  ms <- sort(unique(ss$m))
  rg2 <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0L) c(0, 1) else range(x) }

  ## (a) ratio of symbolic to full-data RMSE of xi against m at B = 10 and B = 50
  ymax <- rg2(ss$rmse_ratio_vs_full)[2]
  graphics::plot(NA, xlim = range(ms), ylim = c(0, max(2, 1.1 * ymax)), log = "x",
                 xlab = "Sample size m (log scale)",
                 ylab = "RMSE(S) / RMSE(F) for xi", main = "")
  graphics::abline(h = 1, lty = 3, col = "grey45")
  for (b in c(10L, 50L)) for (rg in CONFIG$regimes) {
    d <- ss[ss$B == b & ss$regime == rg, ]
    if (nrow(d) == 0L) next
    d <- d[order(d$m), ]
    graphics::lines(d$m, d$rmse_ratio_vs_full, type = "o",
                    pch = ifelse(b == 10L, 16, 17), lty = match(rg, CONFIG$regimes),
                    col = REGIME_COL[[rg]], lwd = 1.7, cex = 1.1)
  }
  graphics::legend("topright", legend = c(paste(CONFIG$regimes, "tail"), "B = 10", "B = 50"),
                   lty = c(1, 2, 3, NA, NA), pch = c(NA, NA, NA, 16, 17),
                   col = c(REGIME_COL[CONFIG$regimes], "grey25", "grey25"),
                   bty = "n", cex = 0.8, lwd = 1.7)
  panel_tag("(a)")

  ## (b) coverage of the 90% symbolic interval
  graphics::plot(NA, xlim = range(ms), ylim = c(0.6, 1.0), log = "x",
                 xlab = "Sample size m (log scale)", ylab = "Coverage of 90% interval", main = "")
  graphics::abline(h = 0.90, lty = 3, col = "grey45")
  for (rg in CONFIG$regimes) {
    d <- ss[ss$regime == rg, ]
    if (nrow(d) == 0L) next
    d <- stats::aggregate(cover90 ~ m, data = d, FUN = function(v) mean(v, na.rm = TRUE))
    d <- d[order(d$m), ]
    graphics::lines(d$m, d$cover90, type = "o", pch = 16, lwd = 1.7, col = REGIME_COL[[rg]])
  }
  graphics::legend("bottomright", legend = CONFIG$regimes, col = REGIME_COL[CONFIG$regimes],
                   lty = 1, pch = 16, bty = "n", cex = 0.85)
  panel_tag("(b)")

  ## (c) cost of one symbolic fit against m: flat, because the cost is O(B)
  ylim_c <- rg2(c(ss$runtime, s$runtime[s$method == "F"]))
  ylim_c <- c(max(1e-3, ylim_c[1] * 0.8), ylim_c[2] * 1.6)
  graphics::plot(NA, xlim = range(ms), ylim = ylim_c, log = "xy",
                 xlab = "Sample size m (log scale)", ylab = "Cost of one fit (s)", main = "")
  for (rg in CONFIG$regimes) {
    d <- ss[ss$regime == rg, ]
    if (nrow(d) == 0L) next
    d <- stats::aggregate(runtime ~ m, data = d, FUN = safe_median); d <- d[order(d$m), ]
    graphics::lines(d$m, d$runtime, type = "o", pch = 16, lwd = 1.7, col = REGIME_COL[[rg]])
  }
  dF <- stats::aggregate(runtime ~ m, data = s[s$method == "F", ], FUN = safe_median)
  dF <- dF[order(dF$m), ]
  graphics::lines(dF$m, dF$runtime, type = "o", pch = 4, lwd = 1.7, col = "grey20")
  graphics::legend("topleft", legend = c(CONFIG$regimes, "Full-data Laplace"),
                   col = c(REGIME_COL[CONFIG$regimes], "grey20"), lty = 1,
                   pch = c(16, 16, 16, 4), bty = "n", cex = 0.85)
  panel_tag("(c)")

  ## (d) cost against resolution B, symbolic versus full-data Laplace
  graphics::plot(NA, xlim = range(ss$B), ylim = ylim_c, log = "y",
                 xlab = "Resolution B", ylab = "Cost of one fit (s)", main = "")
  for (fam in CONFIG$fams) {
    d <- ss[ss$fam == fam, ]
    if (nrow(d) == 0L) next
    d <- stats::aggregate(runtime ~ B, data = d, FUN = safe_median)
    graphics::lines(d$B, d$runtime, type = "o", pch = 16, lwd = 1.7, col = FAM_COL[[fam]])
  }
  graphics::abline(h = safe_median(s$runtime[s$method == "F"]), lty = 2, col = "grey20")
  graphics::legend("topleft", legend = c("Symbolic, GEV", "Symbolic, GPD",
                                         "Full-data Laplace (median)"),
                   col = c(FAM_COL[["gev"]], FAM_COL[["gpd"]], "grey20"),
                   lty = c(1, 1, 2), pch = c(16, 16, NA), bty = "n", cex = 0.85)
  panel_tag("(d)")
  grDevices::dev.off()
  prog("[FIG] fig_sim_main.pdf")
}
## =============================================================================
## 24. FINAL application driver (supersedes the definition in section 21).
##     Self-contained on purpose: every stage is wrapped so that one failing
##     site, window or fit cannot abort the whole run, and nothing is silently
##     imputed - a stage that fails is reported and its cell is left out.
## =============================================================================
run_application <- function(loc_table = LOCATIONS, B = CONFIG$app_B,
                            n_iter = CONFIG$n_iter_app, n_burn = CONFIG$n_burn_app,
                            n_chains = CONFIG$n_chains, n_cores = CONFIG$n_cores,
                            seed = CONFIG$seed, wide_support = FALSE,
                            do_pot = TRUE, do_windows = TRUE) {
  provenance <- list(); desc_rows <- list(); gev_rows <- list(); pot_rows <- list()
  win_rows <- list(); trace_store <- list(); window_traj <- list()
  xl <- if (wide_support) CONFIG$xi_lower_wide else CONFIG$xi_lower

  for (i in seq_len(nrow(loc_table))) {
    loc <- loc_table[i, ]
    prog("[APP] %s ...", loc$location)

    ## ---- retrieval: NASA POWER, or an explicitly requested synthetic series --
    df_raw <- tryCatch({
      if (CONFIG$allow_surrogate) {
        tryCatch(get_power(loc$location, loc$latitude, loc$longitude,
                           CONFIG$app_start, CONFIG$app_end),
          error = function(e) {
            prog("[APP] API unavailable (%s): SYNTHETIC series used", conditionMessage(e))
            list(data = surrogate_power(loc$location, loc$latitude, loc$longitude,
                                        CONFIG$app_start, CONFIG$app_end),
                 url = NA_character_, source = "SYNTHETIC", retrieved = NA_character_)
          })
      } else {
        get_power(loc$location, loc$latitude, loc$longitude,
                  CONFIG$app_start, CONFIG$app_end)
      }
    }, error = function(e) {
      prog("[APP] %s ABORTED: %s", loc$location, conditionMessage(e)); NULL
    })
    if (is.null(df_raw)) next

    d <- df_raw$data
    d$THI <- thi(d$T2M_MAX, d$RH2M)
    d <- d[is.finite(d$THI), ]
    if (nrow(d) < 400L) { prog("[APP] %s skipped: only %d valid days", loc$location, nrow(d)); next }
    provenance[[loc$location]] <- c(url = df_raw$url %||% "NA", source = df_raw$source,
                                    retrieved = df_raw$retrieved %||% NA_character_,
                                    n_days = as.character(nrow(d)))

    thr <- unname(stats::quantile(d$THI, CONFIG$app_threshold_p, type = 7))
    exc <- d$THI[d$THI > thr] - thr
    mm  <- monthly_maxima(d, "THI")
    desc_rows[[loc$location]] <- data.frame(
      location = loc$location, latitude = loc$latitude, longitude = loc$longitude,
      climate = loc$climate, n_days = nrow(d),
      mean_thi = mean(d$THI), median_thi = stats::median(d$THI),
      max_thi = max(d$THI), q95_thi = thr, n_excess = length(exc),
      n_blocks = nrow(mm), stringsAsFactors = FALSE)

    ## ---- monthly-maximum GEV: F, S and the endpoint-free control EF ---------
    xb <- mm$value
    g <- tryCatch(fit_gev_stationary(xb, B = B, n_iter = n_iter, n_burn = n_burn,
                                     n_chains = n_chains, seed = seed + 101L * i,
                                     xi_lower = xl),
                  error = function(e) {
                    prog("[APP] GEV fit failed for %s: %s", loc$location, conditionMessage(e)); NULL
                  })
    if (!is.null(g)) {
      gev_rows[[loc$location]] <- data.frame(
        location = loc$location, n_blocks = length(xb),
        f_xi = g$F$xi, f_se_xi = g$F$se_xi, f_endpoint = g$F$endpoint, f_z100 = g$F$z100,
        f_ci_lo = if (is.finite(g$F$se_xi)) g$F$xi - stats::qnorm(0.95) * g$F$se_xi else NA_real_,
        f_ci_hi = if (is.finite(g$F$se_xi)) g$F$xi + stats::qnorm(0.95) * g$F$se_xi else NA_real_,
        s_xi_mean = g$s$xi[["mean"]], s_xi_lo = g$s$xi[["q05"]], s_xi_hi = g$s$xi[["q95"]],
        s_endpoint_mean = g$s$endpoint[["mean"]], s_z100_mean = g$s$z100[["mean"]],
        s_z10_mean = g$s$z10[["mean"]],
        s_pr_at_lower = g$s$pr_at_lower, s_pr_at_upper = g$s$pr_at_upper,
        s_pr_heavy = g$s$pr_heavy, rhat_xi = g$diag$rhat_xi, ess_xi = g$diag$ess_xi,
        acceptance = g$S$acceptance,
        ef_xi_mean = g$e$xi[["mean"]], ef_xi_lo = g$e$xi[["q05"]], ef_xi_hi = g$e$xi[["q95"]],
        ef_endpoint_mean = g$e$endpoint[["mean"]], ef_z100_mean = g$e$z100[["mean"]],
        ef_runtime = g$EF$elapsed, s_runtime = g$S$elapsed, f_runtime = g$F$elapsed,
        stored_numbers = g$sym$stored_numbers, wide_support = wide_support,
        stringsAsFactors = FALSE)
      if (identical(loc$location, CONFIG$app_trace_site)) {
        trace_store[[loc$location]] <- list(fit = g$S, sum = g$s, diag = g$diag,
                                            xi_lower = CONFIG$xi_lower,
                                            xi_upper = CONFIG$xi_upper)
      }
    }

    ## ---- daily threshold excesses (GPD): full data and symbolic -------------
    if (do_pot) {
      p <- tryCatch(fit_gpd_excesses(exc, B = B, n_iter = n_iter, n_burn = n_burn,
                                     n_chains = n_chains, seed = seed + 211L * i,
                                     xi_lower = xl),
                    error = function(e) {
                      prog("[APP] GPD fit failed for %s: %s", loc$location, conditionMessage(e)); NULL
                    })
      if (!is.null(p)) {
        pot_rows[[loc$location]] <- data.frame(
          location = loc$location, n_excess = length(exc), threshold = thr,
          f_xi = p$F$xi, f_se_xi = p$F$se_xi,
          s_xi_mean = p$s$xi[["mean"]], s_xi_lo = p$s$xi[["q05"]], s_xi_hi = p$s$xi[["q95"]],
          s_endpoint_mean = p$s_endpoint_mean, s_pr_infinite = p$s$pr_heavy,
          ef_xi_mean = p$e$xi[["mean"]], rhat_xi = p$diag$rhat_xi, ess_xi = p$diag$ess_xi,
          s_runtime = p$S$elapsed, f_runtime = p$F$elapsed, stored_numbers = B + 2L,
          wide_support = wide_support, stringsAsFactors = FALSE)
      }
    }

    ## ---- overlapping moving windows: refit F and S in every window ----------
    if (do_windows) {
      wins <- make_windows(xb)
      if (length(wins$mats) > 0L) {
        one_win <- function(wv, wseed) {
          tryCatch({
            st <- standardise(wv, "gev"); sy <- build_symbol(wv, B)
            Fw <- laplace_full(wv, "gev", st)
            Sw <- run_symbolic_fit(wv, "gev", B, st, sym = sy,
                                   n_iter = CONFIG$n_iter, n_burn = CONFIG$n_burn,
                                   n_chains = CONFIG$app_window_chains, seed = wseed)
            sw <- summarise_draws(Sw$draws, "gev", st)
            data.frame(
              xi_F = Fw$xi,
              xi_F_lo = if (is.finite(Fw$se_xi)) Fw$xi - stats::qnorm(0.95) * Fw$se_xi else NA_real_,
              xi_F_hi = if (is.finite(Fw$se_xi)) Fw$xi + stats::qnorm(0.95) * Fw$se_xi else NA_real_,
              ci_ok = is.finite(Fw$se_xi),
              xi_S = sw$xi[["q50"]], xi_S_mean = sw$xi[["mean"]],
              xi_S_lo = sw$xi[["q05"]], xi_S_hi = sw$xi[["q95"]],
              z100_F = Fw$z100, z100_S = sw$z100[["median"]],
              z100_S_lo = sw$z100[["q05"]], z100_S_hi = sw$z100[["q95"]],
              stringsAsFactors = FALSE)
          }, error = function(e) {
            prog("[APP] window fit failed (%s, seed %d): %s", loc$location, wseed,
                 conditionMessage(e)); NULL
          })
        }
        wres <- vector("list", length(wins$mats))
        for (wj in seq_along(wins$mats)) {
          wres[[wj]] <- one_win(wins$mats[[wj]], seed + 5L + 17L * wj)
        }
        keep <- !vapply(wres, is.null, logical(1L))
        prog("[APP] %s: %d of %d window fits succeeded", loc$location, sum(keep), length(keep))
        if (any(keep)) {
          traj <- do.call(rbind, wres[keep])
          traj$location <- loc$location
          traj$window_id <- seq_len(nrow(traj))
          traj$window_start <- mm$date[wins$starts][keep]
          window_traj[[loc$location]] <- traj
          utils::write.csv(traj, file.path(CONFIG$out_dir,
            sprintf("windows_%s.csv", tolower(gsub("[^a-zA-Z0-9]+", "_", loc$location)))),
            row.names = FALSE)
          tr_xi <- neff_trend(traj$xi_S_mean, traj$window_start)
          tr_z  <- neff_trend(traj$z100_S, traj$window_start)
          ok <- traj$ci_ok & is.finite(traj$xi_F_lo) & is.finite(traj$xi_F_hi)
          frac_in <- if (any(ok)) mean(traj$xi_S_mean[ok] >= traj$xi_F_lo[ok] &
                                         traj$xi_S_mean[ok] <= traj$xi_F_hi[ok]) else NA_real_
          win_rows[[loc$location]] <- data.frame(
            location = loc$location, n_windows = nrow(traj), rho1 = tr_xi$rho1,
            n_eff = tr_xi$n_eff, slope_xi = tr_xi$slope, se_slope_xi = tr_xi$se,
            slope_z100 = tr_z$slope, frac_S_in_F = frac_in, stringsAsFactors = FALSE)
        }
      }
    }
  }

  mk <- function(rows) if (length(rows) > 0L) do.call(rbind, rows) else NULL
  list(desc = mk(desc_rows), gev = mk(gev_rows), pot = mk(pot_rows),
       win = mk(win_rows), traj = window_traj, trace = trace_store,
       provenance = provenance)
}

run_all_analysis(mode = Sys.getenv("ENDPOINT_MODE", unset = "main"))
