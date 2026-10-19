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

# ... [KEEP ALL REMAINING FUNCTIONS EXACTLY AS THEY ARE IN YOUR ORIGINAL SCRIPT] ...
# (Sections 2 through 25 remain unchanged, as they contain no confidential data)
