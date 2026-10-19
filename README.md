# Endpoint-Augmented Symbolic Likelihoods for Extreme-Value Inference

This repository contains the complete, reproducible R replication pipeline, dataset acquisition routines, Monte Carlo simulation engines, and figure generation code for the manuscript:

> **"Endpoint-augmented symbolic likelihoods for extreme-value inference"**  
> *Ahmad Hakiim Jamaluddin, Muhammad Aslam Mohd Safari, Farid Zamani Che Rose, Muhammad Jaffri Mohd Nasir, and Hani Syahida Zulkafli.*

---

## Overview

Extreme-value analyses of massive, continuously updated archives (e.g., environmental reanalysis products) face severe computational and memory constraints. This project introduces an **endpoint-augmented symbolic likelihood** framework that compresses a data stream into its sample minimum, sample maximum, and a coarse histogram of the interior, while retaining the exact endpoint densities. 

Unlike traditional grouped-data likelihoods that discard the extrema, this construction ensures that the finite upper endpoint remains identifiable even at minimal resolution. The framework is implemented in pure Base R, requiring no external compiled dependencies, making it highly portable and suitable for streaming or distributed inference.

Key methodological features of this repository include:

* **Endpoint-Augmented Symbolic Likelihood:** A genuine likelihood of the aggregated data (symbol), where the cost of a likelihood evaluation is $O(B)$ and entirely independent of the original sample size $n$.
* **96-Configuration Simulation Grid:** Evaluates the GEV and GPD families across 3 shape regimes (bounded, exponential, heavy), 4 sample sizes ($m \in \{500, 5000, 50000, 500000\}$), and 4 resolutions ($B \in \{5, 10, 20, 50\}$), comparing Full-data MLE, Symbolic (S), and Endpoint-free (EF) estimators.
* **Empirical NASA POWER Benchmark:** Automated API ingestion of 45 years (1981–2025) of daily maximum temperature and relative humidity across 5 contrasting climate zones (Singapore, Yakutsk, Death Valley, Melbourne, Chicago) to compute and analyze the daily Temperature-Humidity Index (THI).
* **Overlapping Moving-Window Analysis:** A robust protocol for assessing trends in extreme heat stress, complete with effective sample size corrections for overlapping window autocorrelation.
* **Strict Reproducibility:** A single master script generates every table (`.tex`) and figure (`.pdf`) in the manuscript, with built-in caching and synthetic data fallbacks for offline testing.

---

## Repository Structure

```text
├── endpoint_augmented_extremes.R    # Master execution script (Base R only)
├── README.md                        # Repository documentation and replication guide
├── .gitignore                       # Git ignore rules for R projects and large outputs
├── outputs/                         # Auto-generated output directory
│   ├── cache/                       # Local CSV cache for downloaded NASA POWER raw data
│   ├── tab_sim_gev.tex              # Table 1: GEV simulation results
│   ├── tab_sim_gpd.tex              # Table 2: GPD simulation results
│   ├── tab_endpoint.tex             # Table 3: Finite-endpoint recovery
│   ├── tab_app_desc.tex             # Table 4: Application site descriptions
│   ├── tab_app_gev.tex              # Table 5: Monthly-maximum GEV fits
│   ├── tab_app_pot.tex              # Table 6: Threshold-excess (GPD) fits
│   ├── tab_app_win.tex              # Table 7: Overlapping moving-window analysis
│   ├── sim_summary.csv              # Detailed simulation summary
│   ├── app_summary.csv              # Application summary
│   └── results.rds                  # Saved R objects for diagnostic verification
└── figures/                         # Auto-generated high-resolution PDF figures
    ├── fig_sim_main.pdf             # Figure 1: Simulation summary
    ├── fig_endpoint.pdf             # Figure 2: Bounded-support comparison
    ├── fig_app_summary.pdf          # Figure 3: Application summary
    ├── fig_pot.pdf                  # Figure 4: Threshold-excess analysis
    ├── fig_real_main.pdf            # Figure 5: Moving-window analysis
    └── fig_app_trace.pdf            # Figure 6: MCMC trace plots and diagnostics
