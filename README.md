# Bayesian Inference with Adaptive Histogram Resolution for Symbolic Likelihoods

This repository contains the complete, production-ready R replication pipeline, dataset acquisition routines, Monte Carlo simulation engines, and figure generation code for the manuscript:

> **"Bayesian Inference with Adaptive Histogram Resolution for Symbolic Likelihoods: A Block-Update MCMC Framework with an Application to Moving-Window Climate Series"**  
> *Ahmad Hakiim Jamaluddin, Hani Syahida Zulkafli, Muhammad Jaffri Mohd Nasir, Farid Zamani Che Rose, and Muhammad Aslam Mohd Safari.*

---

## Overview

Repeated Bayesian inference over overlapping time series windows is computationally intensive. This project introduces a novel adaptive symbolic likelihood framework that treats the histogram resolution ($B$) as an unknown discrete parameter, estimated jointly with continuous model parameters via a block-update Markov chain Monte Carlo (MCMC) algorithm. This compresses raw sample evaluations while allowing the data to dynamically select the optimal amount of symbolic compression per window.

Key methodological features of this repository include:

* **Block-Update MCMC Engine:** Combines adaptive Random-Walk Metropolis (or HMC/NUTS) for continuous parameters with a discrete Metropolis-Hastings step for the histogram resolution $B$, avoiding transdimensional reversible-jump machinery.
* **Factorial Monte Carlo Simulation Grid:** Evaluates 4 data-generating scenarios (Gaussian, Heavy-tailed Student-$t$, Skewed mixture, Window-shift mixture) across 3 sample sizes ($n \in \{500, 2000, 10000\}$) with 100 replicates, comparing Full-data, Fixed-coarse ($B=5$), Fixed-fine ($B=30$), and Adaptive methods.
* **Empirical NASA POWER Benchmark:** Automated API ingestion of 30 years (1994–2023) of daily meteorological data across 5 contrasting climate zones (Subarctic, Tropical Rainforest, Desert Hot, Oceanic Temperate, Humid Continental) to compute moving-window Temperature-Humidity Index (THI) series.
* **HPC-Optimized Bootstrap:** Includes a self-bootstrapping Conda environment setup for seamless, dependency-managed execution on SLURM/LiCO high-performance computing clusters.

---

## Repository Structure

```text
├── run_adaptive_symbolic.R      # Master production execution script (HPC-ready)
├── README.md                    # Repository documentation and replication guide
├── .gitignore                   # Git ignore rules for R projects and large outputs
└── results_adaptive_symbolic/   # Output directory (auto-generated upon execution)
    ├── cache/                   # Local CSV cache for downloaded NASA POWER raw data
    ├── raw/                     # Processed THI series per location
    ├── tables/                  # Publication-ready CSV summary and detailed tables
    ├── figures/                 # High-resolution (300 DPI) PDF figures
    ├── run_config.json          # Execution parameter JSON record
    ├── session_info.txt         # R session environment metadata
    └── sample_results.rds       # Saved MCMC chain draws for diagnostic verification
