# File guide

For the conceptual overview see [`README.md`](README.md).

| File | What it does |
|---|---|
| `CodeA.py` | **Annual-balance estimator.** Sets the model parameters and maps (mean, std) to lognormal (μ, σ²). `expected_min_severity(M)` estimates E[min(X, M)] by Monte Carlo. `simulate_ruin(alpha=..., M=...)` prices the retained book and, for 50,000 replications, draws the yearly claims, applies the treaty and checks ruin on the year-end surplus. `_demo()` runs quota-share for α = 0.1, …, 1.0 and excess-of-loss for M = 500, …, 10,000, prints the results and saves `quota_ruin_plot_jupyter.png` and `excess_ruin_plot_jupyter.png`. |
| `CodeB.py` | **Pathwise estimator.** `simulate_ruin_pathwise(T, alpha or M, ...)` validates the treaty, computes the premium rate on the retained book and, for each path, simulates exponential inter-arrival times, accrues premium continuously and stops as soon as the surplus goes below zero. The main block runs the baseline (α = 1), quota-share α = 0.6 and excess-of-loss M = 2,000. |
