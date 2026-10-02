# CodeA.py
# Annual-balance Monte Carlo estimator of one-year ruin probability
# Author: (Iacopo Cardosi Carrara)
# Description:
#   Computes the one-year (or generic T) ruin probability under the Cramér–Lundberg model
#   using an annual-balance Monte Carlo estimator: it draws the yearly claim count (Poisson),
#   i.i.d. severities (lognormal), applies the selected treaty to retained losses, and
#   checks ruin on the year-end surplus only (U_end = u + premium - sum(retained losses)).
#   Supports proportional (quota-share, alpha) and non-proportional (excess-of-loss, M) treaties.
#   Note: the year-end balance approach tends to underestimate finite-horizon ruin because it
#   ignores intra-year timing; see CodeB.py for the pathwise estimator.
# Dependencies: numpy, matplotlib


import numpy as np
import matplotlib.pyplot as plt

# Set random seed for reproducibility
SEED = 12345
np.random.seed(SEED)
print(f"Random seed set to {SEED}")

# Model parameters
LAMBDA_CLAIMS = 10
MEAN_SEVERITY = 1000.0
STD_SEVERITY = 2000.0
LOADING_INSURER = 0.30
LOADING_REINSURER = 0.15
INITIAL_SURPLUS = 100000.0
N_SIMULATIONS = 50000

# Derived lognormal parameters
SIGMA2 = np.log(1 + (STD_SEVERITY**2) / (MEAN_SEVERITY**2))
MU = np.log(MEAN_SEVERITY) - 0.5 * SIGMA2
E_X = MEAN_SEVERITY

def expected_min_severity(M, sample_size=100000):
    samples = np.random.lognormal(mean=MU, sigma=np.sqrt(SIGMA2), size=sample_size)
    return np.mean(np.minimum(samples, M))

def simulate_ruin(alpha=None, M=None):
    use_quota_share = alpha is not None
    if use_quota_share:
        expected_insured = alpha * E_X
        expected_reinsured = (1 - alpha) * E_X
    else:
        expected_insured = expected_min_severity(M)
        expected_reinsured = E_X - expected_insured
    premium_insurer = (1 + LOADING_INSURER) * LAMBDA_CLAIMS * expected_insured
    premium_reinsurer = (1 + LOADING_REINSURER) * LAMBDA_CLAIMS * expected_reinsured
    ruin_count = 0
    for _ in range(N_SIMULATIONS):
        n_claims = np.random.poisson(LAMBDA_CLAIMS)
        if n_claims > 0:
            severities = np.random.lognormal(mean=MU, sigma=np.sqrt(SIGMA2), size=n_claims)
            if use_quota_share:
                insured_losses = alpha * severities
            else:
                insured_losses = np.minimum(severities, M)
            total_loss = insured_losses.sum()
        else:
            total_loss = 0.0
        final_surplus = INITIAL_SURPLUS + premium_insurer - total_loss
        if final_surplus < 0:
            ruin_count += 1
    return ruin_count / N_SIMULATIONS

def _demo():
    alphas = np.linspace(0.1, 1.0, 10)
    ruin_alpha = [simulate_ruin(alpha=a) for a in alphas]
    retentions = [500, 1000, 2000, 5000, 10000]
    ruin_M = [simulate_ruin(M=m) for m in retentions]
    print("Ruin probability without reinsurance (alpha=1.0): {:.6f}".format(simulate_ruin(alpha=1.0)))
    for a, rp in zip(alphas, ruin_alpha):
        print(f"alpha={a:.2f}: {rp:.6f}")
    for m, rp in zip(retentions, ruin_M):
        print(f"M={m}: {rp:.6f}")
    # Optional: save plots
    plt.figure()
    plt.plot(alphas, ruin_alpha, marker='o')
    plt.title('Ruin probability vs quota-share retention')
    plt.xlabel('alpha')
    plt.ylabel('Ruin probability')
    plt.grid(True)
    plt.savefig('quota_ruin_plot_jupyter.png')
    plt.close()
    plt.figure()
    plt.plot(retentions, ruin_M, marker='o')
    plt.title('Ruin probability vs excess-of-loss retention')
    plt.xlabel('M')
    plt.ylabel('Ruin probability')
    plt.grid(True)
    plt.savefig('excess_ruin_plot_jupyter.png')
    plt.close()

if __name__ == '__main__':
    _demo()
