# CodeB.py
# Pathwise finite-horizon ruin probability
# Author: (Iacopo Cardosi Carrara)
# Description:
#   Computes the one-year (or generic T) ruin probability under the Cramér–Lundberg model
#   using a pathwise Monte Carlo estimator that simulates claim arrival times and accrues
#   premium continuously. Supports quota-share (alpha) and excess-of-loss (M) treaties.
# Dependencies: numpy

import numpy as np

# Set random seed for reproducibility
SEED = 12345
np.random.seed(SEED)
print(f"Random seed set to {SEED}")

def simulate_ruin_pathwise(
    T: float = 1.0,
    alpha: float | None = None,
    M: float | None = None,
    n_paths: int = 50_000,
    lambda_claims: float = 10.0,
    mean_sev: float = 1000.0,
    std_sev: float = 2000.0,
    loading: float = 0.30,
    initial_surplus: float = 100_000.0,
    seed: int | None = None,
    inner_mc_size: int = 100_000,
) -> float:
    """
    Estimate finite-horizon ruin probability on [0, T] by simulating claim arrival times
    (Poisson with rate lambda_claims) and accruing premium continuously at rate
    c = (1 + loading) * lambda_claims * E[retained loss per claim].

    Treaty:
      - Quota-share: specify alpha in (0, 1]; retained loss = alpha * X.
      - Excess-of-loss: specify M > 0; retained loss = min(X, M).

    Exactly one of alpha or M must be provided.

    Args:
        T: Time horizon (years). Default 1.0.
        alpha: Quota-share retention in (0, 1]; set to 1.0 for baseline (no reinsurance).
        M: Excess-of-loss priority (> 0).
        n_paths: Number of Monte Carlo paths.
        lambda_claims: Poisson intensity (expected number of claims per year).
        mean_sev: Target mean of severity X (lognormal).
        std_sev: Target standard deviation of severity X (lognormal).
        loading: Safety loading for pricing the retained book.
        initial_surplus: Initial capital u.
        seed: Optional RNG seed for reproducibility.
        inner_mc_size: Sample size for estimating E[min(X, M)] under XL.

    Returns:
        Estimated finite-horizon ruin probability on [0, T].

    Raises:
        ValueError: If treaty specification is invalid.
    """
    # Treaty validation: exactly one of alpha or M must be set
    if (alpha is None) == (M is None):
        raise ValueError("Specify exactly one treaty: either alpha (quota-share) or M (excess-of-loss).")

    if alpha is not None and not (0.0 < alpha <= 1.0):
        raise ValueError("alpha must be in (0, 1].")

    if M is not None and not (M > 0.0):
        raise ValueError("M must be > 0 for excess-of-loss.")

    # Optional seed for reproducibility
    if seed is not None:
        np.random.seed(seed)

    # Map (mean, std) to (mu, sigma^2) for Lognormal
    if mean_sev <= 0 or std_sev <= 0:
        raise ValueError("mean_sev and std_sev must be positive.")
    sigma2 = float(np.log(1.0 + (std_sev**2) / (mean_sev**2)))
    mu = float(np.log(mean_sev) - 0.5 * sigma2)

    # Retained expected severity and loss transform
    if alpha is not None:
        EX_retained = alpha * mean_sev
        def transform(x: float) -> float:
            return alpha * x
    else:
        # Estimate E[min(X, M)] via inner Monte Carlo (can be cached for repeated M)
        samples = np.random.lognormal(mean=mu, sigma=np.sqrt(sigma2), size=inner_mc_size)
        EX_retained = float(np.minimum(samples, M).mean())
        def transform(x: float) -> float:
            return float(np.minimum(x, M))

    # Premium rate on retained business
    c = (1.0 + loading) * lambda_claims * EX_retained

    ruin = 0
    for _ in range(n_paths):
        t = 0.0
        U = float(initial_surplus)
        while True:
            # Next inter-arrival
            dt = float(np.random.exponential(1.0 / lambda_claims))
            if t + dt >= T:
                # Accrue premium until T and stop (no more claims in horizon)
                U += c * (T - t)
                break

            # Move to claim time and accrue premium
            t += dt
            U += c * dt

            # Draw severity and apply treaty
            x = float(np.random.lognormal(mean=mu, sigma=np.sqrt(sigma2)))
            loss = transform(x)
            U -= loss

            if U < 0.0:
                ruin += 1
                break

    return ruin / n_paths


if __name__ == "__main__":
    # Minimal smoke tests (small n_paths for speed). Increase n_paths for stable estimates.
    print("Pathwise ruin (baseline, no reinsurance):",
          simulate_ruin_pathwise(T=1.0, alpha=1.0, n_paths=50000, seed=12345))
    print("Pathwise ruin (quota-share alpha=0.6):",
          simulate_ruin_pathwise(T=1.0, alpha=0.6, n_paths=50000, seed=12345))
    print("Pathwise ruin (excess-of-loss M=2000):",
          simulate_ruin_pathwise(T=1.0, M=2000, n_paths=50000, seed=12345))
