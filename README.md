# Ruin Theory with Reinsurance — a Monte Carlo simulation

Python code that estimates the **one-year ruin probability** of an insurer under the **Cramér-Lundberg model** with heavy-tailed (lognormal) claims. It compares the baseline without reinsurance with two reinsurance treaties: **quota-share** and **excess-of-loss**.

The code accompanies the project work *Ruin Theory with Reinsurance — A Monte Carlo Simulation* (I. Cardosi Carrara, Sant'Anna School of Advanced Studies, course in Mathematical Methods for Insurance).

---

## 1. Risk model

**Surplus process (Cramér-Lundberg):**

U(t) = u + c·t − Σᵢ₌₁^N(t) Xᵢ

- u: initial surplus;
- N(t) ~ Poisson(λt): number of claims up to time t;
- Xᵢ: i.i.d. claim sizes, independent of N(t);
- c = (1 + θ) λ E[X]: premium rate with safety loading θ.

**Finite-horizon ruin probability:**

ψ_T(u) = P( inf₀≤t≤T U(t) < 0 ),  estimated here for T = 1 year.

**Claim severity.** Claims are lognormal to capture heavy tails. Given a target mean m and standard deviation s:

σ² = ln(1 + s²/m²),  μ = ln m − σ²/2

## 2. Reinsurance

| Treaty | Retained loss per claim | Ceded loss | Insurer's premium rate |
|---|---|---|---|
| No reinsurance | X | 0 | (1+θ) λ E[X] |
| Quota-share, retention α ∈ (0, 1] | α X | (1−α) X | (1+θ) α λ E[X] |
| Excess-of-loss, priority M | min(X, M) | (X − M)⁺ | (1+θ) λ E[min(X, M)] |

For excess-of-loss, E[min(X, M)] is estimated numerically with an inner Monte Carlo sample. Code A also computes the reinsurer's premium on the ceded part (loading θ_R), but only the insurer's premium enters the surplus.

## 3. Two Monte Carlo estimators

**Code A — annual balance.** For each replication:
1. draw the yearly number of claims N ~ Poisson(λ) and the severities;
2. apply the treaty and sum the retained losses;
3. compute the year-end surplus U_end = u + premium − retained losses;
4. count a ruin if U_end < 0.

It is fast, but it only looks at the end of the year: a large claim early in the year that pushes the surplus below zero only temporarily is not detected.

**Code B — pathwise.** For each path:
1. simulate the claim arrival times (exponential inter-arrival times with rate λ);
2. accrue the premium continuously between arrivals;
3. subtract each retained claim at its arrival time;
4. count a ruin as soon as U(t) < 0 during the year.

Because it detects intra-year crossings, it captures the timing of claims, and one typically finds ψ̂_path ≥ ψ̂_annual.

## 4. Parameters

| Parameter | Value |
|---|---|
| Claim intensity λ | 10 claims per year |
| Severity mean m / std s | 1,000 / 2,000 (lognormal) |
| Insurer loading θ | 0.30 |
| Reinsurer loading θ_R | 0.15 |
| Initial surplus u | 100,000 |
| Replications / paths | 50,000 |
| Random seed | 12345 |

## 5. Results (from the project report)

**Code A (annual balance):**
- the baseline estimate is of order 10⁻⁴;
- under quota-share the estimated ruin probability grows with the retention α, from 0 at low retentions to 1.4·10⁻⁴ at α = 1;
- under excess-of-loss no ruin is observed for any priority M ∈ {500, 1,000, 2,000, 5,000, 10,000}.

**Code B (pathwise):**

| Scenario | ψ̂ | 95% CI |
|---|---|---|
| No reinsurance | 0.000100 | [0.000012, 0.000188] |
| Quota-share, α = 0.6 | 0 | [0, 0.000060] |
| Excess-of-loss, M = 2,000 | 0 | [0, 0.000060] |

When no ruin is observed, the 95% upper bound 3/R is reported (rule of three). Otherwise the interval is ψ̂ ± 1.96·√(ψ̂(1−ψ̂)/R).

**Takeaways**:
- With these inputs one-year ruin is a rare event, of order 10⁻⁴.
- The pathwise baseline (1.0·10⁻⁴) is above the annual-balance one (0.8·10⁻⁴), because intra-year dips are captured.
- Both treaties reduce ruin. Larger replications or stressed parameters would be needed to separate them sharply.

## 6. Possible extensions

- **Over-dispersed claim counts**: Negative Binomial instead of Poisson.
- **Alternative heavy tails**: Pareto or Log-Gamma severities.
- **Variance reduction** for rare events: importance sampling, control variates, stratification.
- **Retention optimization**: choose α or M to minimize total cost subject to a ruin or VaR/TVaR constraint.

## 7. How to run

Requirements: Python ≥ 3.10, `numpy` and `matplotlib` (Code A only).

```bash
python CodeA.py   # annual-balance estimates for several alpha and M; saves two plots
python CodeB.py   # pathwise estimates for the baseline, alpha = 0.6 and M = 2000
```

A description of every file is in [`CODE.md`](CODE.md).

## Main references
- Cramér-Lundberg model and reinsurance treaties: course material, Mathematical Methods for Insurance, Sant'Anna School of Advanced Studies.
