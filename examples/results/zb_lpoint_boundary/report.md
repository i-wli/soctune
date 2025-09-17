# Periodic SOC study

Scales multiply the shell-resolved Jha–Heine constants.

The reference source is GPAW 25.7.0 PW/PBE post-SCF SOC (provisional).

## Valence residuals

Values and signed residuals (xTB − reference) are in meV.

| Sample | Observable | Reference | σ | L-1 | L0 | L0 valence-only |
| --- | --- | --- | --- | --- | --- | --- |
| GaAs | l_vb | 212.142 | 0.069 | 10.986 | 14.440 | 1.736 |
| GaP | l_vb | 69.030 | 0.021 | -10.051 | 0.641 | -2.727 |

## Conduction residuals

Values and signed residuals (xTB − reference) are in meV.

| Sample | Observable | Reference | σ | L-1 | L0 | L0 valence-only |
| --- | --- | --- | --- | --- | --- | --- |
| GaAs | l_cb | 88.580 | 0.075 | -58.602 | -49.093 | -58.003 |
| GaP | l_cb | 68.120 | 0.053 | -45.166 | -34.580 | -43.699 |

## Band-gap diagnostic

ΔE_g is the signed SOC-induced change in the band gap on the DFTB+ SCC mesh, in meV.

| Sample | L-1 | L0 | L0 valence-only |
| --- | --- | --- | --- |
| GaAs | -104.656 | -103.559 | -100.270 |
| GaP | -25.154 | -27.523 | -28.571 |

## Fit scales and identifiability

The fit assumes a relative model error of 10%. Reduced χ² uses the data residuals and subtracts one degree of freedom per fitted scale.

The effect column is the largest estimated splitting change across the full scale range. Identifiability requires a change greater than 3σ for at least one eligible observable. Standard errors include the Gaussian prior in the Jacobian.

| Fit | Sample | Shell | Scale | SE | Identifiable | At bound | Max effect (meV) | Reduced χ² |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| L0 | all | Ga p | 1.5000 | 0.1786 | True | True | 20.969 | 56.958 |
| L0 | all | As p | 0.9731 | 0.0961 | True | False | 206.380 | 56.958 |
| L0 | all | P p | 1.0529 | 0.1467 | True | False | 42.207 | 56.958 |
| L0 valence-only | all | Ga p | 1.0561 | 0.2364 | True | False | 17.964 | — |
| L0 valence-only | all | As p | 0.9503 | 0.0967 | True | False | 206.380 | — |
| L0 valence-only | all | P p | 1.1510 | 0.1520 | True | False | 42.207 | — |
