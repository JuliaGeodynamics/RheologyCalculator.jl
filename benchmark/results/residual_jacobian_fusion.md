# Benchmark report: residual_jacobian_fusion

Julia 1.13.0; ForwardDiff 1.4.6; AMD Ryzen 7 7800X3D 8-Core Processor           ; threads=1; rounds=3.

Times are medians of round medians, in nanoseconds. Negative changes mean faster.
A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.

| Fixture | N | Residual ns | Jacobian ns | Solve ns | Warm solve ns | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| scalar_linear | 1 | 2.65 | 2.45 | 17.75 | 17.82 | 1 |
| scalar_maxwell | 1 | 2.66 | 2.82 | 39.66 | 39.66 | 1 |
| scalar_powerlaw | 1 | 17.18 | 33.73 | 461.90 | 113.64 | 6 |
| scalar_creep | 1 | 39.46 | 59.18 | 281.48 | 181.48 | 2 |
| volumetric_maxwell | 2 | 3.10 | 18.52 | 50.00 | 50.00 | 1 |
| nested | 3 | 3.16 | 27.27 | 71.29 | 72.28 | 1 |
| siblings | 4 | 3.11 | 37.50 | 137.04 | 136.54 | 1 |
| kelvin_voigt | 2 | 3.26 | 18.31 | 52.99 | 52.99 | 1 |
| plastic | 2 | 3.30 | 21.19 | 55.38 | 55.38 | 1 |
| dilatant | 3 | 19.15 | 34.69 | 80.25 | 80.25 | 1 |

Interleaved full-solve comparison against the frozen original solver:

| Fixture | Original cold ns | Current cold ns | Cold change | Warm change | Current cold bytes |
| --- | ---: | ---: | --- | --- | ---: |
| scalar_linear | 26.04 | 17.75 | -31.8% (faster) | -31.5% (faster) | 0 |
| scalar_maxwell | 45.14 | 39.66 | -12.2% (faster) | -11.3% (faster) | 0 |
| scalar_powerlaw | 590.48 | 461.90 | -21.8% (faster) | -13.8% (faster) | 0 |
| scalar_creep | 362.96 | 281.48 | -22.4% (faster) | -18.3% (faster) | 0 |
| volumetric_maxwell | 54.35 | 50.00 | -8.0% (faster) | -8.0% (faster) | 0 |
| nested | 74.51 | 71.29 | -4.3% (inconclusive) | -3.9% (inconclusive) | 0 |
| siblings | 140.38 | 137.04 | -2.4% (inconclusive) | -1.7% (inconclusive) | 0 |
| kelvin_voigt | 52.24 | 52.99 | +1.4% (inconclusive) | +1.5% (inconclusive) | 0 |
| plastic | 58.46 | 55.38 | -5.3% (faster) | -5.3% (faster) | 0 |
| dilatant | 102.53 | 80.25 | -21.7% (faster) | -21.7% (faster) | 0 |

Preprocessing experiments below compare Jacobian kernels only, excluding setup.

| Fixture | Cached callable | Cached config | Chunk-1 config | Prepared equations | Config setup ns / bytes |
| --- | --- | --- | --- | --- | ---: |
| scalar_linear | +8.9% (slower) | +8.9% (slower) | +0.0% (inconclusive) | +100.1% (slower) | 4.97 / 32 |
| scalar_maxwell | +0.0% (inconclusive) | -7.1% (faster) | -7.1% (faster) | +582.4% (slower) | 5.19 / 32 |
| scalar_powerlaw | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +1.8% (inconclusive) | 5.15 / 32 |
| scalar_creep | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +10.3% (slower) | 5.49 / 32 |
| volumetric_maxwell | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +2.2% (inconclusive) | 7.41 / 64 |
| nested | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +125.5% (slower) | 9.76 / 112 |
| siblings | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +2000.0% (slower) | 6.98 / 176 |
| kelvin_voigt | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +105.1% (slower) | 5.96 / 64 |
| plastic | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +117.5% (slower) | 5.79 / 64 |
| dilatant | -4.8% (inconclusive) | -4.8% (inconclusive) | -4.8% (inconclusive) | +69.2% (slower) | 8.19 / 112 |

Residual+Jacobian fusion experiment (item 3): a single ForwardDiff.jacobian!/DiffResults
call against the solver's current separate primal-then-Jacobian calls, and against a
Jacobian-only call, all measured within this run.

| Fixture | Separate (residual+jacobian) ns | Fused ns | Fused vs separate | Jacobian-only ns | Fused vs jacobian-only | Fused bytes |
| --- | ---: | ---: | --- | ---: | --- | ---: |
| scalar_linear | 2.86 | 2.85 | -0.7% (inconclusive) | 2.45 | +16.0% (slower) | 0 |
| scalar_maxwell | 4.03 | 4.02 | -0.3% (inconclusive) | 2.82 | +42.3% (slower) | 0 |
| scalar_powerlaw | 55.13 | 38.03 | -31.0% (faster) | 33.73 | +12.8% (slower) | 0 |
| scalar_creep | 96.36 | 59.32 | -38.4% (faster) | 59.18 | +0.2% (inconclusive) | 0 |
| volumetric_maxwell | 13.94 | 14.70 | +5.5% (inconclusive) | 18.52 | -20.6% (inconclusive) | 0 |
| nested | 36.68 | 34.21 | -6.7% (faster) | 27.27 | +25.4% (slower) | 0 |
| siblings | 48.02 | 40.25 | -16.2% (faster) | 37.50 | +7.3% (slower) | 0 |
| kelvin_voigt | 18.01 | 14.17 | -21.4% (faster) | 18.31 | -22.6% (faster) | 0 |
| plastic | 16.83 | 17.22 | +2.3% (inconclusive) | 21.19 | -18.7% (faster) | 0 |
| dilatant | 33.68 | 33.47 | -0.6% (inconclusive) | 34.69 | -3.5% (inconclusive) | 0 |

## Allocations

Kernels with nonzero median allocated bytes (all others are zero):

- dilatant / config_setup: 112 bytes
- dilatant / prepared_equations: 592 bytes
- kelvin_voigt / config_setup: 64 bytes
- kelvin_voigt / prepared_equations: 176 bytes
- nested / config_setup: 112 bytes
- nested / prepared_equations: 448 bytes
- plastic / config_setup: 64 bytes
- plastic / prepared_equations: 416 bytes
- scalar_creep / config_setup: 32 bytes
- scalar_linear / config_setup: 32 bytes
- scalar_maxwell / config_setup: 32 bytes
- scalar_powerlaw / config_setup: 32 bytes
- siblings / config_setup: 176 bytes
- siblings / prepared_equations: 3184 bytes
- volumetric_maxwell / config_setup: 64 bytes
