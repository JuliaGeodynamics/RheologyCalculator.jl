# Benchmark report: residual_jacobian_fusion_integrated

Julia 1.13.0; ForwardDiff 1.4.6; AMD Ryzen 7 7800X3D 8-Core Processor           ; threads=1; rounds=3.

Times are medians of round medians, in nanoseconds. Negative changes mean faster.
A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.

| Fixture | N | Residual ns | Jacobian ns | Solve ns | Warm solve ns | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| scalar_linear | 1 | 2.67 | 2.68 | 18.26 | 18.40 | 1 |
| scalar_maxwell | 1 | 2.68 | 2.64 | 22.18 | 22.18 | 1 |
| scalar_powerlaw | 1 | 17.51 | 34.34 | 463.64 | 100.00 | 6 |
| scalar_creep | 1 | 39.45 | 59.38 | 229.17 | 129.17 | 2 |
| volumetric_maxwell | 2 | 4.20 | 21.90 | 79.57 | 79.57 | 1 |
| nested | 3 | 5.04 | 33.90 | 92.11 | 92.11 | 1 |
| siblings | 4 | 3.12 | 39.13 | 144.00 | 144.00 | 1 |
| kelvin_voigt | 2 | 3.11 | 21.46 | 44.72 | 44.72 | 1 |
| plastic | 2 | 3.32 | 19.84 | 53.49 | 53.49 | 1 |
| dilatant | 3 | 19.36 | 36.28 | 73.17 | 73.17 | 1 |

Interleaved full-solve comparison against the frozen original solver:

| Fixture | Original cold ns | Current cold ns | Cold change | Warm change | Current cold bytes |
| --- | ---: | ---: | --- | --- | ---: |
| scalar_linear | 26.30 | 18.26 | -30.6% (faster) | -29.4% (faster) | 0 |
| scalar_maxwell | 29.93 | 22.18 | -25.9% (faster) | -25.9% (faster) | 0 |
| scalar_powerlaw | 600.00 | 463.64 | -22.7% (faster) | -25.8% (faster) | 0 |
| scalar_creep | 360.87 | 229.17 | -36.5% (faster) | -41.1% (faster) | 0 |
| volumetric_maxwell | 83.87 | 79.57 | -5.1% (faster) | -5.1% (faster) | 0 |
| nested | 103.95 | 92.11 | -11.4% (inconclusive) | -11.4% (inconclusive) | 0 |
| siblings | 150.00 | 144.00 | -4.0% (inconclusive) | -4.0% (inconclusive) | 0 |
| kelvin_voigt | 47.20 | 44.72 | -5.3% (faster) | -5.3% (faster) | 0 |
| plastic | 61.54 | 53.49 | -13.1% (faster) | -12.0% (faster) | 0 |
| dilatant | 106.10 | 73.17 | -31.0% (faster) | -31.0% (faster) | 0 |

Preprocessing experiments below compare Jacobian kernels only, excluding setup.

| Fixture | Cached callable | Cached config | Chunk-1 config | Prepared equations | Config setup ns / bytes |
| --- | --- | --- | --- | --- | ---: |
| scalar_linear | -8.5% (faster) | +0.0% (inconclusive) | +0.0% (inconclusive) | +106.1% (slower) | 5.27 / 32 |
| scalar_maxwell | +9.4% (slower) | +0.0% (inconclusive) | +0.0% (inconclusive) | +384.0% (slower) | 5.22 / 32 |
| scalar_powerlaw | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.5% (inconclusive) | 5.26 / 32 |
| scalar_creep | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +10.5% (slower) | 5.26 / 32 |
| volumetric_maxwell | +0.2% (inconclusive) | +0.0% (inconclusive) | +0.2% (inconclusive) | +8.6% (slower) | 8.97 / 64 |
| nested | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +200.0% (slower) | 11.05 / 112 |
| siblings | -11.1% (inconclusive) | -11.1% (inconclusive) | -11.1% (inconclusive) | +1995.6% (slower) | 6.79 / 176 |
| kelvin_voigt | -20.8% (faster) | -20.8% (faster) | -20.8% (faster) | +87.4% (slower) | 6.19 / 64 |
| plastic | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +158.9% (slower) | 6.38 / 64 |
| dilatant | -3.0% (inconclusive) | -3.0% (inconclusive) | -3.0% (inconclusive) | +117.3% (slower) | 8.58 / 112 |

Residual+Jacobian fusion experiment (item 3): a single ForwardDiff.jacobian!/DiffResults
call against the solver's current separate primal-then-Jacobian calls, and against a
Jacobian-only call, all measured within this run.

| Fixture | Separate (residual+jacobian) ns | Fused ns | Fused vs separate | Jacobian-only ns | Fused vs jacobian-only | Fused bytes |
| --- | ---: | ---: | --- | ---: | --- | ---: |
| scalar_linear | 2.68 | 2.68 | -0.1% (inconclusive) | 2.68 | +0.2% (inconclusive) | 0 |
| scalar_maxwell | 4.10 | 4.10 | -0.2% (inconclusive) | 2.64 | +55.1% (slower) | 0 |
| scalar_powerlaw | 56.40 | 38.82 | -31.2% (faster) | 34.34 | +13.1% (slower) | 0 |
| scalar_creep | 97.35 | 59.48 | -38.9% (faster) | 59.38 | +0.2% (inconclusive) | 0 |
| volumetric_maxwell | 23.34 | 24.26 | +4.0% (inconclusive) | 21.90 | +10.8% (slower) | 0 |
| nested | 46.50 | 42.67 | -8.2% (inconclusive) | 33.90 | +25.9% (slower) | 0 |
| siblings | 41.61 | 42.20 | +1.4% (inconclusive) | 39.13 | +7.8% (slower) | 0 |
| kelvin_voigt | 25.07 | 22.84 | -8.9% (faster) | 21.46 | +6.4% (slower) | 0 |
| plastic | 18.22 | 18.98 | +4.2% (inconclusive) | 19.84 | -4.4% (inconclusive) | 0 |
| dilatant | 33.92 | 33.85 | -0.2% (inconclusive) | 36.28 | -6.7% (faster) | 0 |

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
