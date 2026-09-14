# Benchmark report: residual_jacobian_fusion_integrated_repeat

Julia 1.13.0; ForwardDiff 1.4.6; AMD Ryzen 7 7800X3D 8-Core Processor           ; threads=1; rounds=3.

Times are medians of round medians, in nanoseconds. Negative changes mean faster.
A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.

| Fixture | N | Residual ns | Jacobian ns | Solve ns | Warm solve ns | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| scalar_linear | 1 | 2.49 | 2.54 | 18.62 | 18.97 | 1 |
| scalar_maxwell | 1 | 2.71 | 3.20 | 39.73 | 39.73 | 1 |
| scalar_powerlaw | 1 | 25.86 | 35.15 | 475.00 | 104.55 | 6 |
| scalar_creep | 1 | 56.74 | 81.69 | 331.82 | 190.91 | 2 |
| volumetric_maxwell | 2 | 3.09 | 13.58 | 62.61 | 62.61 | 1 |
| nested | 3 | 3.02 | 27.54 | 64.76 | 64.76 | 1 |
| siblings | 4 | 3.19 | 35.29 | 141.18 | 141.18 | 1 |
| kelvin_voigt | 2 | 3.12 | 17.17 | 48.95 | 49.31 | 1 |
| plastic | 2 | 3.33 | 21.26 | 51.52 | 51.91 | 1 |
| dilatant | 3 | 19.21 | 35.00 | 72.29 | 72.84 | 1 |

Interleaved full-solve comparison against the frozen original solver:

| Fixture | Original cold ns | Current cold ns | Cold change | Warm change | Current cold bytes |
| --- | ---: | ---: | --- | --- | ---: |
| scalar_linear | 26.90 | 18.62 | -30.8% (inconclusive) | -29.5% (inconclusive) | 0 |
| scalar_maxwell | 45.09 | 39.73 | -11.9% (inconclusive) | -11.9% (inconclusive) | 0 |
| scalar_powerlaw | 620.00 | 475.00 | -23.4% (faster) | -25.8% (faster) | 0 |
| scalar_creep | 513.64 | 331.82 | -35.4% (faster) | -38.2% (faster) | 0 |
| volumetric_maxwell | 65.22 | 62.61 | -4.0% (inconclusive) | -4.0% (inconclusive) | 0 |
| nested | 76.19 | 64.76 | -15.0% (faster) | -15.0% (faster) | 0 |
| siblings | 145.10 | 141.18 | -2.7% (inconclusive) | -2.7% (inconclusive) | 0 |
| kelvin_voigt | 53.85 | 48.95 | -9.1% (faster) | -8.4% (faster) | 0 |
| plastic | 60.61 | 51.52 | -15.0% (faster) | -14.4% (faster) | 0 |
| dilatant | 104.88 | 72.29 | -31.1% (faster) | -31.3% (faster) | 0 |

Preprocessing experiments below compare Jacobian kernels only, excluding setup.

| Fixture | Cached callable | Cached config | Chunk-1 config | Prepared equations | Config setup ns / bytes |
| --- | --- | --- | --- | --- | ---: |
| scalar_linear | +0.0% (inconclusive) | +7.9% (inconclusive) | -7.9% (inconclusive) | +102.6% (slower) | 7.30 / 32 |
| scalar_maxwell | +6.5% (inconclusive) | +9.7% (inconclusive) | +0.0% (inconclusive) | +361.3% (slower) | 6.23 / 32 |
| scalar_powerlaw | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +1.7% (inconclusive) | 5.46 / 32 |
| scalar_creep | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +24.2% (slower) | 7.80 / 32 |
| volumetric_maxwell | -0.0% (inconclusive) | -0.0% (inconclusive) | -0.0% (inconclusive) | +6.9% (slower) | 6.46 / 64 |
| nested | +1.1% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +123.0% (slower) | 9.89 / 112 |
| siblings | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +2316.7% (slower) | 6.90 / 176 |
| kelvin_voigt | -2.1% (inconclusive) | -2.5% (inconclusive) | -2.5% (inconclusive) | +129.2% (slower) | 6.49 / 64 |
| plastic | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +137.2% (slower) | 6.32 / 64 |
| dilatant | -1.6% (inconclusive) | -1.6% (inconclusive) | -1.6% (inconclusive) | +97.6% (slower) | 9.29 / 112 |

Residual+Jacobian fusion experiment (item 3): a single ForwardDiff.jacobian!/DiffResults
call against the solver's current separate primal-then-Jacobian calls, and against a
Jacobian-only call, all measured within this run.

| Fixture | Separate (residual+jacobian) ns | Fused ns | Fused vs separate | Jacobian-only ns | Fused vs jacobian-only | Fused bytes |
| --- | ---: | ---: | --- | ---: | --- | ---: |
| scalar_linear | 2.70 | 2.70 | +0.1% (inconclusive) | 2.54 | +6.7% (inconclusive) | 0 |
| scalar_maxwell | 4.44 | 4.43 | -0.3% (inconclusive) | 3.20 | +38.3% (slower) | 0 |
| scalar_powerlaw | 73.71 | 51.30 | -30.4% (faster) | 35.15 | +45.9% (slower) | 0 |
| scalar_creep | 140.28 | 83.02 | -40.8% (faster) | 81.69 | +1.6% (inconclusive) | 0 |
| volumetric_maxwell | 14.06 | 14.36 | +2.1% (inconclusive) | 13.58 | +5.7% (slower) | 0 |
| nested | 38.16 | 35.50 | -7.0% (faster) | 27.54 | +28.9% (slower) | 0 |
| siblings | 48.75 | 35.65 | -26.9% (faster) | 35.29 | +1.0% (inconclusive) | 0 |
| kelvin_voigt | 18.24 | 14.22 | -22.1% (faster) | 17.17 | -17.2% (faster) | 0 |
| plastic | 17.20 | 18.17 | +5.7% (slower) | 21.26 | -14.6% (faster) | 0 |
| dilatant | 34.00 | 33.45 | -1.6% (inconclusive) | 35.00 | -4.4% (inconclusive) | 0 |

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
