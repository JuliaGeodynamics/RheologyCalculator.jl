# Benchmark report: line_search_reuse_inline

Julia 1.13.0; ForwardDiff 1.4.6; AMD Ryzen 7 7800X3D 8-Core Processor           ; threads=1; rounds=3.

Times are medians of round medians, in nanoseconds. Negative changes mean faster.
A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.

| Fixture | N | Residual ns | Jacobian ns | Solve ns | Warm solve ns | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| scalar_linear | 1 | 2.46 | 2.45 | 17.49 | 17.49 | 1 |
| scalar_maxwell | 1 | 2.87 | 2.83 | 21.45 | 21.45 | 1 |
| scalar_powerlaw | 1 | 17.30 | 33.71 | 460.87 | 113.64 | 6 |
| scalar_creep | 1 | 38.62 | 58.59 | 271.43 | 171.43 | 2 |
| volumetric_maxwell | 2 | 3.07 | 12.65 | 48.57 | 48.57 | 1 |
| nested | 3 | 2.97 | 27.03 | 78.95 | 78.95 | 1 |
| siblings | 4 | 3.10 | 41.67 | 138.89 | 137.74 | 1 |
| kelvin_voigt | 2 | 3.12 | 11.36 | 46.20 | 45.86 | 1 |
| plastic | 2 | 3.32 | 16.67 | 50.00 | 49.65 | 1 |
| dilatant | 3 | 19.01 | 34.65 | 86.67 | 86.67 | 1 |

Interleaved full-solve comparison against the frozen original solver:

| Fixture | Original cold ns | Current cold ns | Cold change | Warm change | Current cold bytes |
| --- | ---: | ---: | --- | --- | ---: |
| scalar_linear | 25.36 | 17.49 | -31.0% (faster) | -31.0% (faster) | 0 |
| scalar_maxwell | 29.45 | 21.45 | -27.2% (faster) | -27.2% (faster) | 0 |
| scalar_powerlaw | 586.96 | 460.87 | -21.5% (faster) | -13.8% (faster) | 0 |
| scalar_creep | 367.86 | 271.43 | -26.2% (faster) | -25.0% (faster) | 0 |
| volumetric_maxwell | 52.86 | 48.57 | -8.1% (faster) | -8.1% (faster) | 0 |
| nested | 75.00 | 78.95 | +5.3% (slower) | +4.9% (inconclusive) | 0 |
| siblings | 141.51 | 138.89 | -1.9% (inconclusive) | -2.4% (inconclusive) | 0 |
| kelvin_voigt | 46.20 | 46.20 | +0.0% (inconclusive) | -1.4% (inconclusive) | 0 |
| plastic | 53.90 | 50.00 | -7.2% (faster) | -7.9% (faster) | 0 |
| dilatant | 109.33 | 86.67 | -20.7% (faster) | -21.7% (faster) | 0 |

Preprocessing experiments below compare Jacobian kernels only, excluding setup.

| Fixture | Cached callable | Cached config | Chunk-1 config | Prepared equations | Config setup ns / bytes |
| --- | --- | --- | --- | --- | ---: |
| scalar_linear | +0.0% (inconclusive) | +8.9% (slower) | +8.9% (slower) | +106.7% (slower) | 5.63 / 32 |
| scalar_maxwell | +0.0% (inconclusive) | +0.0% (inconclusive) | -12.1% (faster) | +345.2% (slower) | 5.35 / 32 |
| scalar_powerlaw | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +1.7% (inconclusive) | 5.97 / 32 |
| scalar_creep | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +10.3% (slower) | 5.13 / 32 |
| volumetric_maxwell | +0.8% (inconclusive) | +0.1% (inconclusive) | +0.0% (inconclusive) | +3.8% (inconclusive) | 5.71 / 64 |
| nested | +0.0% (inconclusive) | +0.0% (inconclusive) | +2.3% (inconclusive) | +141.9% (slower) | 8.25 / 112 |
| siblings | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +1900.0% (slower) | 6.35 / 176 |
| kelvin_voigt | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +193.3% (slower) | 5.87 / 64 |
| plastic | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +156.0% (slower) | 5.78 / 64 |
| dilatant | -3.1% (inconclusive) | -3.1% (inconclusive) | -3.1% (inconclusive) | +82.9% (slower) | 8.22 / 112 |

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
