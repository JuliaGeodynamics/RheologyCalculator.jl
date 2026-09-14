# Benchmark report: line_search_reuse

Julia 1.13.0; ForwardDiff 1.4.6; AMD Ryzen 7 7800X3D 8-Core Processor           ; threads=1; rounds=3.

Times are medians of round medians, in nanoseconds. Negative changes mean faster.
A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.

| Fixture | N | Residual ns | Jacobian ns | Solve ns | Warm solve ns | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| scalar_linear | 1 | 2.67 | 2.70 | 26.39 | 26.39 | 1 |
| scalar_maxwell | 1 | 2.92 | 2.98 | 33.85 | 34.38 | 1 |
| scalar_powerlaw | 1 | 17.28 | 34.16 | 547.62 | 119.05 | 6 |
| scalar_creep | 1 | 39.51 | 58.82 | 279.17 | 172.41 | 2 |
| volumetric_maxwell | 2 | 3.06 | 12.84 | 50.71 | 51.43 | 1 |
| nested | 3 | 3.01 | 27.33 | 79.35 | 79.35 | 1 |
| siblings | 4 | 3.13 | 41.67 | 146.67 | 146.67 | 1 |
| kelvin_voigt | 2 | 3.10 | 11.40 | 51.37 | 51.70 | 1 |
| plastic | 2 | 3.31 | 23.21 | 59.35 | 61.79 | 1 |
| dilatant | 3 | 19.28 | 34.48 | 95.65 | 95.65 | 1 |

Interleaved full-solve comparison against the frozen original solver:

| Fixture | Original cold ns | Current cold ns | Cold change | Warm change | Current cold bytes |
| --- | ---: | ---: | --- | --- | ---: |
| scalar_linear | 26.04 | 26.39 | +1.3% (inconclusive) | +1.3% (inconclusive) | 0 |
| scalar_maxwell | 33.85 | 33.85 | +0.0% (inconclusive) | +8.2% (inconclusive) | 0 |
| scalar_powerlaw | 590.48 | 547.62 | -7.3% (faster) | -10.7% (faster) | 0 |
| scalar_creep | 358.33 | 279.17 | -22.1% (faster) | -20.4% (faster) | 0 |
| volumetric_maxwell | 52.86 | 50.71 | -4.1% (inconclusive) | -2.7% (inconclusive) | 0 |
| nested | 76.09 | 79.35 | +4.3% (inconclusive) | +4.3% (inconclusive) | 0 |
| siblings | 143.33 | 146.67 | +2.3% (inconclusive) | +2.3% (inconclusive) | 0 |
| kelvin_voigt | 47.26 | 51.37 | +8.7% (inconclusive) | +9.4% (inconclusive) | 0 |
| plastic | 60.98 | 59.35 | -2.7% (inconclusive) | +15.0% (slower) | 0 |
| dilatant | 104.29 | 95.65 | -8.3% (faster) | -8.3% (faster) | 0 |

Preprocessing experiments below compare Jacobian kernels only, excluding setup.

| Fixture | Cached callable | Cached config | Chunk-1 config | Prepared equations | Config setup ns / bytes |
| --- | --- | --- | --- | --- | ---: |
| scalar_linear | +0.0% (inconclusive) | -8.2% (faster) | -8.2% (faster) | +102.2% (slower) | 5.66 / 32 |
| scalar_maxwell | -7.7% (inconclusive) | -7.7% (inconclusive) | +0.0% (inconclusive) | +347.1% (slower) | 7.49 / 32 |
| scalar_powerlaw | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +1.0% (inconclusive) | 5.35 / 32 |
| scalar_creep | +1.6% (inconclusive) | +1.6% (inconclusive) | +1.6% (inconclusive) | +12.2% (slower) | 5.30 / 32 |
| volumetric_maxwell | -1.3% (inconclusive) | -1.3% (inconclusive) | -1.3% (inconclusive) | +3.9% (inconclusive) | 5.73 / 64 |
| nested | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +130.1% (slower) | 8.21 / 112 |
| siblings | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +1910.0% (slower) | 6.50 / 176 |
| kelvin_voigt | +0.2% (inconclusive) | +0.2% (inconclusive) | +0.2% (inconclusive) | +208.4% (slower) | 6.67 / 64 |
| plastic | +1.7% (inconclusive) | +1.7% (inconclusive) | +1.7% (inconclusive) | +133.3% (slower) | 6.69 / 64 |
| dilatant | -4.0% (inconclusive) | -4.0% (inconclusive) | -4.0% (inconclusive) | +96.1% (slower) | 8.39 / 112 |

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
