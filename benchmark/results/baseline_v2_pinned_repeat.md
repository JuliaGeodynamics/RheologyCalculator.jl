# Benchmark report: baseline_v2_pinned_repeat

Julia 1.13.0; ForwardDiff 1.4.6; AMD Ryzen 7 7800X3D 8-Core Processor           ; threads=1; rounds=3.

Times are medians of round medians, in nanoseconds. Negative changes mean faster.
A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.

| Fixture | N | Residual ns | Jacobian ns | Solve ns | Warm solve ns | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| scalar_linear | 1 | 2.67 | 2.69 | 25.97 | 25.37 | 1 |
| scalar_maxwell | 1 | 2.67 | 2.68 | 29.51 | 29.59 | 1 |
| scalar_powerlaw | 1 | 17.35 | 34.13 | 597.96 | 131.31 | 6 |
| scalar_creep | 1 | 39.74 | 59.77 | 361.33 | 219.49 | 2 |
| volumetric_maxwell | 2 | 3.07 | 12.76 | 53.51 | 53.37 | 1 |
| nested | 3 | 3.01 | 27.86 | 77.01 | 77.02 | 1 |
| siblings | 4 | 3.21 | 38.10 | 143.43 | 144.15 | 1 |
| kelvin_voigt | 2 | 3.09 | 11.46 | 46.33 | 46.46 | 1 |
| plastic | 2 | 3.47 | 16.38 | 53.56 | 53.65 | 1 |
| dilatant | 3 | 19.17 | 40.00 | 103.20 | 103.21 | 1 |

Preprocessing experiments compare Jacobian kernels only, excluding setup. No production solver changes are being compared.

| Fixture | Cached callable | Cached config | Chunk-1 config | Prepared equations | Config setup ns / bytes |
| --- | --- | --- | --- | --- | ---: |
| scalar_linear | -8.5% (faster) | -15.2% (faster) | -8.5% (faster) | +98.5% (slower) | 5.44 / 32 |
| scalar_maxwell | +6.5% (slower) | +0.0% (inconclusive) | +6.5% (slower) | +371.9% (slower) | 5.20 / 32 |
| scalar_powerlaw | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +1.1% (inconclusive) | 5.67 / 32 |
| scalar_creep | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +14.6% (slower) | 5.49 / 32 |
| volumetric_maxwell | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +6.0% (slower) | 5.90 / 64 |
| nested | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +128.4% (slower) | 10.59 / 112 |
| siblings | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +2717.5% (slower) | 6.79 / 176 |
| kelvin_voigt | +1.3% (inconclusive) | +1.3% (inconclusive) | +2.5% (inconclusive) | +193.5% (slower) | 6.40 / 64 |
| plastic | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +185.2% (slower) | 5.82 / 64 |
| dilatant | +3.3% (inconclusive) | +3.3% (inconclusive) | +3.3% (inconclusive) | +109.5% (slower) | 8.31 / 112 |

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
