# Benchmark report: baseline_v2_pinned

Julia 1.13.0; ForwardDiff 1.4.6; AMD Ryzen 7 7800X3D 8-Core Processor           ; threads=1; rounds=3.

Times are medians of round medians, in nanoseconds. Negative changes mean faster.
A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.

| Fixture | N | Residual ns | Jacobian ns | Solve ns | Warm solve ns | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| scalar_linear | 1 | 3.14 | 2.47 | 25.95 | 40.24 | 1 |
| scalar_maxwell | 1 | 4.02 | 4.02 | 46.03 | 46.00 | 1 |
| scalar_powerlaw | 1 | 25.74 | 48.31 | 783.78 | 181.50 | 6 |
| scalar_creep | 1 | 56.56 | 81.94 | 496.61 | 290.20 | 2 |
| volumetric_maxwell | 2 | 3.09 | 12.80 | 53.61 | 53.56 | 1 |
| nested | 3 | 2.99 | 27.56 | 75.68 | 76.06 | 1 |
| siblings | 4 | 3.16 | 36.84 | 140.49 | 140.95 | 1 |
| kelvin_voigt | 2 | 3.12 | 11.59 | 46.63 | 47.01 | 1 |
| plastic | 2 | 3.32 | 16.53 | 53.49 | 53.90 | 1 |
| dilatant | 3 | 20.32 | 33.71 | 106.51 | 106.23 | 1 |

Preprocessing experiments compare Jacobian kernels only, excluding setup. No production solver changes are being compared.

| Fixture | Cached callable | Cached config | Chunk-1 config | Prepared equations | Config setup ns / bytes |
| --- | --- | --- | --- | --- | ---: |
| scalar_linear | +0.0% (inconclusive) | +8.9% (inconclusive) | +0.0% (inconclusive) | +100.0% (slower) | 5.69 / 32 |
| scalar_maxwell | -4.4% (inconclusive) | -4.4% (inconclusive) | -4.4% (inconclusive) | +261.9% (slower) | 7.67 / 32 |
| scalar_powerlaw | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +20.8% (inconclusive) | 7.93 / 32 |
| scalar_creep | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +10.2% (slower) | 7.45 / 32 |
| volumetric_maxwell | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +5.6% (inconclusive) | 5.48 / 64 |
| nested | +0.0% (inconclusive) | +0.3% (inconclusive) | +0.0% (inconclusive) | +128.8% (slower) | 8.31 / 112 |
| siblings | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +2626.6% (slower) | 6.60 / 176 |
| kelvin_voigt | +0.1% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +186.8% (slower) | 6.43 / 64 |
| plastic | +0.0% (inconclusive) | +0.0% (inconclusive) | +0.0% (inconclusive) | +167.2% (slower) | 6.85 / 64 |
| dilatant | +2.8% (inconclusive) | +2.8% (inconclusive) | +2.8% (inconclusive) | +87.2% (slower) | 12.86 / 112 |

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
