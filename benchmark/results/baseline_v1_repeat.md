# Benchmark report: baseline_v1_repeat

Julia 1.13.0; ForwardDiff 1.4.6; AMD Ryzen 7 7800X3D 8-Core Processor           ; threads=1; rounds=3.

Times are medians of round medians, in nanoseconds. Negative changes mean faster.
A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.

| Fixture | N | Residual ns | Jacobian ns | Solve ns | Warm solve ns | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| scalar_linear | 1 | 2.72 | 2.50 | 29.18 | 27.04 | 1 |
| scalar_maxwell | 1 | 3.95 | 3.39 | 46.35 | 32.83 | 1 |
| scalar_powerlaw | 1 | 25.86 | 46.03 | 786.49 | 181.60 | 6 |
| scalar_creep | 1 | 40.59 | 60.96 | 494.74 | 272.80 | 2 |
| volumetric_maxwell | 2 | 4.18 | 22.01 | 67.59 | 67.66 | 1 |
| nested | 3 | 4.83 | 35.29 | 99.31 | 99.32 | 1 |
| siblings | 4 | 4.52 | 54.90 | 180.86 | 184.11 | 1 |
| kelvin_voigt | 2 | 4.18 | 19.58 | 61.52 | 47.65 | 1 |
| plastic | 2 | 3.34 | 23.08 | 62.70 | 62.50 | 1 |
| dilatant | 3 | 19.32 | 33.42 | 104.72 | 104.93 | 1 |

Preprocessing experiments compare Jacobian kernels only, excluding setup. No production solver changes are being compared.

| Fixture | Cached callable | Cached config | Chunk-1 config | Prepared equations | Config setup ns / bytes |
| --- | --- | --- | --- | --- | ---: |
| scalar_linear | +9.4% (slower) | +21.6% (slower) | -0.1% (inconclusive) | +96.0% (slower) | 5.63 / 32 |
| scalar_maxwell | +5.2% (inconclusive) | +0.4% (inconclusive) | -19.2% (inconclusive) | +331.6% (slower) | 7.29 / 32 |
| scalar_powerlaw | +4.0% (inconclusive) | +4.9% (inconclusive) | -24.1% (faster) | +27.5% (inconclusive) | 7.42 / 32 |
| scalar_creep | +5.4% (inconclusive) | -0.7% (inconclusive) | -1.2% (inconclusive) | +30.8% (slower) | 5.69 / 32 |
| volumetric_maxwell | -0.3% (inconclusive) | -1.2% (inconclusive) | +0.6% (inconclusive) | +8.5% (slower) | 8.99 / 64 |
| nested | -0.9% (inconclusive) | +7.2% (inconclusive) | -0.4% (inconclusive) | +73.7% (slower) | 11.07 / 112 |
| siblings | +0.1% (inconclusive) | +0.1% (inconclusive) | +0.1% (inconclusive) | +2318.9% (slower) | 8.83 / 176 |
| kelvin_voigt | -40.3% (inconclusive) | -40.3% (inconclusive) | -31.0% (inconclusive) | +183.6% (slower) | 7.66 / 64 |
| plastic | +1.4% (inconclusive) | +1.2% (inconclusive) | +1.5% (inconclusive) | +120.1% (slower) | 5.81 / 64 |
| dilatant | +1.6% (inconclusive) | +1.4% (inconclusive) | +1.6% (inconclusive) | +65.7% (slower) | 7.96 / 112 |

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
