# Benchmark report: baseline_v1

Julia 1.13.0; ForwardDiff 1.4.6; AMD Ryzen 7 7800X3D 8-Core Processor           ; threads=1; rounds=3.

Times are medians of round medians, in nanoseconds. Negative changes mean faster.
A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.

| Fixture | N | Residual ns | Jacobian ns | Solve ns | Warm solve ns | Iterations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| scalar_linear | 1 | 2.68 | 2.47 | 25.59 | 25.39 | 1 |
| scalar_maxwell | 1 | 2.88 | 2.68 | 29.97 | 30.04 | 1 |
| scalar_powerlaw | 1 | 17.54 | 34.15 | 602.33 | 133.17 | 6 |
| scalar_creep | 1 | 39.44 | 59.51 | 366.67 | 226.23 | 2 |
| volumetric_maxwell | 2 | 3.10 | 17.77 | 57.71 | 58.16 | 1 |
| nested | 3 | 2.94 | 39.12 | 76.99 | 77.60 | 1 |
| siblings | 4 | 4.84 | 50.80 | 175.98 | 178.11 | 1 |
| kelvin_voigt | 2 | 4.96 | 21.47 | 71.06 | 71.07 | 1 |
| plastic | 2 | 3.43 | 24.10 | 54.67 | 55.83 | 1 |
| dilatant | 3 | 19.66 | 43.09 | 131.80 | 119.65 | 1 |

Preprocessing experiments compare Jacobian kernels only, excluding setup. No production solver changes are being compared.

| Fixture | Cached callable | Cached config | Chunk-1 config | Prepared equations | Config setup ns / bytes |
| --- | --- | --- | --- | --- | ---: |
| scalar_linear | +8.4% (slower) | +8.4% (slower) | +0.3% (inconclusive) | +114.7% (slower) | 5.39 / 32 |
| scalar_maxwell | +0.4% (inconclusive) | +0.0% (inconclusive) | -0.3% (inconclusive) | +379.1% (slower) | 4.83 / 32 |
| scalar_powerlaw | -0.4% (inconclusive) | +0.0% (inconclusive) | +0.4% (inconclusive) | +0.7% (inconclusive) | 5.36 / 32 |
| scalar_creep | +0.2% (inconclusive) | -0.4% (inconclusive) | +0.0% (inconclusive) | +11.2% (slower) | 4.78 / 32 |
| volumetric_maxwell | +0.8% (inconclusive) | +0.6% (inconclusive) | +1.1% (inconclusive) | +10.8% (slower) | 5.51 / 64 |
| nested | -0.9% (inconclusive) | -0.7% (inconclusive) | -0.3% (inconclusive) | +50.7% (slower) | 8.40 / 112 |
| siblings | -2.5% (inconclusive) | +1.5% (inconclusive) | +0.9% (inconclusive) | +2557.7% (slower) | 10.17 / 176 |
| kelvin_voigt | +1.3% (inconclusive) | +1.8% (inconclusive) | +0.9% (inconclusive) | +185.6% (slower) | 8.47 / 64 |
| plastic | -32.7% (inconclusive) | -32.6% (inconclusive) | +0.7% (inconclusive) | +178.2% (slower) | 6.51 / 64 |
| dilatant | +0.1% (inconclusive) | -20.6% (inconclusive) | +4.2% (inconclusive) | +97.2% (slower) | 9.90 / 112 |

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
