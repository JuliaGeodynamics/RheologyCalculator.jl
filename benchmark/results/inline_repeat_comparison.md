# Benchmark report: line_search_reuse_inline_repeat

Julia 1.13.0; ForwardDiff 1.4.6; AMD Ryzen 7 7800X3D 8-Core Processor           ; threads=1; rounds=3.

Times are medians of round medians, in nanoseconds. Negative changes mean faster.
A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.

| Fixture | Kernel | Baseline ns | Candidate ns | Change | Bytes before → after |
| --- | --- | ---: | ---: | --- | ---: |
| dilatant | backsolve | 28.36 | 28.30 | -0.2% (inconclusive) | 0 → 0 |
| dilatant | cached_callable | 33.58 | 33.33 | -0.7% (inconclusive) | 0 → 0 |
| dilatant | cached_config | 33.58 | 33.33 | -0.7% (inconclusive) | 0 → 0 |
| dilatant | cached_config_chunk1 | 33.58 | 33.33 | -0.7% (inconclusive) | 0 → 0 |
| dilatant | config_setup | 8.22 | 9.78 | +18.9% (inconclusive) | 112 → 112 |
| dilatant | equations_setup | 5.33 | 5.35 | +0.5% (inconclusive) | 0 → 0 |
| dilatant | jacobian | 34.65 | 34.48 | -0.5% (inconclusive) | 0 → 0 |
| dilatant | prepared_equations | 63.37 | 66.67 | +5.2% (inconclusive) | 592 → 592 |
| dilatant | residual | 19.01 | 30.05 | +58.1% (slower) | 0 → 0 |
| dilatant | residual_jacobian | 33.61 | 33.77 | +0.5% (inconclusive) | 0 → 0 |
| dilatant | solve | 86.67 | 86.49 | -0.2% (inconclusive) | 0 → 0 |
| dilatant | solve_reference | 109.33 | 109.59 | +0.2% (inconclusive) | 0 → 0 |
| dilatant | solve_warm | 86.67 | 87.67 | +1.2% (inconclusive) | 0 → 0 |
| dilatant | solve_warm_reference | 110.67 | 110.96 | +0.3% (inconclusive) | 0 → 0 |
| kelvin_voigt | backsolve | 2.86 | 2.86 | -0.1% (inconclusive) | 0 → 0 |
| kelvin_voigt | cached_callable | 11.36 | 11.51 | +1.3% (inconclusive) | 0 → 0 |
| kelvin_voigt | cached_config | 11.36 | 11.40 | +0.3% (inconclusive) | 0 → 0 |
| kelvin_voigt | cached_config_chunk1 | 11.36 | 11.40 | +0.3% (inconclusive) | 0 → 0 |
| kelvin_voigt | config_setup | 5.87 | 7.58 | +29.2% (inconclusive) | 64 → 64 |
| kelvin_voigt | equations_setup | 2.67 | 2.87 | +7.3% (slower) | 0 → 0 |
| kelvin_voigt | jacobian | 11.36 | 11.24 | -1.1% (inconclusive) | 0 → 0 |
| kelvin_voigt | prepared_equations | 33.33 | 33.99 | +2.0% (inconclusive) | 176 → 176 |
| kelvin_voigt | residual | 3.12 | 3.21 | +2.8% (inconclusive) | 0 → 0 |
| kelvin_voigt | residual_jacobian | 17.50 | 17.48 | -0.1% (inconclusive) | 0 → 0 |
| kelvin_voigt | solve | 46.20 | 46.20 | +0.0% (inconclusive) | 0 → 0 |
| kelvin_voigt | solve_reference | 46.20 | 46.20 | +0.0% (inconclusive) | 0 → 0 |
| kelvin_voigt | solve_warm | 45.86 | 46.20 | +0.7% (inconclusive) | 0 → 0 |
| kelvin_voigt | solve_warm_reference | 46.50 | 46.20 | -0.6% (inconclusive) | 0 → 0 |
| nested | backsolve | 23.21 | 23.68 | +2.0% (inconclusive) | 0 → 0 |
| nested | cached_callable | 27.03 | 27.56 | +2.0% (inconclusive) | 0 → 0 |
| nested | cached_config | 27.03 | 27.56 | +2.0% (inconclusive) | 0 → 0 |
| nested | cached_config_chunk1 | 27.66 | 26.92 | -2.7% (inconclusive) | 0 → 0 |
| nested | config_setup | 8.25 | 8.75 | +6.1% (inconclusive) | 112 → 112 |
| nested | equations_setup | 3.09 | 3.08 | -0.2% (inconclusive) | 0 → 0 |
| nested | jacobian | 27.03 | 27.56 | +2.0% (inconclusive) | 0 → 0 |
| nested | prepared_equations | 65.38 | 64.29 | -1.7% (inconclusive) | 448 → 448 |
| nested | residual | 2.97 | 2.92 | -1.9% (inconclusive) | 0 → 0 |
| nested | residual_jacobian | 30.00 | 29.96 | -0.1% (inconclusive) | 0 → 0 |
| nested | solve | 78.95 | 78.12 | -1.0% (inconclusive) | 0 → 0 |
| nested | solve_reference | 75.00 | 75.34 | +0.5% (inconclusive) | 0 → 0 |
| nested | solve_warm | 78.95 | 79.45 | +0.6% (inconclusive) | 0 → 0 |
| nested | solve_warm_reference | 75.28 | 76.04 | +1.0% (inconclusive) | 0 → 0 |
| plastic | backsolve | 2.86 | 2.86 | -0.1% (inconclusive) | 0 → 0 |
| plastic | cached_callable | 16.67 | 16.31 | -2.1% (inconclusive) | 0 → 0 |
| plastic | cached_config | 16.67 | 16.31 | -2.1% (inconclusive) | 0 → 0 |
| plastic | cached_config_chunk1 | 16.67 | 16.31 | -2.1% (inconclusive) | 0 → 0 |
| plastic | config_setup | 5.78 | 6.15 | +6.4% (inconclusive) | 64 → 64 |
| plastic | equations_setup | 4.35 | 4.41 | +1.4% (inconclusive) | 0 → 0 |
| plastic | jacobian | 16.67 | 16.31 | -2.1% (inconclusive) | 0 → 0 |
| plastic | prepared_equations | 42.67 | 44.08 | +3.3% (inconclusive) | 416 → 416 |
| plastic | residual | 3.32 | 3.13 | -5.5% (inconclusive) | 0 → 0 |
| plastic | residual_jacobian | 20.00 | 23.39 | +16.9% (slower) | 0 → 0 |
| plastic | solve | 50.00 | 50.74 | +1.5% (inconclusive) | 0 → 0 |
| plastic | solve_reference | 53.90 | 53.68 | -0.4% (inconclusive) | 0 → 0 |
| plastic | solve_warm | 49.65 | 50.74 | +2.2% (inconclusive) | 0 → 0 |
| plastic | solve_warm_reference | 53.90 | 53.68 | -0.4% (inconclusive) | 0 → 0 |
| scalar_creep | backsolve | 2.26 | 2.45 | +8.2% (slower) | 0 → 0 |
| scalar_creep | cached_callable | 58.59 | 58.59 | +0.0% (inconclusive) | 0 → 0 |
| scalar_creep | cached_config | 58.59 | 58.59 | +0.0% (inconclusive) | 0 → 0 |
| scalar_creep | cached_config_chunk1 | 58.59 | 58.59 | +0.0% (inconclusive) | 0 → 0 |
| scalar_creep | config_setup | 5.13 | 4.90 | -4.4% (inconclusive) | 32 → 32 |
| scalar_creep | equations_setup | 2.87 | 2.86 | -0.3% (inconclusive) | 0 → 0 |
| scalar_creep | jacobian | 58.59 | 58.95 | +0.6% (inconclusive) | 0 → 0 |
| scalar_creep | prepared_equations | 64.65 | 65.26 | +1.0% (inconclusive) | 0 → 0 |
| scalar_creep | residual | 38.62 | 38.84 | +0.6% (inconclusive) | 0 → 0 |
| scalar_creep | residual_jacobian | 97.69 | 96.33 | -1.4% (inconclusive) | 0 → 0 |
| scalar_creep | solve | 271.43 | 271.43 | +0.0% (inconclusive) | 0 → 0 |
| scalar_creep | solve_reference | 367.86 | 367.86 | +0.0% (inconclusive) | 0 → 0 |
| scalar_creep | solve_warm | 171.43 | 171.43 | +0.0% (inconclusive) | 0 → 0 |
| scalar_creep | solve_warm_reference | 228.57 | 228.57 | +0.0% (inconclusive) | 0 → 0 |
| scalar_linear | backsolve | 2.46 | 2.66 | +7.8% (slower) | 0 → 0 |
| scalar_linear | cached_callable | 2.45 | 2.47 | +0.5% (inconclusive) | 0 → 0 |
| scalar_linear | cached_config | 2.67 | 2.63 | -1.5% (inconclusive) | 0 → 0 |
| scalar_linear | cached_config_chunk1 | 2.67 | 2.47 | -7.7% (faster) | 0 → 0 |
| scalar_linear | config_setup | 5.63 | 5.18 | -8.0% (faster) | 32 → 32 |
| scalar_linear | equations_setup | 2.66 | 2.65 | -0.4% (inconclusive) | 0 → 0 |
| scalar_linear | jacobian | 2.45 | 2.25 | -8.1% (faster) | 0 → 0 |
| scalar_linear | prepared_equations | 5.07 | 5.32 | +4.9% (inconclusive) | 0 → 0 |
| scalar_linear | residual | 2.46 | 2.45 | -0.4% (inconclusive) | 0 → 0 |
| scalar_linear | residual_jacobian | 2.68 | 2.66 | -0.6% (inconclusive) | 0 → 0 |
| scalar_linear | solve | 17.49 | 17.56 | +0.4% (inconclusive) | 0 → 0 |
| scalar_linear | solve_reference | 25.36 | 25.22 | -0.6% (inconclusive) | 0 → 0 |
| scalar_linear | solve_warm | 17.49 | 17.78 | +1.7% (inconclusive) | 0 → 0 |
| scalar_linear | solve_warm_reference | 25.36 | 25.00 | -1.4% (inconclusive) | 0 → 0 |
| scalar_maxwell | backsolve | 2.26 | 2.65 | +17.4% (slower) | 0 → 0 |
| scalar_maxwell | cached_callable | 2.83 | 2.65 | -6.3% (faster) | 0 → 0 |
| scalar_maxwell | cached_config | 2.83 | 2.82 | -0.1% (inconclusive) | 0 → 0 |
| scalar_maxwell | cached_config_chunk1 | 2.48 | 2.48 | -0.1% (inconclusive) | 0 → 0 |
| scalar_maxwell | config_setup | 5.35 | 4.99 | -6.7% (inconclusive) | 32 → 32 |
| scalar_maxwell | equations_setup | 2.67 | 2.65 | -0.8% (inconclusive) | 0 → 0 |
| scalar_maxwell | jacobian | 2.83 | 2.65 | -6.3% (faster) | 0 → 0 |
| scalar_maxwell | prepared_equations | 12.58 | 12.54 | -0.3% (inconclusive) | 0 → 0 |
| scalar_maxwell | residual | 2.87 | 2.86 | -0.4% (inconclusive) | 0 → 0 |
| scalar_maxwell | residual_jacobian | 4.04 | 3.69 | -8.7% (faster) | 0 → 0 |
| scalar_maxwell | solve | 21.45 | 21.51 | +0.2% (inconclusive) | 0 → 0 |
| scalar_maxwell | solve_reference | 29.45 | 29.33 | -0.4% (inconclusive) | 0 → 0 |
| scalar_maxwell | solve_warm | 21.45 | 21.51 | +0.2% (inconclusive) | 0 → 0 |
| scalar_maxwell | solve_warm_reference | 29.45 | 29.33 | -0.4% (inconclusive) | 0 → 0 |
| scalar_powerlaw | backsolve | 2.66 | 2.45 | -7.8% (faster) | 0 → 0 |
| scalar_powerlaw | cached_callable | 33.71 | 33.74 | +0.1% (inconclusive) | 0 → 0 |
| scalar_powerlaw | cached_config | 33.71 | 33.74 | +0.1% (inconclusive) | 0 → 0 |
| scalar_powerlaw | cached_config_chunk1 | 33.71 | 33.74 | +0.1% (inconclusive) | 0 → 0 |
| scalar_powerlaw | config_setup | 5.97 | 5.67 | -5.1% (faster) | 32 → 32 |
| scalar_powerlaw | equations_setup | 2.67 | 2.66 | -0.3% (inconclusive) | 0 → 0 |
| scalar_powerlaw | jacobian | 33.71 | 33.74 | +0.1% (inconclusive) | 0 → 0 |
| scalar_powerlaw | prepared_equations | 34.29 | 34.16 | -0.4% (inconclusive) | 0 → 0 |
| scalar_powerlaw | residual | 17.30 | 17.12 | -1.1% (inconclusive) | 0 → 0 |
| scalar_powerlaw | residual_jacobian | 55.47 | 55.56 | +0.2% (inconclusive) | 0 → 0 |
| scalar_powerlaw | solve | 460.87 | 460.87 | +0.0% (inconclusive) | 0 → 0 |
| scalar_powerlaw | solve_reference | 586.96 | 586.96 | +0.0% (inconclusive) | 0 → 0 |
| scalar_powerlaw | solve_warm | 113.64 | 113.64 | +0.0% (inconclusive) | 0 → 0 |
| scalar_powerlaw | solve_warm_reference | 131.82 | 131.82 | +0.0% (inconclusive) | 0 → 0 |
| siblings | backsolve | 77.51 | 82.85 | +6.9% (slower) | 0 → 0 |
| siblings | cached_callable | 41.67 | 50.00 | +20.0% (slower) | 0 → 0 |
| siblings | cached_config | 41.67 | 50.00 | +20.0% (slower) | 0 → 0 |
| siblings | cached_config_chunk1 | 41.67 | 50.00 | +20.0% (slower) | 0 → 0 |
| siblings | config_setup | 6.35 | 6.79 | +6.9% (inconclusive) | 176 → 176 |
| siblings | equations_setup | 3.44 | 3.44 | +0.2% (inconclusive) | 0 → 0 |
| siblings | jacobian | 41.67 | 50.00 | +20.0% (slower) | 0 → 0 |
| siblings | prepared_equations | 833.33 | 883.33 | +6.0% (slower) | 3184 → 3184 |
| siblings | residual | 3.10 | 3.11 | +0.3% (inconclusive) | 0 → 0 |
| siblings | residual_jacobian | 46.39 | 50.32 | +8.5% (slower) | 0 → 0 |
| siblings | solve | 138.89 | 141.18 | +1.6% (inconclusive) | 0 → 0 |
| siblings | solve_reference | 141.51 | 143.14 | +1.2% (inconclusive) | 0 → 0 |
| siblings | solve_warm | 137.74 | 139.22 | +1.1% (inconclusive) | 0 → 0 |
| siblings | solve_warm_reference | 141.18 | 143.14 | +1.4% (inconclusive) | 0 → 0 |
| volumetric_maxwell | backsolve | 2.87 | 2.86 | -0.6% (inconclusive) | 0 → 0 |
| volumetric_maxwell | cached_callable | 12.75 | 12.26 | -3.8% (inconclusive) | 0 → 0 |
| volumetric_maxwell | cached_config | 12.67 | 12.26 | -3.2% (inconclusive) | 0 → 0 |
| volumetric_maxwell | cached_config_chunk1 | 12.65 | 12.26 | -3.1% (inconclusive) | 0 → 0 |
| volumetric_maxwell | config_setup | 5.71 | 5.85 | +2.4% (inconclusive) | 64 → 64 |
| volumetric_maxwell | equations_setup | 2.88 | 2.86 | -0.5% (inconclusive) | 0 → 0 |
| volumetric_maxwell | jacobian | 12.65 | 12.34 | -2.5% (inconclusive) | 0 → 0 |
| volumetric_maxwell | prepared_equations | 13.14 | 13.11 | -0.2% (inconclusive) | 0 → 0 |
| volumetric_maxwell | residual | 3.07 | 3.06 | -0.3% (inconclusive) | 0 → 0 |
| volumetric_maxwell | residual_jacobian | 28.56 | 13.60 | -52.4% (faster) | 0 → 0 |
| volumetric_maxwell | solve | 48.57 | 48.28 | -0.6% (inconclusive) | 0 → 0 |
| volumetric_maxwell | solve_reference | 52.86 | 52.41 | -0.8% (inconclusive) | 0 → 0 |
| volumetric_maxwell | solve_warm | 48.57 | 48.28 | -0.6% (inconclusive) | 0 → 0 |
| volumetric_maxwell | solve_warm_reference | 52.86 | 52.41 | -0.8% (inconclusive) | 0 → 0 |

| Fixture | Iterations before → after | Residual before → after | Max absolute solution change |
| --- | ---: | ---: | ---: |
| scalar_linear | 1 → 1 | 0.000e+00 → 0.000e+00 | 0.000e+00 |
| scalar_maxwell | 1 → 1 | 0.000e+00 → 0.000e+00 | 0.000e+00 |
| scalar_powerlaw | 6 → 6 | 0.000e+00 → 0.000e+00 | 0.000e+00 |
| scalar_creep | 2 → 2 | 0.000e+00 → 0.000e+00 | 0.000e+00 |
| volumetric_maxwell | 1 → 1 | 1.505e-21 → 1.505e-21 | 0.000e+00 |
| nested | 1 → 1 | 1.455e-17 → 1.455e-17 | 0.000e+00 |
| siblings | 1 → 1 | 4.191e-16 → 4.191e-16 | 0.000e+00 |
| kelvin_voigt | 1 → 1 | 2.910e-17 → 2.910e-17 | 0.000e+00 |
| plastic | 1 → 1 | 2.742e-16 → 2.742e-16 | 0.000e+00 |
| dilatant | 1 → 1 | 2.767e-16 → 2.767e-16 | 0.000e+00 |
