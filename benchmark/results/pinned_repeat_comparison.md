# Benchmark report: baseline_v2_pinned_repeat

Julia 1.13.0; ForwardDiff 1.4.6; AMD Ryzen 7 7800X3D 8-Core Processor           ; threads=1; rounds=3.

Times are medians of round medians, in nanoseconds. Negative changes mean faster.
A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.

| Fixture | Kernel | Baseline ns | Candidate ns | Change | Bytes before → after |
| --- | --- | ---: | ---: | --- | ---: |
| dilatant | backsolve | 37.93 | 28.16 | -25.8% (faster) | 0 → 0 |
| dilatant | cached_callable | 34.65 | 41.33 | +19.3% (slower) | 0 → 0 |
| dilatant | cached_config | 34.65 | 41.33 | +19.3% (slower) | 0 → 0 |
| dilatant | cached_config_chunk1 | 34.65 | 41.33 | +19.3% (slower) | 0 → 0 |
| dilatant | config_setup | 12.86 | 8.31 | -35.4% (faster) | 112 → 112 |
| dilatant | equations_setup | 5.32 | 5.33 | +0.2% (inconclusive) | 0 → 0 |
| dilatant | jacobian | 33.71 | 40.00 | +18.7% (slower) | 0 → 0 |
| dilatant | prepared_equations | 63.11 | 83.78 | +32.7% (slower) | 592 → 592 |
| dilatant | residual | 20.32 | 19.17 | -5.6% (faster) | 0 → 0 |
| dilatant | residual_jacobian | 33.86 | 43.91 | +29.7% (slower) | 0 → 0 |
| dilatant | solve | 106.51 | 103.20 | -3.1% (inconclusive) | 0 → 0 |
| dilatant | solve_warm | 106.23 | 103.21 | -2.8% (inconclusive) | 0 → 0 |
| kelvin_voigt | backsolve | 2.89 | 2.88 | -0.5% (inconclusive) | 0 → 0 |
| kelvin_voigt | cached_callable | 11.59 | 11.61 | +0.1% (inconclusive) | 0 → 0 |
| kelvin_voigt | cached_config | 11.59 | 11.61 | +0.2% (inconclusive) | 0 → 0 |
| kelvin_voigt | cached_config_chunk1 | 11.59 | 11.75 | +1.4% (inconclusive) | 0 → 0 |
| kelvin_voigt | config_setup | 6.43 | 6.40 | -0.5% (inconclusive) | 64 → 64 |
| kelvin_voigt | equations_setup | 2.89 | 2.88 | -0.3% (inconclusive) | 0 → 0 |
| kelvin_voigt | jacobian | 11.59 | 11.46 | -1.1% (inconclusive) | 0 → 0 |
| kelvin_voigt | prepared_equations | 33.23 | 33.63 | +1.2% (inconclusive) | 176 → 176 |
| kelvin_voigt | residual | 3.12 | 3.09 | -1.0% (inconclusive) | 0 → 0 |
| kelvin_voigt | residual_jacobian | 17.97 | 18.08 | +0.6% (inconclusive) | 0 → 0 |
| kelvin_voigt | solve | 46.63 | 46.33 | -0.6% (inconclusive) | 0 → 0 |
| kelvin_voigt | solve_warm | 47.01 | 46.46 | -1.2% (inconclusive) | 0 → 0 |
| nested | backsolve | 23.68 | 23.52 | -0.7% (inconclusive) | 0 → 0 |
| nested | cached_callable | 27.56 | 27.86 | +1.1% (inconclusive) | 0 → 0 |
| nested | cached_config | 27.64 | 27.86 | +0.8% (inconclusive) | 0 → 0 |
| nested | cached_config_chunk1 | 27.56 | 27.86 | +1.1% (inconclusive) | 0 → 0 |
| nested | config_setup | 8.31 | 10.59 | +27.4% (inconclusive) | 112 → 112 |
| nested | equations_setup | 3.08 | 3.09 | +0.3% (inconclusive) | 0 → 0 |
| nested | jacobian | 27.56 | 27.86 | +1.1% (inconclusive) | 0 → 0 |
| nested | prepared_equations | 63.06 | 63.64 | +0.9% (inconclusive) | 448 → 448 |
| nested | residual | 2.99 | 3.01 | +0.6% (inconclusive) | 0 → 0 |
| nested | residual_jacobian | 31.24 | 31.41 | +0.6% (inconclusive) | 0 → 0 |
| nested | solve | 75.68 | 77.01 | +1.8% (inconclusive) | 0 → 0 |
| nested | solve_warm | 76.06 | 77.02 | +1.3% (inconclusive) | 0 → 0 |
| plastic | backsolve | 2.89 | 2.67 | -7.4% (faster) | 0 → 0 |
| plastic | cached_callable | 16.53 | 16.38 | -0.9% (inconclusive) | 0 → 0 |
| plastic | cached_config | 16.53 | 16.38 | -0.9% (inconclusive) | 0 → 0 |
| plastic | cached_config_chunk1 | 16.53 | 16.38 | -0.9% (inconclusive) | 0 → 0 |
| plastic | config_setup | 6.85 | 5.82 | -15.0% (inconclusive) | 64 → 64 |
| plastic | equations_setup | 3.77 | 3.73 | -1.2% (inconclusive) | 0 → 0 |
| plastic | jacobian | 16.53 | 16.38 | -0.9% (inconclusive) | 0 → 0 |
| plastic | prepared_equations | 44.17 | 46.72 | +5.8% (slower) | 416 → 416 |
| plastic | residual | 3.32 | 3.47 | +4.5% (inconclusive) | 0 → 0 |
| plastic | residual_jacobian | 16.95 | 19.47 | +14.8% (slower) | 0 → 0 |
| plastic | solve | 53.49 | 53.56 | +0.1% (inconclusive) | 0 → 0 |
| plastic | solve_warm | 53.90 | 53.65 | -0.5% (inconclusive) | 0 → 0 |
| scalar_creep | backsolve | 3.36 | 2.49 | -25.8% (faster) | 0 → 0 |
| scalar_creep | cached_callable | 81.94 | 59.77 | -27.1% (faster) | 0 → 0 |
| scalar_creep | cached_config | 81.94 | 59.77 | -27.1% (faster) | 0 → 0 |
| scalar_creep | cached_config_chunk1 | 81.94 | 59.77 | -27.1% (faster) | 0 → 0 |
| scalar_creep | config_setup | 7.45 | 5.49 | -26.3% (faster) | 32 → 32 |
| scalar_creep | equations_setup | 4.15 | 2.69 | -35.2% (faster) | 0 → 0 |
| scalar_creep | jacobian | 81.94 | 59.77 | -27.1% (faster) | 0 → 0 |
| scalar_creep | prepared_equations | 90.28 | 68.48 | -24.1% (faster) | 0 → 0 |
| scalar_creep | residual | 56.56 | 39.74 | -29.7% (faster) | 0 → 0 |
| scalar_creep | residual_jacobian | 137.50 | 98.04 | -28.7% (faster) | 0 → 0 |
| scalar_creep | solve | 496.61 | 361.33 | -27.2% (faster) | 0 → 0 |
| scalar_creep | solve_warm | 290.20 | 219.49 | -24.4% (faster) | 0 → 0 |
| scalar_linear | backsolve | 2.51 | 2.51 | -0.2% (inconclusive) | 0 → 0 |
| scalar_linear | cached_callable | 2.47 | 2.46 | -0.2% (inconclusive) | 0 → 0 |
| scalar_linear | cached_config | 2.69 | 2.28 | -15.1% (faster) | 0 → 0 |
| scalar_linear | cached_config_chunk1 | 2.47 | 2.46 | -0.2% (inconclusive) | 0 → 0 |
| scalar_linear | config_setup | 5.69 | 5.44 | -4.4% (inconclusive) | 32 → 32 |
| scalar_linear | equations_setup | 2.69 | 2.67 | -0.8% (inconclusive) | 0 → 0 |
| scalar_linear | jacobian | 2.47 | 2.69 | +9.1% (inconclusive) | 0 → 0 |
| scalar_linear | prepared_equations | 4.94 | 5.35 | +8.3% (inconclusive) | 0 → 0 |
| scalar_linear | residual | 3.14 | 2.67 | -15.0% (inconclusive) | 0 → 0 |
| scalar_linear | residual_jacobian | 2.70 | 2.68 | -0.8% (inconclusive) | 0 → 0 |
| scalar_linear | solve | 25.95 | 25.97 | +0.1% (inconclusive) | 0 → 0 |
| scalar_linear | solve_warm | 40.24 | 25.37 | -37.0% (faster) | 0 → 0 |
| scalar_maxwell | backsolve | 3.07 | 2.67 | -13.0% (faster) | 0 → 0 |
| scalar_maxwell | cached_callable | 3.85 | 2.86 | -25.7% (inconclusive) | 0 → 0 |
| scalar_maxwell | cached_config | 3.85 | 2.68 | -30.3% (inconclusive) | 0 → 0 |
| scalar_maxwell | cached_config_chunk1 | 3.85 | 2.86 | -25.7% (inconclusive) | 0 → 0 |
| scalar_maxwell | config_setup | 7.67 | 5.20 | -32.2% (faster) | 32 → 32 |
| scalar_maxwell | equations_setup | 3.47 | 2.67 | -23.0% (faster) | 0 → 0 |
| scalar_maxwell | jacobian | 4.02 | 2.68 | -33.3% (faster) | 0 → 0 |
| scalar_maxwell | prepared_equations | 14.55 | 12.66 | -13.0% (faster) | 0 → 0 |
| scalar_maxwell | residual | 4.02 | 2.67 | -33.4% (faster) | 0 → 0 |
| scalar_maxwell | residual_jacobian | 4.54 | 4.07 | -10.2% (faster) | 0 → 0 |
| scalar_maxwell | solve | 46.03 | 29.51 | -35.9% (faster) | 0 → 0 |
| scalar_maxwell | solve_warm | 46.00 | 29.59 | -35.7% (faster) | 0 → 0 |
| scalar_powerlaw | backsolve | 3.17 | 2.47 | -21.9% (faster) | 0 → 0 |
| scalar_powerlaw | cached_callable | 48.31 | 34.13 | -29.3% (faster) | 0 → 0 |
| scalar_powerlaw | cached_config | 48.31 | 34.13 | -29.3% (faster) | 0 → 0 |
| scalar_powerlaw | cached_config_chunk1 | 48.31 | 34.13 | -29.3% (faster) | 0 → 0 |
| scalar_powerlaw | config_setup | 7.93 | 5.67 | -28.6% (faster) | 32 → 32 |
| scalar_powerlaw | equations_setup | 3.32 | 2.67 | -19.5% (faster) | 0 → 0 |
| scalar_powerlaw | jacobian | 48.31 | 34.13 | -29.3% (faster) | 0 → 0 |
| scalar_powerlaw | prepared_equations | 58.33 | 34.50 | -40.9% (faster) | 0 → 0 |
| scalar_powerlaw | residual | 25.74 | 17.35 | -32.6% (faster) | 0 → 0 |
| scalar_powerlaw | residual_jacobian | 73.57 | 55.98 | -23.9% (faster) | 0 → 0 |
| scalar_powerlaw | solve | 783.78 | 597.96 | -23.7% (faster) | 0 → 0 |
| scalar_powerlaw | solve_warm | 181.50 | 131.31 | -27.7% (faster) | 0 → 0 |
| siblings | backsolve | 74.33 | 82.53 | +11.0% (slower) | 0 → 0 |
| siblings | cached_callable | 36.84 | 38.10 | +3.4% (inconclusive) | 0 → 0 |
| siblings | cached_config | 36.84 | 38.10 | +3.4% (inconclusive) | 0 → 0 |
| siblings | cached_config_chunk1 | 36.84 | 38.10 | +3.4% (inconclusive) | 0 → 0 |
| siblings | config_setup | 6.60 | 6.79 | +2.9% (inconclusive) | 176 → 176 |
| siblings | equations_setup | 3.55 | 3.41 | -4.0% (inconclusive) | 0 → 0 |
| siblings | jacobian | 36.84 | 38.10 | +3.4% (inconclusive) | 0 → 0 |
| siblings | prepared_equations | 1004.55 | 1073.33 | +6.8% (inconclusive) | 3184 → 3184 |
| siblings | residual | 3.16 | 3.21 | +1.5% (inconclusive) | 0 → 0 |
| siblings | residual_jacobian | 49.72 | 41.32 | -16.9% (faster) | 0 → 0 |
| siblings | solve | 140.49 | 143.43 | +2.1% (inconclusive) | 0 → 0 |
| siblings | solve_warm | 140.95 | 144.15 | +2.3% (inconclusive) | 0 → 0 |
| volumetric_maxwell | backsolve | 2.88 | 2.88 | +0.0% (inconclusive) | 0 → 0 |
| volumetric_maxwell | cached_callable | 12.80 | 12.76 | -0.3% (inconclusive) | 0 → 0 |
| volumetric_maxwell | cached_config | 12.80 | 12.76 | -0.3% (inconclusive) | 0 → 0 |
| volumetric_maxwell | cached_config_chunk1 | 12.80 | 12.76 | -0.3% (inconclusive) | 0 → 0 |
| volumetric_maxwell | config_setup | 5.48 | 5.90 | +7.8% (inconclusive) | 64 → 64 |
| volumetric_maxwell | equations_setup | 2.90 | 2.89 | -0.2% (inconclusive) | 0 → 0 |
| volumetric_maxwell | jacobian | 12.80 | 12.76 | -0.3% (inconclusive) | 0 → 0 |
| volumetric_maxwell | prepared_equations | 13.51 | 13.52 | +0.1% (inconclusive) | 0 → 0 |
| volumetric_maxwell | residual | 3.09 | 3.07 | -0.8% (inconclusive) | 0 → 0 |
| volumetric_maxwell | residual_jacobian | 13.77 | 13.87 | +0.7% (inconclusive) | 0 → 0 |
| volumetric_maxwell | solve | 53.61 | 53.51 | -0.2% (inconclusive) | 0 → 0 |
| volumetric_maxwell | solve_warm | 53.56 | 53.37 | -0.3% (inconclusive) | 0 → 0 |

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
