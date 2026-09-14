# Benchmark report: baseline_v1_repeat

Julia 1.13.0; ForwardDiff 1.4.6; AMD Ryzen 7 7800X3D 8-Core Processor           ; threads=1; rounds=3.

Times are medians of round medians, in nanoseconds. Negative changes mean faster.
A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.

| Fixture | Kernel | Baseline ns | Candidate ns | Change | Bytes before → after |
| --- | --- | ---: | ---: | --- | ---: |
| dilatant | backsolve | 34.85 | 28.93 | -17.0% (faster) | 0 → 0 |
| dilatant | cached_callable | 43.12 | 33.94 | -21.3% (faster) | 0 → 0 |
| dilatant | cached_config | 34.21 | 33.90 | -0.9% (inconclusive) | 0 → 0 |
| dilatant | cached_config_chunk1 | 44.92 | 33.96 | -24.4% (faster) | 0 → 0 |
| dilatant | config_setup | 9.90 | 7.96 | -19.5% (faster) | 112 → 112 |
| dilatant | equations_setup | 5.52 | 5.00 | -9.4% (faster) | 0 → 0 |
| dilatant | jacobian | 43.09 | 33.42 | -22.4% (faster) | 0 → 0 |
| dilatant | prepared_equations | 84.97 | 55.36 | -34.8% (inconclusive) | 592 → 592 |
| dilatant | residual | 19.66 | 19.32 | -1.7% (inconclusive) | 0 → 0 |
| dilatant | residual_jacobian | 35.15 | 34.17 | -2.8% (inconclusive) | 0 → 0 |
| dilatant | solve | 131.80 | 104.72 | -20.5% (faster) | 0 → 0 |
| dilatant | solve_warm | 119.65 | 104.93 | -12.3% (faster) | 0 → 0 |
| kelvin_voigt | backsolve | 3.85 | 2.91 | -24.4% (inconclusive) | 0 → 0 |
| kelvin_voigt | cached_callable | 21.76 | 11.69 | -46.3% (faster) | 0 → 0 |
| kelvin_voigt | cached_config | 21.86 | 11.69 | -46.5% (faster) | 0 → 0 |
| kelvin_voigt | cached_config_chunk1 | 21.65 | 13.51 | -37.6% (faster) | 0 → 0 |
| kelvin_voigt | config_setup | 8.47 | 7.66 | -9.6% (inconclusive) | 64 → 64 |
| kelvin_voigt | equations_setup | 4.21 | 2.97 | -29.5% (faster) | 0 → 0 |
| kelvin_voigt | jacobian | 21.47 | 19.58 | -8.8% (faster) | 0 → 0 |
| kelvin_voigt | prepared_equations | 61.31 | 55.52 | -9.4% (faster) | 176 → 176 |
| kelvin_voigt | residual | 4.96 | 4.18 | -15.7% (inconclusive) | 0 → 0 |
| kelvin_voigt | residual_jacobian | 30.56 | 28.97 | -5.2% (faster) | 0 → 0 |
| kelvin_voigt | solve | 71.06 | 61.52 | -13.4% (faster) | 0 → 0 |
| kelvin_voigt | solve_warm | 71.07 | 47.65 | -32.9% (faster) | 0 → 0 |
| nested | backsolve | 23.45 | 26.76 | +14.1% (inconclusive) | 0 → 0 |
| nested | cached_callable | 38.78 | 34.99 | -9.8% (faster) | 0 → 0 |
| nested | cached_config | 38.84 | 37.83 | -2.6% (inconclusive) | 0 → 0 |
| nested | cached_config_chunk1 | 39.02 | 35.15 | -9.9% (faster) | 0 → 0 |
| nested | config_setup | 8.40 | 11.07 | +31.7% (inconclusive) | 112 → 112 |
| nested | equations_setup | 2.95 | 4.73 | +60.4% (inconclusive) | 0 → 0 |
| nested | jacobian | 39.12 | 35.29 | -9.8% (faster) | 0 → 0 |
| nested | prepared_equations | 58.95 | 61.32 | +4.0% (inconclusive) | 448 → 448 |
| nested | residual | 2.94 | 4.83 | +64.5% (slower) | 0 → 0 |
| nested | residual_jacobian | 37.25 | 29.98 | -19.5% (faster) | 0 → 0 |
| nested | solve | 76.99 | 99.31 | +29.0% (slower) | 0 → 0 |
| nested | solve_warm | 77.60 | 99.32 | +28.0% (inconclusive) | 0 → 0 |
| plastic | backsolve | 2.94 | 2.70 | -8.4% (faster) | 0 → 0 |
| plastic | cached_callable | 16.22 | 23.41 | +44.3% (inconclusive) | 0 → 0 |
| plastic | cached_config | 16.24 | 23.37 | +43.9% (inconclusive) | 0 → 0 |
| plastic | cached_config_chunk1 | 24.26 | 23.43 | -3.4% (inconclusive) | 0 → 0 |
| plastic | config_setup | 6.51 | 5.81 | -10.8% (inconclusive) | 64 → 64 |
| plastic | equations_setup | 3.97 | 3.90 | -1.9% (inconclusive) | 0 → 0 |
| plastic | jacobian | 24.10 | 23.08 | -4.2% (inconclusive) | 0 → 0 |
| plastic | prepared_equations | 67.04 | 50.81 | -24.2% (inconclusive) | 416 → 416 |
| plastic | residual | 3.43 | 3.34 | -2.8% (inconclusive) | 0 → 0 |
| plastic | residual_jacobian | 21.05 | 23.35 | +10.9% (inconclusive) | 0 → 0 |
| plastic | solve | 54.67 | 62.70 | +14.7% (inconclusive) | 0 → 0 |
| plastic | solve_warm | 55.83 | 62.50 | +12.0% (inconclusive) | 0 → 0 |
| scalar_creep | backsolve | 2.28 | 2.49 | +9.3% (slower) | 0 → 0 |
| scalar_creep | cached_callable | 59.62 | 64.27 | +7.8% (slower) | 0 → 0 |
| scalar_creep | cached_config | 59.25 | 60.56 | +2.2% (inconclusive) | 0 → 0 |
| scalar_creep | cached_config_chunk1 | 59.53 | 60.26 | +1.2% (inconclusive) | 0 → 0 |
| scalar_creep | config_setup | 4.78 | 5.69 | +18.9% (slower) | 32 → 32 |
| scalar_creep | equations_setup | 2.92 | 3.30 | +13.2% (inconclusive) | 0 → 0 |
| scalar_creep | jacobian | 59.51 | 60.96 | +2.4% (inconclusive) | 0 → 0 |
| scalar_creep | prepared_equations | 66.17 | 79.76 | +20.5% (inconclusive) | 0 → 0 |
| scalar_creep | residual | 39.44 | 40.59 | +2.9% (inconclusive) | 0 → 0 |
| scalar_creep | residual_jacobian | 98.97 | 103.14 | +4.2% (inconclusive) | 0 → 0 |
| scalar_creep | solve | 366.67 | 494.74 | +34.9% (inconclusive) | 0 → 0 |
| scalar_creep | solve_warm | 226.23 | 272.80 | +20.6% (inconclusive) | 0 → 0 |
| scalar_linear | backsolve | 2.68 | 2.72 | +1.6% (inconclusive) | 0 → 0 |
| scalar_linear | cached_callable | 2.68 | 2.74 | +2.3% (inconclusive) | 0 → 0 |
| scalar_linear | cached_config | 2.68 | 3.05 | +13.8% (slower) | 0 → 0 |
| scalar_linear | cached_config_chunk1 | 2.48 | 2.50 | +1.0% (inconclusive) | 0 → 0 |
| scalar_linear | config_setup | 5.39 | 5.63 | +4.5% (inconclusive) | 32 → 32 |
| scalar_linear | equations_setup | 2.68 | 3.23 | +20.8% (slower) | 0 → 0 |
| scalar_linear | jacobian | 2.47 | 2.50 | +1.4% (inconclusive) | 0 → 0 |
| scalar_linear | prepared_equations | 5.30 | 4.91 | -7.4% (inconclusive) | 0 → 0 |
| scalar_linear | residual | 2.68 | 2.72 | +1.3% (inconclusive) | 0 → 0 |
| scalar_linear | residual_jacobian | 2.68 | 2.73 | +2.1% (inconclusive) | 0 → 0 |
| scalar_linear | solve | 25.59 | 29.18 | +14.0% (slower) | 0 → 0 |
| scalar_linear | solve_warm | 25.39 | 27.04 | +6.5% (slower) | 0 → 0 |
| scalar_maxwell | backsolve | 2.47 | 3.08 | +24.7% (slower) | 0 → 0 |
| scalar_maxwell | cached_callable | 2.69 | 3.57 | +32.7% (slower) | 0 → 0 |
| scalar_maxwell | cached_config | 2.68 | 3.40 | +26.9% (slower) | 0 → 0 |
| scalar_maxwell | cached_config_chunk1 | 2.67 | 2.74 | +2.5% (inconclusive) | 0 → 0 |
| scalar_maxwell | config_setup | 4.83 | 7.29 | +51.0% (slower) | 32 → 32 |
| scalar_maxwell | equations_setup | 2.67 | 3.46 | +29.7% (slower) | 0 → 0 |
| scalar_maxwell | jacobian | 2.68 | 3.39 | +26.5% (slower) | 0 → 0 |
| scalar_maxwell | prepared_equations | 12.84 | 14.64 | +14.0% (slower) | 0 → 0 |
| scalar_maxwell | residual | 2.88 | 3.95 | +37.2% (slower) | 0 → 0 |
| scalar_maxwell | residual_jacobian | 4.06 | 4.34 | +6.9% (slower) | 0 → 0 |
| scalar_maxwell | solve | 29.97 | 46.35 | +54.7% (slower) | 0 → 0 |
| scalar_maxwell | solve_warm | 30.04 | 32.83 | +9.3% (slower) | 0 → 0 |
| scalar_powerlaw | backsolve | 2.68 | 3.35 | +25.1% (slower) | 0 → 0 |
| scalar_powerlaw | cached_callable | 34.01 | 47.85 | +40.7% (slower) | 0 → 0 |
| scalar_powerlaw | cached_config | 34.16 | 48.27 | +41.3% (slower) | 0 → 0 |
| scalar_powerlaw | cached_config_chunk1 | 34.29 | 34.96 | +1.9% (inconclusive) | 0 → 0 |
| scalar_powerlaw | config_setup | 5.36 | 7.42 | +38.5% (slower) | 32 → 32 |
| scalar_powerlaw | equations_setup | 2.69 | 3.35 | +24.3% (slower) | 0 → 0 |
| scalar_powerlaw | jacobian | 34.15 | 46.03 | +34.8% (slower) | 0 → 0 |
| scalar_powerlaw | prepared_equations | 34.39 | 58.66 | +70.6% (slower) | 0 → 0 |
| scalar_powerlaw | residual | 17.54 | 25.86 | +47.4% (slower) | 0 → 0 |
| scalar_powerlaw | residual_jacobian | 56.02 | 73.47 | +31.1% (slower) | 0 → 0 |
| scalar_powerlaw | solve | 602.33 | 786.49 | +30.6% (slower) | 0 → 0 |
| scalar_powerlaw | solve_warm | 133.17 | 181.60 | +36.4% (slower) | 0 → 0 |
| siblings | backsolve | 93.04 | 91.92 | -1.2% (inconclusive) | 0 → 0 |
| siblings | cached_callable | 49.51 | 54.98 | +11.1% (slower) | 0 → 0 |
| siblings | cached_config | 51.55 | 54.98 | +6.6% (slower) | 0 → 0 |
| siblings | cached_config_chunk1 | 51.24 | 54.98 | +7.3% (inconclusive) | 0 → 0 |
| siblings | config_setup | 10.17 | 8.83 | -13.2% (faster) | 176 → 176 |
| siblings | equations_setup | 5.68 | 5.25 | -7.6% (inconclusive) | 0 → 0 |
| siblings | jacobian | 50.80 | 54.90 | +8.1% (slower) | 0 → 0 |
| siblings | prepared_equations | 1350.00 | 1328.00 | -1.6% (inconclusive) | 3184 → 3184 |
| siblings | residual | 4.84 | 4.52 | -6.6% (faster) | 0 → 0 |
| siblings | residual_jacobian | 71.10 | 68.78 | -3.3% (inconclusive) | 0 → 0 |
| siblings | solve | 175.98 | 180.86 | +2.8% (inconclusive) | 0 → 0 |
| siblings | solve_warm | 178.11 | 184.11 | +3.4% (inconclusive) | 0 → 0 |
| volumetric_maxwell | backsolve | 2.91 | 3.91 | +34.3% (slower) | 0 → 0 |
| volumetric_maxwell | cached_callable | 17.91 | 21.95 | +22.6% (slower) | 0 → 0 |
| volumetric_maxwell | cached_config | 17.87 | 21.74 | +21.7% (slower) | 0 → 0 |
| volumetric_maxwell | cached_config_chunk1 | 17.97 | 22.14 | +23.2% (slower) | 0 → 0 |
| volumetric_maxwell | config_setup | 5.51 | 8.99 | +63.2% (slower) | 64 → 64 |
| volumetric_maxwell | equations_setup | 3.31 | 4.28 | +29.3% (slower) | 0 → 0 |
| volumetric_maxwell | jacobian | 17.77 | 22.01 | +23.9% (slower) | 0 → 0 |
| volumetric_maxwell | prepared_equations | 19.68 | 23.89 | +21.4% (slower) | 0 → 0 |
| volumetric_maxwell | residual | 3.10 | 4.18 | +35.0% (slower) | 0 → 0 |
| volumetric_maxwell | residual_jacobian | 18.77 | 28.09 | +49.6% (slower) | 0 → 0 |
| volumetric_maxwell | solve | 57.71 | 67.59 | +17.1% (slower) | 0 → 0 |
| volumetric_maxwell | solve_warm | 58.16 | 67.66 | +16.3% (slower) | 0 → 0 |

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
