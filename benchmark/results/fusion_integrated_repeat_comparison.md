# Benchmark report: residual_jacobian_fusion_integrated_repeat

Julia 1.13.0; ForwardDiff 1.4.6; AMD Ryzen 7 7800X3D 8-Core Processor           ; threads=1; rounds=3.

Times are medians of round medians, in nanoseconds. Negative changes mean faster.
A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.

| Fixture | Kernel | Baseline ns | Candidate ns | Change | Bytes before → after |
| --- | --- | ---: | ---: | --- | ---: |
| dilatant | backsolve | 28.87 | 28.51 | -1.2% (inconclusive) | 0 → 0 |
| dilatant | cached_callable | 35.19 | 34.45 | -2.1% (inconclusive) | 0 → 0 |
| dilatant | cached_config | 35.19 | 34.45 | -2.1% (inconclusive) | 0 → 0 |
| dilatant | cached_config_chunk1 | 35.19 | 34.45 | -2.1% (inconclusive) | 0 → 0 |
| dilatant | config_setup | 8.58 | 9.29 | +8.3% (inconclusive) | 112 → 112 |
| dilatant | equations_setup | 5.34 | 13.27 | +148.4% (slower) | 0 → 0 |
| dilatant | jacobian | 36.28 | 35.00 | -3.5% (inconclusive) | 0 → 0 |
| dilatant | prepared_equations | 78.85 | 69.17 | -12.3% (inconclusive) | 592 → 592 |
| dilatant | residual | 19.36 | 19.21 | -0.8% (inconclusive) | 0 → 0 |
| dilatant | residual_and_jacobian_fused | 33.85 | 33.45 | -1.2% (inconclusive) | 0 → 0 |
| dilatant | residual_jacobian | 33.92 | 34.00 | +0.2% (inconclusive) | 0 → 0 |
| dilatant | solve | 73.17 | 72.29 | -1.2% (inconclusive) | 0 → 0 |
| dilatant | solve_reference | 106.10 | 104.88 | -1.1% (inconclusive) | 0 → 0 |
| dilatant | solve_warm | 73.17 | 72.84 | -0.5% (inconclusive) | 0 → 0 |
| dilatant | solve_warm_reference | 106.10 | 106.10 | +0.0% (inconclusive) | 0 → 0 |
| kelvin_voigt | backsolve | 2.89 | 2.88 | -0.3% (inconclusive) | 0 → 0 |
| kelvin_voigt | cached_callable | 17.00 | 16.80 | -1.2% (inconclusive) | 0 → 0 |
| kelvin_voigt | cached_config | 17.00 | 16.74 | -1.6% (inconclusive) | 0 → 0 |
| kelvin_voigt | cached_config_chunk1 | 17.00 | 16.74 | -1.6% (inconclusive) | 0 → 0 |
| kelvin_voigt | config_setup | 6.19 | 6.49 | +5.0% (inconclusive) | 64 → 64 |
| kelvin_voigt | equations_setup | 2.89 | 2.88 | -0.4% (inconclusive) | 0 → 0 |
| kelvin_voigt | jacobian | 21.46 | 17.17 | -20.0% (faster) | 0 → 0 |
| kelvin_voigt | prepared_equations | 40.21 | 39.34 | -2.1% (inconclusive) | 176 → 176 |
| kelvin_voigt | residual | 3.11 | 3.12 | +0.2% (inconclusive) | 0 → 0 |
| kelvin_voigt | residual_and_jacobian_fused | 22.84 | 14.22 | -37.7% (faster) | 0 → 0 |
| kelvin_voigt | residual_jacobian | 25.07 | 18.24 | -27.2% (faster) | 0 → 0 |
| kelvin_voigt | solve | 44.72 | 48.95 | +9.5% (slower) | 0 → 0 |
| kelvin_voigt | solve_reference | 47.20 | 53.85 | +14.1% (slower) | 0 → 0 |
| kelvin_voigt | solve_warm | 44.72 | 49.31 | +10.3% (slower) | 0 → 0 |
| kelvin_voigt | solve_warm_reference | 47.20 | 53.85 | +14.1% (slower) | 0 → 0 |
| nested | backsolve | 35.32 | 23.29 | -34.1% (faster) | 0 → 0 |
| nested | cached_callable | 33.90 | 27.85 | -17.8% (faster) | 0 → 0 |
| nested | cached_config | 33.90 | 27.54 | -18.8% (faster) | 0 → 0 |
| nested | cached_config_chunk1 | 33.90 | 27.54 | -18.8% (faster) | 0 → 0 |
| nested | config_setup | 11.05 | 9.89 | -10.5% (inconclusive) | 112 → 112 |
| nested | equations_setup | 4.72 | 3.10 | -34.4% (faster) | 0 → 0 |
| nested | jacobian | 33.90 | 27.54 | -18.8% (faster) | 0 → 0 |
| nested | prepared_equations | 101.71 | 61.39 | -39.6% (faster) | 448 → 448 |
| nested | residual | 5.04 | 3.02 | -40.0% (inconclusive) | 0 → 0 |
| nested | residual_and_jacobian_fused | 42.67 | 35.50 | -16.8% (inconclusive) | 0 → 0 |
| nested | residual_jacobian | 46.50 | 38.16 | -17.9% (faster) | 0 → 0 |
| nested | solve | 92.11 | 64.76 | -29.7% (faster) | 0 → 0 |
| nested | solve_reference | 103.95 | 76.19 | -26.7% (faster) | 0 → 0 |
| nested | solve_warm | 92.11 | 64.76 | -29.7% (faster) | 0 → 0 |
| nested | solve_warm_reference | 103.95 | 76.19 | -26.7% (faster) | 0 → 0 |
| plastic | backsolve | 2.89 | 2.68 | -7.3% (faster) | 0 → 0 |
| plastic | cached_callable | 19.85 | 21.26 | +7.1% (slower) | 0 → 0 |
| plastic | cached_config | 19.85 | 21.26 | +7.1% (slower) | 0 → 0 |
| plastic | cached_config_chunk1 | 19.85 | 21.26 | +7.1% (slower) | 0 → 0 |
| plastic | config_setup | 6.38 | 6.32 | -0.8% (inconclusive) | 64 → 64 |
| plastic | equations_setup | 4.51 | 3.67 | -18.5% (faster) | 0 → 0 |
| plastic | jacobian | 19.84 | 21.26 | +7.2% (slower) | 0 → 0 |
| plastic | prepared_equations | 51.37 | 50.43 | -1.8% (inconclusive) | 416 → 416 |
| plastic | residual | 3.32 | 3.33 | +0.3% (inconclusive) | 0 → 0 |
| plastic | residual_and_jacobian_fused | 18.98 | 18.17 | -4.3% (inconclusive) | 0 → 0 |
| plastic | residual_jacobian | 18.22 | 17.20 | -5.6% (faster) | 0 → 0 |
| plastic | solve | 53.49 | 51.52 | -3.7% (inconclusive) | 0 → 0 |
| plastic | solve_reference | 61.54 | 60.61 | -1.5% (inconclusive) | 0 → 0 |
| plastic | solve_warm | 53.49 | 51.91 | -3.0% (inconclusive) | 0 → 0 |
| plastic | solve_warm_reference | 60.77 | 60.61 | -0.3% (inconclusive) | 0 → 0 |
| scalar_creep | backsolve | 2.46 | 2.94 | +19.4% (slower) | 0 → 0 |
| scalar_creep | cached_callable | 59.38 | 81.69 | +37.6% (slower) | 0 → 0 |
| scalar_creep | cached_config | 59.38 | 81.69 | +37.6% (slower) | 0 → 0 |
| scalar_creep | cached_config_chunk1 | 59.38 | 81.69 | +37.6% (slower) | 0 → 0 |
| scalar_creep | config_setup | 5.26 | 7.80 | +48.5% (inconclusive) | 32 → 32 |
| scalar_creep | equations_setup | 2.68 | 4.27 | +59.3% (slower) | 0 → 0 |
| scalar_creep | jacobian | 59.38 | 81.69 | +37.6% (slower) | 0 → 0 |
| scalar_creep | prepared_equations | 65.62 | 101.43 | +54.6% (slower) | 0 → 0 |
| scalar_creep | residual | 39.45 | 56.74 | +43.8% (slower) | 0 → 0 |
| scalar_creep | residual_and_jacobian_fused | 59.48 | 83.02 | +39.6% (inconclusive) | 0 → 0 |
| scalar_creep | residual_jacobian | 97.35 | 140.28 | +44.1% (slower) | 0 → 0 |
| scalar_creep | solve | 229.17 | 331.82 | +44.8% (slower) | 0 → 0 |
| scalar_creep | solve_reference | 360.87 | 513.64 | +42.3% (slower) | 0 → 0 |
| scalar_creep | solve_warm | 129.17 | 190.91 | +47.8% (slower) | 0 → 0 |
| scalar_creep | solve_warm_reference | 219.23 | 309.09 | +41.0% (slower) | 0 → 0 |
| scalar_linear | backsolve | 2.47 | 2.49 | +0.9% (inconclusive) | 0 → 0 |
| scalar_linear | cached_callable | 2.45 | 2.54 | +3.5% (inconclusive) | 0 → 0 |
| scalar_linear | cached_config | 2.68 | 2.74 | +2.2% (inconclusive) | 0 → 0 |
| scalar_linear | cached_config_chunk1 | 2.68 | 2.33 | -12.8% (inconclusive) | 0 → 0 |
| scalar_linear | config_setup | 5.27 | 7.30 | +38.4% (inconclusive) | 32 → 32 |
| scalar_linear | equations_setup | 2.68 | 3.32 | +23.7% (slower) | 0 → 0 |
| scalar_linear | jacobian | 2.68 | 2.54 | -5.3% (inconclusive) | 0 → 0 |
| scalar_linear | prepared_equations | 5.52 | 5.14 | -6.9% (inconclusive) | 0 → 0 |
| scalar_linear | residual | 2.67 | 2.49 | -7.0% (faster) | 0 → 0 |
| scalar_linear | residual_and_jacobian_fused | 2.68 | 2.70 | +0.8% (inconclusive) | 0 → 0 |
| scalar_linear | residual_jacobian | 2.68 | 2.70 | +0.7% (inconclusive) | 0 → 0 |
| scalar_linear | solve | 18.26 | 18.62 | +2.0% (inconclusive) | 0 → 0 |
| scalar_linear | solve_reference | 26.30 | 26.90 | +2.3% (inconclusive) | 0 → 0 |
| scalar_linear | solve_warm | 18.40 | 18.97 | +3.1% (inconclusive) | 0 → 0 |
| scalar_linear | solve_warm_reference | 26.05 | 26.90 | +3.3% (inconclusive) | 0 → 0 |
| scalar_maxwell | backsolve | 2.67 | 2.90 | +8.7% (inconclusive) | 0 → 0 |
| scalar_maxwell | cached_callable | 2.89 | 3.41 | +18.1% (inconclusive) | 0 → 0 |
| scalar_maxwell | cached_config | 2.64 | 3.51 | +33.0% (slower) | 0 → 0 |
| scalar_maxwell | cached_config_chunk1 | 2.64 | 3.20 | +21.3% (inconclusive) | 0 → 0 |
| scalar_maxwell | config_setup | 5.22 | 6.23 | +19.5% (slower) | 32 → 32 |
| scalar_maxwell | equations_setup | 2.67 | 2.72 | +1.7% (inconclusive) | 0 → 0 |
| scalar_maxwell | jacobian | 2.64 | 3.20 | +21.3% (inconclusive) | 0 → 0 |
| scalar_maxwell | prepared_equations | 12.78 | 14.77 | +15.6% (slower) | 0 → 0 |
| scalar_maxwell | residual | 2.68 | 2.71 | +1.1% (inconclusive) | 0 → 0 |
| scalar_maxwell | residual_and_jacobian_fused | 4.10 | 4.43 | +8.1% (slower) | 0 → 0 |
| scalar_maxwell | residual_jacobian | 4.10 | 4.44 | +8.2% (slower) | 0 → 0 |
| scalar_maxwell | solve | 22.18 | 39.73 | +79.1% (slower) | 0 → 0 |
| scalar_maxwell | solve_reference | 29.93 | 45.09 | +50.7% (slower) | 0 → 0 |
| scalar_maxwell | solve_warm | 22.18 | 39.73 | +79.1% (slower) | 0 → 0 |
| scalar_maxwell | solve_warm_reference | 29.93 | 45.09 | +50.7% (slower) | 0 → 0 |
| scalar_powerlaw | backsolve | 2.48 | 2.71 | +9.4% (slower) | 0 → 0 |
| scalar_powerlaw | cached_callable | 34.34 | 35.15 | +2.4% (inconclusive) | 0 → 0 |
| scalar_powerlaw | cached_config | 34.34 | 35.15 | +2.4% (inconclusive) | 0 → 0 |
| scalar_powerlaw | cached_config_chunk1 | 34.34 | 35.15 | +2.4% (inconclusive) | 0 → 0 |
| scalar_powerlaw | config_setup | 5.26 | 5.46 | +3.7% (inconclusive) | 32 → 32 |
| scalar_powerlaw | equations_setup | 2.67 | 3.18 | +19.0% (slower) | 0 → 0 |
| scalar_powerlaw | jacobian | 34.34 | 35.15 | +2.4% (inconclusive) | 0 → 0 |
| scalar_powerlaw | prepared_equations | 34.50 | 35.76 | +3.6% (inconclusive) | 0 → 0 |
| scalar_powerlaw | residual | 17.51 | 25.86 | +47.7% (slower) | 0 → 0 |
| scalar_powerlaw | residual_and_jacobian_fused | 38.82 | 51.30 | +32.1% (slower) | 0 → 0 |
| scalar_powerlaw | residual_jacobian | 56.40 | 73.71 | +30.7% (slower) | 0 → 0 |
| scalar_powerlaw | solve | 463.64 | 475.00 | +2.5% (inconclusive) | 0 → 0 |
| scalar_powerlaw | solve_reference | 600.00 | 620.00 | +3.3% (inconclusive) | 0 → 0 |
| scalar_powerlaw | solve_warm | 100.00 | 104.55 | +4.5% (inconclusive) | 0 → 0 |
| scalar_powerlaw | solve_warm_reference | 134.78 | 140.91 | +4.5% (inconclusive) | 0 → 0 |
| siblings | backsolve | 75.58 | 75.26 | -0.4% (inconclusive) | 0 → 0 |
| siblings | cached_callable | 34.78 | 35.29 | +1.5% (inconclusive) | 0 → 0 |
| siblings | cached_config | 34.78 | 35.29 | +1.5% (inconclusive) | 0 → 0 |
| siblings | cached_config_chunk1 | 34.78 | 35.29 | +1.5% (inconclusive) | 0 → 0 |
| siblings | config_setup | 6.79 | 6.90 | +1.6% (inconclusive) | 176 → 176 |
| siblings | equations_setup | 3.43 | 3.53 | +3.1% (inconclusive) | 0 → 0 |
| siblings | jacobian | 39.13 | 35.29 | -9.8% (faster) | 0 → 0 |
| siblings | prepared_equations | 820.00 | 852.94 | +4.0% (inconclusive) | 3184 → 3184 |
| siblings | residual | 3.12 | 3.19 | +2.1% (inconclusive) | 0 → 0 |
| siblings | residual_and_jacobian_fused | 42.20 | 35.65 | -15.5% (faster) | 0 → 0 |
| siblings | residual_jacobian | 41.61 | 48.75 | +17.2% (slower) | 0 → 0 |
| siblings | solve | 144.00 | 141.18 | -2.0% (inconclusive) | 0 → 0 |
| siblings | solve_reference | 150.00 | 145.10 | -3.3% (inconclusive) | 0 → 0 |
| siblings | solve_warm | 144.00 | 141.18 | -2.0% (inconclusive) | 0 → 0 |
| siblings | solve_warm_reference | 150.00 | 145.10 | -3.3% (inconclusive) | 0 → 0 |
| volumetric_maxwell | backsolve | 3.69 | 2.89 | -21.7% (faster) | 0 → 0 |
| volumetric_maxwell | cached_callable | 21.93 | 13.58 | -38.1% (faster) | 0 → 0 |
| volumetric_maxwell | cached_config | 21.90 | 13.58 | -38.0% (faster) | 0 → 0 |
| volumetric_maxwell | cached_config_chunk1 | 21.93 | 13.58 | -38.1% (faster) | 0 → 0 |
| volumetric_maxwell | config_setup | 8.97 | 6.46 | -28.0% (faster) | 64 → 64 |
| volumetric_maxwell | equations_setup | 4.21 | 2.89 | -31.4% (faster) | 0 → 0 |
| volumetric_maxwell | jacobian | 21.90 | 13.58 | -38.0% (faster) | 0 → 0 |
| volumetric_maxwell | prepared_equations | 23.79 | 14.52 | -39.0% (faster) | 0 → 0 |
| volumetric_maxwell | residual | 4.20 | 3.09 | -26.4% (faster) | 0 → 0 |
| volumetric_maxwell | residual_and_jacobian_fused | 24.26 | 14.36 | -40.8% (faster) | 0 → 0 |
| volumetric_maxwell | residual_jacobian | 23.34 | 14.06 | -39.7% (faster) | 0 → 0 |
| volumetric_maxwell | solve | 79.57 | 62.61 | -21.3% (faster) | 0 → 0 |
| volumetric_maxwell | solve_reference | 83.87 | 65.22 | -22.2% (faster) | 0 → 0 |
| volumetric_maxwell | solve_warm | 79.57 | 62.61 | -21.3% (faster) | 0 → 0 |
| volumetric_maxwell | solve_warm_reference | 83.87 | 65.22 | -22.2% (faster) | 0 → 0 |

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
