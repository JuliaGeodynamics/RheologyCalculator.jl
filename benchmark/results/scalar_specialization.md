# Scalar specialization probe (plan item 4)

Julia 1.13.0; ForwardDiff 1.4.6; threads=1.

| Fixture | jacobian() ns | scalar derivative ns | fused scalar ns | derivative vs jacobian() | fused vs jacobian() |
| --- | ---: | ---: | ---: | --- | --- |
| scalar_linear | 2.50 | 2.71 | 2.92 | +8.5% | +17.1% |
| scalar_maxwell | 2.74 | 2.96 | 4.19 | +7.9% | +52.6% |
| scalar_powerlaw | 38.26 | 38.26 | 38.26 | +0.0% | +0.0% |
| scalar_creep | 60.90 | 60.26 | 62.18 | -1.1% | +2.1% |
