# Elastic strain-rate correction

When a composite model contains springs, its constitutive equations include
backstress history terms: the stress carried by each spring at the *previous*
time step, ``\boldsymbol{\tau}^o``. The solver keeps the element state
functions free of that history. Instead, the history enters through one tensor,
the **history tensor** ``\mathbf{H}``, and the global deviatoric equation of
[`solve`](@ref) is posed on the second invariant of the **assembled strain
rate**

```math
\mathbf{E} = \boldsymbol{\varepsilon} + \mathbf{H} .
```

``\mathbf{E}`` is assembled as a tensor before the invariant is taken, so
backstresses that are not coaxial with ``\boldsymbol{\varepsilon}``, or that
oppose it after a load reversal, enter with their orientation and sign. For a
scalar strain rate every quantity is signed. [`effective_strain_rate_correction`](@ref)
returns ``\mathbf{H}``. This page derives it, starting from the simplest case.
The derivations in full are in `docs/derivations/tensor_reduction.typ` and
`docs/derivations/noncoaxial_reduction.typ`.

All springs use a backward-Euler update:

```math
\boldsymbol{\tau}^{\mathrm{el}} = \boldsymbol{\tau}^o + 2G\,\Delta t\,\boldsymbol{\varepsilon}^{\mathrm{el}}
```

so the elastic strain rate is ``(\boldsymbol{\tau}^{\mathrm{el}} - \boldsymbol{\tau}^o)/(2G\Delta t)``.

---

## Maxwell body

**Model:** `SeriesModel(η, Elasticity(G))`.

In a series model strain rates add:

```math
\boldsymbol{\varepsilon} = \frac{\boldsymbol{\tau}}{2\eta}
  + \frac{\boldsymbol{\tau} - \boldsymbol{\tau}^o}{2G\Delta t}
```

Moving the backstress term to the left-hand side so that the right-hand side
depends only on ``\boldsymbol{\tau}``:

```math
\underbrace{\boldsymbol{\varepsilon} + \frac{\boldsymbol{\tau}^o}{2G\Delta t}}
           _{\boldsymbol{\varepsilon}^{\mathrm{eff}}}
= \boldsymbol{\tau}
  \underbrace{\left(\frac{1}{2\eta} + \frac{1}{2G\Delta t}\right)}
             _{1/(2\eta^{\mathrm{eff}})}
```

The Maxwell effective viscosity is the inverse of the sum of the inverses of
``\eta`` and ``G\Delta t``:

```math
\eta^{\mathrm{eff}}_M = \frac{1}{1/\eta + 1/(G\Delta t)}
```

The strain-rate correction is ``\boldsymbol{\tau}^o/(2G\Delta t)``: the
backstress divided by ``2G\Delta t``, the effective viscosity of the spring
times two.

---

## Kelvin-Voigt body

**Model:** `SeriesModel(η₁, ParallelModel(Elasticity(G), η₂))`.

The outer series splits the total strain rate:

```math
\boldsymbol{\varepsilon} = \frac{\boldsymbol{\tau}}{2\eta_1} + \boldsymbol{\varepsilon}^p
```

Inside the parallel (Kelvin-Voigt) branch stresses add:

```math
\boldsymbol{\tau} = 2\eta_2\,\boldsymbol{\varepsilon}^p
  + \left(2G\Delta t\,\boldsymbol{\varepsilon}^p + \boldsymbol{\tau}^o\right)
= (2\eta_2 + 2G\Delta t)\,\boldsymbol{\varepsilon}^p + \boldsymbol{\tau}^o
```

Defining the Kelvin-Voigt effective viscosity ``\eta^{\mathrm{eff}}_{\mathrm{KV}} = \eta_2 + G\Delta t``
and solving for ``\boldsymbol{\varepsilon}^p``:

```math
\boldsymbol{\varepsilon}^p = \frac{\boldsymbol{\tau} - \boldsymbol{\tau}^o}{2\eta^{\mathrm{eff}}_{\mathrm{KV}}}
```

Substituting back into the series equation and isolating the effective strain
rate:

```math
\boldsymbol{\varepsilon}^{\mathrm{eff}}
= \boldsymbol{\varepsilon} + \frac{\boldsymbol{\tau}^o}{2\eta^{\mathrm{eff}}_{\mathrm{KV}}}
= \boldsymbol{\tau}\!\left(\frac{1}{2\eta_1} + \frac{1}{2\eta^{\mathrm{eff}}_{\mathrm{KV}}}\right)
```

The correction now scales with ``\eta^{\mathrm{eff}}_{\mathrm{KV}}`` rather than
``G\Delta t`` alone, because the elastic backstress must be *distributed across
the stiffness of the whole parallel branch*.  For ``N`` parallel elements the
KV effective viscosity generalizes to a sum:

```math
\eta^{\mathrm{eff}}_{\mathrm{KV}} = \sum_{i=1}^N \eta_i
```

---

## Mixed Kelvin-Voigt and Maxwell body

**Model:** `SeriesModel(η₁, ParallelModel(η₂, SeriesModel(η₃, Elasticity(G))))`.

The inner `SeriesModel(η₃, Elasticity(G))` is itself a Maxwell element with
its own effective viscosity ``\eta^{\mathrm{eff}}_M = 1/(1/\eta_3 + 1/(G\Delta t))``.
Its strain rate equals the branch strain rate ``\boldsymbol{\varepsilon}^p`` of
the outer parallel block:

```math
\boldsymbol{\varepsilon}^p = \frac{\boldsymbol{\tau}^{\mathrm{Max}}}{2\eta^{\mathrm{eff}}_M}
  - \frac{\boldsymbol{\tau}^o}{2G\Delta t}
\quad\Longrightarrow\quad
\boldsymbol{\tau}^{\mathrm{Max}} = 2\eta^{\mathrm{eff}}_M
  \!\left(\boldsymbol{\varepsilon}^p + \frac{\boldsymbol{\tau}^o}{2G\Delta t}\right)
```

The parallel stress balance is:

```math
\boldsymbol{\tau} = 2\eta_2\,\boldsymbol{\varepsilon}^p + \boldsymbol{\tau}^{\mathrm{Max}}
```

Substituting ``\boldsymbol{\tau}^{\mathrm{Max}}`` and collecting
``\boldsymbol{\varepsilon}^p``:

```math
\boldsymbol{\tau}
= 2(\eta_2 + \eta^{\mathrm{eff}}_M)\,\boldsymbol{\varepsilon}^p
  + \frac{\eta^{\mathrm{eff}}_M}{G\Delta t}\,\boldsymbol{\tau}^o
```

```math
\boldsymbol{\varepsilon}^p
= \frac{1}{2\eta^{\mathrm{eff}}_{\mathrm{KV}}}
  \!\left(\boldsymbol{\tau} - \eta^{\star}_M\,\boldsymbol{\tau}^o\right),
\qquad
\eta^{\mathrm{eff}}_{\mathrm{KV}} = \eta_2 + \eta^{\mathrm{eff}}_M,
\qquad
\eta^{\star}_M = \frac{\eta^{\mathrm{eff}}_M}{G\Delta t}
```

Substituting back into the outer series equation, the effective strain rate is:

```math
\boldsymbol{\varepsilon}^{\mathrm{eff}}
= \boldsymbol{\varepsilon} + \frac{\eta^{\star}_M\,\boldsymbol{\tau}^o}{2\eta^{\mathrm{eff}}_{\mathrm{KV}}}
```

The weighting factor ``\eta^{\star}_M \in (0, 1)`` attenuates the backstress:

- **Soft spring** (``G \to 0``): ``\eta^{\star}_M \to 1`` — full backstress, elastic history dominates.
- **Stiff spring** (``G \to \infty``): ``\eta^{\star}_M \to 0`` — backstress vanishes, element behaves as a pure dashpot.

### Example

Under a constant deviatoric strain rate ``\dot{\varepsilon}`` starting from zero stress, this
model obeys a linear first-order ODE whose exact solution is

```math
\tau(t)
= \tau_\infty - \bigl(\tau_\infty - \tau_{\mathrm{init}}\bigr)\,e^{-t/t_{\mathrm{relax}}}
```

with

```math
\tau_{\mathrm{init}} = \frac{2\eta_1\eta_2}{\eta_1+\eta_2}\,\dot{\varepsilon}
\qquad\text{(initial stress — spring unloaded)}
```

```math
\tau_\infty = \frac{2\eta_1(\eta_2+\eta_3)}{\eta_1+\eta_2+\eta_3}\,\dot{\varepsilon}
\qquad\text{(long-time equilibrium)}
```

```math
t_{\mathrm{relax}} = \frac{(\eta_1+\eta_2)\,\eta_3}{G\,(\eta_1+\eta_2+\eta_3)}
\qquad\text{(relaxation timescale)}
```

At ``t=0`` the spring carries no backstress so only the viscous elements ``\eta_1`` and
``\eta_2`` resist deformation; as ``t \to \infty`` the inner Maxwell dashpot reaches its
steady-state load and ``\tau`` approaches the higher value ``\tau_\infty``.

The full runnable example with convergence plot is in `examples/Maxwell_KV_Maxwell.jl`.
The essential time-stepping structure below uses the bundled material models and their
2D tensor helpers:

```julia
using RheologyCalculator
import RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
using RheologyCalculator.RheologyModels

η1 = LinearViscosity(1e22)
η2 = LinearViscosity(1e21)
η3 = LinearViscosity(1e21)
el = IncompressibleElasticity(1e10)

c = SeriesModel(η1, ParallelModel(η2, SeriesModel(η3, el)))

εII    = 1.0e-14
vars   = RheologyModels.vars_2D(εII)
args   = (; τ = 2.0e7, P = 0.0)
others = (; dt = 1.0e9, τ0 = (RheologyModels.zero_stress_tensor_2D(),), P0 = (0.0,))

x = initial_guess_x(c, vars, args, others)

τ_char = 2η1.η * (η2.η + η3.η) / (η1.η + η2.η + η3.η) * εII
xnorm  = normalisation_x(c, τ_char, εII)

τ_e = (RheologyModels.zero_stress_tensor_2D(),)
P_e = (0.0,)
t   = 0.0
for _ in 1:1_000
    others = (; dt = 1.0e9, τ0 = τ_e, P0 = P_e)
    x      = solve(c, x, vars, others; xnorm0 = xnorm).x
    # the spring sits inside the block: its stress is recovered from the
    # whole solver vector, not from the outer stress x[1]
    τ_e    = elastic_stress_history(c, x, vars.ε, τ_e, others)
    t     += 1.0e9
end
```

`x_keys(c)` returns `(:τ, :ε, :τ)` for this composite: the outer stress, the
block strain rate, and the stress of the inner Maxwell sub-branch. The last two
are *shifted* values (see [Unknowns inside a block](@ref)), so the spring stress
is recovered by [`elastic_stress_history`](@ref) from the whole solver vector
and the applied strain rate.

---

## Generalized parallel block

The cases above are instances of one pattern. Consider a `SeriesModel` whose
branches are `ParallelModel` blocks, each holding direct elements and `SeriesModel`
sub-branches of leaf elements. Every element counts with its effective viscosity:
``\eta`` for a dashpot, ``G\Delta t`` for a spring. For block ``b``:

```math
\eta_{M,j} = \Bigl(\sum_{\text{elements of } j} \frac{1}{\eta}\Bigr)^{-1},
\qquad
\eta_{\mathrm{KV}} = \sum_{\text{direct elements}} \eta + \sum_j \eta_{M,j},
```

```math
\eta^{\star}_i =
\begin{cases}
  1 & \text{direct spring of the block} \\[4pt]
  \dfrac{\eta_{M,j}}{G_i\Delta t} & \text{spring } i \text{ of sub-branch } j
\end{cases}
```

where ``G_i`` is the modulus of spring ``i`` itself: a sub-branch with several
springs weighs each with its own stiffness. The block contributes

```math
\frac{1}{2\eta_{\mathrm{KV}}} \sum_{\text{springs } i \text{ of } b} \eta^{\star}_i\,\boldsymbol{\tau}^o_i .
```

A block without springs has no entries in `τ0` and contributes nothing. The
history tensor is the sum of the outer-spring and block terms:

```math
\mathbf{H}
= \underbrace{\sum_{k \in \text{outer springs}} \frac{\boldsymbol{\tau}^o_k}{2G_k\Delta t}}_{\text{Maxwell}}
  + \underbrace{\sum_{b \in \text{blocks}}
      \frac{1}{2\eta^{(b)}_{\mathrm{KV}}}
      \sum_{i \in b} \eta^{\star}_{i}\,\boldsymbol{\tau}^o_i}_{\text{Kelvin-Voigt / generalized Maxwell}}
```

and the global deviatoric equation reads (outer leaf strain rates at ``\tau_{II}``)
``+ \sum_b \tilde{x}_b = E_{II}``, never ``\varepsilon_{II}`` plus a sum of invariants.
`τ0` holds one entry per spring in the order of the composite: outer springs
first, then each block's direct springs followed by the springs of its
sub-branches.

### Supported topologies

The formulas hold for a `SeriesModel` of leaf elements and `ParallelModel`
blocks whose sub-branches hold leaf elements only. A block that carries elastic
history and has a composite nested inside a sub-branch throws an error. The
elements are assumed isotropic, and the inelastic strain rates of the outer
leafs coaxial with ``\boldsymbol{\tau}``. `τ0` is taken in the current frame;
rotating it is up to the caller.

### Unknowns inside a block

The spring state functions are history-free, so inside a block that carries
history the solver unknowns are **shifted** values: the block strain rate
``\tilde{x}_b`` and the sub-branch stresses ``\tilde{x}_j`` of the history-free
problem with the same outer stress. With ``\mathbf{n} = \mathbf{E}/E_{II}``
(``\mathbf{n} = 1`` for a scalar strain rate), ``\boldsymbol{\tau} = \tau_{II}\mathbf{n}``
and the physical block strain rate and sub-branch stresses are

```math
\dot{\boldsymbol{\varepsilon}}_b = \tilde{x}_b\,\mathbf{n} - \frac{\mathbf{W}_b}{2\eta_{\mathrm{KV}}},
\qquad
\boldsymbol{\tau}_j = 2\eta_{M,j}\,\dot{\boldsymbol{\varepsilon}}_b + \mathbf{S}_j,
```

with ``\mathbf{W}_b = \sum_{i \in b} \eta^{\star}_i \boldsymbol{\tau}^o_i`` and
``\mathbf{S}_j`` its part owned by the springs of sub-branch ``j``. Each spring
then carries ``\boldsymbol{\tau}`` (outer spring),
``\boldsymbol{\tau}^o_s + 2G_s\Delta t\,\dot{\boldsymbol{\varepsilon}}_b`` (direct
spring of block ``b``), or ``\boldsymbol{\tau}_j`` (spring of sub-branch ``j``).
These tensors are what [`elastic_stress_history`](@ref) returns and what
[`compute_stress_elastic`](@ref) reduces to invariants;
[`component_partition`](@ref) evaluates the elements of a block at
``\|\dot{\boldsymbol{\varepsilon}}_b\|`` and ``\|\boldsymbol{\tau}_j\|``.
A spring inside a block is in general not coaxial with ``\boldsymbol{\tau}``,
so these functions need the applied strain rate.

### Elements whose viscosity depends on the state

``\eta_{\mathrm{KV}}``, ``\eta_{M,j}`` and ``\eta^{\star}`` must be evaluated at the
physical state of each element, which depends on those same viscosities. For a
block that carries history, the solver closes the loop with extra unknowns:

| Unknown | Key in [`x_keys`](@ref) | Added for |
|---|---|---|
| ``e_b = \|\dot{\boldsymbol{\varepsilon}}_b\|`` | `:εII_block` | a block with such an element among its direct leafs |
| ``t_j = \|\boldsymbol{\tau}_j\|`` | `:τII_subbranch` | each sub-branch with such an element |

Each such element is evaluated through its own state function at ``|e_b|`` or
``|t_j|`` as a secant viscosity: ``\eta = \tau(e_b)/(2e_b)`` for a direct
element, ``\eta = t_j/(2\dot\varepsilon(t_j))`` for an element of a
sub-branch. The graph equations, ``\mathbf{H}`` and the closure equations share
these coefficients, so the converged solution satisfies the tensor equations
exactly. The extra unknowns are the first local equations of their block. For a
scalar strain rate they are signed, like every other unknown.

Which elements need them is decided by [`viscosity_depends_on_state`](@ref).
It defaults to `true`, which is always correct and only costs the extra
unknowns; `LinearViscosity` and the springs return `false`. A user rheology
whose viscosity is independent of ``\boldsymbol{\tau}`` and
``\boldsymbol{\varepsilon}`` may return `false`; returning `false` for a
nonlinear law gives wrong results. Plastic elements in a block that carries
history are evaluated at the shifted block strain rate.

---

## Implementation

``\mathbf{H}`` is computed by [`effective_strain_rate_correction`](@ref) (in
`src/post_processing/strain_rate_correction.jl`). Inside [`compute_residual`](@ref)
it is evaluated at the current solver vector, so it is part of the
differentiated residual. The closure unknowns are in
`src/equation_system/closure_unknowns.jl`.

| Function | Role |
|---|---|
| `effective_strain_rate_correction(leafs, (), ε, τ0, others)` | Maxwell term of the outer springs |
| `_kv_corrections(branches, ε, τ0, others, offset, εb, phys)` | block terms of all `ParallelModel` branches |
| `_block_state(leafs, subs, τ0, args, el_idx_start, phys, seed)` | ``\eta_{\mathrm{KV}}``, ``\mathbf{W}``, ``\eta_{M,j}`` and ``\mathbf{S}_j`` of one block |
| `_η_KV(leafs, subs, args, phys)` | ``\eta_{\mathrm{KV}}`` of one block |
| `_η_eff_maxwell(leafs, args, t, el)` | ``\eta_{M,j}``, the inverse sum over a sub-branch |
| `_secant_viscosity(fn, r, args, s, n)` | secant viscosity of an element at physical argument `s` |
| `closure_unknowns(c, fn)` | which closure unknowns a block adds |
| `add_closure_residuals(residual, eqs, x, E, others)` | residuals of the closure unknowns |
| `elastic_stress_history(c, x, ε, τ0, others)` | spring stress tensors after a solve |

All τ0 index arithmetic and ``\eta^{\star}`` factors are resolved at
compile time via `@generated` functions, so the emitted code is a flat
sequence of multiply-adds with no runtime dispatch or branching.

### Load reversal and change of direction

The test suite (`test/test_elastic_history_tensor.jl`) runs time loops like the
one below for Maxwell, Kelvin-Voigt, Burgers, mixed, two-block and
multi-spring composites, and for power laws inside and outside blocks. It
compares the outer stress and every spring tensor at every step with a solver
that keeps all unknowns as full tensors, under constant loading, reversal, a
change from pure to simple shear, and 3D loading.

```julia
using RheologyCalculator
using RheologyCalculator.RheologyModels

# Burgers body and a Kelvin-Voigt block with a power law
burgers = SeriesModel(LinearViscosity(1e21), IncompressibleElasticity(1e10),
                      ParallelModel(LinearViscosity(5e20), IncompressibleElasticity(2.5e10)))
kv_powerlaw = SeriesModel(LinearViscosity(1e21),
                          ParallelModel(PowerLawViscosity(2.5e34, 3), IncompressibleElasticity(1e10)))

εII, dt = 1e-14, 1e10
pure_shear   = tensor_strain_rate_2D(εII)
simple_shear = tensor_strain_rate_2D(εII; direction = (0.0, 0.0, 1.0))
load(n) = n ≤ 6 ? pure_shear : n ≤ 12 ? map(-, pure_shear) : simple_shear

for (c, nsprings) in ((burgers, 2), (kv_powerlaw, 1))
    τ0 = ntuple(_ -> zero_stress_tensor_2D(), nsprings)
    xnorm = normalisation_x(c, 1e7, εII)
    x = nothing
    for n in 1:18
        vars = (; ε = load(n), θ = 0.0)
        others = (; dt, τ0)
        x0 = x === nothing ? initial_guess_x(c, vars, (; τ = 1e7), others) : x
        x = solve(c, x0, vars, others; xnorm0 = xnorm).x
        τ0 = elastic_stress_history(c, x, vars.ε, τ0, others)
    end
end
```

### Worked example

```julia
using RheologyCalculator
using RheologyCalculator.RheologyModels

η1 = LinearViscosity(1e22)   # outer series dashpot
η2 = LinearViscosity(1e21)   # parallel dashpot (KV branch)
η3 = LinearViscosity(5e20)   # inner Maxwell dashpot
el = IncompressibleElasticity(30e9)

# SeriesModel(η1, ParallelModel(η2, SeriesModel(η3, el)))
c = SeriesModel(η1, ParallelModel(η2, SeriesModel(η3, el)))

dt   = 1e10
τ0   = (1e6,)          # backstress for the one elastic element
ε    = (1e-14,)
others = (; dt, τ0, P0 = (0.0,))

# Analytical expected value
G_dt    = el.G * dt
η_eff_M = 1 / (1/η3.η + 1/G_dt)   # Maxwell eff. viscosity of inner branch
η_star  = η_eff_M / G_dt           # attenuation weight
η_KV    = η2.η + η_eff_M           # KV eff. viscosity of parallel block
correction = η_star * τ0[1] / (2η_KV)

# Computed value
import RheologyCalculator: effective_strain_rate_correction
cor = effective_strain_rate_correction(c, ε, τ0, others)

@assert isapprox(only(cor), correction; rtol = 1e-10)
```
