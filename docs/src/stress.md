# Stress-time curve of a Burger's material

The goal of this section is to better understand how to solve a rheology model in a time loop. We will use a simple Burgers model. We will then see how we can warm-start the Newton iteration from the previous solution to speed-up the convergence. Finally, we will see how to compute the Jacobian at the converged solution, which is useful for computing derivatives through the local solve.

As in the previous section of the tutorial, in addition to the package `RheologyCalculator`, the `RheologyCalculator.RheologyModels` submodule is needed to have access to the material elements but also helpers, such as [`compute_stress_elastic`](@ref RheologyCalculator.compute_stress_elastic)
and [`compute_pressure_elastic`](@ref RheologyCalculator.compute_pressure_elastic):

```julia
using RheologyCalculator
using RheologyCalculator.RheologyModels
```

The Burgers model is composed of a series of a Maxwell element and a Kelvin–Voigt element. The Maxwell element consists of a spring and a dashpot in series, while the Kelvin–Voigt element consists of a spring and a dashpot in parallel.

We start by defining the material properties
```julia
η1, η2, η3 = 5e19, 1e20, 1e21
damper1    = LinearViscosity(η1)
damper2    = LinearViscosity(η2)
damper3    = LinearViscosity(η3)
G, K1, K2  = 10e9, 46.67e9, 30e9
spring1    = Elasticity(G, K1)
spring2    = Elasticity(G, K2)
```

and rheology model

```julia-repl
julia> p = ParallelModel(damper2, spring1)
|--⟦▪̲̅▫̲̅▫̲̅▫̲̅¹--|
|--/\/\/¹--|

julia> c = SeriesModel(damper3, spring2, p)
--⟦▪̲̅▫̲̅▫̲̅▫̲̅¹----/\/\/¹--|--⟦▪̲̅▫̲̅▫̲̅▫̲̅²--|
                    |--/\/\/²--|
```

Next, we define the input variables, which are the strain rate `ε` and the volumetric strain rate `θ`
```julia
vars = (; ε = 1.0e-15, θ = 1.0e-20)  # input variables (constant)
```
and the values of the unknowns from which the initial guess is estimated:
```julia
args = (; τ = 1.0e3, P = 1.0e6)
```
The auxiliary variables are collected in `others`. These values are passed to
state functions but are not differentiated by the Newton solve:

```julia
others = (; dt = 1.0e9, τ0 = (0.0, 0.0), P0 = (0.0, 0.0))
```

and we build the solution vector `x` that contains the initial guess for the variables we want to solve for:
```julia
x   = initial_guess_x(c, vars, args, others)
```

The entries of `x` differ by twenty orders of magnitude — stresses of order
``10^5~\mathrm{Pa}`` alongside strain rates of order ``10^{-15}~\mathrm{s}^{-1}`` — so the residual
norm has to be scaled entry by entry before the requested tolerance means
anything. [`normalisation_x`](@ref) builds that scaling from a characteristic
stress and strain rate, and it is passed to [`solve`](@ref) as `xnorm0`. Without
it the iteration stalls and `solve` raises a `NonConvergenceError`:

```julia
τ_char = 2 * η3 * vars.ε
xnorm  = normalisation_x(c, τ_char, vars.ε)
```

Now we are ready to compute the time evolution of the stress tensor, with some aid from a helper function
```julia
function stress_time(c, vars, args, xnorm; ntime = 200, dt = 1.0e8)
    τ     = zeros(ntime)
    t_v   = zeros(ntime)
    τ_e   = (0.0, 0.0)
    P_e   = (0.0, 0.0)
    t     = 0.0
    for i in 2:ntime
        # non-differentiable variables needed to evaluate the state functions
        others = (; dt = dt, τ0 = τ_e, P0 = P_e)
        # initial guess and solution of this time step
        x      = initial_guess_x(c, vars, args, others)
        sol    = solve(c, x, vars, others; xnorm0 = xnorm)
        # Post-process the results
        τ_e    = compute_stress_elastic(c, sol, others)   # elastic stress
        P_e    = compute_pressure_elastic(c, sol, others) # elastic pressure
        # Store the results
        τ[i]   = τ_e[1]
        t     += others.dt
        t_v[i] = t
    end
    return t_v, τ
end

t_v, τ = stress_time(c, vars, args, xnorm; ntime = 25, dt = 1.0e9);
```

[`solve`](@ref) returns an `RCSolution`, which stores the solved static vector in
`sol.x`, along with the iteration count and final residual; positional indexing
such as `sol[1]` also works. [`inspect`](@ref) says what the four entries are:

```jldoctest
julia> using RheologyCalculator, RheologyCalculator.RheologyModels

julia> c = SeriesModel(LinearViscosity(1e21), Elasticity(10e9, 30e9),
                       ParallelModel(LinearViscosity(1e20), Elasticity(10e9, 46.67e9)));

julia> inspect(c)
4-element ModelInspection:
  index  var  equation                        scope   elements
      1  τ    compute_strain_rate             global  LinearViscosity 1, Elasticity 1
      2  ε    compute_stress                  branch  LinearViscosity 2, Elasticity 2
      3  P    compute_volumetric_strain_rate  global  LinearViscosity 1, Elasticity 1
      4  θ    compute_pressure                branch  LinearViscosity 2, Elasticity 2
```

We can finally compute the analytical solution of the stress time-evolution, and compare it against our results

```julia
using GLMakie

function simulate_series_Burgers_model(E1, η1, E2, η2, ε̇, t_max, dt)
    N    = Int(div(t_max, dt)) + 1
    t    = range(0, step = dt, length = N)
    σ    = zeros(N) # Stress
    ε_KV = zeros(N) # Strain in Kelvin–Voigt element
    σ_KV = zeros(N)
    for i in 2:N
        # Previous values
        ε_KV_prev = ε_KV[i - 1]
        σ_prev    = σ[i - 1]
        dεKVdt    = (σ_prev - E2 * ε_KV_prev) / η2  # Kelvin–Voigt strain rate
        ε_KV[i]   = ε_KV_prev + dt * dεKVdt         # Update ε_KV
        dσdt      = E1 * (ε̇ - σ_prev / η1 - dεKVdt) # Stress rate from Maxwell element
        σ[i]      = σ_prev + dt * dσdt              # Update stress
        σ_KV[i]   = E2 * ε_KV[i] + η2 * dεKVdt      # Calculate σ_KV explicitly at this step
    end
    return t, σ
end

η1 = 2 * c.leafs[1].η
η2 = 2 * c.branches[1].leafs[1].η
G1 = 2 * c.leafs[2].G
G2 = 2 * c.branches[1].leafs[2].G

t_anal, τ_anal = simulate_series_Burgers_model(G1, η1, G2, η2, vars.ε, t_v[end], (t_v[2] - t_v[1]) / 10);

# make figure
SecYear = 3600 * 24 * 365.25
fig     = Figure(fontsize = 30, size = (800, 600))
ax      = Axis(fig[1, 1], title = "Burgers model", xlabel = "t [kyr]", ylabel = L"\tau [MPa]")

lines!(ax, t_anal / SecYear / 1.0e3, τ_anal / 1.0e6, label = "analytical", linewidth = 5, color = :black)
scatter!(ax, t_v / SecYear / 1.0e3, τ / 1.0e6, label = "numerical", color = :red, markersize = 15)

axislegend(ax, position = :rb)
ax.xlabel = L"t [kyr]"
ax.ylabel = L"\tau [MPa]"
display(fig)
```

![](https://raw.githubusercontent.com/albert-de-montserrat/RheologyCalculator.jl/main/docs/assets/Burgers_model.png)

## Warm-starting from the previous solution

[`initial_guess_x`](@ref) builds the starting point of the Newton iteration from
the model itself: each element is evaluated on its own and the results are
combined, harmonically for the stress and pressure of elements in series and as a
plain sum for the strain rates of parallel branches. Plastic elements do not
contribute.

The Burgers model is linear and converges in a single iteration from this
estimate. Adding a Drucker–Prager element in series gives a model that yields
within a few steps. The loop below loads it at a constant strain rate, builds a
new initial guess at every step, and counts the Newton iterations:

```jldoctest warmstart
julia> using RheologyCalculator, RheologyCalculator.RheologyModels
       damper2, damper3 = LinearViscosity(1e20), LinearViscosity(1e21)
       spring1, spring2 = Elasticity(10e9, 46.67e9), Elasticity(10e9, 30e9)
       vars = (; ε = 1.0e-15, θ = 1.0e-20);

julia> cp = SeriesModel(damper3, spring2, ParallelModel(damper2, spring1),
                        DruckerPrager(1.0e5, 30.0, 0.0));

julia> xnorm = normalisation_x(cp, 2 * 1.0e21 * vars.ε, vars.ε);

julia> function loading(c, vars, xnorm, nsteps)
           τ_e, P_e, sol, iterations = (0.0, 0.0), (0.0, 0.0), nothing, 0
           for _ in 1:nsteps
               others = (; dt = 1.0e9, τ0 = τ_e, P0 = P_e)
               x = initial_guess_x(c, vars, (; τ = 1.0e3, P = 1.0e6), others)
               sol = solve(c, x, vars, others; xnorm0 = xnorm)
               iterations += sol.iterations
               τ_e = compute_stress_elastic(c, sol, others)
               P_e = compute_pressure_elastic(c, sol, others)
           end
           return round(sol[1]; sigdigits = 6), iterations
       end;

julia> loading(cp, vars, xnorm, 25)
(86598.8, 372)
```

The stress reaches the yield stress, but 25 steps take 372 Newton iterations,
about 15 per step. Once the material yields, the estimate, which ignores the
plastic element, is far from the solution, while the solution of the previous
step is close to the new one. [`solve`](@ref) accepts the [`RCSolution`](@ref)
it returned as the starting point of the next solve, which reuses every entry,
including the branch unknowns. Only the first step needs an initial guess:

```jldoctest warmstart
julia> function loading_warm(c, vars, xnorm, nsteps)
           τ_e, P_e, sol, iterations = (0.0, 0.0), (0.0, 0.0), nothing, 0
           for _ in 1:nsteps
               others = (; dt = 1.0e9, τ0 = τ_e, P0 = P_e)
               x = isnothing(sol) ? initial_guess_x(c, vars, (; τ = 1.0e3, P = 1.0e6), others) : sol
               sol = solve(c, x, vars, others; xnorm0 = xnorm)
               iterations += sol.iterations
               τ_e = compute_stress_elastic(c, sol, others)
               P_e = compute_pressure_elastic(c, sol, others)
           end
           return round(sol[1]; sigdigits = 6), iterations
       end;

julia> loading_warm(cp, vars, xnorm, 25)
(86598.8, 52)
```

The result is the same, but the warm-started loop needs 52 iterations instead of
372, about two per step. The Burgers loop above can start each step from the
previous solution in the same way; being linear, it converges in one iteration
per step from either start.

## Obtaining the Jacobian at the converged solution

Derivatives of the solution with respect to an input or a material parameter,
for sensitivities or gradient-based fitting for instance, can be obtained without
differentiating through the Newton iterations. Since the residual vanishes at the
converged solution, the implicit function theorem gives, for a parameter ``p``
entering the residual,
``\mathrm{d}\mathbf{x}/\mathrm{d}p = -\mathbf{J}^{-1}\,\partial \mathbf{r}/\partial p``,
where ``\mathbf{J} = \partial \mathbf{r} / \partial \mathbf{x}`` is the
Jacobian at that solution.

It is not one of the Jacobians formed during the
Newton iteration, which are all taken before the last update, so it costs one
extra evaluation. [`solve`](@ref) skips it; [`solve_with_jacobian`](@ref) returns
it in the `jacobian` field of the [`RCSolution`](@ref):

```jldoctest warmstart
julia> c = SeriesModel(damper3, spring2, ParallelModel(damper2, spring1));

julia> others = (; dt = 1.0e9, τ0 = (0.0, 0.0), P0 = (0.0, 0.0));

julia> x0 = initial_guess_x(c, vars, (; τ = 1.0e3, P = 1.0e6), others);

julia> xnorm = normalisation_x(c, 2 * 1.0e21 * vars.ε, vars.ε);

julia> sol = solve_with_jacobian(c, x0, vars, others; xnorm0 = xnorm);

julia> sol.jacobian
4×4 StaticArraysCore.SMatrix{4, 4, Float64, 16} with indices SOneTo(4)×SOneTo(4):
  5.05e-20   1.0      0.0           0.0
 -1.0        2.2e20   0.0           0.0
  0.0        0.0     -3.33333e-20   1.0
 -0.0       -0.0     -1.0          -4.667e19

julia> solve(c, x0, vars, others; xnorm0 = xnorm).jacobian === nothing
true
```

The rows and columns of `sol.jacobian` follow the order given by
[`inspect`](@ref). Plain `solve` leaves the field `nothing`, and `x`,
`iterations`, and `residual` are identical in both cases. For the assembled
material tangent, use [`tangent`](@ref) instead.
