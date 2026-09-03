# Sparsity detection for a coupled Stokes system

For a large nonlinear system, forming a numerical Jacobian just to discover
which entries can be nonzero is expensive. `SparseConnectivityTracer.jl`
instead propagates dependency sets: it answers *which residual entries depend
on which unknowns?* without evaluating derivatives.

This tutorial builds the sparsity pattern of a collocated two-dimensional
thermo-mechanical Stokes residual. The local constitutive relation is supplied
by `RheologyCalculator.jl` at every grid point.

In addition to the `RheologyCalculator` dependency, you will need `SparseConnectivityTracer` to run this example:

```julia
using Pkg
Pkg.add([ "SparseConnectivityTracer"])

using RheologyCalculator
using RheologyCalculator.RheologyModels
using SparseConnectivityTracer
```

## What is being traced?

At every grid point, the global unknown vector contains horizontal velocity,
vertical velocity, pressure, and temperature:

```math
u = [\operatorname{vec}(V_x);\ \operatorname{vec}(V_y);\
     \operatorname{vec}(P);\ \operatorname{vec}(T)].
```

The first three residual blocks are momentum and mass conservation. The fourth
is steady thermal diffusion. Temperature enters the local constitutive solve
through `others.T`, so the resulting stress can depend on both the local strain
rate and temperature. For this tutorial, the local constitutive relation is a simple diffusion creep model. But keep in mind that any local model can be used, including viscoelasticity, plasticity, or a combination of multiple mechanisms and that the performance of the sparsity detection is independent of the complexity of the local model.

```julia
const diffusion = DiffusionCreep(1, 0.0, 0.0, 1.0e-15, 200e3, 0.0, 8.314)
const model = SeriesModel(diffusion)

const dt   = 1.0e10
const τ0ij = ((0.0, 0.0, 0.0),)
const P0   = (0.0,)
const T_ref = 1_573.15

function local_stress(exx, eyy, exy, temperature)
    vars = (; ε = (exx, eyy, exy), θ = zero(exx))
    args = (; τ = 2.0e6, P = 1.0e6)
    others = (; T = temperature, P = 1.0e6, f = 1.0, d = 1.0e-3,
        dt, τ0 = τ0ij, P0)

    x = initial_guess_x(model, vars, args, others)
    τII = solve(model, x, vars, others, verbose = false)[1]

    εII = sqrt(0.5 * (exx^2 + eyy^2 + (-exx - eyy)^2) + exy^2)
    η = τII / (2 * εII + eps(Float64))
    return 2 .* vars.ε .* η
end
```

## Assemble the residual

The following residual uses centred differences in the interior and Dirichlet
conditions on the boundary. This is a simple non-optimised implementation to illustrate the idea. `local_stress` is called once per grid point.

```julia
function stokes_residual(u::AbstractVector{T}, nx, ny, Δx, Δy) where {T}
    N = nx * ny
    Vx = reshape(view(u, 1:N), nx, ny)
    Vy = reshape(view(u, (N + 1):(2N)), nx, ny)
    P  = reshape(view(u, (2N + 1):(3N)), nx, ny)
    temperature = reshape(view(u, (3N + 1):(4N)), nx, ny)

    Exx = Matrix{T}(undef, nx, ny)
    Eyy = Matrix{T}(undef, nx, ny)
    Exy = Matrix{T}(undef, nx, ny)
    for j in 1:ny, i in 1:nx
        ip, im = min(i + 1, nx), max(i - 1, 1)
        jp, jm = min(j + 1, ny), max(j - 1, 1)
        Exx[i, j] = (Vx[ip, j] - Vx[im, j]) / ((ip - im) * Δx)
        Eyy[i, j] = (Vy[i, jp] - Vy[i, jm]) / ((jp - jm) * Δy)
        Exy[i, j] = 0.5 * ((Vx[i, jp] - Vx[i, jm]) / ((jp - jm) * Δy) +
                           (Vy[ip, j] - Vy[im, j]) / ((ip - im) * Δx))
    end

    Txx = Matrix{T}(undef, nx, ny)
    Tyy = Matrix{T}(undef, nx, ny)
    Txy = Matrix{T}(undef, nx, ny)
    for j in 1:ny, i in 1:nx
        Txx[i, j], Tyy[i, j], Txy[i, j] = local_stress(
            Exx[i, j], Eyy[i, j], Exy[i, j], temperature[i, j],
        )
    end

    Rx = Matrix{T}(undef, nx, ny)
    Ry = Matrix{T}(undef, nx, ny)
    Rp = Matrix{T}(undef, nx, ny)
    RT = Matrix{T}(undef, nx, ny)
    for j in 1:ny, i in 1:nx
        if i == 1 || i == nx || j == 1 || j == ny
            Rx[i, j] = Vx[i, j]
            Ry[i, j] = Vy[i, j]
            Rp[i, j] = P[i, j]
            RT[i, j] = temperature[i, j] - T_ref
        else
            Rx[i, j] = (Txx[i + 1, j] - Txx[i - 1, j]) / (2Δx) +
                       (Txy[i, j + 1] - Txy[i, j - 1]) / (2Δy) -
                       (P[i + 1, j] - P[i - 1, j]) / (2Δx)
            Ry[i, j] = (Tyy[i, j + 1] - Tyy[i, j - 1]) / (2Δy) +
                       (Txy[i + 1, j] - Txy[i - 1, j]) / (2Δx) -
                       (P[i, j + 1] - P[i, j - 1]) / (2Δy)
            Rp[i, j] = (Vx[i + 1, j] - Vx[i - 1, j]) / (2Δx) +
                       (Vy[i, j + 1] - Vy[i, j - 1]) / (2Δy)
            RT[i, j] = (temperature[i + 1, j] - 2temperature[i, j] + temperature[i - 1, j]) / Δx^2 +
                       (temperature[i, j + 1] - 2temperature[i, j] + temperature[i, j - 1]) / Δy^2
        end
    end

    return vcat(vec(Rx), vec(Ry), vec(Rp), vec(RT))
end
```

## Detect the Jacobian pattern

Use the out-of-place `jacobian_sparsity(f, u, detector)` interface because
`stokes_residual` returns a new residual vector. The numerical values in `u0`
only establish array sizes here; the tracer records dependencies, not numerical
derivatives.

```julia
function main(n)
    nx, ny = n, n
    Δx, Δy = 1 / (nx - 1), 1 / (ny - 1)

    u0 = zeros(4 * nx * ny)
    u0[(3 * nx * ny + 1):end] .= T_ref
    residual(u) = stokes_residual(u, nx, ny, Δx, Δy)

    detector = TracerSparsityDetector(; gradient_pattern_type = Set{Int})
    t = @elapsed pattern = jacobian_sparsity(residual, u0, detector)
    println("t = $(round(t, digits=2)) seconds")
    display(pattern)
    return pattern
end

pattern = main(300);
```

For `n = 300`, this returns, on my machine:

```
t = 1.34 seconds
360000×360000 SparseArrays.SparseMatrixCSC{Bool, Int64} with 5066612 stored entries:
⎡⠻⣦⡘⢦⡀⠳⣄⠙⢦⡀⎤
⎢⠲⣌⠻⣦⡙⢦⡈⠳⣄⠙⎥
⎢⢤⡈⠳⣌⠛⢄⠙⠂⠈⠓⎥
⎢⠀⠙⠦⠈⠳⠀⠑⣤⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⎦
```

We can see the sparsity pattern is block-structured, with each block corresponding to the four residuals (momentum x, momentum y, mass conservation, and heat equation) and their dependencies on the four unknowns (velocity x, velocity y, pressure, and temperature).


## Why the local solve works under a tracer

Global tracers carry dependency information but no primal numerical value.
The local Newton iteration therefore cannot run normally: it needs numerical
Jacobians and a linear solve. The `SparseConnectivityTracer` extension instead
short-circuits `solve` and returns a conservative dependency pattern for every
local unknown. It collects traced values from `x`, `vars`, and `others`, which
is why `temperature` in `others.T` contributes to the global pattern.

This local approximation can contain extra entries when the true local solve
has conditional or decoupled dependencies, but it does not omit dependencies.
That conservative property is what a sparse Jacobian preallocation or coloring
workflow needs.

## Choosing the dependency-set representation

The default tracer in `SparseConnectivityTracer` uses `BitSet` to represent dependency sets. For a very large grid with small local stencils, this needs to be changed to `Set{Int}` to avoid excessive memory usage.
The `BitSet` representation is efficient for small problems, but it scales with the total number of unknowns, which can be prohibitive for large systems. In contrast,
`Set{Int}` uses less memory because it stores only the active input
indices and scales linearly with the number of dependencies.
