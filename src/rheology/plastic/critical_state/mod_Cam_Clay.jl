import ForwardDiff: ForwardDiff
import ..RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
import ..RheologyCalculator: series_state_functions, parallel_state_functions, _isvolumetric
import ..RheologyCalculator: compute_strain_rate, compute_stress, compute_pressure, compute_volumetric_strain_rate, compute_volumetric_plastic_strain_rate
import ..RheologyCalculator: compute_plastic_strain_rate, compute_lambda

# ModCamClay ------------------------------------------------------
"""
    ModCamClay{T} <: AbstractPlasticity

Represents a Modified Cam-Clay model, see de Souza Neto book (p. 404)

# Fields
- `M::T`    : The flattening factor for yield (1: circle, < 1: ellipse) = 0.9
- `N::T`    : The flattening factor for potential (1: circle, < 1: ellipse), if N = M, then it's associated
- `r::T`    : The radius along pressure axis = 2e8
- `β::T`    : The asymmetry for compaction cap (1:symmetric, <1 asymmetric)  = 0.1
- `Pt::T`   : The Tensile strength
- `η_vp::T` : The regularisation viscosity
"""
struct ModCamClay{T} <: AbstractCapPlasticity
    M::T        # flattening factor for yield (1: circle, < 1: ellipse) = 0.9
    N::T        # flattening factor for potential (1: circle, < 1: ellipse)
    r::T        # radius along pressure axis = 2e8
    β::T        # asymmetry for compaction cap (1:symmetric, <1 asymmetric)  = 0.1
    Pt::T       # Tensile strength
    η_vp::T     # regularisation viscosity
end

function ModCamClay(; M=0.9, N=0.9, r=1e8, β=1.0, Pt=-1e7, η_vp=1e20) 
    return ModCamClay(M, N, r, β, Pt, η_vp)
end

#ModCamClay(args...) = ModCamClay(promote(args...)...)

function compute_F(p::ModCamClay, τII, P)
    (; M, r, β, Pt) = p

    b = P < Pt + r ? one(β) : β

    F  = 1/b *(P - Pt - r)^2  + (τII / M)^2 - r^2 

    # Note that viscoplastic regularisation is taken into account in the residual function
    return F
end

function compute_Q(p::ModCamClay, τII, P) 

    # These parameters are required to compute the constant in the plastic flow
    # potential. Note that this constant does not matter apart when plotting,
    # as we only need derivates of Q in general 
    (; N, r, β, Pt) = p

    b = P < Pt + r ? one(β) : β

    Q  = 1/b *(P - Pt - r)^2  + (τII / N)^2 - r^2 

    return Q
end 

