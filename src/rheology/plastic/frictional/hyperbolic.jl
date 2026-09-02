# This implements a hyperbolic (smoothed-corner) approximation to the
# Drucker-Prager yield surface, as proposed in:
# Abbo, A. J., and Sloan, S. W.: A smooth hyperbolic approximation to the
#   Mohr-Coulomb yield criterion, Computers & Structures, 54(3), 427-441, 1995.

# Hyperbolic ------------------------------------------------------
"""
    Hyperbolic{T} <: AbstractCapPlasticity

Represents a Drucker-Prager plasticity model with a smooth hyperbolic
approximation near the tensile cutoff, avoiding the non-differentiable corner
of the classical yield surface, as described in Abbo & Sloan (1995),
Computers & Structures.

# Fields
- `C::T`: The cohesion parameter.
- `ϕ::T`: The friction angle (in degrees).
- `ψ::T`: The dilatancy angle (in degrees).
- `η_vp::T`: The Duvaut-Lions regeularisation viscosity for the plasticity model.
- `Pt::T`: The tensile strength (should be < 0).
"""
struct Hyperbolic{T} <: AbstractCapPlasticity
    C::T
    ϕ::T        # in degrees for now
    ψ::T        # in degrees for now
    η_vp::T     # regularisation viscosity
    Pt::T       # Tensile strength

    # computational parameters (precomputed, to speed up later calculations)
    sinϕ::T     # Friction angle
    cosϕ::T     # Friction angle
    sinΨ::T     # Dilation angle 
    cosΨ::T     # Dilation angle
end

function Hyperbolic(; C=10e6, ϕ=30.0, ψ=0.0, η_vp=1e20, Pt=-1e5) 
    sinϕ = sind(ϕ) # Friction angle
    cosϕ = cosd(ϕ) # Friction angle
    sinΨ = sind(ψ) # Dilation angle
    cosΨ = cosd(ψ) # Dilation angle
    return Hyperbolic(C, ϕ, ψ, η_vp, Pt, sinϕ, cosϕ, sinΨ, cosΨ)
end

#Hyperbolic(args...) = Hyperbolic(promote(args...)...)

function compute_F(r::Hyperbolic, τII, P)
    cosϕ, sinϕ, C, Pt = r.cosϕ, r.sinϕ, r.C, r.Pt

    F  = sqrt(τII^2 + (C * cosϕ + Pt*sinϕ)^2) - (P * sinϕ + C * cosϕ)

    # Note that viscoplastic regularisation is taken into account in the residual function
    return F #*(F>-1e-8) 
end

function compute_Q(r::Hyperbolic, τII, P) 

    # These parameters are required to compute the constant in the plastic flow
    # potential. Note that this constant does not matter apart when plotting,
    # as we only need derivates of Q in general 
    cosΨ, sinΨ, C, Pt = r.cosΨ, r.sinΨ, r.C, r.Pt

    Q  =  sqrt(τII^2 + (C * cosΨ + Pt*sinΨ)^2) - (P * sinΨ + C * cosΨ)

    return Q
end 

