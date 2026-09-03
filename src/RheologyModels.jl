"""
    RheologyModels

Bundled viscous, elastic, and plastic material models for use with
`RheologyCalculator`.
"""
module RheologyModels

using ..RheologyCalculator
using StaticArrays, LinearAlgebra
import ForwardDiff: ForwardDiff

# postprocessing helpers used by the tensor helpers below and by downstream code
import ..RheologyCalculator: compute_stress_elastic, compute_pressure_elastic
export compute_stress_elastic, compute_pressure_elastic

# Basic viscous/elastic/plastic building blocks, one struct per file, grouped by
# category (and, where the literature draws a clear line, by sub-family).
# Advanced / application-specific models (built from these blocks) sit alongside
# their category and are intentionally NOT exported.

# viscous: Newtonian (linear stress-strain-rate) vs non-Newtonian (rate-dependent).
# Rate-and-state friction lives here too: despite the filename, it has no yield
# surface or plastic multiplier -- compute_strain_rate is a smooth, invertible
# function of stress at every stress level, exactly like the other non-Newtonian
# (sinh-type) laws, so per its own AbstractViscosity supertype it belongs here.
include("rheology/viscous/Newtonian/linear_viscosity.jl")
include("rheology/viscous/Newtonian/bulk_viscosity.jl")
include("rheology/viscous/Newtonian/diffusion_creep.jl")
include("rheology/viscous/nonNewtonian/power_law_viscosity.jl")
include("rheology/viscous/nonNewtonian/LTP_viscosity.jl")
include("rheology/viscous/nonNewtonian/dislocation_creep.jl")
include("rheology/viscous/nonNewtonian/rate_state_hypo_plastic.jl")

include("rheology/elastic/elasticity.jl")
include("rheology/elastic/bulk_elasticity.jl")
include("rheology/elastic/incompressible_elasticity.jl")

# plastic: frictional/Coulomb-type family (cohesion + friction + dilation angle;
# open cone yield surface, with cap/smoothed-corner variants) vs critical-state
# family (M/N + tensile/compaction pressure parametrization; closed elliptical
# yield surface, e.g. Cam-Clay and its extensions).
include("rheology/plastic/frictional/Drucker_Prager.jl")
include("rheology/plastic/frictional/Drucker_Prager_cap.jl")
include("rheology/plastic/frictional/hyperbolic.jl")
include("rheology/plastic/critical_state/mod_Cam_Clay.jl")
include("rheology/plastic/critical_state/Golchin.jl")

include("utils/tensor_helpers.jl")

export LinearViscosity, PowerLawViscosity
export Elasticity, BulkElasticity, IncompressibleElasticity, BulkViscosity
export LTPViscosity, DruckerPrager, DiffusionCreep, DislocationCreep

end
