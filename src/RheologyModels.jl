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

# Names the element files below extend. Importing them once here, rather than
# per file, is what makes an `include`d method extend the core generic: a file
# that defines `compute_viscosity_series` without the import would instead
# create a function of that name local to this module, leaving the core generic
# on its fallback for that element.
import ..RheologyCalculator: series_state_functions, parallel_state_functions,
    _isvolumetric, isvolumetric,
    compute_strain_rate, compute_stress, compute_pressure, compute_volumetric_strain_rate,
    compute_plastic_strain_rate, compute_plastic_stress, compute_volumetric_plastic_strain_rate,
    compute_lambda, compute_lambda_parallel,
    compute_viscosity, compute_viscosity_series, compute_viscosity_parallel

# used by the tensor helpers
import ..RheologyCalculator: second_invariant

# Basic viscous/elastic/plastic building blocks, one struct per file, grouped by
# category (and, where the literature draws a clear line, by sub-family).
# Advanced / application-specific models (built from these blocks) sit alongside
# their category and are intentionally NOT exported. This applies to the
# material models only; the tensor helpers at the end of this file are exported,
# since every example and test needs them to build a `vars` tuple.

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
include("rheology/plastic/cap_plasticity.jl")
include("rheology/plastic/frictional/Drucker_Prager.jl")
include("rheology/plastic/frictional/Drucker_Prager_cap.jl")
include("rheology/plastic/frictional/hyperbolic.jl")
include("rheology/plastic/critical_state/mod_Cam_Clay.jl")
include("rheology/plastic/critical_state/Golchin.jl")

include("utils/tensor_helpers.jl")

export LinearViscosity, PowerLawViscosity
export Elasticity, BulkElasticity, IncompressibleElasticity, BulkViscosity
export LTPViscosity, DruckerPrager, DiffusionCreep, DislocationCreep

# Tensor helpers: flat-tuple deviatoric tensors in Voigt order, used to build
# the `vars` and `τ0` arguments of `solve`.
export second_invariant_2D, second_invariant_3D
export tensor_strain_rate_2D, tensor_strain_rate_3D
export vars_2D, vars_3D
export zero_stress_tensor_2D, zero_stress_tensor_3D
export stress_tensor_from_invariant_2D, stress_tensor_from_invariant_3D
export elastic_stress_history_2D, elastic_stress_history_3D

end
