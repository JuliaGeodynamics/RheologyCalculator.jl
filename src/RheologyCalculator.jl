"""
    RheologyCalculator

Build and solve local rheological models assembled from viscous, elastic, and
plastic elements.

The package core provides the composite containers (`SeriesModel`,
`ParallelModel`), equation generation, solver utilities, and the state-function
interface that concrete rheologies extend. Bundled material models are available
through the `RheologyModels` submodule and can be used as templates for
application-specific material laws.
"""
module RheologyCalculator

using StaticArrays, LinearAlgebra
import ForwardDiff: ForwardDiff

import Base.IteratorsMD.flatten

include("core/rheology_types.jl")
export AbstractViscosity, AbstractPlasticity, AbstractElasticity

include("core/state_functions.jl")

include("core/composite.jl")
export CompositeModel, SeriesModel, ParallelModel

include("core/kwargs.jl")

include("equation_system/equations.jl")
export generate_equations

include("core/others.jl")

include("post_processing/post_calculations.jl")

include("equation_system/initial_guess.jl")
export initial_guess_x, x_keys

include("equation_system/normalize_x.jl")
export normalisation_x

include("equation_system/solver.jl")
export solve, NonConvergenceError

include("post_processing/strain_rate_correction.jl")
export effective_strain_rate_correction

include("display/print_rheology.jl")

# Concrete material catalogue (elements + advanced models) as a submodule.
# Accessed via `RheologyCalculator.RheologyModels` or
# `using RheologyCalculator.RheologyModels`.
include("RheologyModels.jl")

end # module RheologyCalculator
