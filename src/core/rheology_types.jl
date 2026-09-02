"""
    AbstractRheology

Root type for all rheological elements. Concrete rheologies extend this type and
provide `series_state_functions`, `parallel_state_functions`, and state-function
methods such as `compute_strain_rate` or `compute_stress`.
"""
abstract type AbstractRheology end

"""
    AbstractPlasticity <: AbstractRheology

Supertype for plastic yield or flow-rule elements.
"""
abstract type AbstractPlasticity <: AbstractRheology end # in case we need spacilization at some point

"""
    AbstractCapPlasticity <: AbstractPlasticity

Supertype for plastic elements whose yield surface closes in pressure and whose
flow rule is differentiated from a potential. A subtype supplies its struct,
`compute_F(r, τII, P)`, and `compute_Q(r, τII, P)`; the state functions, the
Duvaut-Lions regularization, and the flow rule are inherited from
`src/rheology/plastic/cap_plasticity.jl`.
"""
abstract type AbstractCapPlasticity <: AbstractPlasticity end

"""
    AbstractElasticity <: AbstractRheology

Supertype for elastic elements. Elastic rheologies may consume history fields
such as `τ0`, `P0`, and `dt`.
"""
abstract type AbstractElasticity <: AbstractRheology end # in case we need spacilization at some point

"""
    AbstractViscosity <: AbstractRheology

Supertype for viscous creep or viscosity elements.
"""
abstract type AbstractViscosity <: AbstractRheology end # in case we need spacilization at some point

## METHODS FOR SERIES MODELS
"""
    series_state_functions(r)
    series_state_functions(r, num)

Return the state functions used when `r` participates in a `SeriesModel`.
Concrete rheologies should specialize this method. The `num` method additionally
returns equation and element numbering metadata used during equation generation.
"""
function series_state_functions(r::NTuple{N, AbstractRheology}, num::MVector{N, Int}) where {N}
    statefuns = (series_state_functions(first(r))..., series_state_functions(Base.tail(r))...)
    len = ntuple(i -> length(series_state_functions(r[i])), N)
    statenum = ntuple(i -> val(i, len, num), Val(sum(len)))
    stateelements = ntuple(i -> val_element(i, len), Val(sum(len)))

    return statefuns, statenum, stateelements
end

@inline series_state_functions(r::NTuple{N, AbstractRheology}) where {N} = flatmaptuple(series_state_functions, r)

# Fallbacks
@inline series_state_functions(::AbstractRheology) = error("Rheology not defined")
@inline parallel_state_functions(::AbstractRheology) = error("Rheology not defined")

"""
    parallel_state_functions(r)
    parallel_state_functions(r, num)

Return the state functions used when `r` participates in a `ParallelModel`.
Concrete rheologies should specialize this method. The `num` method additionally
returns equation and element numbering metadata used during equation generation.
"""
@inline parallel_state_functions(r::NTuple{N, AbstractRheology}) where {N} = flatmaptuple(parallel_state_functions, r)
@inline parallel_state_functions(::Tuple{}) = ()

function parallel_state_functions(r::NTuple{N, AbstractRheology}, num::MVector{N, Int}) where {N}
    statefuns = (parallel_state_functions(first(r))..., parallel_state_functions(Base.tail(r))...)

    len = ntuple(i -> length(parallel_state_functions(r[i])), N)
    statenum = ntuple(i -> val(i, len, num), Val(sum(len)))
    stateelements = ntuple(i -> val_element(i, len), Val(sum(len)))

    return statefuns, statenum, stateelements
end

"""
    flatten_repeated_functions(funs)

Return `funs` with duplicate function objects removed while preserving the first
occurrence order. This keeps equation generation from emitting repeated global
state equations.
"""
# Must stay `@generated`: the membership test needs the prefix `funs[1:(i-1)]`
# as a literal at specialisation time. Expressed as a recursion carrying the
# entries seen so far, `f ∈ seen` is deferred to runtime over a tuple of
# distinct function types, which stops inferring and allocates.
@generated function flatten_repeated_functions(funs::NTuple{N, Any}) where {N}
    return quote
        @inline
        f = Base.@ntuple $N i -> i == 1 ? (funs[1],) : (funs[i] ∉ funs[1:(i - 1)] ? (funs[i],) : ())
        Base.IteratorsMD.flatten(f)
    end
end
