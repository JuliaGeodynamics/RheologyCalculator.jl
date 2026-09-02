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

# does not allocate:
# @inline series_state_functions(r::NTuple{N, AbstractRheology}) where {N} = series_state_functions(first(r))..., series_state_functions(Base.tail(r))...
@generated function series_state_functions(r::NTuple{N, AbstractRheology}) where {N} 
    return quote
        @inline
        f = Base.@ntuple $N i -> series_state_functions(r[i]) 
        Base.IteratorsMD.flatten(f)
    end
end

# @inline series_state_functions(::Tuple{}) = ()

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
# @inline parallel_state_functions(r::NTuple{N, AbstractRheology}) where {N} = parallel_state_functions(first(r))..., parallel_state_functions(Base.tail(r))...
@generated function parallel_state_functions(r::NTuple{N, AbstractRheology}) where {N} 
    return quote
        @inline
        f = Base.@ntuple $N i -> parallel_state_functions(r[i]) 
        Base.IteratorsMD.flatten(f)
    end
end
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
@generated function flatten_repeated_functions(funs::NTuple{N, Any}) where {N}
    return quote
        @inline
        f = Base.@ntuple $N i -> i == 1 ? (funs[1],) : (funs[i] ∉ funs[1:(i - 1)] ? (funs[i],) : ())
        Base.IteratorsMD.flatten(f)
    end
end
