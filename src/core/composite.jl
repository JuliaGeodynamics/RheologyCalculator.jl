"""
    AbstractCompositeModel

Root type for composite rheology containers. A composite separates direct
`leafs` from nested `branches` so equation generation can traverse mixed series
and parallel networks.
"""
abstract type AbstractCompositeModel  end

@inline series_state_functions(::AbstractCompositeModel) = ()
@inline parallel_state_functions(::AbstractCompositeModel) = ()

"""
    SeriesModel(elements...)

Build a series composite. Direct rheology elements are stored in `leafs`; nested
`ParallelModel`s are stored in `branches` and traversed as branch equations.

# Example
```julia
c = SeriesModel(LinearViscosity(1e22), IncompressibleElasticity(1e10))
```
"""
struct SeriesModel{L, B} <: AbstractCompositeModel # not 100% about the subtyping here, lets see
    leafs::L     # horizontal stacking
    branches::B  # vertical stacking

    function SeriesModel(c::Vararg{Any, N}) where {N}
        leafs = series_leafs(c)
        branches = series_branches(c)
        return new{typeof(leafs), typeof(branches)}(leafs, branches)
    end
end


# Apply `f` to each element of a tuple and splice the results together. The
# unrolled form keeps the element types inferable, which the residual assembly
# depends on: `f` is a compile-time constant at every call site below.
@generated function _flatmap_state_functions(f::F, funs::NTuple{N, Any}) where {F, N}
    return quote
        @inline
        t = Base.@ntuple $N i -> f(funs[i])
        Base.IteratorsMD.flatten(t)
    end
end

for fun in (:compute_strain_rate, :compute_volumetric_strain_rate)
    @eval @inline _local_series_state_functions(::typeof($fun)) = ()
    @eval @inline _global_series_state_functions(fn::typeof($fun)) = (fn,)
end

@inline _local_series_state_functions(fn::F) where {F <: Function} = (fn,)

@inline local_series_state_functions(funs::NTuple{N, Any}) where {N} = _flatmap_state_functions(_local_series_state_functions, funs)

@inline _global_series_state_functions(::F) where {F <: Function} = ()

@inline global_series_state_functions(funs::NTuple{N, Any}) where {N} = _flatmap_state_functions(_global_series_state_functions, funs)

"""
    ParallelModel(elements...)

Build a parallel composite. Direct rheology elements are stored in `leafs`;
nested `SeriesModel`s are stored in `branches`.

# Example
```julia
c = ParallelModel(LinearViscosity(1e22), LinearViscosity(1e21))
```
"""
struct ParallelModel{L, B} <: AbstractCompositeModel # not 100% about the subtyping here, lets see
    leafs::L     # horizontal stacking
    branches::B  # vertical stacking

    function ParallelModel(c::Vararg{Any, N}) where {N}
        leafs = parallel_leafs(c)
        branches = parallel_branches(c)
        return new{typeof(leafs), typeof(branches)}(leafs, branches)
    end
end

# Split a composite's constructor arguments into direct rheology elements
# (`leafs`) and nested composites (`branches`). What counts as a branch is the
# opposite composition: a SeriesModel branches on ParallelModel and vice versa.
for (prefix, Branch) in ((:series, :ParallelModel), (:parallel, :SeriesModel))
    leafs, branches = Symbol(prefix, :_leafs), Symbol(prefix, :_branches)
    @eval begin
        @inline $leafs(c::NTuple{N, AbstractRheology}) where {N} = c
        @inline $leafs(c::AbstractRheology) = (c,)
        @inline $leafs(::$Branch) = ()
        @inline $leafs(::Tuple{}) = ()
        @inline $leafs(c::NTuple{N, Any}) where {N} = $leafs(first(c))..., $leafs(Base.tail(c))...

        @inline $branches(::NTuple{N, AbstractRheology}) where {N} = ()
        @inline $branches(::AbstractRheology) = ()
        @inline $branches(c::$Branch) = (c,)
        @inline $branches(::Tuple{}) = ()
        @inline $branches(c::NTuple{N, Any}) where {N} = $branches(first(c))..., $branches(Base.tail(c))...
    end
end

Base.size(c::Union{SeriesModel, ParallelModel}) = length(c.leafs), length(c.branches)

for fun in (:compute_stress, :compute_pressure)
    @eval _local_parallel_state_functions(::typeof($fun)) = ()
    @eval @inline _global_parallel_state_functions(fn::typeof($fun)) = (fn,)
end
@inline _local_parallel_state_functions(fn::F) where {F <: Function} = (fn,)

@inline local_parallel_state_functions(funs::NTuple{N, Any}) where {N} = _flatmap_state_functions(_local_parallel_state_functions, funs)

@inline _global_parallel_state_functions(::F) where {F <: Function} = ()

@inline global_parallel_state_functions(funs::NTuple{N, Any}) where {N} = _flatmap_state_functions(_global_parallel_state_functions, funs)

# @inline series_state_functions(c::NTuple{N, ParallelModel}) where {N} = series_state_functions(first(c))..., series_state_functions(Base.tail(c))...
@inline series_state_functions(funs::NTuple{N, Any}) where {N} = _flatmap_state_functions(series_state_functions, funs)
@inline series_state_functions(::Tuple{}) = (compute_strain_rate,)

# @inline series_state_functions(c::ParallelModel)                      = flatten_repeated_functions(parallel_state_functions(c.leafs))
@inline series_state_functions(c::ParallelModel) = flatten_repeated_functions(series_state_functions(c.leafs))


# #######################################################################
# # DEAL FIRST WITH THE SERIES PART
# #######################################################################

"""
    global_series_functions(c::SeriesModel)

Return the unique global state functions that a series composite must satisfy.
For a pure deviatoric series model this is typically `compute_strain_rate`; for
volumetric models it may also include `compute_volumetric_strain_rate`.
"""
@inline function global_series_functions(c::SeriesModel)
    fn_leafs = series_state_functions(c.leafs) |> flatten_repeated_functions |> global_series_state_functions
    fn_branches = series_state_functions(c.branches) |> flatten_repeated_functions |> global_series_state_functions
    return (fn_leafs..., fn_branches...) |> flatten_repeated_functions
end

# simplify working with it
Base.getindex(c::SeriesModel, i::Int) = c.leafs[i]
Base.getindex(c::ParallelModel, i::Int) = c.leafs[i]
