using RheologyCalculator, StaticArrays, ForwardDiff, Chairmarks, DiffResults
using RheologyCalculator.RheologyModels
using Dates, TOML, SHA, Random, Statistics, LinearAlgebra, InteractiveUtils

const RC = RheologyCalculator
include("reference_solver.jl")

function fixture(
        name, c; vars = (; ε = 1.0e-15),
        args = (; τ = 1.0e2, P = 1.0e6), others = (;)
    )
    # The generated guess can already solve simple laws exactly. Deliberately
    # perturb it for cold-start timings; warm starts are measured separately.
    x = initial_guess_x(c, vars, args, others) .* 0.7
    xnorm = normalisation_x(c, 1.0e6, RC.second_invariant_value(vars.ε))
    corrected = merge(
        vars, (;
            ε = RC.second_invariant_value(
                vars.ε .+ RC._direct_leaf_elastic_correction(c, vars.ε, others)
            ),
        )
    )
    return (; name, c, x, vars, corrected, others, xnorm)
end

function fixtures()
    v = LinearViscosity(5.0e19)
    v2 = LinearViscosity(1.0e20)
    e = IncompressibleElasticity(1.0e10)
    history = (; dt = 1.0e10, τ0 = (0.0,), P0 = (0.0,))
    nested = ParallelModel(SeriesModel(v, v2), v2)
    return (
        fixture("scalar_linear", SeriesModel(v)),
        fixture("scalar_maxwell", SeriesModel(v, e); others = history),
        fixture(
            "scalar_powerlaw", SeriesModel(PowerLawViscosity(5.0e19, 3.5));
            args = (; τ = 10.0)
        ),
        fixture(
            "scalar_creep", SeriesModel(
                DiffusionCreep(1, 0, 0, 1.0e-15, 1.0e5, 1.0e-6, 8.314),
                DislocationCreep(3.5, 0, 1.0e-30, 2.0e5, 1.0e-6, 8.314)
            );
            args = (; τ = 1.0e4), others = (; T = 1000.0, P = 1.0e8, f = 1.0, d = 1.0)
        ),
        fixture(
            "volumetric_maxwell", SeriesModel(Elasticity(1.0e10, 1.0e12), v);
            vars = (; ε = 1.0e-15, θ = 1.0e-20), others = history
        ),
        fixture("nested", SeriesModel(v, nested)),
        fixture("siblings", SeriesModel(v, nested, ParallelModel(v, v2))),
        fixture(
            "kelvin_voigt", SeriesModel(e, ParallelModel(e, v));
            others = (; dt = 1.0e10, τ0 = (0.0, 0.0), P0 = (0.0, 0.0))
        ),
        fixture(
            "plastic", SeriesModel(LinearViscosity(1.0e22), e, DruckerPrager(1.0e6, 30.0, 0.0));
            vars = (; ε = 1.0e-14), args = (; τ = 1.0e6), others = merge(history, (; P = 1.0e6))
        ),
        fixture(
            "dilatant", SeriesModel(
                LinearViscosity(1.0e22), Elasticity(1.0e10, 1.0e12),
                DruckerPrager(1.0e6, 30.0, 30.0)
            );
            vars = (; ε = 1.0e-14, θ = 0.0), args = (; τ = 1.0e6, P = 1.0e6), others = history
        ),
    )
end

# Experiment only: exactly the current residual pipeline, with equations supplied
# by preprocessing. Keep parity checks below when the production residual changes.
function prepared_residual(c, eqs, x, vars, others)
    args = RC.generate_args_template(eqs, x, others)
    r = RC.evaluate_state_functions(eqs, args, others)
    r = RC.add_children(r, x, eqs)
    r = RC.subtract_parent(r, x, eqs, vars)
    if c isa RC.SeriesModel
        r = RC.subtract_elastic_correction(c, eqs, r, x, others)
    end
    return SVector(r)
end

# Experiment only: item 3 of the performance plan. Extract both the residual
# and the Jacobian from a single ForwardDiff dual pass via the supported
# jacobian!/DiffResults API, instead of the two separate calls the solver
# currently makes for its pre-loop residual and its first-iteration Jacobian.
# The immutable DiffResults container matches the package's static, isbits,
# allocation-free convention; storage is initialized from `x`'s shape only,
# without evaluating the residual to discover it (both are SVector{N}, so the
# system is square and DiffResults.JacobianResult(x) already assumes that).
function residual_and_jacobian(f, x::SVector)
    result = ForwardDiff.jacobian!(DiffResults.JacobianResult(x), f, x)
    return DiffResults.value(result), DiffResults.jacobian(result)
end

function prepare(q)
    f = y -> compute_residual(q.c, y, q.corrected, q.others)
    eqs = generate_equations(q.c)
    prepared = y -> prepared_residual(q.c, eqs, y, q.corrected, q.others)
    cfg = ForwardDiff.JacobianConfig(f, q.x, ForwardDiff.Chunk{length(q.x)}())
    cfg1 = ForwardDiff.JacobianConfig(f, q.x, ForwardDiff.Chunk{1}())
    sol = solve(q.c, q.x, q.vars, q.others; xnorm0 = q.xnorm)
    r = f(q.x)
    J = jacobian(q.c, q.x, q.corrected, q.others)
    return (; q, f, eqs, prepared, cfg, cfg1, sol, r, J)
end

kernel(::Val{:residual}, s) = s.f(s.q.x)
kernel(::Val{:jacobian}, s) = jacobian(s.q.c, s.q.x, s.q.corrected, s.q.others)
kernel(::Val{:residual_jacobian}, s) = (s.f(s.q.x), kernel(Val(:jacobian), s))
kernel(::Val{:residual_and_jacobian_fused}, s) = residual_and_jacobian(s.f, s.q.x)
kernel(::Val{:cached_callable}, s) = ForwardDiff.jacobian(s.f, s.q.x)
kernel(::Val{:cached_config}, s) = ForwardDiff.jacobian(s.f, s.q.x, s.cfg)
kernel(::Val{:cached_config_chunk1}, s) = ForwardDiff.jacobian(s.f, s.q.x, s.cfg1)
kernel(::Val{:prepared_equations}, s) = ForwardDiff.jacobian(s.prepared, s.q.x)
kernel(::Val{:config_setup}, s) = ForwardDiff.JacobianConfig(s.f, s.q.x, ForwardDiff.Chunk{length(s.q.x)}())
kernel(::Val{:equations_setup}, s) = generate_equations(s.q.c)
kernel(::Val{:backsolve}, s) = RC.backsolve(s.J, s.r)
kernel(::Val{:solve}, s) = solve(s.q.c, s.q.x, s.q.vars, s.q.others; xnorm0 = s.q.xnorm)
kernel(::Val{:solve_warm}, s) = solve(s.q.c, s.sol.x, s.q.vars, s.q.others; xnorm0 = s.q.xnorm)
kernel(::Val{:solve_reference}, s) = SolverReference.solve_reference(s.q.c, s.q.x, s.q.vars, s.q.others; xnorm0 = s.q.xnorm)
kernel(::Val{:solve_warm_reference}, s) = SolverReference.solve_reference(s.q.c, s.sol.x, s.q.vars, s.q.others; xnorm0 = s.q.xnorm)

function measure(kind::Val, input, seconds)
    # Ref loads keep fixed numeric inputs opaque to constant folding. Setup and
    # result boxing are outside the timed kernel; Chairmarks consumes the result.
    run = ref -> kernel(kind, ref[])
    run(input)
    b = @be $input $run seconds = seconds
    return summarize_benchmark(b)
end

function summarize_benchmark(b)
    times = [v.time * 1.0e9 for v in b.samples]
    return Dict{String, Any}(
        "median_ns" => median(times), "min_ns" => minimum(times),
        "p25_ns" => quantile(times, 0.25), "p75_ns" => quantile(times, 0.75),
        "median_bytes" => median([v.bytes for v in b.samples]),
        "median_allocs" => median([v.allocs for v in b.samples]),
        "samples" => length(times)
    )
end

function measure_jacobians(input, seconds)
    a = ref -> kernel(Val(:jacobian), ref[])
    b = ref -> kernel(Val(:cached_callable), ref[])
    c = ref -> kernel(Val(:cached_config), ref[])
    d = ref -> kernel(Val(:cached_config_chunk1), ref[])
    e = ref -> kernel(Val(:prepared_equations), ref[])
    # Chairmarks randomizes variant order within each sample. The shared timing
    # window is less susceptible to frequency drift than sequential benchmarks.
    results = @be $input $a, $b, $c, $d, $e seconds = 5seconds
    return map(summarize_benchmark, results)
end

function set_affinity()
    mask = parse(UInt, get(ENV, "RC_BENCH_AFFINITY", "0"))
    if mask != 0
        Sys.iswindows() || error("RC_BENCH_AFFINITY currently supports Windows only; use taskset on Linux")
        process = ccall((:GetCurrentProcess, "kernel32"), Ptr{Cvoid}, ())
        ok = ccall((:SetProcessAffinityMask, "kernel32"), Int32, (Ptr{Cvoid}, UInt), process, mask)
        ok != 0 || error("Could not set process affinity")
    end
    return string(mask)
end

function measure_solves(input, seconds)
    a = ref -> kernel(Val(:solve_reference), ref[])
    b = ref -> kernel(Val(:solve), ref[])
    c = ref -> kernel(Val(:solve_warm_reference), ref[])
    d = ref -> kernel(Val(:solve_warm), ref[])
    results = @be $input $a, $b, $c, $d seconds = 4seconds
    return map(summarize_benchmark, results)
end

function validate(s)
    for x in (s.q.x, s.sol.x, s.q.x .* 1.01)
        reference = ForwardDiff.jacobian(s.f, x)
        @assert s.prepared(x) == s.f(x)
        @assert ForwardDiff.jacobian(s.f, x, s.cfg) == reference
        @assert ForwardDiff.jacobian(s.f, x, s.cfg1) == reference
        @assert ForwardDiff.jacobian(s.prepared, x) == reference
        r_fused, J_fused = residual_and_jacobian(s.f, x)
        @assert r_fused == s.f(x)
        @assert J_fused == reference
        @assert isbitstype(typeof(r_fused))
        @assert isbitstype(typeof(J_fused))
    end
    @assert (@allocated residual_and_jacobian(s.f, s.q.x)) == 0
    @assert isfinite(s.sol.residual)
    @assert isbitstype(typeof(s.J))
    @assert isbitstype(typeof(s.sol))
    for (old, new) in ((:solve_reference, :solve), (:solve_warm_reference, :solve_warm))
        a, b = kernel(Val(old), s), kernel(Val(new), s)
        @assert isequal(a.x, b.x)
        @assert isequal(a.residual, b.residual)
        @assert a.iterations == b.iterations
    end
    # Outside timed regions: check whether configs change AD pass counts.
    calls = Ref(0)
    counted = x -> begin
        calls[] += 1
        s.f(x)
    end
    for width in (1, length(s.q.x))
        cfg = ForwardDiff.JacobianConfig(counted, s.q.x, ForwardDiff.Chunk{width}())
        calls[] = 0
        @assert ForwardDiff.jacobian(counted, s.q.x, cfg) == s.J
        @assert calls[] == 1
    end
    return nothing
end

function source_hash()
    files = sort([joinpath(dir, f) for (dir, _, fs) in walkdir(joinpath(@__DIR__, "..", "src")) for f in fs if endswith(f, ".jl")])
    return bytes2hex(sha256(join([relpath(f, @__DIR__) * "\n" * read(f, String) for f in files], "\n")))
end

function main()
    label = isempty(ARGS) ? "baseline" : ARGS[1]
    occursin(r"^[A-Za-z0-9_-]+$", label) || error("Use letters, digits, underscores or hyphens for the label")
    seconds = parse(Float64, get(ENV, "RC_BENCH_SECONDS", "0.15"))
    rounds = parse(Int, get(ENV, "RC_BENCH_ROUNDS", "3"))
    seconds > 0 && rounds > 0 || error("Duration and rounds must be positive")
    affinity = set_affinity()
    BLAS.set_num_threads(1)
    meta = Dict(
        "label" => label, "timestamp_utc" => string(now(UTC)),
        "julia" => string(VERSION), "forwarddiff" => string(pkgversion(ForwardDiff)),
        "staticarrays" => string(pkgversion(StaticArrays)), "chairmarks" => string(pkgversion(Chairmarks)),
        "diffresults" => string(pkgversion(DiffResults)),
        "cpu" => Sys.cpu_info()[1].model, "machine" => Sys.MACHINE,
        "threads" => Threads.nthreads(), "blas_threads" => BLAS.get_num_threads(), "affinity_mask" => affinity,
        "seconds" => seconds, "rounds" => rounds, "source_sha256" => source_hash(),
        "harness_sha256" => bytes2hex(sha256(vcat(read(@__FILE__), read(joinpath(@__DIR__, "reference_solver.jl"))))),
        "manifest_sha256" => bytes2hex(sha256(read(joinpath(@__DIR__, "Manifest.toml"))))
    )
    records = Dict{String, Any}[]
    diagnostics = Dict{String, Any}[]
    kinds = (
        :residual, :jacobian, :residual_jacobian, :residual_and_jacobian_fused,
        :cached_callable, :cached_config,
        :cached_config_chunk1, :prepared_equations, :config_setup, :equations_setup,
        :backsolve, :solve, :solve_warm, :solve_reference, :solve_warm_reference,
    )
    jacobian_kinds = (:jacobian, :cached_callable, :cached_config, :cached_config_chunk1, :prepared_equations)
    solve_kinds = (:solve_reference, :solve, :solve_warm_reference, :solve_warm)
    groups = (:jacobian_group, :solve_group, filter(k -> !(k in (jacobian_kinds..., solve_kinds...)), kinds)...)
    rng = MersenneTwister(20260914)
    for q in fixtures()
        println("Preparing ", q.name, " (N=", length(q.x), ")"); flush(stdout)
        # This includes compilation of preparation and the first solve, not import time.
        cold = @timed prepare(q)
        s = cold.value
        validate(s)
        warm = kernel(Val(:solve_warm), s)
        push!(
            diagnostics, Dict(
                "fixture" => q.name, "n" => length(q.x),
                "iterations" => s.sol.iterations, "residual" => s.sol.residual,
                "solution" => collect(s.sol.x), "warm_iterations" => warm.iterations,
                "first_prepare_s" => cold.time, "config_isbits" => isbitstype(typeof(s.cfg)),
                "equations_isbits" => isbitstype(typeof(s.eqs)),
                "config_method" => string(which(ForwardDiff.jacobian, (typeof(s.f), typeof(q.x), typeof(s.cfg))))
            )
        )
        input = Ref(s)
        for kind in kinds
            kernel(Val(kind), s)
        end
        for round in 1:rounds, group in shuffle(rng, collect(groups))
            measured = group === :jacobian_group ?
                zip(jacobian_kinds, measure_jacobians(input, seconds)) :
                group === :solve_group ? zip(solve_kinds, measure_solves(input, seconds)) :
                ((group, measure(Val(group), input, seconds)),)
            for (kind, row) in measured
                merge!(row, Dict("fixture" => q.name, "n" => length(q.x), "kernel" => string(kind), "round" => round))
                push!(records, row)
            end
        end
        println("  validated; iterations=", s.sol.iterations, ", residual=", s.sol.residual); flush(stdout)
    end
    outdir = joinpath(@__DIR__, "results")
    mkpath(outdir)
    path = joinpath(outdir, Dates.format(now(UTC), "yyyymmddTHHMMSSsss") * "_" * label * ".toml")
    ispath(path) && error("Refusing to overwrite $path")
    open(path, "w") do io
        TOML.print(io, Dict("metadata" => meta, "diagnostics" => diagnostics, "measurements" => records); sorted = true)
    end
    return println("Results: ", path)
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
