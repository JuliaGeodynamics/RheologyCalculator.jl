# Plan item 4: aggressive one-unknown specialization. Standalone probe (not
# part of the jacobian.jl kernel matrix, which applies every kernel uniformly
# across fixtures of any N): compares the current `jacobian`/`compute_residual`
# path against a hand-rolled scalar derivative and a fused scalar value+
# derivative, for the four N=1 fixtures only. Reuses `fixtures()`/`prepare()`
# from jacobian.jl; does not call its `main()`.
using RheologyCalculator, StaticArrays, ForwardDiff, Chairmarks
using RheologyCalculator.RheologyModels
using Statistics, Printf

const RC = RheologyCalculator
include("jacobian.jl")

# Scalar derivative only, per the plan's comparison target.
function scalar_derivative(q)
    f1 = t -> compute_residual(q.c, SVector(t), q.corrected, q.others)[1]
    return ForwardDiff.derivative(f1, q.x[1])
end

# Fused scalar value+derivative from one Dual evaluation, via ForwardDiff's own
# public Dual/Tag/value/partials primitives (what `ForwardDiff.derivative`
# does internally, minus discarding the value) rather than DiffResults, which
# item 3 found unsafe for nested AD when its storage type is derived from the
# wrong side of the call.
function scalar_fused(q)
    f1 = t -> compute_residual(q.c, SVector(t), q.corrected, q.others)[1]
    x = q.x[1]
    T = typeof(ForwardDiff.Tag(f1, typeof(x)))
    dual_x = ForwardDiff.Dual{T}(x, one(x))
    y = f1(dual_x)
    return ForwardDiff.value(y), ForwardDiff.partials(y, 1)
end

function validate_scalar(q)
    reference = RC.jacobian(q.c, q.x, q.corrected, q.others)[1]
    d = scalar_derivative(q)
    @assert d == reference "scalar derivative mismatch: $d vs $reference"
    v, d2 = scalar_fused(q)
    @assert v == compute_residual(q.c, q.x, q.corrected, q.others)[1] "fused value mismatch"
    @assert d2 == reference "fused derivative mismatch: $d2 vs $reference"
    return nothing
end

function main()
    seconds = parse(Float64, get(ENV, "RC_BENCH_SECONDS", "0.15"))
    scalar_fixtures = filter(q -> length(q.x) == 1, collect(fixtures()))
    lines = String[]
    push!(lines, "# Scalar specialization probe (plan item 4)\n")
    push!(lines, "Julia $(VERSION); ForwardDiff $(pkgversion(ForwardDiff)); threads=$(Threads.nthreads()).\n")
    push!(lines, "| Fixture | jacobian() ns | scalar derivative ns | fused scalar ns | derivative vs jacobian() | fused vs jacobian() |")
    push!(lines, "| --- | ---: | ---: | ---: | --- | --- |")
    for q0 in scalar_fixtures
        s = prepare(q0)
        validate_scalar(s.q)
        a = ref -> RC.jacobian(ref[].q.c, ref[].q.x, ref[].q.corrected, ref[].q.others)
        b = ref -> scalar_derivative(ref[].q)
        c = ref -> scalar_fused(ref[].q)
        input = Ref(s)
        results = @be $input $a, $b, $c seconds = 5seconds
        rows = map(summarize_benchmark, results)
        ja, sd, sf = rows[1]["median_ns"], rows[2]["median_ns"], rows[3]["median_ns"]
        describe(base, cand) = @sprintf("%+.1f%%", 100 * (cand / base - 1))
        push!(
            lines, @sprintf(
                "| %s | %.2f | %.2f | %.2f | %s | %s |", s.q.name, ja, sd, sf,
                describe(ja, sd), describe(ja, sf)
            )
        )
        println(s.q.name, ": jacobian=", ja, "ns  scalar_derivative=", sd, "ns  fused=", sf, "ns")
    end
    outpath = joinpath(@__DIR__, "results", "scalar_specialization.md")
    write(outpath, join(lines, "\n") * "\n")
    return println("Wrote ", outpath)
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
