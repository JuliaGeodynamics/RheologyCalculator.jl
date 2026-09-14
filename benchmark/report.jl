using TOML, Statistics, Printf

function groups(data)
    result = Dict{Tuple{String, String}, Vector{Dict{String, Any}}}()
    for row in data["measurements"]
        push!(get!(result, (row["fixture"], row["kernel"]), Dict{String, Any}[]), row)
    end
    return result
end

middle(rows, key = "median_ns") = median([row[key] for row in rows])

function describe_change(a, b)
    pct = 100 * (middle(b) / middle(a) - 1)
    arange = extrema([r["median_ns"] for r in a])
    brange = extrema([r["median_ns"] for r in b])
    overlap = max(arange[1], brange[1]) <= min(arange[2], brange[2])
    status = abs(pct) < 5 || overlap ? "inconclusive" : pct < 0 ? "faster" : "slower"
    return @sprintf("%+.1f%% (%s)", pct, status)
end

function main()
    length(ARGS) in (1, 2) || error("Usage: report.jl BASELINE.toml [CANDIDATE.toml]")
    base = TOML.parsefile(ARGS[1])
    candidate = length(ARGS) == 2 ? TOML.parsefile(ARGS[2]) : base
    a, b = groups(base), groups(candidate)
    m = candidate["metadata"]
    println("# Benchmark report: ", m["label"], "\n")
    println(
        "Julia ", m["julia"], "; ForwardDiff ", m["forwarddiff"], "; ", m["cpu"],
        "; threads=", m["threads"], "; rounds=", m["rounds"], ".\n"
    )
    println("Times are medians of round medians, in nanoseconds. Negative changes mean faster.")
    println("A change is marked inconclusive below 5% or when round-median ranges overlap; this is a screening rule, not a statistical confidence interval.\n")
    return if length(ARGS) == 2
        for key in ("julia", "forwarddiff", "staticarrays", "chairmarks", "cpu", "machine", "threads", "blas_threads", "affinity_mask", "seconds", "rounds", "harness_sha256", "manifest_sha256")
            base["metadata"][key] == m[key] || error("Incomparable metadata: $key differs; establish a new baseline")
        end
        Set(keys(a)) == Set(keys(b)) || error("Fixture/kernel coverage differs")
        println("| Fixture | Kernel | Baseline ns | Candidate ns | Change | Bytes before → after |")
        println("| --- | --- | ---: | ---: | --- | ---: |")
        for key in sort(collect(keys(a)))
            @printf(
                "| %s | %s | %.2f | %.2f | %s | %.0f → %.0f |\n", key...,
                middle(a[key]), middle(b[key]), describe_change(a[key], b[key]),
                middle(a[key], "median_bytes"), middle(b[key], "median_bytes")
            )
        end
        println("\n| Fixture | Iterations before → after | Residual before → after | Max absolute solution change |")
        println("| --- | ---: | ---: | ---: |")
        before = Dict(d["fixture"] => d for d in base["diagnostics"])
        for d in candidate["diagnostics"]
            old = before[d["fixture"]]
            @printf(
                "| %s | %d → %d | %.3e → %.3e | %.3e |\n", d["fixture"],
                old["iterations"], d["iterations"], old["residual"], d["residual"],
                maximum(abs.(old["solution"] .- d["solution"]))
            )
        end
    else
        println("| Fixture | N | Residual ns | Jacobian ns | Solve ns | Warm solve ns | Iterations |")
        println("| --- | ---: | ---: | ---: | ---: | ---: | ---: |")
        for d in candidate["diagnostics"]
            name = d["fixture"]
            @printf(
                "| %s | %d | %.2f | %.2f | %.2f | %.2f | %d |\n", name, d["n"],
                (middle(b[(name, k)]) for k in ("residual", "jacobian", "solve", "solve_warm"))..., d["iterations"]
            )
        end
        if haskey(b, (candidate["diagnostics"][1]["fixture"], "solve_reference"))
            println("\nInterleaved full-solve comparison against the frozen original solver:\n")
            println("| Fixture | Original cold ns | Current cold ns | Cold change | Warm change | Current cold bytes |")
            println("| --- | ---: | ---: | --- | --- | ---: |")
            for d in candidate["diagnostics"]
                name = d["fixture"]
                old, new = b[(name, "solve_reference")], b[(name, "solve")]
                @printf(
                    "| %s | %.2f | %.2f | %s | %s | %.0f |\n", name,
                    middle(old), middle(new), describe_change(old, new),
                    describe_change(b[(name, "solve_warm_reference")], b[(name, "solve_warm")]),
                    middle(new, "median_bytes")
                )
            end
        end
        println("\nPreprocessing experiments below compare Jacobian kernels only, excluding setup.\n")
        println("| Fixture | Cached callable | Cached config | Chunk-1 config | Prepared equations | Config setup ns / bytes |")
        println("| --- | --- | --- | --- | --- | ---: |")
        for d in candidate["diagnostics"]
            name = d["fixture"]
            baseline = b[(name, "jacobian")]
            setup = b[(name, "config_setup")]
            @printf(
                "| %s | %s | %s | %s | %s | %.2f / %.0f |\n", name,
                (
                    describe_change(baseline, b[(name, k)]) for k in
                        ("cached_callable", "cached_config", "cached_config_chunk1", "prepared_equations")
                )...,
                middle(setup), middle(setup, "median_bytes")
            )
        end
        if haskey(b, (candidate["diagnostics"][1]["fixture"], "residual_and_jacobian_fused"))
            println("\nResidual+Jacobian fusion experiment (item 3): a single ForwardDiff.jacobian!/DiffResults")
            println("call against the solver's current separate primal-then-Jacobian calls, and against a")
            println("Jacobian-only call, all measured within this run.\n")
            println("| Fixture | Separate (residual+jacobian) ns | Fused ns | Fused vs separate | Jacobian-only ns | Fused vs jacobian-only | Fused bytes |")
            println("| --- | ---: | ---: | --- | ---: | --- | ---: |")
            for d in candidate["diagnostics"]
                name = d["fixture"]
                separate, fused, jac_only = b[(name, "residual_jacobian")], b[(name, "residual_and_jacobian_fused")], b[(name, "jacobian")]
                @printf(
                    "| %s | %.2f | %.2f | %s | %.2f | %s | %.0f |\n", name,
                    middle(separate), middle(fused), describe_change(separate, fused),
                    middle(jac_only), describe_change(jac_only, fused),
                    middle(fused, "median_bytes")
                )
            end
        end
        println("\n## Allocations\n")
        nonzero = [(key, rows) for (key, rows) in b if middle(rows, "median_bytes") > 0]
        if isempty(nonzero)
            println("All measured kernels have zero median allocated bytes.")
        else
            println("Kernels with nonzero median allocated bytes (all others are zero):\n")
            for (key, rows) in sort(nonzero; by = first)
                @printf("- %s / %s: %.0f bytes\n", key..., middle(rows, "median_bytes"))
            end
        end
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
