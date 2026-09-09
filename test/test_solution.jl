@testset "RCSolution" begin
    ε, dt = 1.0e-14, 1.0e10

    # `:τ` twice: once for the series composite, once for the Maxwell branch of
    # the parallel element
    c = SeriesModel(
        LinearViscosity(1.0e18),
        ParallelModel(
            SeriesModel(LinearViscosity(1.0e19), IncompressibleElasticity(1.0e10)),
            LinearViscosity(1.0e18),
        ),
    )
    vars = (; ε)
    others = (; dt, τ0 = (0.0,))
    x0 = initial_guess_x(c, vars, (; τ = 1.0), others)
    sol = solve(c, x0, vars, others)

    @testset "runs on a GPU" begin
        # A solution holds numbers only, so that it can be built and read inside
        # a GPU kernel.
        @test isbitstype(typeof(sol))
    end

    @testset "vector interface" begin
        @test sol isa AbstractVector{Float64}
        @test length(sol) == length(inspect(c))
        @test sol.x isa typeof(x0)
        @test collect(sol) == collect(sol.x)
        @test sol[1] == sol.x[1]
        @test propertynames(sol) == (:x, :iterations, :residual)
        # a solution can seed the next solve
        @test solve(c, sol, vars, others) ≈ sol
    end

    @testset "convergence diagnostics" begin
        @test sol.iterations ≥ 1
        @test isfinite(sol.residual)
        @test sol.residual ≤ 1.0e-6
    end

    @testset "inspect describes the entries" begin
        entries = inspect(c)
        @test entries isa AbstractVector
        @test [e.var for e in entries] == collect(x_keys(c)) == [:τ, :ε, :τ]
        @test [e.equation for e in entries] ==
            [:compute_strain_rate, :compute_stress, :compute_strain_rate]

        # the global `:τ` is the stress of the composite as a whole; the other
        # belongs to the Maxwell sub-branch
        @test [e.isglobal for e in entries] == [true, false, false]
        @test findfirst(e -> e.var === :τ && e.isglobal, entries) == stress_index(c)
        @test entries[3].elements ==
            (:LinearViscosity => 3, :IncompressibleElasticity => 1)

        # a Kelvin chain repeats `:ε`; the elements each equation spans are what
        # separate the two branches
        kelvin = SeriesModel(
            LinearViscosity(1.0e18),
            ParallelModel(LinearViscosity(1.0e19), IncompressibleElasticity(1.0e10)),
            ParallelModel(LinearViscosity(1.0e20), IncompressibleElasticity(1.0e10)),
        )
        branches = inspect(kelvin)
        @test [e.var for e in branches] == [:τ, :ε, :ε]
        @test branches[2].elements == (:LinearViscosity => 2, :IncompressibleElasticity => 1)
        @test branches[3].elements == (:LinearViscosity => 3, :IncompressibleElasticity => 2)
    end

    @testset "inspect displays as a table" begin
        rendered = sprint(show, MIME"text/plain"(), inspect(c))
        lines = split(rendered, '\n')

        @test lines[1] == "3-element ModelInspection:"
        @test lines[2] == "  index  var  equation             scope   elements"
        @test lines[5] ==
            "      3  τ    compute_strain_rate  branch  LinearViscosity 3, IncompressibleElasticity 1"
        # the last column is left unpadded
        @test !any(endswith(line, " ") for line in lines)
    end

    @testset "inference" begin
        @test (@inferred solve(c, x0, vars, others)) isa RCSolution
    end
end
