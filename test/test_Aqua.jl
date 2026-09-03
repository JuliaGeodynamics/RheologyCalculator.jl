if VERSION ≤ v"1.12.3"
    using Aqua
    const RCM = RheologyCalculator.RheologyModels

    @testset "Project extras" begin
        @test Aqua.test_project_extras(RheologyCalculator).value
    end

    @testset "Undefined exports" begin
        @test Aqua.test_undefined_exports(RheologyCalculator).value
        @test Aqua.test_undefined_exports(RCM).value
    end

    @testset "Compats" begin
        @test !Aqua.test_deps_compat(
            RheologyCalculator;
            check_julia = true,
            check_extras = false,
        ).anynonpass
    end

    @testset "Persistent tasks" begin
        Aqua.test_persistent_tasks(RheologyCalculator)
    end

    @testset "Ambiguities" begin
        @test Aqua.test_ambiguities([RheologyCalculator, RCM]).value
    end

    @testset "Piracy" begin
        @test Aqua.test_piracies(RheologyCalculator).value
    end
end
