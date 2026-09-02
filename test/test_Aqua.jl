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

@testset "Stale dependencies" begin
    Aqua.test_stale_deps(RheologyCalculator)
end

@testset "Persistent tasks" begin
    Aqua.test_persistent_tasks(RheologyCalculator)
end

@testset "Ambiguities" begin
    # RheologyModels is a submodule, which Aqua's ambiguity check does not
    # accept as a target; methods defined there extend generics owned by
    # RheologyCalculator and are covered through it.
    @test Aqua.test_ambiguities(RheologyCalculator).value
end

@testset "Piracy" begin
    @test Aqua.test_piracies(RheologyCalculator).value
end
