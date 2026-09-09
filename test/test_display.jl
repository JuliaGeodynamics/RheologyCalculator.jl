import RheologyCalculator: print_rheology_matrix, remove_colors_string

@testset "display" begin
    viscous1 = LinearViscosity(1.0e22)
    viscous2 = LinearViscosity(1.0e21)
    viscous3 = LinearViscosity(5.0e20)
    elastic = IncompressibleElasticity(30.0e9)
    elastic2 = IncompressibleElasticity(20.0e9)
    volumetric = Elasticity(1.0e10, 1.0e12)
    bulkvisc = BulkViscosity(1.0e18)
    plastic = DruckerPrager(1.0e6, 30.0, 0.0)

    composites = (
        SeriesModel(viscous1, viscous2),
        SeriesModel(elastic, viscous1),
        ParallelModel(viscous1, viscous2),
        ParallelModel(elastic, viscous1),
        SeriesModel(elastic, ParallelModel(elastic2, viscous1)),
        SeriesModel(viscous1, ParallelModel(SeriesModel(viscous2, elastic), viscous2)),
        SeriesModel(viscous1, ParallelModel(viscous1, viscous2), ParallelModel(viscous2, viscous3)),
        SeriesModel(viscous1, elastic, plastic),
        SeriesModel(volumetric, bulkvisc),
        ParallelModel(SeriesModel(viscous1, elastic), SeriesModel(viscous2, elastic2)),
        SeriesModel(viscous1, ParallelModel(elastic, SeriesModel(viscous2, elastic2))),
        # more elements than the single-digit superscripts cover
        SeriesModel(ntuple(i -> LinearViscosity(1.0e20 * i), 12)...),
    )

    @testset "every composite renders" begin
        for c in composites
            str = repr(c)
            @test !isempty(str)
            # one line per row of the rendered block, each carrying element art
            @test all(!isempty, split(rstrip(str, '\n'), '\n'))
        end
    end

    @testset "series and parallel read differently" begin
        # a parallel block is bracketed and stacks vertically; a series block
        # does not bracket and runs horizontally
        series = remove_colors_string(repr(SeriesModel(viscous1, viscous2)))
        parallel = remove_colors_string(repr(ParallelModel(viscous1, viscous2)))
        @test count(==('\n'), series) == 1
        @test count(==('\n'), parallel) == 2
        @test !occursin('|', series)
        @test occursin('|', parallel)
    end

    @testset "element numbering" begin
        # superscripts number each element type independently, and switch to two
        # digits once any counter passes 9
        @test occursin('¹', remove_colors_string(repr(SeriesModel(viscous1, elastic))))
        @test occursin('²', remove_colors_string(repr(SeriesModel(viscous1, viscous2))))
        many = remove_colors_string(repr(SeriesModel(ntuple(i -> LinearViscosity(1.0e20 * i), 12)...)))
        @test occursin("¹⁰", many)
    end

    @testset "print_rheology_matrix returns lines" begin
        for c in composites
            block = print_rheology_matrix(c)
            @test block isa AbstractVector
            @test !isempty(block)
            @test all(s -> s isa AbstractString, block)
        end
    end

    @testset "composites larger than a screen render" begin
        # The layout buffer is sized from the element count, so there is no
        # fixed ceiling on how many elements a composite may draw.
        wide = SeriesModel(ntuple(i -> LinearViscosity(1.0e20 * i), 45)...)
        @test !isempty(repr(wide))
        deep = SeriesModel(ntuple(i -> ParallelModel(LinearViscosity(1.0e20), IncompressibleElasticity(1.0e10)), 15)...)
        @test !isempty(repr(deep))
    end
end
