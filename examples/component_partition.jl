# Strain-rate partitioning of a visco-elasto-plastic point across six orders of
# magnitude of imposed strain rate.
#
# Below yield the stress rises linearly with the imposed rate and the split between
# the viscous dashpot and the elastic spring is fixed by the ratio η/(G*dt). Once
# the Drucker-Prager yield stress is reached the stress saturates, and every further
# increment of imposed strain rate is taken up by plastic flow. `component_partition`
# is what makes that handover visible: the stress alone cannot show it, because the
# stress stops changing exactly where the interesting behaviour starts.
#
# Regenerates docs/assets/component_partition.png

using RheologyCalculator
using RheologyCalculator.RheologyModels
using GLMakie
using LaTeXStrings

function component_partition_figure()
    η = 1.0e19                                  # Pa s
    G = 1.0e10                                  # Pa
    dt = 1.0e10                                 # s
    C, ϕ = 1.0e6, 30.0                          # Pa, degrees

    c = SeriesModel(LinearViscosity(η), IncompressibleElasticity(G), DruckerPrager(C, ϕ, 0.0))
    others = (; dt, τ0 = (0.0,))

    ε̇ = 10.0 .^ range(-16, -11; length = 120)
    τ = similar(ε̇)
    fviscous, felastic, fplastic = similar(ε̇), similar(ε̇), similar(ε̇)

    for (i, e) in enumerate(ε̇)
        vars = (; ε = e)
        sol = solve(c, initial_guess_x(c, vars, (; τ = 1.0e6), others), vars, others)
        p = component_partition(c, sol, vars, others)

        τ[i] = only(p.viscous_τ)                # series: every component shares it
        v = only(p.viscous_ε)
        el = isempty(p.elastic_ε) ? 0.0 : only(p.elastic_ε)
        pl = isempty(p.plastic_ε) ? 0.0 : only(p.plastic_ε)
        total = v + el + pl
        fviscous[i], felastic[i], fplastic[i] = v / total, el / total, pl / total
    end

    fig = Figure(size = (900, 380))

    xlims = (first(ε̇), last(ε̇))            # no padding beyond the swept range

    ax1 = Axis(
        fig[1, 1];
        xscale = log10, yscale = log10,
        limits = (xlims, nothing),
        xlabel = L"Imposed $\dot{\varepsilon}_{II}$ [s$^{-1}$]",
        ylabel = L"$\tau_{II}$ [MPa]",
        title = "Stress saturates at yield",
    )
    lines!(ax1, ε̇, τ ./ 1.0e6; color = :black, linewidth = 2)
    hlines!(ax1, [C * cosd(ϕ) / 1.0e6]; color = :red, linestyle = :dash)
    # anchored just under the yield line: top-aligned so the text hangs below it
    text!(
        ax1, ε̇[8], C * cosd(ϕ) / 1.0e6 * 0.93; text = L"C\cos\phi",
        color = :red, align = (:left, :top)
    )

    ax2 = Axis(
        fig[1, 2];
        xscale = log10,
        xlabel = L"Imposed $\dot{\varepsilon}_{II}$ [s$^{-1}$]",
        ylabel = "Fraction of strain rate",
        title = "…and the partition shows where it goes",
        limits = (xlims, (0, 1)),
        # the bands are the content here; gridlines only show through them
        xgridvisible = false, ygridvisible = false,
    )
    # stacked bands: viscous at the bottom, then elastic, then plastic
    band!(ax2, ε̇, zeros(length(ε̇)), fviscous; color = (:steelblue, 0.85), label = "viscous")
    band!(ax2, ε̇, fviscous, fviscous .+ felastic; color = (:seagreen, 0.85), label = "elastic")
    band!(ax2, ε̇, fviscous .+ felastic, ones(length(ε̇)); color = (:indianred, 0.85), label = "plastic")
    axislegend(
        ax2;
        position = :lb,
        framevisible = true,
        framecolor = :black,
        backgroundcolor = :white,
        padding = (8, 8, 6, 6),
    )

    save(joinpath(@__DIR__, "..", "docs", "assets", "component_partition.png"), fig)
    return fig
end

component_partition_figure()
