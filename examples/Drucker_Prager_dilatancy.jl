# Dilatancy in a closed system.
#
# A visco-elasto-plastic column is sheared at a constant deviatoric strain rate
# while its total volumetric strain rate is held at zero. A Drucker-Prager
# element with a dilatancy angle ψ > 0 expands as it yields, and with the volume
# fixed that expansion has to be taken up elastically. The pressure therefore
# climbs, and because the Drucker-Prager yield stress grows with pressure the
# material keeps hardening instead of settling on a plateau.
#
# The run checks the three identities that make the result meaningful:
#
#   flow rule       θ_pl = λ sinψ
#   mass balance    θ_el + θ_pl = θ  (the prescribed total, here zero)
#   yield           τ = C cosϕ + P sinϕ
#
# With ψ = 0 the pressure never moves and the stress settles on the classical
# C cosϕ plateau, which is the reference curve the dilatant runs depart from.

using RheologyCalculator
using RheologyCalculator.RheologyModels
using RheologyCalculator.RheologyModels: second_invariant_2D, tensor_strain_rate_2D,
    zero_stress_tensor_2D, elastic_stress_history_2D, compute_F, compute_Q

using GLMakie
using Printf

const C = 10.0e6      # cohesion [Pa]
const ϕ = 30.0        # friction angle [degrees]
const G = 10.0e9      # shear modulus [Pa]
const K = 20.0e9      # bulk modulus [Pa]
const η = 1.0e22      # shear viscosity [Pa s]
const ε̇ = 1.0e-14     # applied deviatoric strain rate [1/s]

dilatant_column(ψ) = SeriesModel(LinearViscosity(η), Elasticity(G, K), DruckerPrager(C, ϕ, ψ))

"""
    dilatancy_history(ψ; ntime, dt)

Shear a dilatant column at fixed deviatoric strain rate and zero total
volumetric strain rate. Returns the time series of stress, pressure and
volumetric plastic strain rate, together with the largest violation of each of
the three identities above.
"""
function dilatancy_history(ψ; ntime = 1500, dt = 1.0e8)
    c = dilatant_column(ψ)
    vars = (; ε = tensor_strain_rate_2D(ε̇), θ = 0.0)
    others = (; dt = dt, τ0 = (zero_stress_tensor_2D(),), P0 = (0.0,))

    x = initial_guess_x(c, vars, (; τ = 0.0, P = 0.0, λ = 0.0), others)
    xnorm = normalisation_x(c, C, second_invariant_2D(vars.ε))

    t_v, τ_v, P_v, θ_v = ntuple(_ -> zeros(ntime), 4)
    τ_e, P_e = (zero_stress_tensor_2D(),), (0.0,)
    t = 0.0
    err_flow = err_mass = err_yield = 0.0

    for i in 2:ntime
        P_prev = P_e[1]
        others = (; dt = dt, τ0 = τ_e, P0 = P_e)
        x = solve(c, x, vars, others; xnorm0 = xnorm)
        s = NamedTuple{x_keys(c)}(Tuple(x))

        θ_pl = volumetric_plastic_strain_rate(c, x, others)
        # The elastic volumetric strain rate carries the opposite sign to the
        # pressure increment, so a dilating element drives the pressure up.
        θ_el = -(s.P - P_prev) / (K * dt)

        err_flow = max(err_flow, abs(θ_pl - s.λ * sind(ψ)))
        err_mass = max(err_mass, abs(θ_el + θ_pl - vars.θ))
        # Only meaningful while the element is yielding, i.e. λ > 0.
        s.λ > 0 && (err_yield = max(err_yield, abs(s.τ - C * cosd(ϕ) - s.P * sind(ϕ))))

        t += dt
        t_v[i], τ_v[i], P_v[i], θ_v[i] = t, s.τ, s.P, θ_pl
        τ_e = elastic_stress_history_2D(c, s.τ, vars.ε, τ_e, others)
        P_e = compute_pressure_elastic(c, x, others)
    end

    return (; t = t_v, τ = τ_v, P = P_v, θ_pl = θ_v, err_flow, err_mass, err_yield)
end

ψ_all = (0.0, 10.0, 20.0, 30.0)
runs = map(dilatancy_history, ψ_all)

println("Drucker-Prager dilatancy, series column")
println("  ψ [°]   τ_final [MPa]   P_final [MPa]   θ_pl [1/s]     flow rule    mass balance   yield")
for (ψ, r) in zip(ψ_all, runs)
    @printf(
        "  %5.1f   %13.4f   %13.4f   %10.3e   %10.2e   %12.2e   %8.2e\n",
        ψ, r.τ[end] / 1.0e6, r.P[end] / 1.0e6, r.θ_pl[end], r.err_flow, r.err_mass, r.err_yield
    )
end

let
    # ψ = 0 must reproduce the pressure-independent plateau exactly, and every
    # dilatant run must both harden and stay on the yield surface.
    @assert runs[1].P[end] == 0
    @assert runs[1].τ[end] ≈ C * cosd(ϕ)
    for (ψ, r) in zip(ψ_all, runs)
        @assert r.err_flow < 1.0e-24
        @assert r.err_mass < 1.0e-24
        @assert r.err_yield < 1.0e-6
        ψ > 0 && @assert r.P[end] > 0 && r.τ[end] > runs[1].τ[end]
    end
    # Steeper dilatancy hardens harder.
    @assert issorted(map(r -> r.τ[end], runs))
end

let
    # The same flow rule has to hold when the plastic element sits in a parallel
    # branch, where the branch carries its own pressure and volumetric strain rate.
    ψ = 30.0
    c = SeriesModel(
        LinearViscosity(η), Elasticity(G, K),
        ParallelModel(LinearViscosity(1.0e21), DruckerPrager(C, ϕ, ψ))
    )
    vars = (; ε = tensor_strain_rate_2D(ε̇), θ = 0.0)
    others = (; dt = 1.0e11, τ0 = (zero_stress_tensor_2D(),), P0 = (0.0,))
    x = initial_guess_x(c, vars, (; τ = 0.0, P = 0.0, λ = 0.0), others)
    x = solve(c, x, vars, others; xnorm0 = normalisation_x(c, C, second_invariant_2D(vars.ε)))
    s = NamedTuple{x_keys(c)}(Tuple(x))

    θ_pl = volumetric_plastic_strain_rate(c, x, others)
    println("\nparallel branch (ψ = $(ψ)°): ", x_keys(c))
    @printf(
        "  λ = %.6e   θ_pl = %.6e   θ_branch = %.6e   P = %.4f MPa   P_pl = %.4f MPa\n",
        s.λ, θ_pl, s.θ, s.P / 1.0e6, s.P_pl / 1.0e6
    )

    @assert θ_pl ≈ s.λ * sind(ψ)
    # Every element of a parallel branch shares its volumetric strain rate, and
    # the branch viscosity carries no volumetric stress, so it also shares P.
    @assert s.θ ≈ θ_pl
    @assert s.P_pl ≈ s.P
    @assert s.λ > 0
end

let
    # The stress path is the clearest view of what dilatancy does. Loading is
    # elastic and volume-preserving until the path reaches the yield envelope
    # F = 0; from there it cannot leave the envelope, so as the plastic dilation
    # raises the pressure the path travels along it and the stress keeps rising.
    # Flow is non-associated: the plastic potential Q has slope sinψ against the
    # envelope's sinϕ, and the gap between the two is the dilation.
    P_grid = range(-2.0e6, 22.0e6, length = 400)
    τ_grid = range(0.0, 22.0e6, length = 400)
    elements = map(ψ -> last(dilatant_column(ψ).leafs), ψ_all)
    F_grid = [compute_F(first(elements), τ, P) for P in P_grid, τ in τ_grid]
    # The potential contour through the point where the elastic path first meets
    # the envelope, one per run. Its slope is sinψ against the envelope's sinϕ,
    # and that difference is the dilation; they coincide when ψ = ϕ.
    Q_grids = map(dp -> [compute_Q(dp, τ, P) for P in P_grid, τ in τ_grid], elements)
    Q_level = C * cosd(ϕ)

    function figure()
        SecYear = 3600 * 24 * 365.25
        fig = Figure(fontsize = 20, size = (1000, 900))
        ax1 = Axis(fig[1, 1], title = "Deviatoric stress", xlabel = L"$t$ [kyr]", ylabel = L"$\tau_{II}$ [MPa]")
        ax2 = Axis(fig[1, 2], title = "Pressure", xlabel = L"$t$ [kyr]", ylabel = L"$P$ [MPa]")
        ax3 = Axis(fig[2, 1], title = "Volumetric plastic strain rate", xlabel = L"$t$ [kyr]", ylabel = L"$\dot{\theta}_{pl}$ [1/s]")
        ax4 = Axis(fig[2, 2], title = "Stress paths", xlabel = L"$P$ [MPa]", ylabel = L"$\tau_{II}$ [MPa]")
        colors = (:black, :steelblue, :darkorange, :crimson)

        # All four elements share C and ϕ, so one envelope serves every run. The
        # ψ = 30° potential lies on top of it: that run is associated flow.
        contour!(ax4, P_grid / 1.0e6, τ_grid / 1.0e6, F_grid, levels = [0.01], color = :black, linewidth = 4)
        for (Q_grid, col) in zip(Q_grids, colors)
            contour!(ax4, P_grid / 1.0e6, τ_grid / 1.0e6, Q_grid, levels = [Q_level], color = col, linestyle = :dash)
        end

        for (ψ, r, col) in zip(ψ_all, runs, colors)
            t_kyr = r.t / SecYear / 1.0e3
            label = L"\psi = %$(Int(ψ))°"
            lines!(ax1, t_kyr, r.τ / 1.0e6, color = col, linewidth = 3, label = label)
            lines!(ax2, t_kyr, r.P / 1.0e6, color = col, linewidth = 3)
            lines!(ax3, t_kyr, r.θ_pl, color = col, linewidth = 3)
        end

        # Every path rides the same envelope, so draw the longest first and
        # shrink the marker as they get shorter to keep all four visible.
        for (i, (r, col)) in enumerate(reverse(collect(zip(runs, colors))))
            scatter!(ax4, r.P[2:30:end] / 1.0e6, r.τ[2:30:end] / 1.0e6, color = col, markersize = 17 - 3i)
        end
        hlines!(ax1, C * cosd(ϕ) / 1.0e6, color = :grey, linestyle = :dash)

        axislegend(ax1, position = :lt)
        limits!(ax4, -2, 22, 0, 22)
        linkxaxes!(ax1, ax2, ax3)

        save(joinpath(@__DIR__, "..", "docs", "assets", "DruckerPrager_dilatancy.png"), fig)
        return display(fig)
    end
    with_theme(figure, theme_latexfonts())
end
