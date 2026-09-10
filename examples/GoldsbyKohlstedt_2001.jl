using RheologyCalculator
using RheologyCalculator.RheologyModels
using GLMakie

const R = 8.314

diffusion = GoldsbyKohlstedtDiffusion(R, 9.1e-4, 59.4e3, 1.0e-4, 59.4e3, 1.97e-5, 1.0e-9)
dislocation = DislocationCreep(4, 0, 4.0e4 * 10.0^(-24), 64.0e3, 0.0, R)
gbs = GoldsbyKohlstedtCreep(1.7, 1.4, 3.9e-3 * 10.0^(-10.2), 49.0e3, R)
basal = GoldsbyKohlstedtCreep(2.4, 0, 5.5e7 * 10.0^(-14.4), 60.0e3, R)

rates = 10.0 .^ range(-16, -4, length = 180)

function flow_curve(ice, rates, T, d)
    others = (; T, d = (d, d, d, d), f = 1.0)
    x = initial_guess_x(ice, (; ε = first(rates)), (; τ = 1.0e5, ε = 1.0e-12), others)
    stress = similar(rates)
    for (i, ε) in pairs(rates)
        sol = solve(ice, x, (; ε), others; verbose = false)
        stress[i] = sol[1]
        x = sol
    end
    stress
end

ice = SeriesModel(ParallelModel(diffusion, dislocation, SeriesModel(gbs, basal)))
temperatures = (230.0, 240.0, 250.0, 260.0)
grain_sizes = (1.0e-5, 1.0e-4, 1.0e-3)
fig = Figure(size = (1200, 900))
for (j, d) in pairs(grain_sizes), (i, T) in pairs(temperatures)
    ax = Axis(fig[i, j], xscale = log10, yscale = log10,
        xlabel = "strain rate [s⁻¹]", ylabel = "stress [Pa]",
        title = "T = $(T) K, d = $(d) m")
    lines!(ax, rates, flow_curve(ice, rates, T, d), linewidth = 2)
end
Label(fig[0, :], "Goldsby–Kohlstedt (2001) Figure 7 attempt — native RC solve", fontsize = 24)
display(fig)
