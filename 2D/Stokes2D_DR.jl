# Initialisation
using Printf, Statistics, LinearAlgebra
# import Plots as pt
using GLMakie

using ParallelStencil
@init_parallel_stencil(Threads, Float64, 2)
# Macros
@views    av(A) = 0.25 * (A[1:(end - 1), 1:(end - 1)] .+ A[2:end, 1:(end - 1)] .+ A[1:(end - 1), 2:end] .+ A[2:end, 2:end])
@views av_xa(A) = 0.5 * (A[1:(end - 1), :] .+ A[2:end, :])
@views av_ya(A) = 0.5 * (A[:, 1:(end - 1)] .+ A[:, 2:end])
@views av4_harm(A) = 1.0 ./ (0.25 .* (1.0 ./ A[1:(end - 1), 1:(end - 1)] .+ 1.0 ./ A[2:end, 1:(end - 1)] .+ 1.0 ./ A[1:(end - 1), 2:end] .+ 1.0 ./ A[2:end, 2:end]))

# can be replaced by AD
function Gershgorin_Stokes2D_SchurComplement(ηc, ηv, γ, Δx, Δy, ncx, ncy)

    ηN = ones(ncx - 1, ncy)
    ηS = ones(ncx - 1, ncy)
    ηN[:, 1:(end - 1)] .= ηv[2:(end - 1), 2:(end - 1)]
    ηS[:, 2:(end - 0)] .= ηv[2:(end - 1), 2:(end - 1)]
    ηW = ηc[1:(end - 1), :]
    ηE = ηc[2:(end - 0), :]
    ebW = γ[1:(end - 1), :]
    ebE = γ[2:(end - 0), :]
    Cxx = ones(ncx - 1, ncy)
    Cxy = ones(ncx - 1, ncy)
    @. Cxx = abs.(ηN ./ Δy .^ 2) + abs.(ηS ./ Δy .^ 2) + abs.(ebE ./ Δx .^ 2 + (4 // 3) * ηE ./ Δx .^ 2) + abs.(ebW ./ Δx .^ 2 + (4 // 3) * ηW ./ Δx .^ 2) + abs.(-(-ηN ./ Δy - ηS ./ Δy) ./ Δy + (ebE ./ Δx + ebW ./ Δx) ./ Δx + ((4 // 3) * ηE ./ Δx + (4 // 3) * ηW ./ Δx) ./ Δx)
    @. Cxy = abs.(ebE ./ (Δx .* Δy) - 2 // 3 * ηE ./ (Δx .* Δy) + ηN ./ (Δx .* Δy)) + abs.(ebE ./ (Δx .* Δy) - 2 // 3 * ηE ./ (Δx .* Δy) + ηS ./ (Δx .* Δy)) + abs.(ebW ./ (Δx .* Δy) + ηN ./ (Δx .* Δy) - 2 // 3 * ηW ./ (Δx .* Δy)) + abs.(ebW ./ (Δx .* Δy) + ηS ./ (Δx .* Δy) - 2 // 3 * ηW ./ (Δx .* Δy))

    Dx = ones(ncx - 1, ncy)
    @. Dx .= -(-ηN ./ Δy - ηS ./ Δy) ./ Δy + (ebE ./ Δx + ebW ./ Δx) ./ Δx + ((4 // 3) * ηE ./ Δx + (4 // 3) * ηW ./ Δx) ./ Δx

    ηE = ones(ncx, ncy - 1)
    ηW = ones(ncx, ncy - 1)
    ηE[1:(end - 1), :] .= ηv[2:(end - 1), 2:(end - 1)]
    ηW[2:(end - 0), :] .= ηv[2:(end - 1), 2:(end - 1)]
    ηS = ηc[:, 1:(end - 1)]
    ηN = ηc[:, 2:(end - 0)]
    ebS = γ[:, 1:(end - 1)]
    ebN = γ[:, 2:(end - 0)]
    Cyy = ones(ncx, ncy - 1)
    Cyx = ones(ncx, ncy - 1)
    @. Cyy = abs.(ηE ./ Δx .^ 2) + abs.(ηW ./ Δx .^ 2) + abs.(ebN ./ Δy .^ 2 + (4 // 3) * ηN ./ Δy .^ 2) + abs.(ebS ./ Δy .^ 2 + (4 // 3) * ηS ./ Δy .^ 2) + abs.((ebN ./ Δy + ebS ./ Δy) ./ Δy + ((4 // 3) * ηN ./ Δy + (4 // 3) * ηS ./ Δy) ./ Δy - (-ηE ./ Δx - ηW ./ Δx) ./ Δx)
    @. Cyx = abs.(ebN ./ (Δx .* Δy) + ηE ./ (Δx .* Δy) - 2 // 3 * ηN ./ (Δx .* Δy)) + abs.(ebN ./ (Δx .* Δy) - 2 // 3 * ηN ./ (Δx .* Δy) + ηW ./ (Δx .* Δy)) + abs.(ebS ./ (Δx .* Δy) + ηE ./ (Δx .* Δy) - 2 // 3 * ηS ./ (Δx .* Δy)) + abs.(ebS ./ (Δx .* Δy) - 2 // 3 * ηS ./ (Δx .* Δy) + ηW ./ (Δx .* Δy))

    Dy = ones(ncx, ncy - 1)
    @. Dy .= (ebN ./ Δy + ebS ./ Δy) ./ Δy + ((4 // 3) * ηN ./ Δy + (4 // 3) * ηS ./ Δy) ./ Δy - (-ηE ./ Δx - ηW ./ Δx) ./ Δx

    λmaxVx = 1.0 ./ Dx .* (Cxx .+ Cxy)
    λmaxVy = 1.0 ./ Dy .* (Cyx .+ Cyy)

    return Dx, Dy, λmaxVx, λmaxVy
end

# 2D Stokes routine
@views function Stokes2D_VEP(n)
    sc = (σ = 1.0e0, t = 1.0e0, L = 1.0e0)
    # sc = (σ=1e6, t=1e10, L=1e3)
    # Physics
    Lx, Ly = 500.0e3 / sc.L, 500.0e3 / sc.L # domain size
    radi = 0.1e3 / sc.L             # inclusion radius
    η0 = 1.0e21 / sc.σ / sc.t         # viscous viscosity
    ηi = 1.0e23 / sc.σ / sc.t         # block viscosity
    δρ = 1000                   # density contrast
    ρ0 = 3.2e3 * 0              # viscous density
    ρi = 3.2e3 * 0 + δρ         # block density
    G = Inf / sc.σ
    C = Inf / sc.σ
    Δt = 10.0e3 * 365.25 * 24 * 3600 / sc.t # time step
    εbg = 0 * sc.t                         # background strain-rate
    comp = false
    K = Inf / sc.σ
    ϕ = 35.0
    ψ = 5.0
    ηvp = 2.0e20 / sc.σ / sc.t
    g = 9.81 * 1
    # Rectangular density anomaly
    xc_anomaly = 250.0e3         # origin of thermal anomaly
    yc_anomaly = -(Ly - 400.0e3) # origin of thermal anomaly
    r_anomaly = 50.0e3          # radius of perturbation
    # Numerics
    ncx, ncy = n * 32, n * 32 # numerical grid resolution
    nt = 30          # time steps
    ϵ = 1.0e-6         # tolerance
    iterMax = 20000        # max number of iters
    nout = 100          # check frequency
    c_fact = 0.5          # damping factor
    dτ_local = true         # helps a little bit sometimes, sometimes not!
    CFL_v = 0.99         # CFL: can't make it larger
    γfact = 20           # penalty: multiplier to the arithmetic mean of η
    rel_drop = 1.0e-2         # relative drop of velocity residual per PH iteration
    λ̇rel = 1.075        # overrelaxation helps!
    # Preprocessing
    Δx, Δy = Lx / ncx, Ly / ncy
    # Array initialisation
    Pt = zeros(ncx, ncy)
    Pt0 = zeros(ncx, ncy)
    Ptv = zeros(ncx + 1, ncy + 1)
    ΔPψ = zeros(ncx, ncy)
    ∇V = zeros(ncx, ncy)
    Vx = zeros(ncx + 1, ncy + 2)
    Vy = zeros(ncx + 2, ncy + 1)
    dVx = zeros(ncx - 1, ncy)
    dVy = zeros(ncx, ncy - 1)
    EIIc = zeros(ncx, ncy)
    EIIv = zeros(ncx + 1, ncy + 1)
    Exx = zeros(ncx, ncy)
    Eyy = zeros(ncx, ncy)
    Exy = zeros(ncx + 1, ncy + 1)
    Exxv = zeros(ncx + 1, ncy + 1)
    Eyyv = zeros(ncx + 1, ncy + 1)
    Exyc = zeros(ncx, ncy)
    TIIc = zeros(ncx, ncy)
    TIIv = zeros(ncx + 1, ncy + 1)
    Txx = zeros(ncx, ncy)
    Tyy = zeros(ncx, ncy)
    Txy = zeros(ncx + 1, ncy + 1)
    Txxv = zeros(ncx + 1, ncy + 1)
    Tyyv = zeros(ncx + 1, ncy + 1)
    Txyc = zeros(ncx, ncy)
    Txx0 = zeros(ncx, ncy)
    Tyy0 = zeros(ncx, ncy)
    Txy0 = zeros(ncx + 1, ncy + 1)
    Txxv0 = zeros(ncx + 1, ncy + 1)
    Tyyv0 = zeros(ncx + 1, ncy + 1)
    Txy0c = zeros(ncx, ncy)
    Fc = zeros(ncx, ncy)
    Fv = zeros(ncx + 1, ncy + 1)
    λ̇c = zeros(ncx, ncy)
    λ̇v = zeros(ncx + 1, ncy + 1)
    λ̇_true_c = zeros(ncx, ncy)
    λ̇_true_v = zeros(ncx + 1, ncy + 1)
    Rx = zeros(ncx - 1, ncy)
    Ry = zeros(ncx, ncy - 1)
    Rp = zeros(ncx, ncy)
    Rx0 = zeros(ncx - 1, ncy)
    Ry0 = zeros(ncx, ncy - 1)
    dVxdτ = zeros(ncx - 1, ncy)
    dVydτ = zeros(ncx, ncy - 1)
    βVx = zeros(ncx - 1, ncy)  # this disappears is dτ is not local
    βVy = zeros(ncx, ncy - 1)  # this disappears is dτ is not local
    cVx = zeros(ncx - 1, ncy)  # this disappears is dτ is not local
    cVy = zeros(ncx, ncy - 1)  # this disappears is dτ is not local
    αVx = zeros(ncx - 1, ncy)  # this disappears is dτ is not local
    αVy = zeros(ncx, ncy - 1)  # this disappears is dτ is not local
    ρg = zeros(ncx, ncy)
    ηb = zeros(ncx, ncy)
    ηc = zeros(ncx, ncy)
    ηv = zeros(ncx + 1, ncy + 1)
    ηve_c = zeros(ncx, ncy)
    ηve_v = zeros(ncx + 1, ncy + 1)
    ηvp_c = zeros(ncx, ncy)
    ηvp_v = zeros(ncx + 1, ncy + 1)
    ηvep_c = zeros(ncx, ncy)
    ηvep_v = zeros(ncx + 1, ncy + 1)
    ηc_sharp = zeros(ncx, ncy)
    ηv_sharp = zeros(ncx + 1, ncy + 1)
    P_num = zeros(ncx, ncy)
    # Initialisation
    xce, yce = LinRange(-Lx / 2 - Δx / 2, Lx / 2 + Δx / 2, ncx + 2), LinRange(-Ly - Δy / 2, 0.0e0, ncy + 2)
    xc, yc = LinRange(-Lx / 2 + Δx / 2, Lx / 2 - Δx / 2, ncx), LinRange(-Ly + Δy / 2, 0.0e0, ncy)
    xv, yv = LinRange(-Lx / 2, Lx / 2, ncx + 1), LinRange(-Ly, 0.0e0, ncy + 1)
    # Rectangular density anomaly
    xc_anomaly = 0.0e0             # origin of thermal anomaly
    yc_anomaly = -(Ly - 400.0e3) # origin of thermal anomaly
    r_anomaly = 50.0e3          # radius of perturbation
    isinblock_c = [(x^2 ≤ r_anomaly^2) & ((y - yc_anomaly)^2 ≤ r_anomaly^2) for x in xc, y in yc]
    isinblock_v = [(x^2 ≤ r_anomaly^2) & ((y - yc_anomaly)^2 ≤ r_anomaly^2) for x in xv, y in yv]
    # Block buoyancy
    ρg .= ρ0 .* g
    ρg[isinblock_c] .= ρi .* g
    ρg_vy = av_ya(ρg)
    Pt .= reverse(cumsum(reverse((ρg) .* Δy, dims = 2), dims = 2), dims = 2)
    Pt0 .= Pt
    # Block viscosity
    # ηv_sharp              .= η0
    # ηc_sharp              .= η0
    # ηv_sharp[isinblock_v] .= ηi
    # ηc_sharp[isinblock_c] .= ηi
    ηv .= η0
    ηc .= η0
    ηv[isinblock_v] .= ηi
    ηc[isinblock_c] .= ηi
    # Harmonic averaging mimicking PIC interpolation
    # ηc                  .= av4_harm(ηv_sharp)
    # ηv[2:end-1,2:end-1] .= av4_harm(ηc_sharp)
    ηv[1, :] .= ηv[2, :]; ηv[end, :] .= ηv[end - 1, :]
    ηv[:, 1] .= ηv[:, 2]; ηv[:, end] .= ηv[:, end - 1]
    # Effective viscosity
    ηve_c .= (1 ./ ηc .+ 1 ./ (G * Δt)) .^ -1
    ηve_v .= (1 ./ ηv .+ 1 ./ (G * Δt)) .^ -1
    # Bulk viscosity
    ηb .= K .* Δt
    # Select γ
    γi = γfact * mean(ηc) .* ones(size(ηc))
    # (Pseudo-)compressibility
    γ_eff = zeros(size(ηb))
    if comp
        γ_num = γi .* ones(size(ηb)) * 1.0
        γ_phy = ηb
        γ_eff = ((γ_phy .* γ_num) ./ (γ_phy .+ γ_num))
    else
        γ_eff .= γi
        γ_eff .= γ_eff
        ηb .= γ_eff
    end
    # Optimal pseudo-time steps - can be replaced by AD
    Dx, Dy, λmaxVx, λmaxVy = Gershgorin_Stokes2D_SchurComplement(ηve_c, ηve_v, γ_eff, Δx, Δy, ncx, ncy)
    # Select dτ
    if dτ_local
        dτVx = 2.0 ./ sqrt.(λmaxVx) * CFL_v
        dτVy = 2.0 ./ sqrt.(λmaxVy) * CFL_v
    else
        dτVx = 2.0 ./ sqrt.(maximum(λmaxVx)) * CFL_v
        dτVy = 2.0 ./ sqrt.(maximum(λmaxVy)) * CFL_v
    end
    βVx .= 2 .* dτVx ./ (2 .+ cVx .* dτVx)
    βVy .= 2 .* dτVy ./ (2 .+ cVy .* dτVy)
    αVx .= (2 .- cVx .* dτVx) ./ (2 .+ cVx .* dτVx)
    αVy .= (2 .- cVy .* dτVy) ./ (2 .+ cVy .* dτVy)
    # Initial condition
    # Vx     .=   εbg.*xv .+    0*yce'
    # Vy     .=     0*xce .- εbg.*yv'
    # Vx[2:end-1,:] .= 0 # ensure non zero initial pressure residual
    # Vy[:,2:end-1] .= 0 # ensure non zero initial pressure residual
    # Time
    Tii_evo = zeros(nt)
    it_evo = zeros(nt)
    itg = 0
    err_evo_it, err_evo_V, err_evo_P = zeros(iterMax), zeros(iterMax), zeros(iterMax)
    for it in 1:1
        Txx0 .= Txx; Tyy0 .= Tyy; Txy0 .= Txy; Txy0c .= Txyc;  Txxv0 .= Txxv; Tyyv0 .= Tyyv; Pt0 .= Pt
        # Iteration loop
        errVx0 = 1.0;  errVy0 = 1.0;  errPt0 = 1.0
        errVx00 = 1.0;  errVy00 = 1.0
        iter = 0; err = 2 * ϵ; err_evo_it .= 0.0; err_evo_V .= 0.0; err_evo_P .= 0.0
        @time for itPH in 1:1000
            # Boundaries
            Vx[:, 1] .= Vx[:, 2]; Vx[:, end] .= Vx[:, end - 1]
            Vy[1, :] .= Vy[2, :]; Vy[end, :] .= Vy[end - 1, :]
            # Divergence
            ∇V .= (Vx[2:end, 2:(end - 1)] .- Vx[1:(end - 1), 2:(end - 1)]) ./ Δx .+ (Vy[2:(end - 1), 2:end] .- Vy[2:(end - 1), 1:(end - 1)]) ./ Δy
            # Pressure on vertices
            Ptv[2:(end - 1), 2:(end - 1)] .= av(Pt)
            Ptv[1, :] .= Ptv[2, :]; Ptv[end, :] .= Ptv[end - 1, :]
            Ptv[:, 1] .= Ptv[:, 2]; Ptv[:, end] .= Ptv[:, end - 1]
            # Deviatoric strain rate
            Exx .= (Vx[2:end, 2:(end - 1)] .- Vx[1:(end - 1), 2:(end - 1)]) ./ Δx .- 1.0 / 3.0 .* ∇V
            Eyy .= (Vy[2:(end - 1), 2:end] .- Vy[2:(end - 1), 1:(end - 1)]) ./ Δy .- 1.0 / 3.0 .* ∇V
            Exy .= 0.5 .* ((Vx[:, 2:end] .- Vx[:, 1:(end - 1)]) ./ Δy .+ (Vy[2:end, :] .- Vy[1:(end - 1), :]) ./ Δx)
            Exxv[2:(end - 1), 2:(end - 1)] .= av(Exx)
            Eyyv[2:(end - 1), 2:(end - 1)] .= av(Eyy)
            Exyc .= av(Exy)
            EIIc .= sqrt.(0.5 .* ((Exx .+ Txx0 ./ (2 * G * Δt)) .^ 2 .+ (Eyy .+ Tyy0 ./ (2 * G * Δt)) .^ 2 .+ (.-(Exx .+ Txx0 ./ (2 * G * Δt)) .- (Eyy .+ Tyy0 ./ (2 * G * Δt))) .^ 2) .+ (Exyc .+ Txy0c ./ (2 * G * Δt)) .^ 2)
            EIIv .= sqrt.(0.5 .* ((Exxv .+ Txxv0 ./ (2 * G * Δt)) .^ 2 .+ (Eyyv .+ Tyyv0 ./ (2 * G * Δt)) .^ 2 .+ (.-(Exxv .+ Txxv0 ./ (2 * G * Δt)) .- (Eyyv .+ Tyyv0 ./ (2 * G * Δt))) .^ 2) .+ (Exy .+ Txy0 ./ (2 * G * Δt)) .^ 2)
            # Deviatoric stress
            Txx .= 2.0 .* ηve_c .* (Exx .+ Txx0 ./ (2 * G * Δt))
            Tyy .= 2.0 .* ηve_c .* (Eyy .+ Tyy0 ./ (2 * G * Δt))
            Txy .= 2.0 .* ηve_v .* (Exy .+ Txy0 ./ (2 * G * Δt))
            Txxv .= 2.0 .* ηve_v .* (Exxv .+ Txxv0 ./ (2 * G * Δt))
            Tyyv .= 2.0 .* ηve_v .* (Eyyv .+ Tyyv0 ./ (2 * G * Δt))
            Txyc .= 2.0 .* ηve_c .* (Exyc .+ Txy0c ./ (2 * G * Δt))
            TIIc .= sqrt.(0.5 .* (Txx .^ 2 .+ Tyy .^ 2 .+ (.-Txx .- Tyy) .^ 2) .+ Txyc .^ 2)
            TIIv .= sqrt.(0.5 .* (Txxv .^ 2 .+ Tyyv .^ 2 .+ (.-Txxv .- Tyyv) .^ 2) .+ Txy .^ 2)
            # # Plasticity
            # λ̇c            .= 0.
            # λ̇v            .= 0.
            # Fc            .= TIIc .- C.*cosd(ϕ) .- Pt .*sind(ϕ)
            # Fv            .= TIIv .- C.*cosd(ϕ) .- Ptv.*sind(ϕ)
            # λ̇c[Fc.>0]     .= Fc[Fc.>0]./(ηve_c[Fc.>0] .+ ηvp .+ K.*Δt.*sind(ϕ).*sind.(ψ))
            # λ̇v[Fv.>0]     .= Fv[Fv.>0]./(ηve_v[Fv.>0] .+ ηvp .+ K.*Δt.*sind(ϕ).*sind.(ψ))
            # ηvep_c        .= ηve_c
            # ηvep_v        .= ηve_v
            # ηvp_c .= (TIIc.-λ̇c.*ηve_c) ./ (2 .* EIIc)
            # ηvp_v .= (TIIv.-λ̇v.*ηve_v) ./ (2 .* EIIv)
            # ηvep_c[Fc.>0] .= ηvp_c[Fc.>0]
            # ηvep_v[Fv.>0] .= ηvp_v[Fv.>0]
            # Txx   .= 2.0.*ηvep_c.*(Exx .+ Txx0./(2*G*Δt))
            # Tyy   .= 2.0.*ηvep_c.*(Eyy .+ Tyy0./(2*G*Δt))
            # Txy   .= 2.0.*ηvep_v.*(Exy .+ Txy0./(2*G*Δt))
            # Txxv  .= 2.0.*ηvep_v.*(Exxv .+ Txxv0./(2*G*Δt))
            # Tyyv  .= 2.0.*ηvep_v.*(Eyyv .+ Tyyv0./(2*G*Δt))
            # Txyc  .= 2.0.*ηvep_c.*(Exyc .+ Txy0c./(2*G*Δt))
            # ΔPψ   .= λ̇c.*sind(ψ).*K.*Δt
            # Check
            TIIc .= sqrt.(0.5 .* (Txx .^ 2 .+ Tyy .^ 2 .+ (.-Txx .- Tyy) .^ 2) .+ Txyc .^ 2)
            TIIv .= sqrt.(0.5 .* (Txxv .^ 2 .+ Tyyv .^ 2 .+ (.-Txxv .- Tyyv) .^ 2) .+ Txy .^ 2)
            Fc .= TIIc .- C .* cosd(ϕ) .- (Pt .+ λ̇c .* sind(ψ) .* K .* Δt) .* sind(ϕ) .- ηvp .* λ̇c
            Fv .= TIIv .- C .* cosd(ϕ) .- (Ptv .+ λ̇v .* sind(ψ) .* K .* Δt) .* sind(ϕ) .- ηvp .* λ̇v
            # Residuals
            Rx .= .-(Pt[2:end, :] .- Pt[1:(end - 1), :]) ./ Δx .- (ΔPψ[2:end, :] .- ΔPψ[1:(end - 1), :]) ./ Δx .+ (Txx[2:end, :] .- Txx[1:(end - 1), :]) ./ Δx .+ (Txy[2:(end - 1), 2:end] .- Txy[2:(end - 1), 1:(end - 1)]) ./ Δy
            Ry .= .-(Pt[:, 2:end] .- Pt[:, 1:(end - 1)]) ./ Δy .- (ΔPψ[:, 2:end] .- ΔPψ[:, 1:(end - 1)]) ./ Δy .+ (Tyy[:, 2:end] .- Tyy[:, 1:(end - 1)]) ./ Δy .+ (Txy[2:end, 2:(end - 1)] .- Txy[1:(end - 1), 2:(end - 1)]) ./ Δx .- ρg_vy
            Rp .= .-∇V .- comp * (Pt .- Pt0) ./ ηb
            # Residual check
            errVy = norm(Ry) / length(Ry)
            errVx = norm(Rx) / length(Rx)
            errPt = norm(Rp) / length(Rp)
            if itPH == 1
                errVx0 = errVx; errVy0 = errVy; errPt0 = errPt
            end
            err = maximum([min(errVx / errVx0, errVx), min(errVy / errVy0, errVy)]) #, min(errPt/errPt0, errPt)
            @printf("itPH = %02d iter = %06d iter/nx = %03d, err = %1.3e norm[Rx=%1.3e, Ry=%1.3e, Rp=%1.3e] \n", itPH, iter, iter / ncx, err, min(errVx / errVx0, errVx), min(errVy / errVy0, errVy), min(errPt / errPt0, errPt))
            if (err < ϵ)
                break
            end
            # Set tolerance of velocity solve proportional to residual
            ϵ_vel = err * rel_drop
            itPT = 0.0
            while (err > ϵ_vel && itPT <= iterMax)
                iter += 1
                itPT += 1
                itg += 1
                # Pseudo-old dudes
                Rx0 .= Rx
                Ry0 .= Ry
                # Boundaries
                Vx[:, 1] .= Vx[:, 2]; Vx[:, end] .= Vx[:, end - 1]
                Vy[1, :] .= Vy[2, :]; Vy[end, :] .= Vy[end - 1, :]
                # Divergence
                ∇V .= (Vx[2:end, 2:(end - 1)] .- Vx[1:(end - 1), 2:(end - 1)]) ./ Δx .+ (Vy[2:(end - 1), 2:end] .- Vy[2:(end - 1), 1:(end - 1)]) ./ Δy
                Rp .= .-∇V .- comp * (Pt .- Pt0) ./ ηb
                # Deviatoric strain rate
                Exx .= (Vx[2:end, 2:(end - 1)] .- Vx[1:(end - 1), 2:(end - 1)]) ./ Δx .- 1.0 / 3.0 .* ∇V
                Eyy .= (Vy[2:(end - 1), 2:end] .- Vy[2:(end - 1), 1:(end - 1)]) ./ Δy .- 1.0 / 3.0 .* ∇V
                Exy .= 0.5 .* ((Vx[:, 2:end] .- Vx[:, 1:(end - 1)]) ./ Δy .+ (Vy[2:end, :] .- Vy[1:(end - 1), :]) ./ Δx)
                Exxv[2:(end - 1), 2:(end - 1)] .= av(Exx)
                Eyyv[2:(end - 1), 2:(end - 1)] .= av(Eyy)
                Exyc .= av(Exy)
                EIIc .= sqrt.(0.5 .* ((Exx .+ Txx0 ./ (2 * G * Δt)) .^ 2 .+ (Eyy .+ Tyy0 ./ (2 * G * Δt)) .^ 2 .+ (.-(Exx .+ Txx0 ./ (2 * G * Δt)) .- (Eyy .+ Tyy0 ./ (2 * G * Δt))) .^ 2) .+ (Exyc .+ Txy0c ./ (2 * G * Δt)) .^ 2)
                EIIv .= sqrt.(0.5 .* ((Exxv .+ Txxv0 ./ (2 * G * Δt)) .^ 2 .+ (Eyyv .+ Tyyv0 ./ (2 * G * Δt)) .^ 2 .+ (.-(Exxv .+ Txxv0 ./ (2 * G * Δt)) .- (Eyyv .+ Tyyv0 ./ (2 * G * Δt))) .^ 2) .+ (Exy .+ Txy0 ./ (2 * G * Δt)) .^ 2)
                # Deviatoric stress
                Txx .= 2.0 .* ηve_c .* (Exx .+ Txx0 ./ (2 * G * Δt))
                Tyy .= 2.0 .* ηve_c .* (Eyy .+ Tyy0 ./ (2 * G * Δt))
                Txy .= 2.0 .* ηve_v .* (Exy .+ Txy0 ./ (2 * G * Δt))
                Txxv .= 2.0 .* ηve_v .* (Exxv .+ Txxv0 ./ (2 * G * Δt))
                Tyyv .= 2.0 .* ηve_v .* (Eyyv .+ Tyyv0 ./ (2 * G * Δt))
                Txyc .= 2.0 .* ηve_c .* (Exyc .+ Txy0c ./ (2 * G * Δt))
                TIIc .= sqrt.(0.5 .* (Txx .^ 2 .+ Tyy .^ 2 .+ (.-Txx .- Tyy) .^ 2) .+ Txyc .^ 2)
                TIIv .= sqrt.(0.5 .* (Txxv .^ 2 .+ Tyyv .^ 2 .+ (.-Txxv .- Tyyv) .^ 2) .+ Txy .^ 2)
                # # Plasticity
                # Fc              .= TIIc .- C.*cosd(ϕ) .- Pt .*sind(ϕ)
                # Fv              .= TIIv .- C.*cosd(ϕ) .- Ptv.*sind(ϕ)
                # λ̇_true_c        .= 0.
                # λ̇_true_v        .= 0.
                # λ̇_true_c[Fc.>0] .= Fc[Fc.>0]./(ηve_c[Fc.>0] .+ ηvp .+ K.*Δt.*sind(ϕ).*sind.(ψ))
                # λ̇_true_v[Fv.>0] .= Fv[Fv.>0]./(ηve_v[Fv.>0] .+ ηvp .+ K.*Δt.*sind(ϕ).*sind.(ψ))
                # λ̇c              .= λ̇rel*λ̇_true_c .+ (1-λ̇rel).*λ̇c
                # λ̇v              .= λ̇rel*λ̇_true_v .+ (1-λ̇rel).*λ̇v
                # ηvep_c .= ηve_c
                # ηvep_v .= ηve_v
                # ηvp_c  .= (TIIc.-λ̇c.*ηve_c) ./ (2 .* EIIc)
                # ηvp_v  .= (TIIv.-λ̇v.*ηve_v) ./ (2 .* EIIv)
                # ηvep_c[Fc.>0] .= ηvp_c[Fc.>0]
                # ηvep_v[Fv.>0] .= ηvp_v[Fv.>0]
                # Txx    .= 2.0.*ηvep_c.*(Exx  .+ Txx0 ./(2*G*Δt))
                # Tyy    .= 2.0.*ηvep_c.*(Eyy  .+ Tyy0 ./(2*G*Δt))
                # Txy    .= 2.0.*ηvep_v.*(Exy  .+ Txy0 ./(2*G*Δt))
                # ΔPψ    .= λ̇c.*sind(ψ).*K.*Δt
                # Residuals
                P_num .= γ_eff .* Rp
                Rx .= (1.0 ./ Dx) .* (.-(P_num[2:end, :] .- P_num[1:(end - 1), :]) ./ Δx .- (Pt[2:end, :] .- Pt[1:(end - 1), :]) ./ Δx .- (ΔPψ[2:end, :] .- ΔPψ[1:(end - 1), :]) ./ Δx .+ (Txx[2:end, :] .- Txx[1:(end - 1), :]) ./ Δx .+ (Txy[2:(end - 1), 2:end] .- Txy[2:(end - 1), 1:(end - 1)]) ./ Δy)
                Ry .= (1.0 ./ Dy) .* (.-(P_num[:, 2:end] .- P_num[:, 1:(end - 1)]) ./ Δy .- (Pt[:, 2:end] .- Pt[:, 1:(end - 1)]) ./ Δy .- (ΔPψ[:, 2:end] .- ΔPψ[:, 1:(end - 1)]) ./ Δy .+ (Tyy[:, 2:end] .- Tyy[:, 1:(end - 1)]) ./ Δy .+ (Txy[2:end, 2:(end - 1)] .- Txy[1:(end - 1), 2:(end - 1)]) ./ Δx .- ρg_vy)
                # Damping-pong
                dVxdτ .= αVx .* dVxdτ .+ Rx
                dVydτ .= αVy .* dVydτ .+ Ry
                # PT updates
                Vx[2:(end - 1), 2:(end - 1)] .+= dVxdτ .* βVx .* dτVx
                Vy[2:(end - 1), 2:(end - 1)] .+= dVydτ .* βVy .* dτVy
                # Residual check
                if mod(iter, nout) == 0
                    errVx = norm(Dx .* Rx) / length(Ry)
                    errVy = norm(Dy .* Ry) / length(Rx)
                    if iter == nout
                        errVx00 = errVx; errVy00 = errVy
                    end
                    err = maximum([errVx ./ errVx00, errVy ./ errVy00])
                    # err_evo_V[iter] = errVx/errVx00; err_evo_P[iter] = errPt/errPt0; err_evo_it[iter] =  iter
                    dVx .= dVxdτ .* βVx .* dτVx
                    dVy .= dVydτ .* βVy .* dτVy
                    # @printf("it = %d, iter = %d, err = %1.3e norm[Rx=%1.3e, Ry=%1.3e] \n", it, iter, err, norm_Rx, norm_Ry)
                    # λminV  = abs.((sum(dVx.*(Rx .- Rx0))) + abs.((sum(dVy.*(Ry .- Ry0))) )/ ( sum(dVx.*dVx)) + sum(dVy.*dVy) )
                    λminV = abs(sum(dVx .* (Rx .- Rx0)) + sum(dVy .* (Ry .- Ry0))) / (sum(dVx .* dVx) .+ sum(dVy .* dVy))
                    cVx .= 2 * sqrt.(λminV) * c_fact
                    cVy .= 2 * sqrt.(λminV) * c_fact
                    # Optimal pseudo-time steps - can be replaced by AD
                    Dx, Dy, λmaxVx, λmaxVy = Gershgorin_Stokes2D_SchurComplement(ηve_c, ηve_v, γ_eff, Δx, Δy, ncx, ncy)
                    # Select dτ
                    if dτ_local
                        dτVx = 2.0 ./ sqrt.(λmaxVx) * CFL_v
                        dτVy = 2.0 ./ sqrt.(λmaxVy) * CFL_v
                    else
                        dτVx = 2.0 ./ sqrt.(maximum(λmaxVx)) * CFL_v
                        dτVy = 2.0 ./ sqrt.(maximum(λmaxVy)) * CFL_v
                    end
                    βVx .= 2 .* dτVx ./ (2 .+ cVx .* dτVx)
                    βVy .= 2 .* dτVy ./ (2 .+ cVy .* dτVy)
                    αVx .= (2 .- cVx .* dτVx) ./ (2 .+ cVx .* dτVx)
                    αVy .= (2 .- cVy .* dτVy) ./ (2 .+ cVy .* dτVy)
                end
            end
            Pt .+= γ_eff .* Rp
        end

        Tii_evo[it] = maximum(TIIc)
        it_evo[it] = iter / ncx

        # Plotting
        EIIc .= sqrt.(0.5 .* ((Exx) .^ 2 .+ (Eyy) .^ 2 .+ (.-(Exx) .- (Eyy)) .^ 2) .+ (Exyc) .^ 2)
        p1 = pt.heatmap(Vy', aspect_ratio = 1)
        p2 = pt.heatmap(∇V', aspect_ratio = 1)
        display(pt.plot(p1, p2))
        @show iter / ncx
        @show itg

    end
    n = length(ηc)
    @show η_h = 1.0 / sum(1.0 / n ./ ηc)
    @show η_g = exp(sum(1.0 / n * log.(ηc)))
    @show η_a = mean(ηc)
    @show extrema(∇V)

    return Vx, Vy, ∇V, Pt
end

n = 4
Vx, Vy, ∇V, Pt = Stokes2D_VEP(n)

f1, ax, h1 = heatmap(Vx)
f2, ax, h2 = heatmap(Vy)
f3, ax, h3 = heatmap(log10.(abs.(∇V)))

Colorbar(f1[1, 2], h1); f1
Colorbar(f2[1, 2], h2); f2
Colorbar(f3[1, 2], h3); f3
