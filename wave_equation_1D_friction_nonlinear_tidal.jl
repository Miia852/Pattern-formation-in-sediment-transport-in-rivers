using LinearAlgebra
using SparseArrays
using Plots
using DifferentialEquations
using Statistics
using BenchmarkTools
using Printf

g = 9.80665    # [m/s²]
H = 0.1        # [m]
LeftX = 0
RightX = 1
tStart = 0
tEnd = 500

L = RightX - LeftX  # Domain length
Δx = 1/100          # Spatial spacing
Nx = Int64(L/Δx)    # Number of sub-intervals in space domain
μ₀ = 0.001          # Friction parameter -> r = μ₀*u²
A = 5               # Amplitude of tidal forcing
ω = 5               # Frequency of tidal forcing

x = [j for j in LeftX:Δx:RightX]   # include boundary points

function SystemMatrix(Nx)                           # A
    k = [1.0 for i in 1:Nx]                         # k=1 and k=-1 diagonal array
    A = SparseArrays.spdiagm(-1 => k, 1 => -k)

    A[1, end] = 1
    A[end, 1] = -1

    return A
end

println("Friction parameter amplitude: ", μ₀)
println("Amplitude of the tidal force: ", A)
println("Frequency of the tidal force: ", ω)

Aₓ = SystemMatrix(Nx)

### Using DifferentialEquations.jl ###
ζ⁰ = sin.(π .* x / L)    # ζ(x, t) at t = tStart
u⁰ = zeros(Nx+1, 1)      # dζ/dt at t = tStart
tspan = (tStart, tEnd)

z⁰ = vcat(ζ⁰, u⁰)

function RHS!(dz, z, p, t)
    ζ = z[1:Nx+1]
    u = z[Nx+2:end]
    F = (A * sin(ω*t)) .* ones(Nx+1, 1)
    dz[1:Nx+1] .= -H .* (Aₓ * u) 
    dz[Nx+2:end] .= F .- g .* (Aₓ * ζ) .- (μ₀/H) .* u.^3
end

problem = ODEProblem(RHS!, z⁰, tspan)

solution = solve(problem, Tsit5())

t_values = solution.t
println("Number of nodes in time: ", length(t_values))

function return_extrema(sol, slice_start, slice_end)
    sol_max = maximum(sol.u[1][slice_start:slice_end])
    sol_min = minimum(sol.u[1][slice_start:slice_end])
    
    for i in 2:length(sol.t)
        max = maximum(sol.u[i][slice_start:slice_end])
        min = minimum(sol.u[i][slice_start:slice_end]) 
        if max > sol_max
            sol_max = max
        end
        if min < sol_min
            sol_min = min
        end
    end
    
    return (sol_min, sol_max)
end

ζ_range = return_extrema(solution, 1, Nx+1)

animation = @animate for i in 1:length(t_values)
    ζₚ = solution.u[i][1:Nx+1]
    formatted_t = @sprintf("%.8f", solution.t[i])
    
    plot(x, ζₚ, ylims = ζ_range, xlabel = "x", ylabel = "ζ", title = "Time: $formatted_t [s]", legend = false)
end

current_directory = @__DIR__
gif_path = joinpath(current_directory, "animations", "wave_equation_1D_friction_nonlinear.gif")
gif(animation, gif_path, fps = 60)

Δt_array = collect(t_values[i+1]-t_values[i] for i in 1:length(t_values)-1)

println("Mean time step size: ", mean(Δt_array))
p = scatter(1:length(t_values)-1, Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", marker = :circle, markersize = 2, markercolor = :blue, markerstrokecolor = :blue, seriescolor = :blue, legend = false)
p1 = plot(1:length(t_values)-1, Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", marker = :circle, markersize = 2, markercolor = :blue, markerstrokecolor = :blue, seriescolor = :blue, legend = false)
p2 = bar(Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", legend = false)

savefig(p1,joinpath(current_directory, "time_steps_line_nonlinear.png"))
savefig(p2, joinpath(current_directory, "time_steps_bar_nonlinear.png"))
savefig(p, joinpath(current_directory, "time_steps_scatter_nonlinear.png"))
