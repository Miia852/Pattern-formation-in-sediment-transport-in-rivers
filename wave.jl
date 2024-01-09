using LinearAlgebra
using SparseArrays
using Plots
using DifferentialEquations
using Statistics
using BenchmarkTools
using Printf
using StaticArrays
using Sundials
using Trapz

g = 9.80665    # [m/s²]
LeftX = 0
RightX = 1
tStart = 0
tEnd = 100

L = RightX - LeftX  # Domain length
Δx = 1/200         # Spatial spacing
Nx = Int64(L/Δx)    # Number of sub-intervals in space domain
μ = 0.5             # Friction parameter
A = 5               # Amplitude of tidal forcing   -> F = A*sin(ω*t)
ω = 5               # Frequency of tidal forcing   -> F = A*sin(ω*t)
x = [j for j in LeftX:Δx:RightX][1:end-1]   # include boundary points
H = -0.1 .* (atan.(30 .* (x .- 0.3)) .- atan.(20 .* (x .-0.5))).+2
algorithms = [
    # Tsit5(),         # Tsitouras 5/4 Runge-Kutta method
    # DP5(),           # Dormand-Prince 5(4) Runge-Kutta method
    RK4(),           # Classic 4th order Runge-Kutta method
    Vern7(),         # Vernon 7/8 Runge-Kutta method
    # CVODE_BDF(),     # CVODE's BDF method  
    # Midpoint(),      # Midpoint method (2nd order)  ## time-consuming
    BS3(),           # Bogacki-Shampine 3(2) Runge-Kutta method
    # CVODE_Adams(),   # CVODE's Adams method
]


function SystemMatrix(Nx, Δx)                       # A
    k = [1/(2*Δx) for i in 1:Nx-1]                    # k=1 and k=-1 diagonal array
    A = SparseArrays.spdiagm(-1 => -k, 1 => k)

    A[1, end] = -1/(2*Δx)
    A[end, 1] = 1/(2*Δx)

    return A
end

println("Friction parameter: ", μ)
println("Amplitude of the tidal force: ", A)
println("Frequency of the tidal force: ", ω)

Aₓ = SystemMatrix(Nx, Δx)

### Using DifferentialEquations.jl ###
# ζ⁰ = sin.(π .* x / L)    # ζ(x, t) at t = tStart
# ζ⁰ = sin.((2*π) .* x / L)    # ζ(x, t) at t = tStart?
ζ⁰ = sin.(π .* x / L)
u⁰ = zeros(Nx, 1)        # dζ/dt at t = tStart
tspan = (tStart, tEnd)

z⁰ = MVector(ζ⁰..., u⁰...)         # combine ζ and u for the ODEProblem

function RHS!(dz, z, p, t)
    Nx, A, ω, H, g, μ, Aₓ = p
    ζ = @view z[1:Nx]
    u = @view z[Nx+1:end]
    F = (A * cos(ω*t)) .* ones(Nx, 1)
    dz[1:Nx] .= -Aₓ * (H .* u)                            # dζ/dt      
    dz[Nx+1:end] .= F .- g .* (Aₓ * ζ) .- μ .* (u ./ H)   # du/dt
end

p = (; Nx, A, ω, H, g, μ, Aₓ)
problem = ODEProblem(RHS!, z⁰, tspan, p)

### This is for testing which algorithm could lead to convergence
# function run(solver)
#     println("Starting with " * string(solver))

#     solution = solve(problem, solver)

#     uₚ = solution.u[end][Nx+1:end]
#     ζₚ = solution.u[end][1:Nx]

#     plot(x, ζₚ)
#     savefig(joinpath(@__DIR__, "plots_t_end=50", split(string(solver), r"[\({]")[1]*"__zeta.png")) 

#     plot(x, uₚ, ylims = (-1.1, 1.1))
#     savefig(joinpath(@__DIR__, "plots_t_end=50", split(string(solver), r"[\({]")[1]*"__u.png"))  
#     println("Finished with " * string(solver))
# end
# for solver in algorithms
#     try
#         run(solver)
#     catch e
#         println("Error: $e")
#     end
    
# end

solution = solve(problem, RK4())

t_values = solution.t
println("Number of nodes in time: ", length(t_values))
println("Final time value: ", t_values[end])

## calculate the area under the curve
ζₚ = solution.u[end][1:Nx]
area_0 = trapz(x, ζ⁰)
area_end = trapz(x, ζₚ)
diff = area_0-area_end
println("Initial area: $area_0")
println("Area at the end: $area_end")
println("Difference: $diff")


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

ζ_range = return_extrema(solution, 1, Nx)
u_range = return_extrema(solution, Nx+1, 2*Nx)

println("ζ extrema: ", ζ_range)
Δt_array = collect(t_values[i+1]-t_values[i] for i in 1:length(t_values)-1)
println("Mean time step size: ", mean(Δt_array))

step = 5

animation = @animate for i in 1:step:length(t_values)
    uₚ = solution.u[i][Nx+1:end]
    formatted_t = @sprintf("%.8f", solution.t[i])
    
    plot(x, uₚ, ylims = u_range, xlabel = "x", ylabel = "ζ", title = "Time: $formatted_t [s]", legend = false)
end

## calculating framesPerSecond to make the animation consistent with real time
frames = floor(length(collect(1:step:length(t_values)))/tEnd)

gif(animation, "test_u.gif", fps = frames)

animation1 = @animate for i in 1:step:length(t_values)
    ζₚ = solution.u[i][1:Nx]
    formatted_t = @sprintf("%.8f", solution.t[i])
    
    plot(x, ζₚ, ylims = ζ_range, xlabel = "x", ylabel = "ζ", title = "Time: $formatted_t [s]", legend = false)
end

gif(animation1, "test_zeta.gif", fps = frames)

p = scatter(1:length(t_values)-1, Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", marker = :circle, markersize = 2, markercolor = :blue, markerstrokecolor = :blue, seriescolor = :blue, legend = false)
p1 = plot(1:length(t_values)-1, Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", marker = :circle, markersize = 2, markercolor = :blue, markerstrokecolor = :blue, seriescolor = :blue, legend = false)
p2 = bar(Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", legend = false)

savefig(p1, "test_time_steps_line.png")
savefig(p2, "test_time_steps_bar.png")
savefig(p, "test_time_steps_scatter.png")
