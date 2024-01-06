using LinearAlgebra
using SparseArrays
using Plots
using DifferentialEquations
using BenchmarkTools
using Printf

g = 9.80665 # [m/s^2]
H = 0.1 # [m]
LeftX = 0
RightX = 1
tStart = 0
tEnd = 1

T = tEnd - tStart   # Simulation time
L = RightX - LeftX  # Domain length
c = sqrt(g*H)       # wave propagation speed
Δt = 1/200          # time step
Δx = 1/100          # spatial spacing
r = c/Δx            # Courant number adjusted
Nt = Int64(T/Δt)    # number of sub-intervals in time domain
Nx = Int64(L/Δx)    # number of sub-intervals in space domain
μ = 0.5             # friction parameter

x = [j for j in LeftX:Δx:RightX]   # include boundary points

function Laplacian1D(Nx, hx)    # A
    k = [1.0 for i in 1:Nx]                                # k=1 and k=-1 diagonal array
    A = Array(Tridiagonal(k, [-2.0 for i in 1:Nx+1], k))   # excluding 1/Δx^2

    # # A[1, 1:end] .= 0
    # # A[end, 1:end] .= 0

    A[1, end] = 1
    A[end, 1] = 1

    return A

end

println("Wave propagation speed: ", c)
println("Courant number: ", r*Δt)
println("Friction parameter: ", μ)

A = Laplacian1D(Nx, Δx)

### Using DifferentialEquations.jl ###
ζ⁰ = sin.(π .* x / L)    # ζ(x, t) at t = tStart
u⁰ = zeros(Nx+1, 1)      # dζ/dt at t = tStart
tspan = (tStart, tEnd)

function RHS!(ddu, du, u, p, t)
    ddu .= r^2 * A * u .- (μ/H) .* du 
end

problem = SecondOrderODEProblem(RHS!, u⁰, ζ⁰, tspan)

# solution = solve(problem, Euler(), dt = Δt, save_everystep=true)
solution = solve(problem, Tsit5())
# solution = solve(problem, save_everystep=false)

# solutions = solution.u[end]
# duEnd = solutions[1:Nx+1]
# uEnd = solutions[Nx+2:end]

# plot(x, uEnd)

t_values = solution.t
println("Number of nodes in time: ", length(t_values))

animation = @animate for i in 1:length(t_values)
    u_k = solution.u[i][Nx+2:end]
    formatted_t = @sprintf("%.8f", solution.t[i])
    
    plot(x, u_k, ylims=(0, 1), title="Time: $formatted_t")
end

gif(animation, "animations/wave_equation_1D_friction.gif", fps = 15)

Δt_array = collect(solution.t[i+1]-solution.t[i] for i in 1:length(solution.t)-1)

println("Mean time step size: ", mean(Δt_array))
# p = scatter(1:length(t_values), Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", markersize = 2)
p1 = plot(1:length(t_values), Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", marker = :circle, markersize = 2, markercolor = :blue, markerstrokecolor = :blue, seriescolor = :blue, legend = false)
p2 = bar(Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", legend = false)

savefig(p1, "time_steps_line.png")
savefig(p2, "time_steps_bar.png")