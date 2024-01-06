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

x = [j for j in LeftX:Δx:RightX]   # include boundary points

function Laplacian1D(Nx, hx)                               # A
    k = [1.0 for i in 1:Nx]                                # k=1 and k=-1 diagonal array
    A = Array(Tridiagonal(k, [-2.0 for i in 1:Nx+1], k))   # excluding 1/Δx^2
    
    # A[1, 1:end] .= 0
    # A[end, 1:end] .= 0

    A[1, end] = 1
    A[end, 1] = 1

    return A
end

println("Wave propagation speed: ", c)
println("Courant number: ", r*Δt)

A = Laplacian1D(Nx, Δx)

### Using DifferentialEquations.jl ###
ζ⁰ = sin.(π .* x / L)    # ζ(x, t) at t = tStart
u⁰ = zeros(Nx+1, 1)      # dζ/dt at t = tStart
tspan = (tStart, tEnd)

function RHS!(ddu, du, u, p, t)
    ddu .= (r*r) .* (A*u)
end

problem = SecondOrderODEProblem(RHS!, u⁰, ζ⁰, tspan)
solution = solve(problem, Tsit5())

t_values = solution.t
println("Number of nodes in time: ", length(t_values))

pyplot()
animation = @animate for i in 1:length(t_values)
    formatted_t = @sprintf("%.3f", solution.t[i])
    
    plot(x, solution.u[i][Nx+2:end], xlabel = "x", ylabel = "ζ", title = "Time: $formatted_t", xlims = (LeftX, RightX), ylims = (0, 1))
end

# Save the animation
gif(animation, "animations/wave_equation_1D_periodic_DE.gif", fps = 15)
