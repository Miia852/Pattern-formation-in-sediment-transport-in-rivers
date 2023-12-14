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
c = 2*sqrt(g*H)     # wave propagation speed
Δt = 1/200          # time step
Δx = 1/100          # spatial spacing
r = c*(Δt/Δx)       # Courant number
Nt = Int64(T/Δt)    # number of sub-intervals in time domain
Nx = Int64(L/Δx)    # number of sub-intervals in space domain

x = [j for j in LeftX:Δx:RightX][2:end-1]   # exclude boundary points
t = 0:Δt:T                                  # Time values array

function Laplacian1D(Nx, hx)    # A
    k = [1.0 for i in 1:Nx-2]                                 # k=1 and k=-1 diagonal array
    return Array(Tridiagonal(k, [-2.0 for i in 1:Nx-1], k))   # excluding 1/Δx^2
end

function time_march(C, un, un_m1)
    # C = 2*I + r^2*A where r = c * (Δt/Δx)
    un_p1 = C*un .- un_m1
    return un_p1
end

println("Wave propagation speed: ", c)
println("Courant number: ", r)

A = Laplacian1D(Nx, Δx)

# Initial conditions
ζ⁰ = sin.(π .* x / L)           # ζ(x, t) at t = tStart so, ζ(x, t)
u⁰ = (π/L).*cos.(π .* x / L)    # dζ/dx at t = tStart

# ζ⁰ = zeros(Nx-1, 1)
# u⁰ = zeros(Nx-1, 1)             # dζ/dx at t = tStart
# ζ⁰[49] = 0.1                    # to form a triangle
# ζ⁰[50] = 0.2                    # to form a triangle
# ζ⁰[51] = 0.1                    # to form a triangle

ζ¹ = ζ⁰ .+ u⁰ .* Δt .+ (0.5*(r/c)^2) .* (A*ζ⁰) # r/c = Δt/Δx => (0.5*(Δt/Δx)^2) .* (A*u)

ζₚ = zeros(Nx - 1, Nt + 1)  # zeta_plot
ζₚ[:, 1] = ζ⁰
ζₚ[:, 2] = ζ¹

ζⁿ⁻¹ = ζ⁰          # zeta_n_minus_one
ζⁿ = ζ¹            # zeta_n
ζⁿ⁺¹ = nothing     # zeta_n_plus_one

Identity = 1* Matrix(I, Nx-1, Nx-1)
C = 2 .* Identity .+ r^2 .* A

for i in range(2, Nt) # from 2 to Nt (including)
    global ζⁿ⁺¹ = time_march(C, ζⁿ, ζⁿ⁻¹)
    global ζₚ[:, i+1] = ζⁿ⁺¹
    global ζⁿ⁻¹ = ζⁿ
    global ζⁿ = ζⁿ⁺¹
    # ζⁿ⁺¹ and ζⁿ of the current step are ζⁿ and ζⁿ⁻¹ of the next step
end

# Get the range of ζ values across all frames
ζ_range = extrema(ζₚ)

# plotly()
# surface(x, t, ζₚ', xlabel = "x", ylabel = "t", zlabel = "ζ", title = "Δt = $Δt Δx = $Δx", xlims = (LeftX, RightX), ylims = (tStart, tEnd), zlims = ζ_range, clim = ζ_range)

pyplot()
animation = @animate for i in 1:Nt+1
    formatted_t = @sprintf("%.3f", t[i])
    
    # surface(x, t[1:i], ζₚ[:, 1:i]', xlabel = "x", ylabel = "t", zlabel = "ζ", title = "Time: $formatted_t", xlims = (LeftX, RightX), ylims = (tStart, tEnd), zlims = ζ_range, clim = ζ_range)
    plot(x, ζₚ[:, i], xlabel = "x", ylabel = "ζ", title = "Time: $formatted_t", xlims = (LeftX, RightX), ylims = ζ_range)
end

# animation = @animate for i in 1:Nx-1
#     formatted_x = @sprintf("%.3f", x[i])

#     surface(x[1:i], t, ζₚ[1:i, :]', xlabel = "x", ylabel = "t", zlabel = "ζ", title = "x: $formatted_x", xlims = (LeftX, RightX), ylims = (tStart, tEnd), zlims = ζ_range, clim = ζ_range)
#     # plot(t, ζₚ[i, :], xlabel = "t", ylabel = "ζ", title = "x: $formatted_x", xlims = (tStart, tEnd), ylims = ζ_range)
# end

# Save the animation
gif(animation, "animations/wave_equation_1D_sin(pix_L)_zetadx.gif", fps = 15)


### Using DifferentialEquations.jl ###
# ζ⁰ = zeros(Nx-1, 1)
# u⁰ = ζ⁰
# tspan = (tStart, tEnd)
# function RHS!(ddu, du, u, p, t)
#     ddu .= -courant*A*u 
# end
# problem = SecondOrderODEProblem(RHS!, u⁰, ζ⁰, tspan)
# solution = solve(problem, Euler(), dt = ht, save_everystep=false)
# solution = solve(problem, Tsit5(), save_everystep=false)
# solution = solve(problem, save_everystep=false)

# solutions = solution.u[end]
# duEnd = solutions[1:Nx-1]
# uEnd = solutions[Nx:end]

# plot(x, uEnd)