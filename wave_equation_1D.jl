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

T = tEnd - tStart
L = RightX - LeftX
c = 2*sqrt(g*H)   # wave propagation speed
Δt = 1/200
Δx = 1/100
r = c*(Δt/Δx)   # Courant number
Nt = Int64(T/Δt)
Nx = Int64(L/Δx)

x = [j for j in LeftX:Δx:RightX][2:end-1] # exclude boundary points
t = 0:Δt:T

function Laplacian1D(Nx, hx)
    return Array(Tridiagonal([1.0 for i in 1:Nx-2],[-2.0 for i in 1:Nx-1],[1.0 for i in 1:Nx-2])) # excluding 1/Δx^2
end

function time_march(C, un, un_m1)
    # C = 2*I + r^2*A where r = c * (Δt/Δx)
    un_p1 = C*un .- un_m1
    return un_p1
end

println("Wave propagation speed: ", c)
println("Courant number: ", r)

A = Laplacian1D(Nx, Δx)
u0 = sin.(π .* x / L)
du0 = (π/L).*cos.(π .* x / L)
# du0 = zeros(Nx-1, 1)
u1 = u0 .+ du0 .* Δt .+ (0.5*r^2) .* (A*u0)

u_plot = zeros(Nx - 1, Nt + 1)
u_plot[:, 1] = u0
u_plot[:, 2] = u1
U0 = u0
U1 = u1
U2 = nothing
Identity = 1* Matrix(I, Nx-1, Nx-1)
C = 2 .* Identity .+ r^2 .* A
for i in range(2, Nt) # from 2 to Nt (including)
    global U2 = time_march(C, U1, U0)
    global u_plot[:, i+1] = U2
    global U0 = U1
    global U1 = U2
end

# plotly()
# surface(x, t, u_plot', xlabel="Position (x)", ylabel="Time (t)", zlabel="Solution", title="Δt = $Δt Δx = $Δx")

pyplot()
animation = @animate for i in 1:Nt+1
    # surface(x, t[1:i], u_plot[:, 1:i]', xlabel="Position (x)", ylabel="Time (t)", zlabel="Solution", title="Δt = $Δt Δx = $Δx")
    t_current = t[i]
    formatted_t = @sprintf("%.3f", t_current)
    plot(x, u_plot[:, i], xlabel="Position (x)", ylabel="Solution", title="Time: $formatted_t", ylims=(-1.5, 1.5))
end

# Save the animation
gif(animation, "wave_equation_1D.gif", fps = 15)


### Using DifferentialEquations.jl ###
# u0 = zeros(Nx-1, 1)
# du0 = u0
# tspan = (tStart, tEnd)
# function RHS!(ddu, du, u, p, t)
#     ddu .= -courant*A*u 
# end
# problem = SecondOrderODEProblem(RHS!, du0, u0, tspan)
# solution = solve(problem, Euler(), dt = ht, save_everystep=false)
# solution = solve(problem, Tsit5(), save_everystep=false)
# solution = solve(problem, save_everystep=false)

# solutions = solution.u[end]
# duEnd = solutions[1:Nx-1]
# uEnd = solutions[Nx:end]

# plot(x, uEnd)