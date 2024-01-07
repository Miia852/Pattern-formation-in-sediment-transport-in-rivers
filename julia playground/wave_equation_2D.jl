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
LeftY = 0
RightY = 1
tStart = 0
tEnd = 1

T = tEnd - tStart
Lₓ = RightX - LeftX
Ly = RightY - LeftY
c = 2*sqrt(g*H)     # wave propagation speed
Δt = 1/200          # time step
Δx = 1/50           # spatial spacing in x
Δy = 1/50           # spatial spacing in y
rx = c*(Δt/Δx)      # Courant number in x
ry = c*(Δt/Δy)      # Courant number in y
Nₜ = Int64(T/Δt)     # number of sub-intervals in time domain  
Nₓ = Int64(Lₓ/Δx)   # number of sub-intervals in x domain 
Ny = Int64(Ly/Δy)   # number of sub-intervals in y domain

t = 0:Δt:T    # Time values array 

# Create grid matrices
x_points = [j for j in LeftX:Δx:RightX][2:end-1]
y_points = [i for i in LeftY:Δy:RightY][2:end-1]
x = zeros((Ny-1, Nₓ-1))
y = zeros((Ny-1, Nₓ-1))
for j in 1:length(y_points)
    x[j, :] = x_points
end
for i in 1:length(x_points)
    y[:, i] = y_points
end

function Laplacian_x(Nₓ, Ny)    # Ax
    matrix_size = (Nₓ-1)*(Ny-1)
    k = [1.0 for i in 1:matrix_size-1]                             # k=1 and k=-1 diagonal array
    return Array(Tridiagonal(k, [-2.0 for i in 1:matrix_size], k)) # excluding 1/Δx^2
end

function Laplacian_y(Nₓ, Ny)    # Ay
    matrix_size = (Nₓ-1)*(Ny-1)
    k0 = [-2.0 for i in 1:matrix_size]                      # main diagonal array
    k = ones(matrix_size-(Nₓ-1))                            # k=1 and k=-1 diagonal array
    result = spdiagm(0 => k0, (-Nₓ+1) => k, (Nₓ-1) => k)    # excluding 1/Δy^2
    return result
end

function time_march(C, uⁿ, uⁿ⁻¹)
    # C = 2*I + rx^2*Aₓ + ry^2*Ay 
    uⁿ⁺¹ = C*uⁿ .- uⁿ⁻¹
    return uⁿ⁺¹
end

println("Wave propagation speed: ", c)
println("Courant number in x: ", rx)
println("Courant number in y: ", ry)

Aₓ = Laplacian_x(Nₓ, Ny)
Ay = Laplacian_y(Nₓ, Ny)

# Initial conditions
# ζ⁰ = sin.(π .* x ./ Lₓ) .+ sin.(π .* y ./ Ly)    # ζ(x, y, t) at t = tStart so, ζ(x, y, 0)
# u⁰ = (π/Lₓ).*cos.(π .* x ./ Lₓ)                  # dζ/dx at t = tStart
# v⁰ = (π/Ly).*cos.(π .* y ./ Ly)                  # dζ/dy at t = tStart
# ζ⁰_xy = zeros((Nₓ-1)*(Ny-1), 1)                  # dζ/dxdy at t = tStart   

# ζ⁰ = sin.(π .* x ./ Lₓ)                       # ζ(x, y, t) at t = tStart so, ζ(x, y, 0)
# ζ⁰ = sin.(π .* y ./ Ly)                       # ζ(x, y, t) at t = tStart so, ζ(x, y, 0)
# u⁰ = zeros((Nₓ-1)*(Ny-1), 1)               
# v⁰ = zeros((Nₓ-1)*(Ny-1), 1)               

ζ⁰ = sin.(π .* x ./ Lₓ) .* sin.(π .* y ./ Ly)    # ζ(x, y, t) at t = tStart so, ζ(x, y, 0)
u⁰ = (π/Lₓ).*cos.(π .* x ./ Lₓ) .* sin.(π .* y ./ Ly)   # dζ/dx at t = tStart
v⁰ = sin.(π .* x ./ Lₓ) .* (π/Ly) .* cos.(π .* y ./ Ly) # dζ/dy at t = tStart
ζ⁰_xy = (π/Lₓ).*cos.(π .* x ./ Lₓ) .* (π/Ly) .* cos.(π .* y ./ Ly) # dζ/dxdy at t = tStart   

ζ⁰ = vec(ζ⁰')         # flatten the matrix to 1D
u⁰ = vec(u⁰')         # flatten the matrix to 1D
v⁰ = vec(v⁰')         # flatten the matrix to 1D
ζ⁰_xy = vec(ζ⁰_xy')   # flatten the matrix to 1D

ζ¹ = ζ⁰ .+ u⁰ .* Δt .+ v⁰ .* Δt .+ 0.5 * ((rx / c)^2 .* (Aₓ * ζ⁰) + (ry / c)^2 .* (Ay * ζ⁰)) .+ ζ⁰_xy .* Δt^2  

ζₚ = zeros((Nₓ - 1)*(Ny - 1), Nₜ + 1)   # zeta_plot
ζₚ[:, 1] = ζ⁰
ζₚ[:, 2] = ζ¹

ζⁿ⁻¹ = ζ⁰          # zeta_n_minus_one
ζⁿ = ζ¹            # zeta_n
ζⁿ⁺¹ = nothing     # zeta_n_plus_one

Identity = 1* Matrix(I, (Nₓ-1)*(Ny-1), (Nₓ-1)*(Ny-1))
C = @. 2 * Identity + rx^2 * Aₓ + ry^2 * Ay

for i in range(2, Nₜ) # from 2 to Nₜ (including)
    global ζⁿ⁺¹ = time_march(C, ζⁿ, ζⁿ⁻¹)
    global ζₚ[:, i+1] = ζⁿ⁺¹
    global ζⁿ⁻¹ = ζⁿ
    global ζⁿ = ζⁿ⁺¹
    # ζⁿ⁺¹ and ζⁿ of the current step are ζⁿ and ζⁿ⁻¹ of the next step, respectively
end

# Get the range of ζ values across all frames
ζ_range = extrema(ζₚ)

# plotly()
# # surface(x, y, reshape(ζₚ[:, end], Nₓ-1, Ny-1)', xlabel="Position (x)", ylabel="Position (y)", zlabel="Solution (ζ)", title="Δt = $Δt Δx = $Δx Δy = $Δy")
# surface(x, y, reshape(ζₚ[:, 1], Nₓ-1, Ny-1)', xlabel="x", ylabel="y", zlabel="ζ", title="Δt = $Δt Δx = $Δx Δy = $Δy", xlims = (LeftX, RightX), ylims = (LeftY, RightY), zlims = ζ_range, clim = ζ_range)

pyplot()
animation = @animate for i in 1:Nₜ+1
    formatted_t = @sprintf("%.3f", t[i])    # always use 3 decimal points
    ζ = reshape(ζₚ[:, i], Nₓ-1, Ny-1)'       # convert 1D array to 2D matrix

    surface(x, y, ζ, xlabel = "x", ylabel = "y", zlabel = "ζ", title = "Time: $formatted_t", xlims = (LeftX, RightX), ylims = (LeftY, RightY), zlims = ζ_range, clim = ζ_range)    
end

# Save the animation
gif(animation, "animations/wave_equation_2D_sin(pix_Lx)tsin(piy_Ly)_zetadx_zetady_0.gif", fps = 15)
