include(joinpath(@__DIR__, "setup.jl"))

## general parameters
g = 9.80665 # [m/s^2]
H = 0.1 # [m]
LeftX = 0
RightX = 1
tStart = 0
tEnd = 10
T = tEnd - tStart   # Simulation time
L = RightX - LeftX  # Domain length
c = sqrt(g*H)       # wave propagation speed
Nt = 200    # number of sub-intervals in time domain
Nx = 100    # number of sub-intervals in space domain
Δt = T/Nt          # time step
Δx = L/Nx          # spatial spacing
r = c/Δx            # Courant number adjusted
x = [j for j in LeftX:Δx:RightX]   # include boundary points

## friction and forcing parameters
μ = 0.5
μ₀ = 0.5            # Friction parameter amplitude -> r = μ₀*u²
A = 5               # Amplitude of tidal forcing   -> F = A*sin(ω*t)
ω = 5               # Frequency of tidal forcing   -> F = A*sin(ω*t)

# simplestModel()
waveEquationFriction()
waveEquationNonlinear()