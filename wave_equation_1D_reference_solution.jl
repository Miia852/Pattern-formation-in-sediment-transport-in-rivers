using Plots

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
Nt = Int64(T/Δt)    # number of sub-intervals in time domain
Nx = Int64(L/Δx)    # number of sub-intervals in space domain

ζₚ = zeros(Nx + 1, Nt + 1)         # zeta_plot
x = [j for j in LeftX:Δx:RightX]   # include boundary points
t = 0:Δt:T                         # Time values array

for j in range(1, length(t))
    for i in range(1, length(x))
        ζₚ[i, j] = sin(π*x[i])*cos(π*t[j])
    end
end

# Get the range of ζ values across all frames
ζ_range = extrema(ζₚ)

plotly()
surface(x, t, ζₚ', xlabel = "x", ylabel = "t", zlabel = "ζ", title = "Δt = $Δt Δx = $Δx", xlims = (LeftX, RightX), ylims = (tStart, tEnd), zlims = ζ_range, clim = ζ_range)
