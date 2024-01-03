using DataFrames, CSV
using Statistics

LeftX = 0
RightX = 1
tStart = 0
tEnd = 1

T = tEnd - tStart   # Simulation time
L = RightX - LeftX  # Domain length
Δt = 1/200          # time step
Δx = 1/100          # spatial spacing
Nt = Int64(T/Δt)    # number of sub-intervals in time domain
Nx = Int64(L/Δx)    # number of sub-intervals in space domain

x = [j for j in LeftX:Δx:RightX]   # include boundary points
t = 0:Δt:T                         # Time values array

numerical = Matrix(CSV.read("wave_1D_numerical_solution.csv", DataFrame))
reference = Matrix(CSV.read("wave_1D_reference_solution.csv", DataFrame))

abs_diff = abs.(numerical .- reference)

println("Maximum difference\t:\t", maximum(abs_diff))
println("Mean difference   \t:\t", mean(abs_diff))

ζ_range = extrema(numerical)

pyplot()
animation = @animate for i in 1:Nt+1
    formatted_t = @sprintf("%.3f", t[i])
    
    plot(x, numerical[i, :], xlabel = "x", ylabel = "ζ", title = "Time: $formatted_t", xlims = (LeftX, RightX), ylims = ζ_range)
    plot!(reference[i, :])
end

# Save the animation
gif(animation, "solution_comparison.gif", fps = 15)
