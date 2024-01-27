using DataFrames, CSV
using Statistics
using Plots
using Printf
using BenchmarkTools
using Printf

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


numerical_path = joinpath(@__DIR__, "output", "csv", "wave_1D_numerical_solution.csv")
reference_path = joinpath(@__DIR__, "output", "csv", "wave_1D_reference_solution.csv")

# Read CSV files
numerical = Matrix(CSV.read(numerical_path, DataFrame))
reference = Matrix(CSV.read(reference_path, DataFrame))

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
output_directory = joinpath(@__DIR__, "output", "images")
gif_filename = "solution_comparison.gif"
gif(animation, joinpath(output_directory, gif_filename), fps = 15)
