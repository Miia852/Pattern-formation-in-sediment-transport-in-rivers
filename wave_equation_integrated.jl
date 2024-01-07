using Pkg

function setup()
    function add_package(pkg::String)
        try
            @eval import $(Symbol(pkg))
            println("$pkg is already installed.")
        catch
            println("Installing $pkg...")
            Pkg.add(pkg)
            @eval using $pkg
        end
    end
    
    # List of packages to add
    packages_to_add = ["LinearAlgebra", "SparseArrays", "Plots", "DifferentialEquations", "BenchmarkTools", "Printf", "Statistics"]
    
    # Add packages
    for pkg in packages_to_add
        add_package(pkg)
    end
end


## Run this function to import the packages if they are not installed yet
setup()

using LinearAlgebra
using SparseArrays
using Plots
using DifferentialEquations
using BenchmarkTools
using Printf
using Statistics

## general parameters
g = 9.80665 # [m/s^2]
H = 0.1 # [m]
LeftX = 0
RightX = 1
tStart = 0
tEnd = 500
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
μ₀ = 0.001          # Friction parameter -> r = μ₀*u²
A = 5               # Amplitude of tidal forcing
ω = 5               # Frequency of tidal forcing

function simplestModel()
    println("Now running the simplest model...")
    include("wave_equation_1D_DE.jl")
    println("Finished running the simplest model...")
    println("---------------------------------------------------------------")
end

function waveEquationFriction()
    println("Now running the model with friction term...")
    include("wave_equation_1D_friction.jl")
    println("Finished running the model with friction term...")
    println("---------------------------------------------------------------")
end

function waveEquationNonlinear()
    println("Now running the nonlinear model...")
    include("wave_equation_1D_friction_nonlinear_tidal.jl")
    println("Finished running the nonlinear model...")
    println("---------------------------------------------------------------")
end

# simplestModel()
waveEquationFriction()
waveEquationNonlinear()


