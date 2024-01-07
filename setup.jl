using Pkg

function setup()
    println("Now installing packages...")

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
    packages_to_add = ["LinearAlgebra", "SparseArrays", "Plots", "DifferentialEquations", "BenchmarkTools", "Printf", "Statistics", "InteractiveUtils"]
    
    # Add packages
    for pkg in packages_to_add
        add_package(pkg)
    end

    println("All packages installed.")
    println("-------------------------------------")
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

function simplestModel()
    println("Now running the simplest model...")
    include(joinpath(@__DIR__, "solvers", "wave_equation_1D_DE.jl"))
    println("Finished running the simplest model...")
    println("---------------------------------------------------------------")
end

function waveEquationFriction()
    println("Now running the model with friction term...")
    include(joinpath(@__DIR__, "solvers", "wave_equation_1D_friction.jl"))
    println("Finished running the model with friction term...")
    println("---------------------------------------------------------------")
end

function waveEquationNonlinear()
    println("Now running the nonlinear model...")
    include(joinpath(@__DIR__, "solvers", "wave_equation_1D_friction_nonlinear_tidal.jl"))
    println("Finished running the nonlinear model...")
    println("---------------------------------------------------------------")
end