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
    packages_to_add = ["LinearAlgebra", "SparseArrays", "Plots", "DifferentialEquations", "BenchmarkTools", "Printf", "Statistics", "InteractiveUtils", "Sundials", "DataFrames", "CSV"]
    
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
using Sundials
using DataFrames
using CSV

include(joinpath(@__DIR__, "solvers", "wave_equation_1D_DE.jl"))
include(joinpath(@__DIR__, "solvers", "wave_equation_1D_friction.jl"))
include(joinpath(@__DIR__, "solvers", "wave_equation_1D_friction_nonlinear_tidal.jl"))
include(joinpath(@__DIR__, "solvers", "algorithm_picker.jl"))

function simplestModel()
    println("Now running the simplest model...")
    try
        run_basic(Tsit5(), true) # make_animation = true
        println("Finished running the simplest model...")
    catch e
        @warn "Error: $e"
    end
    println("---------------------------------------------------------------")
end

function waveEquationFriction()
    println("Now running the model with friction term...")
    try
        run_friction(Tsit5(), true) # make_animation = true
        println("Finished running the model with friction term...")
    catch e
        @warn "Error: $e"
    end
    println("---------------------------------------------------------------")
end

function waveEquationNonlinear()
    println("Now running the nonlinear model...")
    try
        run_nonlinear(Tsit5(), true) # make_animation = true
        println("Finished running the nonlinear model...")
    catch e
        @warn "Error: $e"
    end
    println("---------------------------------------------------------------")
end