function run_nonlinear(solver)
    function SystemMatrix(Nx, Δx)                       # A
        k = [1/(2*Δx) for i in 1:Nx]                    # k=1 and k=-1 diagonal array
        A = SparseArrays.spdiagm(-1 => -k, 1 => k)
    
        A[1, end] = -1/(2*Δx)
        A[end, 1] = 1/(2*Δx)
    
        return A
    end
    
    println("Friction parameter amplitude: ", μ₀)
    println("Amplitude of the tidal force: ", A)
    println("Frequency of the tidal force: ", ω)
    
    Aₓ = SystemMatrix(Nx, Δx)
    
    ### Using DifferentialEquations.jl ###
    ζ⁰ = sin.(π .* x / L)    # ζ(x, t) at t = tStart
    u⁰ = zeros(Nx+1, 1)      # dζ/dt at t = tStart
    tspan = (tStart, tEnd)
    
    z⁰ = vcat(ζ⁰, u⁰)        # combine ζ and u for the ODEProblem
    
    function RHS!(dz, z, p, t)
        Nx, A, ω, H, g, μ₀, Aₓ = p
        ζ = z[1:Nx+1]
        u = z[Nx+2:end]
        F = (A * sin(ω*t)) .* ones(Nx+1, 1)
        dz[1:Nx+1] .= -H .* (Aₓ * u)                            # dζ/dt
        dz[Nx+2:end] .= F .- g .* (Aₓ * ζ) .- (μ₀/H) .* u.^3    # du/dt
    end
    
    p = (; Nx, A, ω, H, g, μ₀, Aₓ)
    problem = ODEProblem(RHS!, z⁰, tspan, p)
    
    solution = solve(problem, solver)
    
    t_values = solution.t
    println("Number of nodes in time: ", length(t_values))
    
    function return_extrema(sol, slice_start, slice_end)
        sol_max = maximum(sol.u[1][slice_start:slice_end])
        sol_min = minimum(sol.u[1][slice_start:slice_end])
        
        for i in 2:length(sol.t)
            max = maximum(sol.u[i][slice_start:slice_end])
            min = minimum(sol.u[i][slice_start:slice_end]) 
            if max > sol_max
                sol_max = max
            end
            if min < sol_min
                sol_min = min
            end
        end
        
        return (sol_min, sol_max)
    end
    
    ζ_range = return_extrema(solution, 1, Nx+1)
    println("ζ extrema: ", ζ_range)
    
    animation = @animate for i in 1:length(t_values)
        ζₚ = solution.u[i][1:Nx+1]
        formatted_t = @sprintf("%.8f", solution.t[i])
        
        plot(x, ζₚ, ylims = ζ_range, xlabel = "x", ylabel = "ζ", title = "Time: $formatted_t [s]", legend = false)
    end
    
    
    gif_path = joinpath(@__DIR__, "..", "output", "animations", "wave_equation_1D_friction_nonlinear.gif")
    gif(animation, gif_path, fps = 60)
    
    Δt_array = collect(t_values[i+1]-t_values[i] for i in 1:length(t_values)-1)
    
    println("Mean time step size: ", mean(Δt_array))
    p = scatter(1:length(t_values)-1, Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", marker = :circle, markersize = 2, markercolor = :blue, markerstrokecolor = :blue, seriescolor = :blue, legend = false)
    p1 = plot(1:length(t_values)-1, Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", marker = :circle, markersize = 2, markercolor = :blue, markerstrokecolor = :blue, seriescolor = :blue, legend = false)
    p2 = bar(Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", legend = false)
    
    imagePath = joinpath(@__DIR__, "..", "output", "images")
    savefig(p1,joinpath(imagePath, "time_steps_line_nonlinear.png"))
    savefig(p2, joinpath(imagePath, "time_steps_bar_nonlinear.png"))
    savefig(p, joinpath(imagePath, "time_steps_scatter_nonlinear.png"))
end
