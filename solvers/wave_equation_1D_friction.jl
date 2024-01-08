function run_friction(solver)
    function Laplacian1D(Nx, hx)                               # A
        k = [1.0 for i in 1:Nx]                                # k=1 and k=-1 diagonal array
        A = Array(Tridiagonal(k, [-2.0 for i in 1:Nx+1], k))   # excluding 1/Δx^2
    
        A[1, end] = 1
        A[end, 1] = 1
    
        return A
    end
    
    println("Wave propagation speed: ", c)
    println("Courant number: ", r*Δt)
    println("Friction parameter: ", μ)
    
    A = Laplacian1D(Nx, Δx)
    
    ### Using DifferentialEquations.jl ###
    ζ⁰ = sin.(π .* x / L)    # ζ(x, t) at t = tStart
    u⁰ = zeros(Nx+1, 1)      # dζ/dt at t = tStart
    tspan = (tStart, tEnd)
    
    function RHS!(ddu, du, u, p, t)
        r, A, μ, H = p
        ddu .= r^2 * A * u .- (μ/H) .* du 
    end
    
    p = (; r, A, μ, H)
    problem = SecondOrderODEProblem(RHS!, u⁰, ζ⁰, tspan, p)
    
    solution = solve(problem, solver)
    
    t_values = solution.t
    println("Number of nodes in time: ", length(t_values))
    
    function return_extrema(sol, Nx)
        sol_max = maximum(sol.u[1][Nx+2:end])
        sol_min = minimum(sol.u[1][Nx+2:end])
        
        for i in 2:length(sol.t)
            max = maximum(sol.u[i][Nx+2:end])
            min = minimum(sol.u[i][Nx+2:end]) 
            if max > sol_max
                sol_max = max
            end
            if min < sol_min
                sol_min = min
            end
        end
        
        return (sol_min, sol_max)
    end
    
    ζ_range = return_extrema(solution, Nx)
    
    animation = @animate for i in 1:length(t_values)
        ζₚ = solution.u[i][Nx+2:end]
        formatted_t = @sprintf("%.8f", solution.t[i])
        
        plot(x, ζₚ, ylims = ζ_range, title = "Time: $formatted_t", legend = false)
    end
    
    gif_path = joinpath(@__DIR__, "..", "output", "animations", "wave_equation_1D_friction.gif")
    gif(animation, gif_path, fps = 15)
    
    Δt_array = collect(t_values[i+1]-t_values[i] for i in 1:length(t_values)-1)
    
    println("Mean time step size: ", mean(Δt_array))
    p = scatter(1:length(t_values)-1, Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", marker = :circle, markersize = 2, markercolor = :blue, markerstrokecolor = :blue, seriescolor = :blue, legend = false)
    p1 = plot(1:length(t_values)-1, Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", marker = :circle, markersize = 2, markercolor = :blue, markerstrokecolor = :blue, seriescolor = :blue, legend = false)
    p2 = bar(Δt_array, xaxis = "Step number", yaxis = "Δt", title = "Time Step Size", legend = false)
    
    savefig(p, joinpath(@__DIR__, "..", "output", "images", "time_steps_scatter.png"))
    savefig(p1, joinpath(@__DIR__, "..", "output", "images", "time_steps_line.png"))
    savefig(p2, joinpath(@__DIR__, "..", "output", "images", "time_steps_bar.png"))
end