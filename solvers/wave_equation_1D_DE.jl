function run(solver)
    function Laplacian1D(Nx, hx)                               # A
        k = [1.0 for i in 1:Nx]                                # k=1 and k=-1 diagonal array
        A = Array(Tridiagonal(k, [-2.0 for i in 1:Nx+1], k))   # excluding 1/Δx^2
        A[1, end] = 1
        A[end, 1] = 1
    
        return A
    end
    
    println("Wave propagation speed: ", c)
    println("Courant number: ", r*Δt)
    
    A = Laplacian1D(Nx, Δx)
    
    ### Using DifferentialEquations.jl ###
    ζ⁰ = sin.(π .* x / L)    # ζ(x, t) at t = tStart
    u⁰ = zeros(Nx+1, 1)      # dζ/dt at t = tStart
    tspan = (tStart, tEnd)
    
    function RHS!(ddu, du, u, p, t)
        r, A = p
        ddu .= (r*r) .* (A*u)
    end
    
    p = (; r, A)
    problem = SecondOrderODEProblem(RHS!, u⁰, ζ⁰, tspan, p)
    solution = solve(problem, solver)
    
    t_values = solution.t
    println("Number of nodes in time: ", length(t_values))
    
    pyplot()
    animation = @animate for i in 1:length(t_values)
        formatted_t = @sprintf("%.3f", solution.t[i])
        
        plot(x, solution.u[i][Nx+2:end], xlabel = "x", ylabel = "ζ", title = "Time: $formatted_t", xlims = (LeftX, RightX), ylims = (0, 1))
    end
    
    # Save the animation
    gif_path = joinpath(@__DIR__, "..", "output", "animations", "wave_equation_1D_periodic_DE.gif")
    gif(animation, gif_path, fps = 15)
    
end