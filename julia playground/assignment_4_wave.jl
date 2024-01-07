using LinearAlgebra
using SparseArrays
using Plots
using DifferentialEquations
using BenchmarkTools
gr()

Nx = 200
Ny = 100
LeftX = 0
LeftY = 0
RightX = 10
RightY = 5
hx = (RightX - LeftX)/Nx 
hy = (RightY - LeftY)/Ny 
tStart = 0
tEnd = 4
# Nt = tEnd*1600

x_points = [j for j in LeftX:hx:RightX][2:end-1]
y_points = [i for i in LeftY:hy:RightY][2:end-1]

x = zeros((Ny-1, Nx-1))
y = zeros((Ny-1, Nx-1))

for j in 1:length(y_points)
    x[j, :] = x_points
end

for i in 1:length(x_points)
    y[:, i] = y_points
end

function sourcefunc(x, y, t)
    f = sin(4*π*t) .* (exp.(-40.0 * (x .- 0.25 * RightX) .^ 2 .- 40.0 * (y .- 0.25 * RightY) .^ 2) +
        exp.(-40.0 * (x .- 0.25 * RightX) .^ 2 .- 40.0 * (y .- 0.75 * RightY) .^ 2) +
        exp.(-40.0 * (x .- 0.75 * RightX) .^ 2 .- 40.0 * (y .- 0.75 * RightY) .^ 2) +
        exp.(-40.0 * (x .- 0.75 * RightX) .^ 2 .- 40.0 * (y .- 0.25 * RightY) .^ 2))
    return f
end

function coeffK(x, y)
    if !(x isa Array)
        if x<(RightX/2) && y<(RightY/2)
            return 0.1
        elseif x<(RightX/2) && y>=(RightY/2)
            return 0.4
        elseif x>=(RightX/2) && y>=(RightY/2)
            return 0.7
        else
            return 1.0
        end
    else
        result = zeros(size(x))
        line_x = findall(x.<(RightX/2))[end][2]
        line_y = findall(y.<(RightY/2))[end][1]
        result[1:line_y, 1:line_x] .= 0.1
        result[line_y+1:end, 1:line_x] .= 0.4
        result[line_y+1:end, line_x+1:end] .= 0.7
        result[1:line_y, line_x+1:end] .= 1.0
    end
    return result
end

function create2DLFVM(Nx, hx, hy, x, y, coeffFun)
    main_arr = ((1 / (hx * hx)) .* (coeffFun(x .- 0.5 * hx, y) .+ coeffFun(x .+ 0.5 * hx, y)) .+ (1 / (hy * hy)) .* (coeffFun(x, y .- 0.5 * hy) .+ coeffFun(x, y .+ 0.5 * hy)))
    main_arr = vec(main_arr')
    k3_arr = -(coeffFun(x, y .+ 0.5 * hy)) ./ (hy * hy)
    k3_arr = vec(k3_arr')
    k3_arr = k3_arr[1:(end - Nx + 1)]
    k1_arr = -(coeffFun(x .+ 0.5 * hx, y)) ./ (hx * hx)
    k1_arr[:, end] .= 0
    k1_arr = vec(k1_arr')
    k1_arr = k1_arr[1:(end - 1)]
    A = SparseArrays.spdiagm(0 => main_arr, -1 => k1_arr, 1 => k1_arr, (-Nx+1) => k3_arr, (Nx-1) => k3_arr)
    return A
end

function RHS(ddu, du, u, p, t)
    f = sourcefunc(x, y, t)
    fLX = vec(f')
    ddu .= -A*u .+ fLX
end

# @btime create2DLFVM(Nx, hx, hy, x, y, coeffK)
# @code_warntype create2DLFVM(Nx, hx, hy, x, y, coeffK)
A = create2DLFVM(Nx, hx, hy, x, y, coeffK)
u0 = zeros((Nx-1)*(Ny-1), 1)
du0 = u0
tspan = (tStart, tEnd)

problem = SecondOrderODEProblem(RHS, du0, u0, tspan)
solution = solve(problem, save_everystep=true)

# solutions = solution.u[end]
# duEnd = solutions[1:((Nx-1)*(Ny-1))]
# uEnd = solutions[((Nx-1)*(Ny-1)+1):end]

# duEnd1 = reshape(duEnd, Nx-1, Ny-1)'
# uEnd1 = reshape(uEnd, Nx-1, Ny-1)'
# heatmap(x_points, y_points, uEnd1, xlim=(LeftX, RightX), ylim=(LeftY, RightY), origin=:lower, c=:viridis, xlabel="X", ylabel="Y", title="uEnd Nx=$Nx Ny=$Ny", colorbar=true)

pyplot()
animation = @animate for i in range(1, size(solution.t)[1])
    formatted_t = @sprintf("%.3f", solution.t[i])
    u_plot = solution.u[i][((Nx-1)*(Ny-1)+1):end]
    u_plot = reshape(u_plot, Nx-1, Ny-1)'

    heatmap(x_points, y_points, u_plot, xlim=(LeftX, RightX), ylim=(LeftY, RightY), origin=:lower, c=:viridis, xlabel="X", ylabel="Y", title="uEnd Nx=$Nx Ny=$Ny", colorbar=true)

end

# Save the animation
gif(animation, "test.gif", fps = 15)
