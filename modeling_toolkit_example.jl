using ModelingToolkit
using MethodOfLines
using DomainSets
using DifferentialEquations
using SciMLBase

# Define the parameters
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

# Define the variables and functions
@parameters t, x, y
@variables u(..)
Dx = Differential(x); Dy = Differential(y); Dt = Differential(t)
Dxx = Dx^2; Dyy = Dy^2; Dtt = Dt^2

α = 40
ω = 4*π

f(x, y, t) = sin(ω*t)(exp(-α * (x - 0.25 * RightX)^2 - α * (y - 0.25 * RightY)^2) +
exp(-α * (x -0.25 * RightX)^2 - α * (y - 0.75 * RightY)^2) +
exp(-α * (x -0.75 * RightX)^2 - α * (y - 0.75 * RightY)^2) +
exp(-α * (x -0.75 * RightX)^2 - α * (y - 0.25 * RightY)^2)) 

@register_symbolic f(x, y, t)

function coeffK(x, y, t)
    if x<(RightX/2) && y<(RightY/2)
        return 0.1
    elseif x<(RightX/2) && y>=(RightY/2)
        return 0.4
    elseif x>=(RightX/2) && y>=(RightY/2)
        return 0.7
    else
        return 1.0
    end
end

# k(x, y, t) = ...
# @register_symbolic k(x, y, t)

eq = Dtt(u(x, y, t)) - Dxx(u(x, y, t)) - Dyy(u(x, y, t)) ~ f(x, y, t)

domain = [x ∈ Interval(LeftX, RightX),
          y ∈ Interval(LeftY, RightY),
          t ∈ Interval(tStart, tEnd)]

ic_bc = [u(x, y, 0) ~ 0.0,
         Dt(u(x, y, 0)) ~ 0.0,
         u(0.0, y, t) ~ 0.0,
         u(x, 0.0, t) ~ 0.0,
         u(RightX, y, t) ~ 0.0,
         u(x, RightY, t) ~ 0.0]

@named sys = PDESystem(eq, ic_bc, domain, [x, y, t], [u(x, y, t)])

discretization = MOLFiniteDifference([x => hx, y => hy], t, approx_order = 2)

prob = discretize(sys, discretization)

# DifferentialEquations.discretize(sys, discretization)

sol = solve(prob, Euler(), save_everystep = false)

grid = get_discrete(sys, discretization)
