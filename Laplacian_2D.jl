using SparseArrays
using BenchmarkTools

Nx = 4
Ny = 4
LeftX = 0
LeftY = 0
RightX = 5
RightY = 5
hx = (RightX - LeftX)/Nx 
hy = (RightY - LeftY)/Ny 

function Laplacian2D(Nx, Ny, hx, hy)
    hx2 = 1/(hx*hx)
    hy2 = 1/(hy*hy)

    main_len = (Nx-1)*(Ny-1)
    k_outer = ones(main_len-(Nx-1)) .* (-hy2)
    k1 = ones(main_len-1) .* (-hx2)
    k0 = ones(Ny-1, Nx-1) .* (2*(hx2+hy2))
    
    k0[1, 2:end-1] .= hx2 + 2*hy2
    k0[end, 2:end-1] .= hx2 + 2*hy2
    
    k0[2:end-1, 1] .= 2*hx2 + hy2
    k0[2:end-1, end] .= 2*hx2 + hy2
    
    k0[1, 1] = hx2 + hy2
    k0[1, end] = hx2 + hy2
    k0[end, 1] = hx2 + hy2
    k0[end, end] = hx2 + hy2
    
    k0 = vec(k0')
    A = spdiagm(0 => k0, -1 => k1, 1 => k1, (-Nx+1) => k_outer, (Nx-1) => k_outer)
    return A
end

@btime Laplacian2D(Nx, Ny, hx, hy)
