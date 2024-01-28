using LinearAlgebra, Plots, LeastSquaresOptim, Optim, SparseArrays


function getBase(H,w)
    array = []

    for i in range(0,H-1)
        function sinus(t)
            return sin.((H-i)*w*t)
        end
        push!(array,sinus)
    end

    function ones(t)
        return 0 .* t .+ 1
    end
    
    push!(array,ones)

    for i in range(0,H-1)
        function cosinus(t)
            return cos.((i+1)*w*t)            
        end
        push!(array,cosinus)
    end
    return array
end

function getBNum(t, w, H)
    Bnum = zeros(2 * H + 1, length(t))
    for i in range(1, 2 * H + 1)
        B = getBase(H,w)
        Bnum[i,:] = B[i](t)
    end
    return Bnum
end

function getPTE(H)
    P = zeros(2*H+1, 2*H+1)
    P[H+1, H+1] = 1
   
    halfDiag = zeros(2*H+1, 2*H+1)
    for i in range(1,H)
        halfDiag[i, i] = 1
    end
    
    function rotateMatrix(matrix)
        return reverse(transpose(matrix), dims=1)     
    end

    P = P - halfDiag * 1im + rotateMatrix(halfDiag) + rotateMatrix(rotateMatrix(halfDiag)) + rotateMatrix(rotateMatrix(rotateMatrix(halfDiag))) * 1im
    return P
end

function diffMat(H, w, n=1)
    vector = []
    for i in range(-H, H)
        push!(vector,i)
    end

    vector = [convert(Int64, x) for x in vector]

    derivativeMatrix = diagm(vector) * w * 1im


    for i in range(0,n-2)
        derivativeMatrix = derivativeMatrix * derivativeMatrix
    end
    return derivativeMatrix
end

function getF(Nx, trunc, h_0, r, g)
    mat = zeros((2 * trunc + 1, Nx))
    mat[trunc, :] += π * (g - 1 / h_0) * cos.(π * x_lin)
    mat[trunc + 2, :] += r / (h_0 * h_0) * cos.(π * x_lin)
    vec = reshape(mat,(Nx * (2 * trunc + 1)))
    return vec
end

function getDifdxMat(h_x, Nx)
    result = -diagm(2 => ones(Nx-2)) + 8 * diagm(1 => ones(Nx-1)) - 8 * diagm(-1 => ones(Nx-1)) + diagm(-2 => ones(Nx-2))
    result[Nx-1, 1] = -1
    result[Nx, 1] = 3
    result[1, 1:5] .= [-25, 48, -36, 16, -3]
    result[2, 1:5] .= [-3, -10, 18, -6, 1]
    result[Nx, Nx-3:Nx] .= [-1, 6, -18, 10]
    result = result / (12 * h_x)
    return sparse(result)
end

function solveHB(r, g, L, h_0, omega_HB, Nx, Nt, trunc)

    function getA(g, r, D_t, D_x, H_arr, H_arr_inv)
        j_0_0 = D_t
        j_0_1 = D_x * H_arr
        j_1_0 = g * D_x
        j_1_1 = D_t + r * H_arr_inv
        return sparse([j_0_0 j_0_1; j_1_0 j_1_1])
    end

    dx = L/Nx

    h_arr = diagm(0 => ones(Nx-0)) * h_0
    h_arr = sparse(h_arr)
    h_arr_inv = diagm(0 => ones(Nx-0)) / h_0
    h_arr_inv = sparse(h_arr_inv)


    P_TE_small = getPTE(trunc)
    P_ET_small = inv(P_TE_small)
    P_TE = kron(sparse(diagm(0 => ones(Nx-0))), P_TE_small)
    P_ET = kron(sparse(diagm(0 => ones(Nx-0))), P_ET_small)
    P_TE_w = kron(diagm(0 => ones(2-0)), P_TE)
    P_ET_w = kron(diagm(0 => ones(2-0)), P_ET)

    println("checkmark 1")
    D_1 = real(P_TE_small * diffMat(trunc, omega_HB, 1) * P_ET_small)

    F_vec = getF(Nx, trunc, h_0, r, g)
    diffdx_mat = getDifdxMat(dx, Nx)

    println("checkmark 2")

    D_t = kron(sparse(diagm(0 => ones(Nx-0))), D_1)
    H_arr = kron(diagm(0 => ones(2 * trunc + 1)), h_arr)
    H_arr_inv = kron(diagm(0 => ones(2 * trunc + 1)), h_arr_inv)
    D_x = kron(diffdx_mat, sparse(diagm(0 => ones(2 * trunc + 1))))

    println("checkmark 3")
    A = getA(g, r, D_t, D_x, H_arr, H_arr_inv)

    A[trunc+1, :] .= 0
    A[trunc+1, trunc+1] = 1

    println("checkmark 4")


    result = A \ vcat(zeros(length(F_vec)), F_vec)


    t_lin = range(0, stop=2 * π / omega_HB, length=Nt)
    B_num = getBNum(t_lin, omega_HB, trunc)

    dzeta_coeff = reshape(result[1:Nx * (2 * trunc + 1)], (2 * trunc + 1, Nx))
    u_coeff = reshape(result[Nx * (2 * trunc + 1)+ 1:end], (2 * trunc + 1, Nx))

    dzeta_num = transpose(dzeta_coeff) * B_num
    u_num = transpose(u_coeff) * B_num

    return dzeta_num, u_num
end

paramaters = Dict("r"=> 0.5, "g"=> 9.81, "L"=> 6, "h_0"=> 0.2, "omega_HB"=> π, "Nx"=> 10000, "Nt"=> 1000, "trunc"=> 10)

dx = paramaters["L"] / paramaters["Nx"]
t_lin = range(0, stop=2, length=paramaters["Nt"])
x_lin = range(0, paramaters["L"] - (paramaters["L"] / paramaters["Nx"]), paramaters["Nx"])


start = time()
dzeta_num, u_num = solveHB(paramaters["r"], paramaters["g"], paramaters["L"], paramaters["h_0"], paramaters["omega_HB"], paramaters["Nx"], paramaters["Nt"], paramaters["trunc"]) 
println("time spent optimizing: ", time() - start, " sec")


errors = []
for ii in 1:length(t_lin)
    y_1 = (1 / paramaters["h_0"]) * cos.(π .* t_lin[ii]) * cos.(π * x_lin)
    y_2 = u_num[:, ii]

    # For zeta:
    # y_1 = sin.(π .* x_lin) * sin.(π .* t_lin[ii])
    # y_2 = dzeta

    push!(errors, sum((y_2 .- y_1) .^ 2))
end

i_max = argmax(errors)


# Define the figure with subplots
fig = plot(layout=(1,3), size=(900, 300))

# Plot error over time
plot!(fig[1], t_lin, errors, xlabel="t", ylabel="error (sum of squared differences)")


# Calculate and plot squared differences
y_1 = (1 / paramaters["h_0"]) .* cos.(π .* x_lin) .* cos.(π .* t_lin[i_max])
y_2 = u_num[:, i_max]

# Plot error over time
plot_errors = plot(t_lin, errors, xlabel="t", ylabel="error (sum of squared differences)")

# Plot squared differences
plot_squared_diff = plot(x_lin, (y_2 .- y_1) .^ 2, xlabel="x", ylabel="error (sum of squared differences)")

# Plot found and reference solutions
plot_solutions = plot(x_lin, y_2, label="found solution")
plot!(x_lin, y_1, ls=:dash, label="reference solution", xlabel="x", ylabel="ζ(x)")
plot!(legend=:topleft)

# Display each plot separately
# display(plot_errors)
# display(plot_squared_diff)
display(plot_solutions) 
savefig("/Users/Kai/Downloads/sol.pdf")
