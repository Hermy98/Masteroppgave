
# function setinitialu_active(N, F, l, m, y, Ly, T, ϕ)
    
#     ux = zeros(2*N+1, m, T)

#     uy = zeros(2*N+1, m, T)
    
#     for i in 1:2*N+1

#         ux[i, :, 1] = F*sin.(((2 * π )/(Ly))*y ) + cos.(ϕ[i, :, 1])

#         uy[i, :, 1] = F*sin.(((2 * π )/(Ly))*y ) + sin.(ϕ[i, :, 1])

#         #ux[i, :, 1] = F*y.^3

#         #uy[i, :, 1] = F*y.^3 .*sin.(((2 * π )/(Ly))*y )

#     end 

#     return ux, uy

# end

function calculateuxuy_active(coefficientmatrix, N, Lx, x, m, ux, uy, ϕ, F_a)

    #uxstep = zeros(2*N+1, m-2)

    #uystep = zeros(2*N+1, m-2)

    for j in 2:m-1

        for i in 0:N 

            ux[:, j] += (coefficientmatrix[j, 1, i+1] * cos.(((2 * π * i)/Lx)*x ) + coefficientmatrix[j, 2, i+1] * sin.(((2 * π * i)/Lx)*x ))


            uy[:, j] += (coefficientmatrix[j, 3, i+1] * cos.(((2 * π * i)/Lx)*x ) + coefficientmatrix[j, 4, i+1] * sin.(((2 * π * i)/Lx)*x ))
        end 
     

    end

    return (1/(2*N+1))*ux + F_a*cos.(ϕ), (1/(2*N+1))*uy + F_a*sin.(ϕ)

end

function initializeelastic(N, F, m, y, Ly, T, ϕ)

    ux, uy = setinitialu(N, F, m, y, Ly, T)

    coefficientmatrix = initcoefficientmatrix(m, N, ux, uy)

    ux[:, :, 1] += cos.(ϕ[:, :, 1])

    uy[:, :, 1] += sin.(ϕ[:, :, 1])

    K1matrix= deepcopy(coefficientmatrix)

    K2matrix = deepcopy(coefficientmatrix)

    K3matrix = deepcopy(coefficientmatrix)

    K4matrix = deepcopy(coefficientmatrix)

    Kmatrix = [K1matrix, K2matrix, K3matrix, K4matrix]

    return ux, uy, coefficientmatrix, Kmatrix


end 

function run_activesystem(N::Int, dy::Float64, Lx::Int, Ly::Int, F::Float64, T::Int, dt::Float64, K::Float64, μ::Float64, λ::Float64, F_a::Float64)

    x, y, ϕ, factors, m = runsystem(N, dy, Lx, Ly, F, T, dt, λ)

    ux, uy, coefficientmatrix, Kmatrix = initializeelastic(N, F, m, y, Ly, T, ϕ)

    for i in 2:T

        coefficientmatrix = rk4(coefficientmatrix, dt, Kmatrix, N, factors, K, μ, y)

        ux[:, :, i], uy[:, :, i] = calculateuxuy_active(coefficientmatrix, N, Lx, x, m, ux[:, :, i], uy[:, :, i], ϕ[:, :, i], F_a)

    end

    return ux, uy, ϕ, x, y

end


