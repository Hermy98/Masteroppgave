using Plots, LinearAlgebra


function get_matrix(N, α, β, r)


    beta_array = β*ones(N)./(2 .*r)

    
    off_diagup = -1*ones(N-2)*(α/2)
    diag = ones(N-1)*(1 + α)
    off_diaglo = -1*ones(N-2)*(α/2)

    off_diaglo .+= beta_array[3:end]

    off_diagup .-= beta_array[2:end-1]


    M = Tridiagonal(off_diaglo, diag, off_diagup)


    return M

end

function initialcondition(N)

    u_1 = zeros(N)

    for i in 1:N
        u_1[i] = 5+((i-1)*(i-100))/2000
    end

    return  u_1

end


function initialcondition2(r)
    return r.^2 .+ 4
end

origoInterpolasjon(u) = (1/3)*(4*u[2] - u[3])

#u_originLaplace(u, d_r) = u[2] - ((d_r^2)/4)*u[1] 



function initialilize(N, J, d_t, D)

    r = range(0, 1, length=N)

    Θ = range(0, 2*π, length=J)

    d_r = r[2] - r[1]

    α = (D*d_t)/(d_r^2)

    β = (D*d_t)/(2*d_r)

    #γ = 4*(d_t^2)/(d_r^2*D^2)

    A = get_matrix(N, α, β, r)

    B = get_matrix(N, -α, -β, r)

    #setting the boundary conditions


    A[1,2] = -α

    B[1,2] = α

    A[end, end] = 1
    A[end, end-1] = 0

    B[end, end] = 1
    B[end, end-1] = 0

    return r, Θ, A, B, d_r


end


function main(N, J, T, d_t, D)

    r, Θ, A, B , d_r = initialilize(N, J, d_t, D)

    u = zeros(N, T)

    u[:, 1] = initialcondition2(r)


    for i in 1:T-1
        u[2:end, i+1] = (A\B)*u[2:end, i]

        u[1, i+1] = origoInterpolasjon(u[:, i+1])
    end

    return u, r, Θ
 

end



u , r, Θ = main(100, 10, 10000, 1E-3, 1)


plot(u[:, 1])
plot!(u[:, 100])




