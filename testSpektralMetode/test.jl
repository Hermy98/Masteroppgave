using Plots, LinearAlgebra


function get_matrix(N, α, β, r)

    β_array = ones(N).*(β./r[1:end])
    
    off_diag = -1*ones(N-2)*(α)
    diag = ones(N-1)*(1 + α)

    return Tridiagonal(off_diag .+ β_array[2:end-1], diag, off_diag .+β_array[3:end])


end

function initialcondition(r)

    return  r.^2 .+4

end

origoInterpolasjon(u) = (1/3)*(4*u[2] + u[3]) 

#u_originLaplace(u, d_r) = u[2] - ((d_r^2)/4)*u[1] 



function initialilize(N, J, d_t, D)

    r = range(0, 1, length=N)

    Θ = range(0, 2*π, length=J)

    d_r = r[2] - r[1]

    α = (D*d_t)/(d_r^2)

    β = (D*d_t)/(2*d_r)

    println(α, " ", β)

    #γ = 4*(d_t^2)/(d_r^2*D^2)

    A = get_matrix(N, α, β, r)

    B = get_matrix(N, -α, -β, r)

    return r, Θ, A, B, d_r


end


function main(N, J, T, d_t, D)

    r, Θ, A, B , d_r = initialilize(N, J, d_t, D)

    u = zeros(N, T)

    u[:, 1] = initialcondition(r)


    for i in 1:T-1
        u[2:end, i+1] = (A\B)*u[2:end, i]

        u[1, i+1] = origoInterpolasjon(u[:, i+1])
    end

    return u, r, Θ
 

end



u , r, Θ = main(100, 10, 10000, 1E-7, 1)


plot(u[:, 1])
plot!(u[:, 20])




