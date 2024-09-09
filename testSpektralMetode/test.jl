using Plots, LinearAlgebra


function get_matrix(N, α, β, r)

    off_diag = -1*ones(N-2)*(α+β)
    diag = ones(N-1)*(1 + α)

    off_diaglo = off_diag
    off_diaglo[1:end] = off_diaglo[1:end]./r[2:N-1]




    return Tridiagonal(off_diaglo, diag, off_diag./r[3:N])


end

function initialcondition(r)

    return r.^2 .+3

end

origoInterpolasjon(u) = (1/3)*(4*u[2] + u[3]) #må tilpasses for funksjoner som ikke er konstante for alle theta

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

    return r, Θ, A, B, d_r


end


function main(N, J, T, d_t, D)

    r, Θ, A, B , d_r = initialilize(N, J, d_t, D)

    u = zeros(N, T)

    u[:, 1] = initialcondition(r)


    A[end, end] *= 5

    B[end, end] *= 5


    for i in 1:T-1
        u[2:end, i+1] = (A\B)*u[2:end, i]

        u[1, i+1] = origoInterpolasjon(u[:, i+1])
    end

    return u, r, Θ
 

end



u , r, Θ = main(10, 10, 10, 1E-7, 1)

print(u)

plot(u[end, :])




