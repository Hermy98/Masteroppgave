using Plots, LinearAlgebra

function euler_matrix(N, α, β, r)


    diag = ones(N-1)*(1 -2*α)

    off_diag = ones(N-2)*(α)

    beta_array = β*ones(N)./r[1:end]


    off_diaglo = off_diag .-β #.- beta_array[2:end-1]

    off_diagup = off_diag  .+ β #.+ beta_array[3:end]

    A = Tridiagonal(off_diaglo, diag, off_diagup)

    A[1, 2] = 2*α

    A[end, end] = 1
    
    A[end, end-1] = 0

    return A



end


function initialcondition(r)

    return  r.^2 .+4

end


origoInterpolasjon(u) = (1/3)*(4*u[2] - u[3])


function initialilize(N, J, d_t, D)


    r = range(0, 1, length=N)

    Θ = range(0, 2*π, length=J)

    d_r = r[2] - r[1]

    α = (D*d_t)/(d_r^2)

    β = (D*d_t)/(2*d_r)


    A = euler_matrix(N, α, β, r)

    return r, Θ, A

end


function run(N, J, T, d_t, D)


    r, Θ, A = initialilize(N, J, d_t, D)

    u = zeros(N, T)

    u[:, 1] = initialcondition(r)


    for i in 1:T-1

        u[2:end, i+1] = A*u[2:end, i]


        u[1, i+1] = origoInterpolasjon(u[:, i+1])

    end


    return u

end

u = run(100, 10, 1000000, 1E-5, 1)

plot(u[:, 1])
plot!(u[:, 500])
