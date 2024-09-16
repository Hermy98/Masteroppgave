include("./euler.jl")


function angularMatrix(N, α, β, γ , r)

    Mr = euler_matrix(N, α, β, r)

    diag = ones(N-1)*(-2*γ)

    off_diag = ones(N-2)*(γ)

    Mtheta = Tridiagonal(off_diag, diag, off_diag)

    return Mr + Mtheta

end


function initialcondition(r, θ)

    return r.^2 .* cos.(θ) .+ 4

end

function initialilize(N, J, d_t, D)


    r = range(0, 1, length=N)

    Θ = range(0, 2*π, length=J)

    d_r = r[2] - r[1]

    d_θ = Θ[2] - Θ[1]

    α = (D*d_t)/(d_r^2)

    β = (D*d_t)/(2*d_r)

    γ = (D*d_t)/(d_θ^2)


    A = angularMatrix(N, α, β, γ, r)

    u = zeros(N, J, T)


    return r, Θ, A, u

end

function runmultidim(N, J, T, d_t, D)

    r, Θ, A, u = initialilize(N, J, d_t, D)

    u[: , :, 1] = initialcondition(r, Θ)


    for t in 1:T-1
        u[2:end,2:end,t+1] = A \ u[2:end,2:end,t]
    end





end