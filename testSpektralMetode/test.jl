

function get_matrix(α, β, r)

    off_diag = -1*ones(N-1)*(α+β)
    diag = ones(N)*(1 + α)

    return Tridiagonal(off_diag.*r[1:N-1], diag, off_diag.*r[2:N])


end

function initialcondition(r)

    return r.^2 .+5 

end


function initialilize(N, J, T, d_t, D)

    r = range(0, 1, length=N)

    Θ = range(0, 2*π, length=J)

    d_r = r[2] - r[1]

    t = StepRange(0, d_t, T)

    α = (D*d_t)/(d_r^2)

    β = (D*d_t)/(2*d_r)

    γ = 4*(d_t^2)/(d_r^2*D^2)

    A = get_matrix(α, β, r)

    B = get_matrix(-α, -β, r)

    return r, Θ, t, A, B, d_t, d_r, α, β, γ


end


function main(N, J, T, d_t, D)

    r, Θ, t, A, B, d_t, d_r, α, β, γ = initialilize(N, J, T, d_t, D)

    u0 = initialcondition(r)

    A[1, 1] = 1
    A[1, 2] = 0
    A[end, end] = 5

    B[1, 1] = 2 + 4*γ(u0[2]- 1)
    B[1, 2] = 0

    
 

end



