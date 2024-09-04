

function get_matrix(α, β)

    off_diag = -1*ones(N-1)*(α-β)
    diag = ones(N)*(1-β)

    return Tridiagonal(off_diag.*r[1:N-1], diag, off_diag.*r[2:N])


end








