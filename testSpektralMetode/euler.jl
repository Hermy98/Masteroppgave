using Plots, LinearAlgebra

function euler_matrix(N, α, β, r)


    diag = ones(N-1)*(1 -2*α)

    off_diag = ones(N-2)*(α)

    beta_array = ones(N).*(β./r[1:end])

    off_diaglo = off_diag .+= beta_array[2:end-1]

    off_diagup = off_diag .-= beta_array[3:end]

    return Tridiagonal(off_diaglo, diag, off_diagup)



end