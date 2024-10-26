using LinearAlgebra, Plots, SparseArrays

function dy2matrix(dy::Float64, m::Int) #l er antal punkter i x-retning, m er antall punkter i y-retning

    diag = ones(m).*((-2)/(dy^2))

    off_diagup = ones(m-1)./(dy^2)
    off_diaglo = ones(m-1)./(dy^2)


    return Tridiagonal(off_diaglo, diag, off_diagup)
end



function dymatrix(dy::Float64, m::Int)

    off_diagup = ones(m-1)./(2*dy)
    off_diaglo = -1. *ones(m-1)./(2*dy)

    return spdiagm(-1 => off_diaglo, 1 => off_diagup)

end


prefactor(L::Int, n::Int) = (2*π*n/L)

