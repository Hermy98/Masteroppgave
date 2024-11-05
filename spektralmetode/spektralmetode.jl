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

function initializederivatives(m::Int, dy::Float64, N::Int)

    doublederiv = dy2matrix(dy, m)

    doublederiv[1,2] = 0

    doublederiv[end, end-1] = 0

    deriv = dymatrix(dy, m)

    deriv[1,2] = 0
    deriv[1, 1] = 1

    deriv[end, end-1] = 0
    deriv[end, end] = 1

    doublederivarivesmatrix = zeros(m, 4, N)

    derivativematrix = zeros(m, 4, N)

    return doublederiv, deriv, doublederivarivesmatrix, derivativematrix

end

function coefficientmatrix(m::Int, N::Int, F::Float64)

    coefficientmatrix = zeros(m, 4,  N)

    coefficientmatrix[1, :, :] .= -F/2

    coefficientmatrix[end, :, :] .= F/2

    for i in 2:m-1

        coefficientmatrix[i, : , :] .= sin((y[i]*π)/(Ly*2))
          
    end

    return coefficientmatrix

end


prefactor(L::Int, n::Int) = (2*π * n)/L


function prefactors(Lx::Int,N::Int)
    
    return [prefactor(Lx, n) for n in 1:N ]

end


function initialzsystem(m::Int, l::Int, dy::Float64, N::Int, Lx::Int, Ly::Int, F::Float64)

    x = LinRange(-Lx/2, Lx/2, l)

    y = LinRange(-Ly/2, Ly/2, m)

    ux = zeros(l, m)

    uy = zeros(l, m)

    factors = prefactors(Lx, N)

    doublederiv, deriv, doublederivarivesmatrix, derivativematrix = initializederivatives(m, dy, N)

    coefficientmatrix = initializecoefficientmatrix(m, N, F)

    Kmatrix = deepcopy(coefficientmatrix)

    return x, y, ux, uy, factors, coefficientmatrix, doublederiv, deriv, doublederivarivesmatrix, derivativematrix, Kmatrix

end


function derivatives(coefficientmatrix, doublederiv, deriv, derivativematrix, doublederivarivesmatrix, N)

    for i in 1:4
        for j in 1:N

            doublederivativesmatrix[:, i, j] = doublederiv * coefficientmatrix[:, i, j]

            derivativesmatrix[:, i, j] = deriv * derivativematrix[:, i, j]
        end
    end

    return derivativematrix, doublederivarivesmatrix

end




function timederivatives(coefficientmatrix, Kmatrix,  factor, K, μ ,dublederiv, deriv, N)


    for i in 1:N

        Kmatrix[2:end-1, 1, i] = (K + (3/2)* μ) * dublederiv * coefficientmatrix[2:end-1, 1, i] + (K +(3/2)*μ) * (-factor[i]) * coefficientmatrix[2:end-1, 1, i] + (K + 1/2*μ) * factor[i] * deriv * coefficientmatrix[2:end-1, 4,i]


        Kmatrix[2:end-1, 2, i] = (K + (3/2)* μ) * dublederiv * coefficientmatrix[2:end-1, 2 ,i] + (K +(3/2)*μ) * (-factor[i]) *coefficientmatrix[2:end-1, 2, i] - (K + 1/2*μ) * factot[i] * deriv * c[2:end-1, 3,i]


        Kmatrix[2:end-1, 3, i] = (K + (3/2)* μ) * dublederiv * c[2:end-1, i] + (K +(3/2)*μ) * (-factor) *c[2:end-1, i] + (K + 1/2*μ) * factor[i] * deriv * b[2:end-1, i]


        Kmatrix[2:end-1, 4, i] = (K + (3/2)* μ) * dublederiv* d[2:end-1, i] + (K +(3/2)*μ) * (-factor) *d[2:end-1, i] - (K + 1/2*μ) *  factro[i] *deriv * a[2:end-1, i]

    end

    return a_dot, b_dot, c_dot, d_dot

end
