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

    #=
    doublederiv[1, 1] = 1

    doublederiv[1, 2] = 0

    doublederiv[end, end-1] = 0

    doublederiv[end, end] = 1
    =#

    #doublederiv[1, 2] *= 2
    #doublederiv[end, end-1] *= 2



    deriv = dymatrix(dy, m)

    
    #deriv[1, 1] = 1

    #deriv[end, end] = 1
    
    
    #deriv[1, 2] = 0

    #deriv[end, end-1] = 0
    


    doublederivarivesmatrix = zeros(m, 4, N+1)

    derivativematrix = zeros(m, 4, N+1)

    return doublederiv, deriv, doublederivarivesmatrix, derivativematrix

end

function derivatives(coefficientmatrix, doublederiv, deriv, derivativematrix, doublederivarivesmatrix, N)

    for i in 1:4

        for j in 1:N+1

            doublederivarivesmatrix[:, i, j] = doublederiv * coefficientmatrix[:, i, j]

            derivativematrix[:, i, j] = deriv * coefficientmatrix[:, i, j]
        end
    end

    return derivativematrix, doublederivarivesmatrix

end
