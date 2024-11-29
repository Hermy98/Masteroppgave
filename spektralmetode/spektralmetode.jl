using LinearAlgebra, CairoMakie, SparseArrays, Integrals

function dy2matrix(dy::Float64, m::Int) #l er antal punkter i x-retning, m er antall punkter i y-retning

    diag = ones(m).*((-2)/(dy^2))

    off_diagup = ones(m-1)./(dy^2)
    off_diaglo = ones(m-1)./(dy^2)


    return Tridiagonal(off_diaglo, diag, off_diagup)
end



function dymatrix(dy::Float64, m::Int)

    off_diagup = ones(m-1)./(2*dy)
    off_diaglo = -1. *ones(m-1)./(2*dy)

    print(off_diaglo)
    print(off_diagup)

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

function initcoefficientmatrix(m::Int, N::Int, Lx::Int, ux::Array{Float64, 3}, uy::Array{Float64, 3}, x, l)

    coefficientmatrix = zeros(m, 4,  N)

    yint = ones(4, l)


    for i in 1:N

        for j in 1:m

            yint[1, :] = ux[:, j, 1].*cos.(((2 * π * i)/Lx)*x )

            yint[2, :] = ux[:, j, 1].*sin.(((2 * π * i)/Lx)*x )


            yint[3, :] = uy[:, j, 1].*cos.(((2 * π * i)/Lx)*x )

            yint[4, :] = uy[:, j, 1].*sin.(((2 * π * i)/Lx)*x )
            

            method = SimpsonsRule()

            problem = SampledIntegralProblem(yint, x)

            coefficientmatrix[j, :, i] = solve(problem, method)

        end
    end


    return coefficientmatrix

end


prefactor(L::Int, n::Int) = (2*π * n)/L


function prefactors(Lx::Int,N::Int)
    
    return [prefactor(Lx, n) for n in 1:N ]

end


function initialzsystem(dy::Float64, l::Int, N::Int, Lx::Int, Ly::Int, F::Float64, T::Int)

    x = LinRange(-Lx/2, Lx/2, l)

    y = -Ly/2:dy:Ly/2

    m = length(y)

    ux, uy = setinitialu(N, F, l, m, y, Ly, T)

    factors = prefactors(Lx, N)

    coefficientmatrix = initcoefficientmatrix(m, N, Lx, ux, uy, x, l)

    K1matrix= deepcopy(coefficientmatrix)

    K2matrix = deepcopy(coefficientmatrix)

    return x, y, m, ux, uy, factors, coefficientmatrix, K1matrix, K2matrix

end


function derivatives(coefficientmatrix, doublederiv, deriv, derivativematrix, doublederivarivesmatrix, N)

    for i in 1:4

        for j in 1:N

            doublederivarivesmatrix[:, i, j] = doublederiv * coefficientmatrix[:, i, j]

            derivativematrix[:, i, j] = deriv * coefficientmatrix[:, i, j]
        end
    end

    return derivativematrix, doublederivarivesmatrix

end




function timederivatives(coefficientmatrix, Kmatrix,  factor, K, μ ,dublederiv, deriv, N)

    for i in 1:N


        Kmatrix[2:end-1, 1, i] = (K + (3/2)* μ) * (dublederiv[2:end-1, 1, i] - (factor[i])^2 * coefficientmatrix[2:end-1, 1, i]) + (K + 1/2*μ) * factor[i] * deriv[2:end-1, 4, i]


        Kmatrix[2:end-1, 2, i] = (K + (3/2)* μ) * (dublederiv[2:end-1, 2, i] - (factor[i])^2 * coefficientmatrix[2:end-1, 2, i]) - (K + 1/2*μ) * factor[i] * deriv[2:end-1, 3, i]


        Kmatrix[2:end-1, 3, i] = (K + (3/2)* μ) * (dublederiv[2:end-1, 3, i] - (factor[i])^2 * coefficientmatrix[2:end-1, 3, i]) + (K + 1/2*μ) * factor[i] * deriv[2:end-1, 2, i]


        Kmatrix[2:end-1, 4, i] = (K + (3/2)* μ) * (dublederiv[2:end-1, 4, i] - (factor[i])^2 * coefficientmatrix[2:end-1, 4, i]) - (K + 1/2*μ) * factor[i] * deriv[2:end-1, 1, i]

    end

    return Kmatrix

end

function calculateuxuy(coefficientmatrix, N, Lx, x, l, m)

    uxstep = zeros(l, m-2)

    uystep = zeros(l, m-2)

    for j in 1:m-2

        for i in 1:N #indekseringsfeil per nå

            uxstep[:, j] += coefficientmatrix[j+1, 1, i] * cos.(((2 * π * i)/Lx)*x ) + coefficientmatrix[j+1, 2, i] * sin.(((2 * π * i)/Lx)*x )


            uystep[:, j] += coefficientmatrix[j+1, 3, i] * cos.(((2 * π * i)/Lx)*x ) + coefficientmatrix[j+1, 4, i] * sin.(((2 * π * i)/Lx)*x )
        end 
     

    end

    #print(uxstep, "\n")

    return uxstep, uystep

end

function setinitialu(N, F, l, m, y, Ly, T)
    
    ux = zeros(l, m, T)

    uy = zeros(l, m, T)
    
    for i in 1:l
        ux[i, :, 1] = F*cos.(((2 * π )/(Ly))*y )

        uy[i, :, 1] = F*sin.(((2 * π )/(Ly))*y )

    end

    ux[:, 1, :] .= ux[:, 1, 1] 

    uy[:, 1, :] .= uy[:, 1, 1]

    ux[:, end, :] .= ux[:, end, 1]

    uy[:, end, :] .= uy[:, end, 1]

    return ux, uy

end

function rk2(coefficientmatrix, dt, K1matrix, K2matrix,  N, factor, K, μ, doublederiv, deriv, doublederivarivesmatrix, derivativematrix)

    derivativematrix, doublederivarivesmatrix = derivatives(coefficientmatrix, doublederiv, deriv, derivativematrix, doublederivarivesmatrix, N)


    K1matrix = timederivatives(coefficientmatrix, K1matrix,  factor, K, μ ,doublederivarivesmatrix, derivativematrix, N)


    derivativematrix, doublederivarivesmatrix = derivatives(coefficientmatrix + dt/2 * K1matrix , doublederiv, deriv, derivativematrix, doublederivarivesmatrix, N)


    K2matrix = timederivatives(coefficientmatrix + dt/2 * K1matrix, K2matrix,  factor, K, μ ,doublederivarivesmatrix, derivativematrix, N)


    coefficientmatrix += dt * K2matrix


    return coefficientmatrix


end

function sumulation(N::Int, dy::Float64, l::Int, Lx::Int, Ly::Int, F::Float64, T::Int, dt::Float64, K::Float64, μ::Float64)

    x, y, m, ux, uy, factors, coefficientmatrix, K1matrix, K2matrix = initialzsystem(dy, l, N, Lx, Ly, F, T)

    doublederiv, deriv, doublederivarivesmatrix, derivativematrix = initializederivatives(m, dy, N)

    for i in 2:T

        coefficientmatrix = rk2(coefficientmatrix, dt, K1matrix, K2matrix, N, factors, K, μ, doublederiv, deriv, doublederivarivesmatrix, derivativematrix)


        ux[:, 2:end-1, i], uy[:, 2:end-1, i] = calculateuxuy(coefficientmatrix, N, Lx, x, l, m)

    end

    return ux, uy, x, y

end

ux, uy, x, y = sumulation(20, 0.025, 10, 1, 2, 2., 120, 0.00001, 5., 1.)

f = Figure(size = (800, 800))

Axis(f[1, 1])

arrows(x, y, ux[:, :, end], uy[:, :, end], arrowsize = 10, lengthscale = 0.1)

