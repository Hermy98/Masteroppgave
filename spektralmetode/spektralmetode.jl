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

    doublederiv[1, 2] *= 2
    doublederiv[end, end-1] *= 2



    deriv = dymatrix(dy, m)

    
    #deriv[1, 1] = 1

    #deriv[end, end] = 1
    
    
    deriv[1, 2] = 0

    deriv[end, end-1] = 0
    


    doublederivarivesmatrix = zeros(m, 4, N)

    derivativematrix = zeros(m, 4, N)

    return doublederiv, deriv, doublederivarivesmatrix, derivativematrix

end

function initcoefficientmatrix(m::Int, N::Int, Lx::Int, ux::Array{Float64, 3}, uy::Array{Float64, 3}, x, l)

    coefficientmatrix = zeros(m, 4,  N)

    yint = ones(4, l)


    for i in 1:N

        for j in 1:m

            yint[1, :] = (1/Lx) * ux[:, j, 1].*cos.(((2 * π * i)/Lx)*x )

            yint[2, :] = (1/Lx) * ux[:, j, 1].*sin.(((2 * π * i)/Lx)*x )


            yint[3, :] = (1/Lx) * uy[:, j, 1].*cos.(((2 * π * i)/Lx)*x )

            yint[4, :] = (1/Lx) * uy[:, j, 1].*sin.(((2 * π * i)/Lx)*x )
            

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

    savecoefficients = zeros(m, 4, N, T)

    K1matrix= deepcopy(coefficientmatrix)

    K2matrix = deepcopy(coefficientmatrix)

    K3matrix = deepcopy(coefficientmatrix)

    K4matrix = deepcopy(coefficientmatrix)

    return x, y, m, ux, uy, factors, coefficientmatrix, K1matrix, K2matrix, K3matrix, K4matrix, savecoefficients

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


       # Kmatrix[:, 1,  i] = (K + (3/2)* μ) * (dublederiv[:, 1, i] - (factor[i])^2 * coefficientmatrix[:, 1, i]) + (K + 1/2*μ) * factor[i] * deriv[:, 4, i]

        Kmatrix[:, 1,  i] =  μ* dublederiv[:, 1, i] - (K + 3/2*μ) * (factor[i])^2 * coefficientmatrix[:, 1, i] + (K + 1/2*μ) * factor[i] * deriv[:, 4, i]

        #Kmatrix[:, 2, i] = (K + (3/2)* μ) * (dublederiv[:, 2, i] - (factor[i])^2 * coefficientmatrix[:, 2, i]) - (K + 1/2*μ) * factor[i] * deriv[:, 3, i]
        Kmatrix[:, 2, i] = μ * dublederiv[:, 2, i] - (K + 3/2*μ) * (factor[i])^2 * coefficientmatrix[:, 2, i]  - (K + 1/2*μ) * factor[i] * deriv[:, 3, i]

        #Kmatrix[:, 3, i] = (K + (3/2)* μ) * (dublederiv[:, 3, i] - (factor[i])^2 * coefficientmatrix[:, 3, i]) + (K + 1/2*μ) * factor[i] * deriv[:, 2, i]

        Kmatrix[:, 3, i] = (K + 3/2 *μ) * dublederiv[:, 3, i] - μ * (factor[i])^2 * coefficientmatrix[:, 3, i] + (K + 1/2*μ) * factor[i] * deriv[:, 2, i]
        #Kmatrix[:, 4, i] = (K + (3/2)* μ) * (dublederiv[:, 4, i] - (factor[i])^2 * coefficientmatrix[:, 4, i]) - (K + 1/2*μ) * factor[i] * deriv[:, 1, i]

        Kmatrix[:, 4, i] = (K + (3/2)* μ) * dublederiv[:, 4, i] - μ * (factor[i])^2 * coefficientmatrix[:, 4, i] - (K + 1/2*μ) * factor[i] * deriv[:, 1, i]
    end

    return Kmatrix

end

function calculateuxuy(coefficientmatrix, N, Lx, x, l, m)

    uxstep = zeros(l, m)

    uystep = zeros(l, m)

    for j in 1:m

        for i in 1:N #indekseringsfeil per nå

            uxstep[:, j] += coefficientmatrix[j, 1, i] * cos.(((2 * π * i)/Lx)*x ) + coefficientmatrix[j, 2, i] * sin.(((2 * π * i)/Lx)*x )


            uystep[:, j] += coefficientmatrix[j, 3, i] * cos.(((2 * π * i)/Lx)*x ) + coefficientmatrix[j, 4, i] * sin.(((2 * π * i)/Lx)*x )
        end 
     

    end


    return uxstep, uystep

end

function setinitialu(N, F, l, m, y, Ly, T)
    
    ux = zeros(l, m, T)

    uy = zeros(l, m, T)
    
    for i in 1:l
        ux[i, :, 1] = F*cos.(((2 * π )/(Ly))*y )

        uy[i, :, 1] = F*sin.(((2 * π )/(Ly))*y )

        #ux[i, :, 1] = F*y.^3

        #uy[i, :, 1] = F*y.^3 .*sin.(((2 * π )/(Ly))*y )

    end 

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


function rk4(coefficientmatrix, dt, K1matrix, K2matrix, K3matrix, K4matrix, N, factor, K, μ, doublederiv, deriv, doublederivarivesmatrix, derivativematrix)

    derivativematrix, doublederivarivesmatrix = derivatives(coefficientmatrix, doublederiv, deriv, derivativematrix, doublederivarivesmatrix, N)


    K1matrix = timederivatives(coefficientmatrix, K1matrix,  factor, K, μ ,doublederivarivesmatrix, derivativematrix, N)


    derivativematrix, doublederivarivesmatrix = derivatives(coefficientmatrix + dt/2 * K1matrix , doublederiv, deriv, derivativematrix, doublederivarivesmatrix, N)


    K2matrix = timederivatives(coefficientmatrix + dt/2 * K1matrix, K2matrix,  factor, K, μ ,doublederivarivesmatrix, derivativematrix, N)


    derivativematrix, doublederivarivesmatrix = derivatives(coefficientmatrix + dt/2 * K2matrix , doublederiv, deriv, derivativematrix, doublederivarivesmatrix, N)


    K3matrix = timederivatives(coefficientmatrix + dt/2 * K2matrix, K3matrix,  factor, K, μ ,doublederivarivesmatrix, derivativematrix, N)


    derivativematrix, doublederivarivesmatrix = derivatives(coefficientmatrix + dt * K3matrix , doublederiv, deriv, derivativematrix, doublederivarivesmatrix, N)


    K4matrix = timederivatives(coefficientmatrix + dt * K3matrix, K4matrix,  factor, K, μ ,doublederivarivesmatrix, derivativematrix, N)


    coefficientmatrix += dt/6 * (K1matrix + 2 * K2matrix + 2 * K3matrix + K4matrix)


    return coefficientmatrix

end

function sumulation(N::Int, dy::Float64, l::Int, Lx::Int, Ly::Int, F::Float64, T::Int, dt::Float64, K::Float64, μ::Float64)

    x, y, m, ux, uy, factors, coefficientmatrix, K1matrix, K2matrix, K3matrix, K4matrix, savecoefficients = initialzsystem(dy, l, N, Lx, Ly, F, T)

    doublederiv, deriv, doublederivarivesmatrix, derivativematrix = initializederivatives(m, dy, N)

    savecoefficients[:, :, :, 1] = coefficientmatrix

    for i in 2:T

        #coefficientmatrix = rk2(coefficientmatrix, dt, K1matrix, K2matrix, N, factors, K, μ, doublederiv, deriv, doublederivarivesmatrix, derivativematrix)


        coefficientmatrix = rk4(coefficientmatrix, dt, K1matrix, K2matrix, K3matrix, K4matrix, N, factors, K, μ, doublederiv, deriv, doublederivarivesmatrix, derivativematrix)

        ux[:, :, i], uy[:, :, i] = calculateuxuy(coefficientmatrix, N, Lx, x, l, m)

        savecoefficients[:, :, :, i] = coefficientmatrix

    end

    return ux, uy, x, y, savecoefficients

end

function sum(savecoefficients, N, T, m)

    coeffient = zeros(m, 4, T)

    for j in 1:T

        for i in 1:N 

            coeffient[:, :, j] += savecoefficients[:, :, i, j]

        end 

        
    end

    return coeffient

end

ux, uy, x, y, savecoefficients = sumulation(20, 0.025, 10, 1, 2, 2., 100, 0.00001, 5., 1.)

f = Figure(size = (800, 800))

Axis(f[1, 1])

arrows(x, y, ux[:, :, 13], uy[:, :, 13], arrowsize = 10, lengthscale = 0.1)

