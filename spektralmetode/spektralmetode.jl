using LinearAlgebra, CairoMakie, SparseArrays, Integrals, FFTW, Dierckx


function initerpolatderivate(coefficientmatrix, y , n)


    spilnes = [Spline1D(y, coefficientmatrix[:, i, n]) for i in 1:4]


    derivativ = [derivative(spilne, y) for spilne in spilnes]

    doublederivativ = [derivative(spilne, y, nu = 2) for spilne in spilnes]

    return  derivativ[1], derivativ[2], derivativ[3], derivativ[4], doublederivativ[1], doublederivativ[2], doublederivativ[3], doublederivativ[4]

end

function initcoefficientmatrix(m::Int, N::Int, Lx::Int, ux::Array{Float64, 3}, uy::Array{Float64, 3}, x, l)

    coefficientmatrix = zeros(m, 4,  N+1)

    for j in 1:m

        ft_x = fft(ux[:, j, 1])

        ft_y = fft(uy[:, j, 1])

        coefficientmatrix[j, 1, :] = real.(ft_x)[1:N+1]

        coefficientmatrix[j, 2, :] = imag.(ft_x)[1:N+1]

        coefficientmatrix[j, 3, :] = real.(ft_y)[1:N+1]

        coefficientmatrix[j, 4, :] = imag.(ft_y)[1:N+1]


    end

    return coefficientmatrix

end


prefactor(L::Int, n::Int) = (2*π * n)/L


function prefactors(Lx::Int,N::Int)
    
    return [prefactor(Lx, n) for n in 0:N+1 ]

end


function initialzsystem(dy::Float64, l::Int, N::Int, Lx::Int, Ly::Int, F::Float64, T::Int)

    x = LinRange(-Lx/2, Lx/2, 2*N+1)

    y = -Ly/2:dy:Ly/2

    m = length(y)

    ux, uy = setinitialu(N, F, l, m, y, Ly, T)

    factors = prefactors(Lx, N+1)

    coefficientmatrix = initcoefficientmatrix(m, N, Lx, ux, uy, x, l)

    savecoefficients = zeros(m, 4, N+1, T)

    K1matrix= deepcopy(coefficientmatrix)

    K2matrix = deepcopy(coefficientmatrix)

    K3matrix = deepcopy(coefficientmatrix)

    K4matrix = deepcopy(coefficientmatrix)

    return x, y, m, ux, uy, factors, coefficientmatrix, K1matrix, K2matrix, K3matrix, K4matrix, savecoefficients

end



function timederivatives(coefficientmatrix, Kmatrix,  factor, K, μ , N, y)

    for i in 1:N+1

        an_deriv, bn_deriv, cn_deriv, dn_deriv, an_doublederiv, bn_doublederiv, cn_doublederiv, dn_doublederiv = initerpolatderivate(coefficientmatrix, y, i)


       #Kmatrix[:, 1,  i] = (K + (3/2)* μ) * (dublederiv[:, 1, i] - (factor[i])^2 * coefficientmatrix[:, 1, i]) + (K + 1/2*μ) * factor[i] * deriv[:, 4, i]

        Kmatrix[2:end-1, 1,  i] =  μ * an_doublederiv[2:end-1] - (K + μ) * (factor[i])^2 * coefficientmatrix[2:end-1, 1, i] + K  * factor[i] * dn_deriv[2:end-1]

        #Kmatrix[:, 2, i] = (K + (3/2)* μ) * (dublederiv[:, 2, i] - (factor[i])^2 * coefficientmatrix[:, 2, i]) - (K + 1/2*μ) * factor[i] * deriv[:, 3, i]
        Kmatrix[2:end-1, 2, i] = μ * bn_doublederiv[2:end-1] - (K + μ) * (factor[i])^2 * coefficientmatrix[2:end-1, 2, i]  - K  * factor[i] * cn_deriv[2:end-1]

        #Kmatrix[:, 3, i] = (K + (3/2)* μ) * (dublederiv[:, 3, i] - (factor[i])^2 * coefficientmatrix[:, 3, i]) + (K + 1/2*μ) * factor[i] * deriv[:, 2, i]

        Kmatrix[2:end-1, 3, i] = (K +  μ) * cn_doublederiv[2:end-1] - μ * (factor[i])^2 * coefficientmatrix[2:end-1, 3, i] + K * factor[i] * bn_deriv[2:end-1]
        #Kmatrix[:, 4, i] = (K + (3/2)* μ) * (dublederiv[:, 4, i] - (factor[i])^2 * coefficientmatrix[:, 4, i]) - (K + 1/2*μ) * factor[i] * deriv[:, 1, i]

        Kmatrix[2:end-1, 4, i] = (K + μ)* dn_doublederiv[2:end-1] - μ * (factor[i])^2 * coefficientmatrix[2:end-1, 4, i] - K * factor[i] * an_deriv[2:end-1] 
    end

    return Kmatrix

end

function calculateuxuy(coefficientmatrix, N, Lx, x, l, m)



    uxstep = zeros(2*N+1, m-2)

    uystep = zeros(2*N+1, m-2)

    for j in 2:m-1

        for i in 0:N 

            uxstep[:, j-1] += (coefficientmatrix[j, 1, i+1] * cos.(((2 * π * i)/Lx)*x ) - coefficientmatrix[j, 2, i+1] * sin.(((2 * π * i)/Lx)*x ))


            uystep[:, j-1] += (coefficientmatrix[j, 3, i+1] * cos.(((2 * π * i)/Lx)*x ) - coefficientmatrix[j, 4, i+1] * sin.(((2 * π * i)/Lx)*x ))
        end 
     

    end


    return (1/(2*N+1))*uxstep, (1/(2*N+1))*uystep

end

function setinitialu(N, F, l, m, y, Ly, T)
    
    ux = zeros(2*N+1, m, T)

    uy = zeros(2*N+1, m, T)
    
    for i in 1:2*N+1
        ux[i, :, 1] = F*sin.(((2 * π )/(Ly))*y )

        uy[i, :, 1] = F*sin.(((2 * π )/(Ly))*y )

        #ux[i, :, 1] = F*y.^3

        #uy[i, :, 1] = F*y.^3 .*sin.(((2 * π )/(Ly))*y )

    end 

    return ux, uy

end

function rk2(coefficientmatrix, dt, K1matrix, K2matrix,  N, factor, K, μ, y)



    K1matrix = timederivatives(coefficientmatrix, K1matrix,  factor, K, μ , N, y)


    K2matrix = timederivatives(coefficientmatrix + dt/2 * K1matrix, K2matrix,  factor, K, μ , N, y)


    coefficientmatrix += dt * K2matrix


    return coefficientmatrix


end


function rk4(coefficientmatrix, dt, K1matrix, K2matrix, K3matrix, K4matrix, N, factor, K, μ, y)


    K1matrix = timederivatives(coefficientmatrix, K1matrix,  factor, K, μ , N, y)


    K2matrix = timederivatives(coefficientmatrix + dt/2 * K1matrix, K2matrix,  factor, K, μ , N, y)


    K3matrix = timederivatives(coefficientmatrix + dt/2 * K2matrix, K3matrix,  factor, K, μ , N, y)


    K4matrix = timederivatives(coefficientmatrix + dt * K3matrix, K4matrix,  factor, K, μ , N, y)


    coefficientmatrix += dt/6 * (K1matrix + 2 * K2matrix + 2 * K3matrix + K4matrix)


    return coefficientmatrix

end

function sumulation(N::Int, dy::Float64, l::Int, Lx::Int, Ly::Int, F::Float64, T::Int, dt::Float64, K::Float64, μ::Float64)

    x, y, m, ux, uy, factors, coefficientmatrix, K1matrix, K2matrix, K3matrix, K4matrix, savecoefficients = initialzsystem(dy, l, N, Lx, Ly, F, T)

    savecoefficients[:, :, :, 1] = coefficientmatrix

    for i in 2:T

        #coefficientmatrix = rk2(coefficientmatrix, dt, K1matrix, K2matrix, N, factors, K, μ, y)

        coefficientmatrix = rk4(coefficientmatrix, dt, K1matrix, K2matrix, K3matrix, K4matrix, N, factors, K, μ,y)

        ux[:, 2:end-1, i], uy[:, 2:end-1, i] = calculateuxuy(coefficientmatrix, N, Lx, x, l, m)

        savecoefficients[:, :, :, i] = coefficientmatrix

    end

    return ux, uy, x, y, savecoefficients

end

function sum(savecoefficients, N, T, m)

    coeffient = zeros(m, 4, T)

    for j in 1:T

        for i in 1:N+1

            coeffient[:, :, j] += savecoefficients[:, :, i, j]

        end 

        
    end

    return coeffient

end




@time u_x, u_y, x, y, savecoefficients = sumulation(20, 0.1, 20, 1, 2, 0.2, 100, 0.00001, 5., 2.)

f = Figure(size = (800, 800))

Axis(f[1, 1])
arrows(x, y, u_x[:, :, 100], u_y[:, :, 100], arrowsize = 10, lengthscale = 0.1)



