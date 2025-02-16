


function polarisationtimederivs(polarisatincoefficients, Dmatrix, factor, λ, y, N)

    for n in 1:N+1

        fn_doublederiv, gn_doublederiv = interpolatderivpolar(polarisatincoefficients, y, n)

        Dmatrix[:, 1, n] = λ*( -(factor[n]^2) * polarisatincoefficients[:, 1, n] + fn_doublederiv)

        Dmatrix[:, 2, n] = λ*( -(factor[n]^2) * polarisatincoefficients[:, 2, n] + gn_doublederiv)


    return Dmatrix

    end

end

# function elastictimederivs(elcoefficientmatrix, Kmatrix, factor, K, μ, y, n)

#     an_deriv, bn_deriv, cn_deriv, dn_deriv, an_doublederiv, bn_doublederiv, cn_doublederiv, dn_doublederiv = initerpolatderivate(elcoefficientmatrix, y, n)

#     Kmatrix[2:end-1, 1,  n] =  μ * an_doublederiv[2:end-1] - (K + μ) * (factor[n])^2 * elcoefficientmatrix[2:end-1, 1, n] + K  * factor[n] * dn_deriv[2:end-1]

#     Kmatrix[2:end-1, 2, n] = μ * bn_doublederiv[2:end-1] - (K + μ) * (factor[n])^2 * elcoefficientmatrix[2:end-1, 2, n]  - K  * factor[n] * cn_deriv[2:end-1]

#     Kmatrix[2:end-1, 3, n] = (K +  μ) * cn_doublederiv[2:end-1] - μ * (factor[n])^2 * elcoefficientmatrix[2:end-1, 3, n] + K * factor[n] * bn_deriv[2:end-1]

#     Kmatrix[2:end-1, 4, n] = (K + μ)* dn_doublederiv[2:end-1] - μ * (factor[n])^2 * elcoefficientmatrix[2:end-1, 4, n] - K * factor[n] * an_deriv[2:end-1] 

#     return Kmatrix


# end

# function runge_kuttastep(elcoefficientmatrix, polarisationtioncoefficients, Dmatrix, Kmatrix, N, factor, K, μ, λ,  y)

#     for i in 1:N+1

#         Dmatrix[:, :, i] = polarisationtimederivs(polarisationtioncoefficients, Dmatrix, factor, λ, y, i)

#         #Kmatrix[:, :, i] = elastictimederivs(elcoefficientmatrix, Kmatrix, factor, K, μ, y, i)


#     end
    
#     return Dmatrix#, Kmatrix

# end

function initploarcoefficients(N, m, ϕ)

    polarisationtioncoefficients = zeros(m, 2, N+1)

    for i in 1:m

        ft_ϕ = fft(ϕ[:, i, 1])

        polarisationtioncoefficients[i, 1, :] = real.(ft_ϕ)[1:N+1]

        polarisationtioncoefficients[i, 2, :] = imag.(ft_ϕ)[1:N+1]


    end

    return polarisationtioncoefficients

end


function interpolatderivpolar(polarisatincoefficients, y, n)

    #println(n)

    spilnes = [Spline1D(y, polarisatincoefficients[:, j, n]) for j in 1:2]


    derivativ = [derivative(spilne, y, nu =2) for spilne in spilnes]

    return derivativ[1], derivativ[2]

end


function initialcondition(N, F, m, y, Ly, T)

    ϕ = zeros(2*N+1, m, T)
    

    for i in 1:2*N+1

        ϕ[i, :, 1] =  rand(1:10) * sin.(((2 * π )/Ly)*y)

        # ϕ[10, 10, 1] = 1000        

    end

    return ϕ

end

function calculate_ϕ(polarisationtioncoefficients, N, Lx, x, m)

    ϕ_step = zeros(2*N+1, m)

    for j in 1:m
        
        for i in 0:N
            ϕ_step[:, j] += polarisationtioncoefficients[j, 1, i+1] * cos.(((2 * π * i)/Lx)*x ) + polarisationtioncoefficients[j, 2, i+1] * sin.(((2 * π * i)/Lx)*x )
        end
    end

    return (1/(2*N+1))*ϕ_step
end


function runge_kutta_polarization(polarisatincoefficients, Dmatrix, N, factor, λ, y, dt)



    Dmatrix[1] = polarisationtimederivs(polarisatincoefficients, Dmatrix[1], factor, λ, y, N)

    Dmatrix[2] = polarisationtimederivs(polarisatincoefficients + dt/2 * Dmatrix[1], Dmatrix[2], factor, λ, y, N)

    Dmatrix[3] = polarisationtimederivs(polarisatincoefficients + dt/2 * Dmatrix[2], Dmatrix[3], factor, λ, y, N)

    Dmatrix[4] = polarisationtimederivs(polarisatincoefficients + dt * Dmatrix[3], Dmatrix[4], factor, λ, y, N)

    polarisatincoefficients += dt/6 * (Dmatrix[1] + 2 * Dmatrix[2] + 2 * Dmatrix[3] + Dmatrix[4])


    return polarisatincoefficients

end



function initalize(dy::Float64, N::Int, Lx::Int, Ly::Int, F::Float64, T::Int)

    x = LinRange(-Lx/2, Lx/2, 2*N+1)

    y = -Ly/2:dy:Ly/2

    m = length(y)

    ϕ = initialcondition(N, F, m, y, Ly, T)

    factors = prefactors(Lx, N+1)

    polarisationtioncoefficients = initploarcoefficients(N, m, ϕ)

    D1matrix = zeros(m, 2, N+1)

    D2matrix = deepcopy(D1matrix)

    D3matrix = deepcopy(D1matrix)

    D4matrix = deepcopy(D1matrix)

    Dmatrix = [D1matrix, D2matrix, D3matrix, D4matrix]

    return x, y, m, ϕ, factors, polarisationtioncoefficients, Dmatrix

end

function runsystem(N::Int, dy::Float64, Lx::Int, Ly::Int, F::Float64, T::Int, dt::Float64, λ::Float64)

    x, y, m, ϕ, factors, polarisationtioncoefficients, Dmatrix = initalize(dy, N, Lx, Ly, F, T)


    for i in 2:T

        polarisationtioncoefficients = runge_kutta_polarization(polarisationtioncoefficients, Dmatrix, N, factors, λ, y, dt)

        ϕ[:, :, i] = calculate_ϕ(polarisationtioncoefficients, N, Lx, x, m)

    end

    return x, y, ϕ, factors, m

end




