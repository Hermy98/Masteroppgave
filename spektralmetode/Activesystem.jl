function initalize_active(dy::Float64, N::Int, Lx::Int, Ly::Int, F::Float64, F_a::Float64, T::Int)

    x = LinRange(-Lx/2, Lx/2, 2*N+1)

    y = -Ly/2:dy:Ly/2

    m = length(y)

    ux, uy, ϕ = initialcondition_active(N, F, F_a, m, y, Ly, T)

    factors = prefactors(Lx, N+1)

    elcoefficientmatrix, polarisationtioncoefficients = initialcoefficients(ux, uy, ϕ, m, N)

    Dmatrix = initialize_dmatrix(m, N)

    Kmatrix = initialze_kmatrix(m, N)

    return x, y, m, ϕ, ux, uy,factors, elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix

end

function initialize_dmatrix(m, N)

    D1matrix = zeros(m, 2, N+1)

    D2matrix = deepcopy(D1matrix)

    D3matrix = deepcopy(D1matrix)

    D4matrix = deepcopy(D1matrix)

    Dmatrix = [D1matrix, D2matrix, D3matrix, D4matrix]


    return Dmatrix

end

function initialze_kmatrix(m, N)

    K1matrix = zeros(m, 4, N+1)

    K2matrix = deepcopy(K1matrix)

    K3matrix = deepcopy(K1matrix)

    K4matrix = deepcopy(K1matrix)

    Kmatrix = [K1matrix, K2matrix, K3matrix, K4matrix]


    return Kmatrix

end


function initialcondition_active(N, F, F_a, m, y, Ly, T)
    
    ux = zeros(2*N+1, m, T)

    uy = zeros(2*N+1, m, T)

    ϕ = zeros(2*N+1, m, T)
    
    for i in 1:2*N+1

        ϕ[i, :, 1] =  sin.(((2 * π) /Ly)*y)

        ux[i, :, 1] = F*sin.(((2 * π )/(Ly))*y ) + F_a * cos.(ϕ[i, :, 1])

        uy[i, :, 1] = F*sin.(((2 * π )/(Ly))*y ) + F_a * sin.(ϕ[i, :, 1])

        #ux[i, :, 1] = F*y.^3

        #uy[i, :, 1] = F*y.^3 .*sin.(((2 * π )/(Ly))*y )

    end 

    return ux, uy, ϕ

end

function initialcoefficients(ux, uy, ϕ, m, N)


    polarisationtioncoefficients = zeros(m, 2, N+1)

    elcoefficientmatrix = zeros(m, 4, N+1)

    for i in 1:m

        ft_ϕ = fft(ϕ[:, i, 1])

        polarisationtioncoefficients[i, 1, :] = real.(ft_ϕ)[1:N+1]

        polarisationtioncoefficients[i, 2, :] = imag.(ft_ϕ)[1:N+1]

        ft_ux = fft(ux[:, i, 1])

        ft_uy = fft(uy[:, i, 1])

        elcoefficientmatrix[i, 1, :] = real.(ft_ux)[1:N+1]

        elcoefficientmatrix[i, 2, :] = imag.(ft_ux)[1:N+1]

        elcoefficientmatrix[i, 3, :] = real.(ft_uy)[1:N+1]

        elcoefficientmatrix[i, 4, :] = imag.(ft_uy)[1:N+1]

    end

    return elcoefficientmatrix,polarisationtioncoefficients

end


function derivativeelastic_active(coefficientmatrix, y)

    spilnes = [Spline1D(y, coefficientmatrix[:, i]) for i in 1:4]

    derivativ = [derivative(spilne, y) for spilne in spilnes]

    doublederivativ = [derivative(spilne, y, nu = 2) for spilne in spilnes]

    return  derivativ[1], derivativ[2], derivativ[3], derivativ[4], doublederivativ[1], doublederivativ[2], doublederivativ[3], doublederivativ[4]

end

function derivativepolar_active(polarisatincoefficients, y)

    #println(n)

    spilnes = [Spline1D(y, polarisatincoefficients[:, j]) for j in 1:2]


    derivativ = [derivative(spline, y, nu =2) for spline in spilnes]

    return derivativ[1], derivativ[2]

end


function polarisationtimederivs_active(polarisatincoefficients, Dmatrix, factor, λ, y)

    fn_doublederiv, gn_doublederiv = derivativepolar_active(polarisatincoefficients, y)

    Dmatrix[2:end-1, 1] = λ*( -(factor^2) * polarisatincoefficients[2:end-1, 1] + fn_doublederiv[2:end-1])

    Dmatrix[2:end-1, 2] = λ*( -(factor^2) * polarisatincoefficients[2:end-1, 2] + gn_doublederiv[2:end-1])

    Dmatrix[1, 1] = λ*( (- factor^2) * polarisatincoefficients[1, 1] + (-2*polarisatincoefficients[1, 1] + 2*polarisatincoefficients[2, 1]))

    Dmatrix[1, 2] = λ*( (- factor^2) * polarisatincoefficients[1, 2] + (-2*polarisatincoefficients[1, 2] + 2* polarisatincoefficients[2, 2]))

    Dmatrix[end, 1] = λ* ((- factor^2) * polarisatincoefficients[end, 1] + (-2*polarisatincoefficients[end, 1] +2* polarisatincoefficients[end-1, 1]))

    Dmatrix[end, 2] = λ* ((- factor^2) * polarisatincoefficients[end, 2] + (-2*polarisatincoefficients[end, 2] +2* polarisatincoefficients[end-1, 2]))

    

    return Dmatrix

end


function elastictimederivs_active(elcoefficientmatrix, Kmatrix, factor, K, μ, y, F_a, active_contribution)

    an_deriv, bn_deriv, cn_deriv, dn_deriv, an_doublederiv, bn_doublederiv, cn_doublederiv, dn_doublederiv = derivativeelastic_active(elcoefficientmatrix, y)

    Kmatrix[2:end-1, 1] =  μ * an_doublederiv[2:end-1] - (K + μ) * (factor)^2 * elcoefficientmatrix[2:end-1, 1] + K  * factor * dn_deriv[2:end-1] + F_a * active_contribution[2:end-1, 1]

    Kmatrix[2:end-1, 2] = μ * bn_doublederiv[2:end-1] - (K + μ) * (factor)^2 * elcoefficientmatrix[2:end-1, 2]  - K  * factor * cn_deriv[2:end-1] + F_a * active_contribution[2:end-1, 2]

    Kmatrix[2:end-1, 3] = (K +  μ) * cn_doublederiv[2:end-1] - μ * (factor)^2 * elcoefficientmatrix[2:end-1, 3] + K * factor * bn_deriv[2:end-1] + F_a * active_contribution[2:end-1, 3]

    Kmatrix[2:end-1, 4] = (K + μ)* dn_doublederiv[2:end-1] - μ * (factor)^2 * elcoefficientmatrix[2:end-1, 4] - K * factor * an_deriv[2:end-1] + F_a * active_contribution[2:end-1, 4]

    Kmatrix[1, 1] = F_a * active_contribution[1, 1]

    Kmatrix[1, 2] = F_a * active_contribution[1, 2]

    Kmatrix[1, 3] = F_a * active_contribution[1, 3]

    Kmatrix[1, 4] = F_a * active_contribution[1, 4]

    Kmatrix[end, 1] = F_a * active_contribution[end, 1]

    Kmatrix[end, 2] = F_a * active_contribution[end, 2]

    Kmatrix[end, 3] = F_a * active_contribution[end, 3]

    Kmatrix[end, 4] = F_a * active_contribution[end, 4]


    return Kmatrix

end

function fouriertransform_active(ϕ, active_contribution, N, m)

    for i in 1:m

        ft_cosϕ = fft(cos.(ϕ[:, i]))

        ft_sinϕ = fft(sin.(ϕ[:, i]))

        active_contribution[i, 1, :] = real.(ft_cosϕ)[1:N+1]

        active_contribution[i, 2, :] = imag.(ft_cosϕ)[1:N+1]

        active_contribution[i, 3, :] = real.(ft_sinϕ)[1:N+1]

        active_contribution[i, 4, :] = imag.(ft_sinϕ)[1:N+1]

    end

    return active_contribution

end




function runge_kuttastep(elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix, N, factor, K, μ, λ,  y, F_a, active_contribution)

    for i in 1:N+1

        Dmatrix[:, :, i] = polarisationtimederivs_active(polarisationtioncoefficients[:, :, i], Dmatrix[:, :, i], factor[i], λ, y)

        Kmatrix[:, :, i] = elastictimederivs_active(elcoefficientmatrix[:, :, i], Kmatrix[:, :, i], factor[i], K, μ, y, F_a, active_contribution[:, :, i])


    end
    
    return  Kmatrix, Dmatrix

end

function rk4_active(elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix, dt, N, factor, K, μ, λ, y, F_a, active_contribution)

    Kmatrix[1], Dmatrix[1] = runge_kuttastep(elcoefficientmatrix, Kmatrix[1], polarisationtioncoefficients, Dmatrix[1], N, factor, K, μ, λ, y, F_a, active_contribution)


    Kmatrix[2], Dmatrix[2] = runge_kuttastep(elcoefficientmatrix + dt/2 * Kmatrix[1], Kmatrix[2],polarisationtioncoefficients + dt/2 * Dmatrix[1], Dmatrix[2], N, factor, K, μ, λ, y, F_a, active_contribution)


    Kmatrix[3], Dmatrix[3] = runge_kuttastep(elcoefficientmatrix + dt/2 * Kmatrix[2], Kmatrix[3],polarisationtioncoefficients + dt/2 * Dmatrix[2], Dmatrix[3], N, factor, K, μ, λ, y, F_a, active_contribution)


    Kmatrix[4], Dmatrix[4] = runge_kuttastep(elcoefficientmatrix + dt * Kmatrix[3], Kmatrix[4],polarisationtioncoefficients + dt * Dmatrix[3], Dmatrix[4], N, factor, K, μ, λ, y, F_a, active_contribution)


    elcoefficientmatrix += dt/6 * (Kmatrix[1] + 2 * Kmatrix[2] + 2 * Kmatrix[3] + Kmatrix[4])


    polarisationtioncoefficients += dt/6 * (Dmatrix[1] + 2 * Dmatrix[2] + 2 * Dmatrix[3] + Dmatrix[4])


    return elcoefficientmatrix, polarisationtioncoefficients

end


function calculateuxuy_active(elcoefficientmatrix, polarisationtioncoefficients, N, Lx, x, m, ux, uy, ϕ)


    for j in 1:m

        for i in 0:N 

            ϕ[:, j] += polarisationtioncoefficients[j, 1, i+1]*cos.(((2 * π * i)/Lx)*x ) + polarisationtioncoefficients[j, 2, i+1]*sin.(((2 * π * i)/Lx)*x )



            ux[:, j] += (elcoefficientmatrix[j, 1, i+1] * cos.(((2 * π * i)/Lx)*x ) + elcoefficientmatrix[j, 2, i+1] * sin.(((2 * π * i)/Lx)*x ))


            uy[:, j] += (elcoefficientmatrix[j, 3, i+1] * cos.(((2 * π * i)/Lx)*x ) + elcoefficientmatrix[j, 4, i+1] * sin.(((2 * π * i)/Lx)*x ))
                
        end 
     

    end

    return (1/(2*N+1))*ux, (1/(2*N+1))*uy, (1/(2*N+1))*ϕ

end



function run_activesystem(N::Int, dy::Float64, Lx::Int, Ly::Int, F::Float64, T::Int, dt::Float64, K::Float64, μ::Float64, λ::Float64, F_a::Float64)

    x, y, m, ϕ, ux, uy, factors, elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix = initalize_active(dy, N, Lx, Ly, F, F_a, T)

    active_contribution = zeros(m, 4, N+1)

    active_contribution = fouriertransform_active(ϕ[:, :, 1], active_contribution, N, m)


    for i in 2:T

        elcoefficientmatrix, polarisationtioncoefficients = rk4_active(elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix, dt, N, factors, K, μ, λ, y, F_a, active_contribution)

        ux[:, :, i], uy[:, :, i], ϕ[:, :, i] = calculateuxuy_active(elcoefficientmatrix, polarisationtioncoefficients, N, Lx, x, m, ux[:, :, i], uy[:, :, i], ϕ[:, :, i])

        active_contribution = fouriertransform_active(ϕ[:, :, i], active_contribution, N, m)

    end

    return ux, uy, ϕ, x, y

end


