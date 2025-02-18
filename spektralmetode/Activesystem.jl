
# function setinitialu_active(N, F, l, m, y, Ly, T, ϕ)
    
#     ux = zeros(2*N+1, m, T)

#     uy = zeros(2*N+1, m, T)
    
#     for i in 1:2*N+1

#         ux[i, :, 1] = F*sin.(((2 * π )/(Ly))*y ) + cos.(ϕ[i, :, 1])

#         uy[i, :, 1] = F*sin.(((2 * π )/(Ly))*y ) + sin.(ϕ[i, :, 1])

#         #ux[i, :, 1] = F*y.^3

#         #uy[i, :, 1] = F*y.^3 .*sin.(((2 * π )/(Ly))*y )

#     end 

#     return ux, uy

# end


function derivativeelastic_active(coefficientmatrix, y)

    spilnes = [Spline1D(y, coefficientmatrix[:, i]) for i in 1:4]

    derivativ = [derivative(spilne, y) for spilne in spilnes]

    doublederivativ = [derivative(spilne, y, nu = 2) for spilne in spilnes]

    return  derivativ[1], derivativ[2], derivativ[3], derivativ[4], doublederivativ[1], doublederivativ[2], doublederivativ[3], doublederivativ[4]

end

function interpolatderivpolar(polarisatincoefficients, y)

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


function elastictimederivs_active(elcoefficientmatrix, Kmatrix, factor, K, μ, y)

    an_deriv, bn_deriv, cn_deriv, dn_deriv, an_doublederiv, bn_doublederiv, cn_doublederiv, dn_doublederiv = derivativeelastic_active(elcoefficientmatrix, y)

    Kmatrix[2:end-1, 1] =  μ * an_doublederiv[2:end-1] - (K + μ) * (factor)^2 * elcoefficientmatrix[2:end-1, 1] + K  * factor * dn_deriv[2:end-1]

    Kmatrix[2:end-1, 2] = μ * bn_doublederiv[2:end-1] - (K + μ) * (factor)^2 * elcoefficientmatrix[2:end-1, 2]  - K  * factor * cn_deriv[2:end-1]

    Kmatrix[2:end-1, 3] = (K +  μ) * cn_doublederiv[2:end-1] - μ * (factor)^2 * elcoefficientmatrix[2:end-1, 3] + K * factor * bn_deriv[2:end-1]

    Kmatrix[2:end-1, 4] = (K + μ)* dn_doublederiv[2:end-1] - μ * (factor)^2 * elcoefficientmatrix[2:end-1, 4] - K * factor * an_deriv[2:end-1] 

    return Kmatrix


end

function runge_kuttastep(elcoefficientmatrix, polarisationtioncoefficients, Dmatrix, Kmatrix, N, factor, K, μ, λ,  y)

    for i in 1:N+1

        Dmatrix[:, :, i] = polarisationtimederivs_active(polarisationtioncoefficients[:, :, i], Dmatrix[:, :, i], factor[i], λ, y)

        Kmatrix[:, :, i] = elastictimederivs_active(elcoefficientmatrix[:, :, i], Kmatrix[:, :, i], factor[i], K, μ, y)


    end
    
    return Dmatrix, Kmatrix

end

function rk4_active(elcoefficientmatrix, polarisationtioncoefficients, Dmatrix, Kmatrix, dt, N, factor, K, μ, λ, y)

    Dmatrix[1], Kmatrix[1] = runge_kuttastep(elcoefficientmatrix, polarisationtioncoefficients, Dmatrix[1], Kmatrix[1], N, factor, K, μ, λ, y)


    Dmatrix[2], Kmatrix[2] = runge_kuttastep(elcoefficientmatrix + dt/2 * Kmatrix[1], polarisationtioncoefficients + dt/2 * Dmatrix[1], Dmatrix[2], Kmatrix[2], N, factor, K, μ, λ, y)


    Dmatrix[3], Kmatrix[3] = runge_kuttastep(elcoefficientmatrix + dt/2 * Kmatrix[2], polarisationtioncoefficients + dt/2 * Dmatrix[2], Dmatrix[3], Kmatrix[3], N, factor, K, μ, λ, y)


    Dmatrix[4], Kmatrix[4] = runge_kuttastep(elcoefficientmatrix + dt * Kmatrix[3], polarisationtioncoefficients + dt * Dmatrix[3], Dmatrix[4], Kmatrix[4], N, factor, K, μ, λ, y)


    elcoefficientmatrix += dt/6 * (Kmatrix[1] + 2 * Kmatrix[2] + 2 * Kmatrix[3] + Kmatrix[4])


    polarisationtioncoefficients += dt/6 * (Dmatrix[1] + 2 * Dmatrix[2] + 2 * Dmatrix[3] + Dmatrix[4])


    return elcoefficientmatrix, polarisationtioncoefficients

end


function calculateuxuy_active(coefficientmatrix, N, Lx, x, m, ux, uy, ϕ, F_a)

    #uxstep = zeros(2*N+1, m-2)

    #uystep = zeros(2*N+1, m-2)

    for j in 2:m-1

        for i in 0:N 

            ux[:, j] += (coefficientmatrix[j, 1, i+1] * cos.(((2 * π * i)/Lx)*x ) + coefficientmatrix[j, 2, i+1] * sin.(((2 * π * i)/Lx)*x ))


            uy[:, j] += (coefficientmatrix[j, 3, i+1] * cos.(((2 * π * i)/Lx)*x ) + coefficientmatrix[j, 4, i+1] * sin.(((2 * π * i)/Lx)*x ))
        end 
     

    end

    return (1/(2*N+1))*ux + F_a*cos.(ϕ), (1/(2*N+1))*uy + F_a*sin.(ϕ)

end

function initializeelastic(N, F, m, y, Ly, T, ϕ)

    ux, uy = setinitialu(N, F, m, y, Ly, T)

    coefficientmatrix = initcoefficientmatrix(m, N, ux, uy)

    ux[:, :, 1] += cos.(ϕ[:, :, 1])

    uy[:, :, 1] += sin.(ϕ[:, :, 1])

    K1matrix= deepcopy(coefficientmatrix)

    K2matrix = deepcopy(coefficientmatrix)

    K3matrix = deepcopy(coefficientmatrix)

    K4matrix = deepcopy(coefficientmatrix)

    Kmatrix = [K1matrix, K2matrix, K3matrix, K4matrix]

    return ux, uy, coefficientmatrix, Kmatrix


end 

function run_activesystem(N::Int, dy::Float64, Lx::Int, Ly::Int, F::Float64, T::Int, dt::Float64, K::Float64, μ::Float64, λ::Float64, F_a::Float64)

    x, y, ϕ, factors, m = runsystem(N, dy, Lx, Ly, F, T, dt, λ)

    ux, uy, coefficientmatrix, Kmatrix = initializeelastic(N, F, m, y, Ly, T, ϕ)

    for i in 2:T

        elcoefficientmatrix, polarisationtioncoefficients = rk4_active(elcoefficientmatrix, polarisationtioncoefficients, Dmatrix, Kmatrix, dt, N, factors, K, μ, λ, y)

        ux[:, :, i], uy[:, :, i] = calculateuxuy_active(coefficientmatrix, N, Lx, x, m, ux[:, :, i], uy[:, :, i], ϕ[:, :, i], F_a)

    end

    return ux, uy, ϕ, x, y

end


