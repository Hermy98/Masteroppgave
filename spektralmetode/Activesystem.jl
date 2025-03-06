function initalize_active(dy::Float64, dx::Float64, N::Int, Lx::Int, Ly::Int, F::Float64, F_a::Float64, T::Int)

    x = -Lx/2:dx:Lx/2

    y = -Ly/2:dy:Ly/2

    l = length(x)

    m = length(y)

    ux, uy, ϕ = initialcondition_active(F, F_a, m, l, y, Ly, T, Lx, x)

    factors = prefactors(Lx, N)

    elcoefficientmatrix, polarisationtioncoefficients = initialcoefficients(ux, uy, ϕ, m, N)

    Dmatrix = initialize_dmatrix(m, N)

    Kmatrix = initialze_kmatrix(m, N)

    return x, y, m, l, ϕ, ux, uy,factors, elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix

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


function initialcondition_active(F, F_a, m, l,y, Ly, T, Lx, x)
    
    ux = zeros(l, m, T)

    uy = zeros(l, m, T)

    ϕ = zeros(l, m, T)
    
    for i in 1:m


        ϕ[:, i, 1] = sin.(((2 * π )/(Lx))*x)

        ux[:, i, 1] = F_a*cos.(ϕ[:, i, 1])

        uy[:, i, 1] = F_a*sin.(ϕ[:, i, 1])

        # ux[:, 1, 1] .= 0

        # ux[:, end, 1] .= 0

        # uy[:, 1, 1] .= 0

        # uy[:, end, 1] .= 0

        #ϕ[:, 1, 1] .= 0

        #ϕ[:, end, 1] .= 0

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

    spilnes = [Spline1D(y, polarisatincoefficients[:, j]) for j in 1:2]

    derivativ = [derivative(spline, y, nu =2) for spline in spilnes]

    return derivativ[1], derivativ[2]

end


function polarisationtimederivs_active(polarisatincoefficients, Dmatrix, nonlinear_active, factor, λ, ξ,y)

    fn_doublederiv, gn_doublederiv = derivativepolar_active(polarisatincoefficients, y)

    Dmatrix[2:end-1, 1] = nonlinear_active[2:end-1, 1]  #+ λ*( -(factor^2) * polarisatincoefficients[2:end-1, 1] + fn_doublederiv[2:end-1])

    Dmatrix[2:end-1, 2] = nonlinear_active[2:end-1, 2]  #+ λ*( -(factor^2) * polarisatincoefficients[2:end-1, 2] + gn_doublederiv[2:end-1])

    Dmatrix[1, :] = nonlinear_active[1, :]  #+ λ*(- factor^2) * polarisatincoefficients[1, :] #+ (-2*polarisatincoefficients[1, :] + 2*polarisatincoefficients[2, :]))

    Dmatrix[end, :] = nonlinear_active[end, :]  #+ λ * (- factor^2) * polarisatincoefficients[end, :] #+ (-2*polarisatincoefficients[end, :] +2* polarisatincoefficients[end-1, :]))

    return Dmatrix

end


function elastictimederivs_active(elcoefficientmatrix, Kmatrix, factor, K, μ, y, F_a, active_contribution)

    an_deriv, bn_deriv, cn_deriv, dn_deriv, an_doublederiv, bn_doublederiv, cn_doublederiv, dn_doublederiv = derivativeelastic_active(elcoefficientmatrix, y)

    Kmatrix[2:end-1, 1] =  μ * an_doublederiv[2:end-1] - (K + (3/2)*μ) * (factor)^2 * elcoefficientmatrix[2:end-1, 1] + (K + (1/2)*μ ) * factor * dn_deriv[2:end-1] + active_contribution[2:end-1, 1]

    Kmatrix[2:end-1, 2] = μ * bn_doublederiv[2:end-1] - (K + (3/2)*μ) * (factor)^2 * elcoefficientmatrix[2:end-1, 2]  - (K + (1/2)*μ ) * factor * cn_deriv[2:end-1] + active_contribution[2:end-1, 2]

    Kmatrix[2:end-1, 3] = (K + (3/2)* μ) * cn_doublederiv[2:end-1] - μ * (factor)^2 * elcoefficientmatrix[2:end-1, 3] + (K + (1/2)*μ )* factor * bn_deriv[2:end-1] + active_contribution[2:end-1, 3]

    Kmatrix[2:end-1, 4] = (K + (3/2)* μ)* dn_doublederiv[2:end-1] - μ * (factor)^2 * elcoefficientmatrix[2:end-1, 4] - (K + (1/2)*μ) * factor * an_deriv[2:end-1] + active_contribution[2:end-1, 4]


    Kmatrix[1, 1] = μ * (elcoefficientmatrix[end, 1] -2*elcoefficientmatrix[1, 1] + elcoefficientmatrix[2, 1]) - (K + (3/2)*μ) * (factor)^2 * elcoefficientmatrix[1,1]  + (K + (1/2)*μ) * factor * (elcoefficientmatrix[end, 4] - elcoefficientmatrix[2, 4])/2  + active_contribution[1, 1]

    Kmatrix[1, 2] = μ * (elcoefficientmatrix[end, 2] -2*elcoefficientmatrix[1, 2] + elcoefficientmatrix[2, 2]) - (K + (3/2)*μ) * (factor)^2 * elcoefficientmatrix[1,2]  - (K + (1/2)*μ) * factor * (elcoefficientmatrix[end, 3] - elcoefficientmatrix[2, 3])/2  + active_contribution[1, 2]

    Kmatrix[1, 3] = μ * (elcoefficientmatrix[end, 3] -2*elcoefficientmatrix[1, 3] + elcoefficientmatrix[2, 3]) - (K + (3/2)*μ) * (factor)^2 * elcoefficientmatrix[1,3]  + (K + (1/2)*μ) * factor * (elcoefficientmatrix[end, 2] - elcoefficientmatrix[2, 2])/2  + active_contribution[1, 3]

    Kmatrix[1, 4] = μ * (elcoefficientmatrix[end, 4] -2*elcoefficientmatrix[1, 4] + elcoefficientmatrix[2, 4]) - (K + (3/2)*μ) * (factor)^2 * elcoefficientmatrix[1,4]  - (K + (1/2)*μ) * factor * (elcoefficientmatrix[end, 1] - elcoefficientmatrix[2, 1])/2  + active_contribution[1, 4]

    Kmatrix[end, 1] = μ * (elcoefficientmatrix[1, 1] -2*elcoefficientmatrix[end, 1] + elcoefficientmatrix[end-1 , 1]) - (K + (3/2)*μ) * (factor)^2 * elcoefficientmatrix[end,1]  + (K + (1/2)*μ) * factor * (elcoefficientmatrix[1, 4] - elcoefficientmatrix[end-1, 4])/2  + active_contribution[end, 1]

    Kmatrix[end, 2] =  μ * (elcoefficientmatrix[1, 2] -2*elcoefficientmatrix[end, 2] + elcoefficientmatrix[end-1 , 2]) - (K + (3/2)*μ) * (factor)^2 * elcoefficientmatrix[end,2]  - (K + (1/2)*μ) * factor * (elcoefficientmatrix[1, 3] - elcoefficientmatrix[end-1, 3])/2  + active_contribution[end, 2]

    Kmatrix[end, 3] =  μ * (elcoefficientmatrix[1, 3] -2*elcoefficientmatrix[end, 3] + elcoefficientmatrix[end-1 , 3]) - (K + (3/2)*μ) * (factor)^2 * elcoefficientmatrix[end,3]  + (K + (1/2)*μ) * factor * (elcoefficientmatrix[1, 2] - elcoefficientmatrix[end-1, 2])/2  + active_contribution[end, 3]

    Kmatrix[end, 4] =  μ * (elcoefficientmatrix[1, 4] -2*elcoefficientmatrix[end, 4] + elcoefficientmatrix[end-1 , 4]) - (K + (3/2)*μ) * (factor)^2 * elcoefficientmatrix[end,4]  - (K + (1/2)*μ) * factor * (elcoefficientmatrix[1, 1] - elcoefficientmatrix[end-1, 1])/2  + active_contribution[end, 4]

    return Kmatrix

end

function velocity_firststep(ux, uy, x, y, μ, K, ϕ, F_a, dt_u)

    ux_periodic = hcat(ux[:, end], ux, ux[:, 1])

    uy_periodic = hcat(uy[:, end], uy, uy[:, 1])

    y_periodic = vcat(y[1] - 1, y, y[end] + 1)

    ux_spline = Spline2D(x, y_periodic, ux_periodic)

    uy_spline = Spline2D(x, y_periodic, uy_periodic)

    ux_xxderiv = derivative(ux_spline, x, y_periodic, nux = 2, nuy = 0)[:, 2:end-1]

    ux_yyderiv = derivative(ux_spline, x, y_periodic, nux = 0, nuy = 2)[:, 2:end-1]

    ux_xyderiv = derivative(ux_spline, x, y_periodic, nux = 1, nuy = 1)[:, 2:end-1]

    uy_xxderiv = derivative(uy_spline, x, y_periodic, nux = 2, nuy = 0)[:, 2:end-1]

    uy_yyderiv = derivative(uy_spline, x, y_periodic, nux = 0, nuy = 2)[:, 2:end-1]

    uy_xyderiv = derivative(uy_spline, x, y_periodic, nux = 1, nuy = 1)[:, 2:end-1]


    dt_u[:, :, 1] = μ * (ux_xxderiv + ux_yyderiv) + (K + (1/2)*μ) * (ux_xxderiv + uy_xyderiv) + F_a*cos.(ϕ)

    dt_u[:, :, 2] = μ * (uy_xxderiv + uy_yyderiv) + (K + (1/2)*μ) * (ux_xyderiv + uy_yyderiv) + F_a*sin.(ϕ)

    #dt_u[:, 1, :] .= 0 #for boundry conditions u = 0

    #dt_u[:, end, :] .= 0 #for boundry conditions u = 0


    return dt_u
    

end

function velocity(ux, uy, dt_u, dt, i)

    dt_u[:, :, 1] = (ux[:, :, i-1] - ux[:, :, i-2])/dt

    dt_u[:, :, 2] = (uy[:, :, i-1] - uy[:, :, i-2])/dt

    #dt_u[:, end, :] .= 0 #for boundry conditions u = 0
    #dt_u[:, 1, :] .= 0 #for boundry conditions u = 0

    return dt_u


end

function ϕ_transform(ϕ, active_contribution, N, F_a)

    ft_cosϕ = fft(F_a * cos.(ϕ))

    ft_sinϕ = fft(F_a * sin.(ϕ))

    active_contribution[1, :] = real.(ft_cosϕ)[1:N+1]

    active_contribution[2, :] = imag.(ft_cosϕ)[1:N+1]

    active_contribution[3, :] = real.(ft_sinϕ)[1:N+1]

    active_contribution[4, :] = imag.(ft_sinϕ)[1:N+1]

    return active_contribution

end

function non_lineartransform(dt_u, ϕ, non_linearactive, N, ξ)

    fftU = fft(-ξ* sin.(ϕ) .* dt_u[:, 1])
    fftU2 = fft(ξ* cos.(ϕ) .* dt_u[:, 2])

    non_linearactive[1, :] = real.(fftU)[1:N+1] + real.(fftU2)[1:N+1]

    non_linearactive[2, :] = imag.(fftU)[1:N+1] + imag.(fftU2)[1:N+1]


    return non_linearactive

end

function fouriertransform_active(ϕ, active_contribution, dt_u, nonlinear_active, N, m, F_a, ξ)

    for i in 1:m

        active_contribution[i , :, :] = ϕ_transform(ϕ[:, i], active_contribution[i, :, :], N, F_a)

        nonlinear_active[i, :, :] = non_lineartransform(dt_u[:, i, :], ϕ[:, i], nonlinear_active[i, :, :], N, ξ)


    end

    return active_contribution, nonlinear_active

end




function runge_kuttastep(elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix, N, factor, K, μ, λ, ξ, y, F_a, active_contribution, nonlinear_active)

    for i in 2:N+1 #zeroth mode does not change in time

        Dmatrix[:, :, i] = polarisationtimederivs_active(polarisationtioncoefficients[:, :, i], Dmatrix[:, :, i], nonlinear_active[:, :, i], factor[i-1], λ, ξ, y)

        Kmatrix[:, :, i] = elastictimederivs_active(elcoefficientmatrix[:, :, i], Kmatrix[:, :, i], factor[i-1], K, μ, y, F_a, active_contribution[:, :, i])


    end
    
    return  Kmatrix, Dmatrix

end

function rk4_active(elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix, dt, N, factor, K, μ, λ, y, F_a, active_contribution, nonlinear_active, ξ)

    Kmatrix[1], Dmatrix[1] = runge_kuttastep(elcoefficientmatrix, Kmatrix[1], polarisationtioncoefficients, Dmatrix[1], N, factor, K, μ, λ, ξ, y, F_a, active_contribution, nonlinear_active)


    Kmatrix[2], Dmatrix[2] = runge_kuttastep(elcoefficientmatrix + dt/2 * Kmatrix[1], Kmatrix[2],polarisationtioncoefficients + dt/2 * Dmatrix[1], Dmatrix[2], N, factor, K, μ, λ, ξ, y, F_a, active_contribution + dt/2 * Kmatrix[1], nonlinear_active + dt/2*Dmatrix[1])


    Kmatrix[3], Dmatrix[3] = runge_kuttastep(elcoefficientmatrix + dt/2 * Kmatrix[2], Kmatrix[3],polarisationtioncoefficients + dt/2 * Dmatrix[2], Dmatrix[3], N, factor, K, μ, λ, ξ, y, F_a, active_contribution + dt/2 * Kmatrix[2], nonlinear_active + dt/2*Dmatrix[2])


    Kmatrix[4], Dmatrix[4] = runge_kuttastep(elcoefficientmatrix + dt * Kmatrix[3], Kmatrix[4],polarisationtioncoefficients + dt * Dmatrix[3], Dmatrix[4], N, factor, K, μ, λ, ξ, y, F_a, active_contribution + dt * Kmatrix[3], nonlinear_active + dt*Dmatrix[3])


    elcoefficientmatrix += dt/6 * (Kmatrix[1] + 2 * Kmatrix[2] + 2 * Kmatrix[3] + Kmatrix[4])


    polarisationtioncoefficients += dt/6 * (Dmatrix[1] + 2 * Dmatrix[2] + 2 * Dmatrix[3] + Dmatrix[4])


    return elcoefficientmatrix, polarisationtioncoefficients

end


function calculateuxuy_active(elcoefficientmatrix, polarisationtioncoefficients, N, Lx, x, m, ux, uy, ϕ, l, T, d_t)

    
    for j in 1:m

        ϕ[:, j] .+= polarisationtioncoefficients[j, 1, 1] *(T*d_t)/2

        ux[:, j] .+= elcoefficientmatrix[j, 1, 1] * (T*d_t)/2

        uy[:, j] .+= elcoefficientmatrix[j, 3, 1] * (T*d_t)/2

        for i in 2:N+1 

            ϕ[:, j] += polarisationtioncoefficients[j, 1, i]*cos.(((2 * π * (i-1))/Lx)*x ) + polarisationtioncoefficients[j, 2, i]*sin.(((2 * π * (i-1))/Lx)*x )

            ux[:, j] += (elcoefficientmatrix[j, 1, i] * cos.(((2 * π * (i-1))/Lx)*x ) + elcoefficientmatrix[j, 2, i] * sin.(((2 * π * (i-1))/Lx)*x ))

            uy[:, j] += (elcoefficientmatrix[j, 3, i] * cos.(((2 * π * (i-1))/Lx)*x ) + elcoefficientmatrix[j, 4, i] * sin.(((2 * π * (i-1))/Lx)*x ))
                
        end 
     

    end

    return (1/N)*ux, (1/N)*uy, (1/N)*ϕ# ux, uy, ϕ 

end



function run_activesystem(N::Int, dy::Float64, dx::Float64, Lx::Int, Ly::Int, F::Float64, T::Int, dt::Float64, K::Float64, μ::Float64, λ::Float64, F_a::Float64, ξ::Float64)

    x, y, m, l, ϕ, ux, uy, factors, elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix = initalize_active(dy, dx, N, Lx, Ly, F, F_a, T)

    ux[:, :, 1], uy[:, :, 1], ϕ[:, :, 1] = calculateuxuy_active(elcoefficientmatrix, polarisationtioncoefficients, N, Lx, x, m, ux[:, :, 1], uy[:, :, 1], ϕ[:, :, 1], l, 1, dt)

    active_contribution = zeros(m, 4, N+1)

    nonlinear_active = zeros(m, 2, N+1)

    dt_u = zeros(l, m, 2)

    dt_u = velocity_firststep(ux[:, :, 1], uy[:, :, 1], x, y, μ, K, ϕ[:, :, 1], F_a, dt_u)

    active_contribution, nonlinear_active = fouriertransform_active(ϕ[:, :, 1], active_contribution, dt_u, nonlinear_active, N, m, F_a, ξ)
    
    elcoefficientmatrix, polarisationtioncoefficients = rk4_active(elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix, dt, N, factors, K, μ, λ, y, F_a, active_contribution, nonlinear_active, ξ)

    ux[:, :, 2], uy[:, :, 2], ϕ[:, :, 2] = calculateuxuy_active(elcoefficientmatrix, polarisationtioncoefficients, N, Lx, x, m, ux[:, :, 2], uy[:, :, 2], ϕ[:, :, 2], l, 2, dt)




    for i in 3:T

        dt_u = velocity_firststep(ux[:, :, i-1], uy[:, :, i-1], x, y, μ, K, ϕ[:, :, i-1], F_a, dt_u)

        #dt_u = velocity(ux, uy, dt_u, dt, i)

        active_contribution, nonlinear_active = fouriertransform_active(ϕ[:, :, i-1], active_contribution, dt_u, nonlinear_active, N, m, F_a, ξ)

        elcoefficientmatrix, polarisationtioncoefficients = rk4_active(elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix, dt, N, factors, K, μ, λ, y, F_a, active_contribution, nonlinear_active, ξ)

        ux[:, :, i], uy[:, :, i], ϕ[:, :, i] = calculateuxuy_active(elcoefficientmatrix, polarisationtioncoefficients, N, Lx, x, m, ux[:, :, i], uy[:, :, i], ϕ[:, :, i], l, i, dt)


    end

    return ux, uy, ϕ, x, y, m

end



