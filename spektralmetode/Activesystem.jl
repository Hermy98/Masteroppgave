function calculateuxuy_active(elcoefficientmatrix::Array{Float64, 3}, polarisationtioncoefficients::Array{Float64, 3}, N::Int, Lx::Int, x::StepRangeLen, m::Int, ux::Matrix{Float64}, uy::Matrix{Float64}, ϕ::Matrix{Float64}, l::Int, T::Int, d_t::Float64)
    # Pre-compute wave numbers and trig functions
    k = [(2π * (i-1))/Lx for i in 1:N+1]
    cos_kx = [cos.(k[i] .* x) for i in 1:N+1]
    sin_kx = [sin.(k[i] .* x) for i in 1:N+1]
    
    # Pre-compute normalization factor
    norm_factor = 1/sqrt(l)
    
    @inbounds for j in 1:m
        # Initialize with first mode (i=1)
        @views begin
            ϕ[:, j] .= polarisationtioncoefficients[j, 1, 1]
            ux[:, j] .= elcoefficientmatrix[j, 1, 1]
            uy[:, j] .= elcoefficientmatrix[j, 3, 1]
        end
        
        # Add higher modes (i≥2)
        @simd for i in 2:N+1
            # Get coefficients once
            ϕ_cos = polarisationtioncoefficients[j, 1, i]
            ϕ_sin = polarisationtioncoefficients[j, 2, i]
            ux_cos = elcoefficientmatrix[j, 1, i]
            ux_sin = elcoefficientmatrix[j, 2, i]
            uy_cos = elcoefficientmatrix[j, 3, i]
            uy_sin = elcoefficientmatrix[j, 4, i]
            
            # Use pre-computed trig functions and broadcasting
            @views begin
                ϕ[:, j] .+= ϕ_cos .* cos_kx[i] .+ ϕ_sin .* sin_kx[i]
                ux[:, j] .+= ux_cos .* cos_kx[i] .+ ux_sin .* sin_kx[i]
                uy[:, j] .+= uy_cos .* cos_kx[i] .+ uy_sin .* sin_kx[i]
            end
        end
    end
    
    # Apply normalization
    return norm_factor * ux, norm_factor * uy, norm_factor * ϕ
end

function initalize_active(Loggingstep::Int, dy::Float64, dx::Float64, N::Int, Lx::Int, Ly::Int, F::Float64, F_a::Float64, T::Int, BC0::Bool)

    x = -Lx/2:dx:Lx/2

    y = -Ly/2:dy:Ly/2

    l = length(x)

    m = length(y)

    ux_temp = zeros(l, m)

    uy_temp = zeros(l, m)

    ϕ_temp = zeros(l, m)

    ux, uy, ϕ = initialcondition_active(Loggingstep, F, F_a, m, l, y, Ly, T, Lx, x, BC0)

    ux_temp = ux[:, :, 1]

    uy_temp = uy[:, :, 1]

    ϕ_temp = ϕ[:, :, 1]

    factors = prefactors(Lx, N)

    elcoefficientmatrix, polarisationtioncoefficients = initialcoefficients(ux, uy, ϕ, m, N, l)

    Dmatrix = initialize_dmatrix(m, N)

    Kmatrix = initialze_kmatrix(m, N)


    return x, y, m, l, ϕ, ux, uy,factors, elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix, ux_temp, uy_temp, ϕ_temp
end

function initialize_dmatrix(m::Int, N::Int)

    D1matrix = zeros(m, 2, N+1)

    D2matrix = deepcopy(D1matrix)

    D3matrix = deepcopy(D1matrix)

    D4matrix = deepcopy(D1matrix)

    Dmatrix = [D1matrix, D2matrix, D3matrix, D4matrix]


    return Dmatrix

end

function initialze_kmatrix(m::Int, N::Int)

    K1matrix = zeros(m, 4, N+1)

    K2matrix = deepcopy(K1matrix)

    K3matrix = deepcopy(K1matrix)

    K4matrix = deepcopy(K1matrix)

    Kmatrix = [K1matrix, K2matrix, K3matrix, K4matrix]


    return Kmatrix

end

function smothinginitialcond(m::Int, l::Int, ϕ::Array{Float64, 3})

    ϕ_temp = similar(ϕ[:, :, 1])

    for i in 2:l-1

        for j in 2:m-1

            ax =  (2*cos.(ϕ[i, j, 1]) + cos.(ϕ[i-1, j, 1]) + cos.(ϕ[i+1, j, 1]) )/4

            
            ay =  (2*sin.(ϕ[i, j, 1]) + sin.(ϕ[i-1, j, 1]) + sin.(ϕ[i+1, j, 1]) )/4

            ϕ_temp[i, j] = atan.(ay, ax)

            
        end

    end

    ϕ[:, :, 1] .= ϕ_temp

    return ϕ

end

function band_polarization(x, Lx)

    if x < Lx/3

        return 0.0

    elseif x < 2Lx/3

        return π/2 

    else

        return 0.0 

    end

end

function band_polarization_smooth(x::Float64, Lx::Int)
    w = Lx/15  # width of transition region
    middle = π/2 * (tanh((x - Lx/3)/w) - tanh((x - 2Lx/3)/w))/2
    return middle
end

function initialcondition_active(Loggingstep::Int, F::Float64, F_a::Float64, m::Int, l::Int,y::StepRangeLen, Ly::Int, T::Int, Lx::Int, x::StepRangeLen, BC0::Bool)
    
    t = T ÷ Loggingstep

    ux = zeros(l, m, t)

    uy = zeros(l, m, t)

    ϕ = zeros(l, m, t)


    for i in 1:m

        ϕ[:, i, 1] .= π/4 #(π/2 *0.9) .+ rand((-0.25, 0.25))

    end

    


    return ux, uy, ϕ

end

function initialcoefficients(ux::Array{Float64, 3}, uy::Array{Float64, 3}, ϕ::Array{Float64, 3}, m::Int, N::Int, l::Int) 

    polarisationtioncoefficients = zeros(m, 2, N+1)
    elcoefficientmatrix = zeros(m, 4, N+1)

    norm_factor = 1/sqrt(l)

    for i in 1:m
        ft_ϕ = norm_factor * fft(ϕ[:, i, 1])
        polarisationtioncoefficients[i, 1, :] = real.(ft_ϕ)[1:N+1]
        polarisationtioncoefficients[i, 2, :] = imag.(ft_ϕ)[1:N+1]

        ft_ux = norm_factor * fft(ux[:, i, 1])
        ft_uy = norm_factor * fft(uy[:, i, 1])
        
        elcoefficientmatrix[i, 1, :] = real.(ft_ux)[1:N+1]
        elcoefficientmatrix[i, 2, :] = imag.(ft_ux)[1:N+1]
        elcoefficientmatrix[i, 3, :] = real.(ft_uy)[1:N+1]
        elcoefficientmatrix[i, 4, :] = imag.(ft_uy)[1:N+1]
    end

    return elcoefficientmatrix,polarisationtioncoefficients

end


# function derivativeelastic_active(coefficientmatrix::Matrix{Float64}, y::StepRangeLen)



#     coefficientmatrix_periodic = vcat(coefficientmatrix[end, :]', coefficientmatrix, coefficientmatrix[1, :]')

#     y_periodic = vcat(y[1] - 1, y, y[end] + 1)

#     spilnes = [Spline1D(y_periodic, coefficientmatrix_periodic[:, i]) for i in 1:4]

#     derivativ = [derivative(spilne, y_periodic) for spilne in spilnes]

#     doublederivativ = [derivative(spilne, y_periodic, nu = 2) for spilne in spilnes]

#     return  derivativ[1][2:end-1], derivativ[2][2:end-1], derivativ[3][2:end-1], derivativ[4][2:end-1], doublederivativ[1][2:end-1], doublederivativ[2][2:end-1], doublederivativ[3][2:end-1], doublederivativ[4][2:end-1]

# end

function derivativeelastic_active(coefficientmatrix::Matrix{Float64}, y::StepRangeLen, derivselastic::Matrix{Float64}, derivdoubleelastic::Matrix{Float64})
 
    
    @views begin
        coefficientmatrix_periodic = vcat(coefficientmatrix[end, :]', coefficientmatrix, coefficientmatrix[1, :]')
        y_periodic = vcat(y[1] - 1, y, y[end] + 1)
        
        @inbounds for i in 1:4
            spline = Spline1D(y_periodic, coefficientmatrix_periodic[:, i])
            derivselastic[:, i] = derivative(spline, y_periodic)[2:end-1]
            derivdoubleelastic[:, i] = derivative(spline, y_periodic, nu=2)[2:end-1]
        end
    end
    
    return derivselastic[:, 1], derivselastic[:, 2], derivselastic[:, 3], derivselastic[:, 4],
           derivdoubleelastic[:, 1], derivdoubleelastic[:, 2], derivdoubleelastic[:, 3], derivdoubleelastic[:, 4]
end

function derivativepolar_active(polarisatincoefficients::Matrix{Float64}, y::StepRangeLen, derivspolar::Matrix{Float64})
    
    @views begin
        @inbounds for j in 1:2
            spline = Spline1D(y, polarisatincoefficients[:, j])
            derivspolar[:, j] = derivative(spline, y, nu=2)
        end
        
        return derivspolar[:, 1], derivspolar[:, 2]
    end
end



# function derivativepolar_active(polarisatincoefficients::Array{Float64, 2}, y::StepRangeLen)

#     spilnes = [Spline1D(y, polarisatincoefficients[:, j]) for j in 1:2]

#     derivativ = [derivative(spline, y, nu =2) for spline in spilnes]

#     return derivativ[1], derivativ[2]

# end

function polarisationtimederivs_periodic(Dmatrix::Array{Float64, 2}, nonlinear_active::Array{Float64, 2})


        Dmatrix = nonlinear_active

    return Dmatrix

end


function polarisationtimederivs_zero(polarisatincoefficients::Array{Float64, 2}, Dmatrix::Array{Float64, 2}, nonlinear_active::Array{Float64, 2}, factor::Float64, λ::Float64, y::StepRangeLen, derivspolar::Array{Float64, 2})

    fn_doublederiv, gn_doublederiv = derivativepolar_active(polarisatincoefficients, y, derivspolar)

    @views begin

        Dmatrix[2:end-1, 1] = nonlinear_active[2:end-1, 1]  + λ*( -(factor^2) * polarisatincoefficients[2:end-1, 1] + fn_doublederiv[2:end-1])

        Dmatrix[2:end-1,  2] = nonlinear_active[2:end-1, 2]  + λ*( -(factor^2) * polarisatincoefficients[2:end-1, 2] + gn_doublederiv[2:end-1])

        Dmatrix[1, :] =  λ*(-( factor^2) * polarisatincoefficients[1, :] + (-2*polarisatincoefficients[1, :] + 2*polarisatincoefficients[2, :]))

        Dmatrix[end, :] =  λ * (-( factor^2) * polarisatincoefficients[end, :] + (-2*polarisatincoefficients[end, :] + 2*polarisatincoefficients[end-1, :]))



    end

    return Dmatrix

end



function elastictimederivs_periodic(elcoefficientmatrix::Matrix{Float64}, Kmatrix::Matrix{Float64}, factor::Float64, K::Float64, μ::Float64, y::StepRangeLen, active_contribution::Matrix{Float64}, derivselastic::Matrix{Float64}, double_derivsel::Matrix{Float64})

    an_deriv, bn_deriv, cn_deriv, dn_deriv, an_doublederiv, bn_doublederiv, cn_doublederiv, dn_doublederiv = derivativeelastic_active(elcoefficientmatrix, y, derivselastic, double_derivsel)

    @views begin

        Kmatrix[:, 1] =  μ * an_doublederiv - (K +  μ) * (factor)^2 * elcoefficientmatrix[:, 1] + K   * factor * dn_deriv + active_contribution[:, 1]

        Kmatrix[:, 2] = μ * bn_doublederiv - (K + μ) * (factor)^2 * elcoefficientmatrix[:, 2]  - K  * factor * cn_deriv + active_contribution[:, 2]

        Kmatrix[:, 3] = (K + μ) * cn_doublederiv - μ * (factor)^2 * elcoefficientmatrix[: ,3] + K  * factor * bn_deriv+ active_contribution[:, 3]

        Kmatrix[:, 4] = (K + μ) * dn_doublederiv - μ * (factor)^2 * elcoefficientmatrix[:, 4] - K  * factor * an_deriv + active_contribution[:, 4]

    end

    return Kmatrix

end



function elastictimederivs_zero(elcoefficientmatrix::Matrix{Float64}, Kmatrix::Matrix{Float64}, factor::Float64, K::Float64, μ::Float64, y::StepRangeLen, active_contribution::Matrix{Float64}, derivselastic::Matrix{Float64}, double_derivsel::Matrix{Float64})

    an_deriv, bn_deriv, cn_deriv, dn_deriv, an_doublederiv, bn_doublederiv, cn_doublederiv, dn_doublederiv = derivativeelastic_active(elcoefficientmatrix, y, derivselastic, double_derivsel)

    @views begin

        Kmatrix[2:end-1, 1] = μ * an_doublederiv[2:end-1] - (K + μ) * (factor)^2 * elcoefficientmatrix[2:end-1, 1] + K * factor * dn_deriv[2:end-1] + active_contribution[2:end-1, 1]

        Kmatrix[2:end-1, 2] = μ * bn_doublederiv[2:end-1] - (K + μ) * (factor)^2 * elcoefficientmatrix[2:end-1, 2] - K  * factor * cn_deriv[2:end-1] + active_contribution[2:end-1, 2]

        Kmatrix[2:end-1, 3] = (K +  μ) * cn_doublederiv[2:end-1] - μ * (factor)^2 * elcoefficientmatrix[2:end-1 ,3] + K  * factor * bn_deriv[2:end-1] + active_contribution[2:end-1, 3]

        Kmatrix[2:end-1, 4] = (K +  μ) * dn_doublederiv[2:end-1] - μ * (factor)^2 * elcoefficientmatrix[2:end-1, 4] - K  * factor * an_deriv[2:end-1] + active_contribution[2:end-1, 4]

    end

    return Kmatrix

end

# function velocity(ux::Array{Float64, 3}, uy::Array{Float64, 3}, dt_u::Array{Float64, 3}, dt::Float64, i::Int64, x::StepRangeLen, y::StepRangeLen, μ::Float64, K::Float64, ϕ::Array{Float64, 3}, F_a::Float64)
#     # Use current state for derivatives
#     ux_periodic = hcat(ux[:, end, i], ux[:, :, i], ux[:, 1, i])
#     uy_periodic = hcat(uy[:, end, i], uy[:, :, i], uy[:, 1, i])
    
#     y_periodic = vcat(y[1] - 1, y, y[end] + 1)
    
#     # Calculate spatial derivatives using splines
#     ux_spline = Spline2D(x, y_periodic, ux_periodic)
#     uy_spline = Spline2D(x, y_periodic, uy_periodic)
    
#     # Second derivatives
#     ux_xxderiv = derivative(ux_spline, x, y_periodic, nux = 2, nuy = 0)[:, 2:end-1]
#     ux_yyderiv = derivative(ux_spline, x, y_periodic, nux = 0, nuy = 2)[:, 2:end-1]
#     uy_xxderiv = derivative(uy_spline, x, y_periodic, nux = 2, nuy = 0)[:, 2:end-1]
#     uy_yyderiv = derivative(uy_spline, x, y_periodic, nux = 0, nuy = 2)[:, 2:end-1]
    
#     # Mixed derivatives for elastic coupling
#     ux_xyderiv = derivative(ux_spline, x, y_periodic, nux = 1, nuy = 1)[:, 2:end-1]
#     uy_xyderiv = derivative(uy_spline, x, y_periodic, nux = 1, nuy = 1)[:, 2:end-1]
    
#     # Full equation including elastic terms
#     dt_u[:, :, 1] = μ * (ux_xxderiv + ux_yyderiv) + (K + (1/2)*μ) * (ux_xxderiv + uy_xyderiv) + F_a*cos.(ϕ[:, :, i])
#     dt_u[:, :, 2] = μ * (uy_xxderiv + uy_yyderiv) + (K + (1/2)*μ) * (ux_xyderiv + uy_yyderiv) + F_a*sin.(ϕ[:, :, i])
    
#     return dt_u
# end


function velocity(ux_temp::Matrix{Float64}, uy_temp::Matrix{Float64}, dt_u::Array{Float64, 3}, x::StepRangeLen, y::StepRangeLen, μ::Float64, K::Float64, ϕ_temp::Matrix{Float64}, F_a::Float64)
    
    
    @views begin
        ux_periodic = hcat(ux_temp[:, end], ux_temp, ux_temp[:, 1])
        uy_periodic = hcat(uy_temp[:, end], uy_temp, uy_temp[:, 1])
    end
    
    y_periodic = vcat(y[1] - 1, y, y[end] + 1)
    
   
    @views begin
        ux_spline = Spline2D(x, y_periodic, ux_periodic)
        uy_spline = Spline2D(x, y_periodic, uy_periodic)
        
   
        ux_xxderiv = similar(ux_temp, size(ux_temp, 1), size(ux_temp, 2))
        ux_yyderiv = similar(ux_xxderiv)
        ux_xyderiv = similar(ux_xxderiv)
        uy_xxderiv = similar(ux_xxderiv)
        uy_yyderiv = similar(ux_xxderiv)
        uy_xyderiv = similar(ux_xxderiv)
        

        ux_xxderiv = derivative(ux_spline, x, y_periodic, nux = 2, nuy = 0)[:, 2:end-1]
        ux_yyderiv = derivative(ux_spline, x, y_periodic, nux = 0, nuy = 2)[:, 2:end-1]
        uy_xxderiv = derivative(uy_spline, x, y_periodic, nux = 2, nuy = 0)[:, 2:end-1]
        uy_yyderiv = derivative(uy_spline, x, y_periodic, nux = 0, nuy = 2)[:, 2:end-1]
        
        ux_xyderiv = derivative(ux_spline, x, y_periodic, nux = 1, nuy = 1)[:, 2:end-1]
        uy_xyderiv = derivative(uy_spline, x, y_periodic, nux = 1, nuy = 1)[:, 2:end-1]
    end
    
    @views begin
        dt_u[:, :, 1] = μ * (ux_xxderiv + ux_yyderiv) + K * (ux_xxderiv + uy_xyderiv) + F_a*cos.(ϕ_temp)
        dt_u[:, :, 2] = μ * (uy_xxderiv + uy_yyderiv) + K * (ux_xyderiv + uy_yyderiv) + F_a*sin.(ϕ_temp)
    end
    
    return dt_u
end


function ϕ_transform(ϕ_temp::Array{Float64, 1}, active_contribution::Array{Float64, 2}, N::Int, F_a::Float64, l::Int, fft_plan::FFTW.FFTWPlan)

    norm_factor = 1/sqrt(l)
    
    ft_cosϕ = norm_factor * fft_plan * (F_a * cos.(ϕ_temp))
    ft_sinϕ = norm_factor * fft_plan * (F_a * sin.(ϕ_temp))

    active_contribution[1, :] = real.(ft_cosϕ)[1:N+1]
    active_contribution[2, :] = imag.(ft_cosϕ)[1:N+1]
    active_contribution[3, :] = real.(ft_sinϕ)[1:N+1]
    active_contribution[4, :] = imag.(ft_sinϕ)[1:N+1]

    return active_contribution
end

function non_lineartransform(dt_u::Array{Float64, 2}, ϕ_temp::Array{Float64, 1}, non_linearactive::Array{Float64, 2}, N::Int, ξ::Float64, l::Int, fft_plan::FFTW.FFTWPlan)
 
    norm_factor = 1/sqrt(l)
    
    ux_term = -ξ * sin.(ϕ_temp) .* dt_u[:, 1]

    uy_term = ξ * cos.(ϕ_temp) .* dt_u[:, 2]

    fftux = norm_factor * fft_plan * (ux_term)

    fftuy = norm_factor * fft_plan * (uy_term)

    non_linearactive[1, :] = real.(fftux)[1:N+1] + real.(fftuy)[1:N+1]
    non_linearactive[2, :] = imag.(fftux)[1:N+1] + imag.(fftuy)[1:N+1]

    return non_linearactive
end

function fouriertransform_active(ϕ_temp::Matrix{Float64}, active_contribution::Array{Float64, 3}, dt_u::Array{Float64, 3}, nonlinear_active::Array{Float64, 3}, N::Int, m::Int, F_a::Float64, ξ::Float64, l::Int, fft_plan::FFTW.FFTWPlan)

    for i in 1:m

        active_contribution[i , :, :] = ϕ_transform(ϕ_temp[:, i],    active_contribution[i, :, :], N, F_a, l, fft_plan)

        nonlinear_active[i, :, :] = non_lineartransform(dt_u[:, i, :], ϕ_temp[:, i], nonlinear_active[i , :, :], N, ξ, l, fft_plan)


    end

    return active_contribution, nonlinear_active

end


function runge_kuttastep(elcoefficientmatrix::Array{Float64, 3}, Kmatrix::Array{Float64, 3}, polarisationtioncoefficients::Array{Float64, 3}, Dmatrix::Array{Float64, 3}, N::Int, factor::Array{Float64, 1}, K::Float64, μ::Float64, λ::Float64, y::StepRangeLen, active_contribution::Array{Float64, 3}, nonlinear_active::Array{Float64, 3}, BC0::Bool, derivspolar::Array{Float64, 2}, derivselastic::Array{Float64, 2}, double_derivsel::Array{Float64, 2})

    if BC0 == true

        for i in 1:N+1

            Dmatrix[:, :, i] = polarisationtimederivs_zero(polarisationtioncoefficients[:, :, i], Dmatrix[:, :, i], nonlinear_active[:, :, i], factor[i], λ, y, derivspolar)

            Kmatrix[:, :, i] = elastictimederivs_zero(elcoefficientmatrix[:, :, i], Kmatrix[:, :, i], factor[i], K, μ, y, active_contribution[:, :, i], derivselastic, double_derivsel)

     end

    elseif BC0 == false
        for i in 1:N+1

            Dmatrix[:, :, i] = polarisationtimederivs_periodic(Dmatrix[:, :, i], nonlinear_active[:, :, i])

            Kmatrix[:, :, i] = elastictimederivs_periodic(elcoefficientmatrix[:, :, i], Kmatrix[:, :, i], factor[i], K, μ, y, active_contribution[:, :, i], derivselastic, double_derivsel)


        end
    end 


    return  Kmatrix, Dmatrix

end


function rk4_active(elcoefficientmatrix::Array{Float64, 3}, Kmatrix::Vector, polarisationtioncoefficients::Array{Float64, 3}, Dmatrix::Vector, dt::Float64, N::Int, factor::Array{Float64, 1}, K::Float64, μ::Float64, λ::Float64, y::StepRangeLen, active_contribution::Array{Float64, 3}, nonlinear_active::Array{Float64, 3}, zero::Bool, derivspolar::Array{Float64, 2}, derivselastic::Array{Float64, 2}, double_derivsel::Array{Float64, 2})

    Kmatrix[1], Dmatrix[1] = runge_kuttastep(elcoefficientmatrix, Kmatrix[1], polarisationtioncoefficients, Dmatrix[1], N, factor, K, μ, λ, y, active_contribution, nonlinear_active, zero, derivspolar, derivselastic, double_derivsel)


    Kmatrix[2], Dmatrix[2] = runge_kuttastep(elcoefficientmatrix + dt/2 * Kmatrix[1], Kmatrix[2],polarisationtioncoefficients + dt/2 * Dmatrix[1], Dmatrix[2], N, factor, K, μ, λ, y, active_contribution + dt/2 * Kmatrix[1], nonlinear_active + dt/2*Dmatrix[1], zero, derivspolar, derivselastic, double_derivsel)


    Kmatrix[3], Dmatrix[3] = runge_kuttastep(elcoefficientmatrix + dt/2 * Kmatrix[2], Kmatrix[3],polarisationtioncoefficients + dt/2 * Dmatrix[2], Dmatrix[3], N, factor, K, μ, λ, y, active_contribution + dt/2 * Kmatrix[2], nonlinear_active + dt/2*Dmatrix[2], zero, derivspolar, derivselastic, double_derivsel)


    Kmatrix[4], Dmatrix[4] = runge_kuttastep(elcoefficientmatrix + dt * Kmatrix[3], Kmatrix[4],polarisationtioncoefficients + dt * Dmatrix[3], Dmatrix[4], N, factor, K, μ, λ, y, active_contribution + dt * Kmatrix[3], nonlinear_active + dt*Dmatrix[3], zero, derivspolar, derivselastic, double_derivsel)


    elcoefficientmatrix += dt/6 * (Kmatrix[1] + 2 * Kmatrix[2] + 2 * Kmatrix[3] + Kmatrix[4])


    polarisationtioncoefficients += dt/6 * (Dmatrix[1] + 2 * Dmatrix[2] + 2 * Dmatrix[3] + Dmatrix[4])


    return elcoefficientmatrix, polarisationtioncoefficients

end


function run_activesystem(Loggingstep::Int, N::Int, dy::Float64, dx::Float64, Lx::Int, Ly::Int, F::Float64, T::Int, dt::Float64, K::Float64, μ::Float64, λ::Float64, F_a::Float64, ξ::Float64, BC0::Bool)

    x, y, m, l, ϕ, ux, uy, factors, elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix, ux_temp, uy_temp, ϕ_temp = initalize_active(Loggingstep, dy, dx, N, Lx, Ly, F, F_a, T, BC0)

    active_contribution = zeros(m, 4, N+1)

    nonlinear_active = zeros(m, 2, N+1)

    elastic_modes = zeros(N+1, T)

    polar_modes = zeros(N+1, T)

    dt_u = zeros(l, m, 2)

    derivselastic = Array{Float64}(undef, m, 4)
    double_derivsel = Array{Float64}(undef, m, 4)

    derivspolar = Array{Float64}(undef, m, 2)

    fft_plan = plan_fft(zeros(l), flags = FFTW.MEASURE)

    for i in 2:T-1

        dt_u = velocity(ux_temp, uy_temp, dt_u, x, y, μ, K, ϕ_temp, F_a)

       

        active_contribution, nonlinear_active = fouriertransform_active(ϕ_temp, active_contribution, dt_u, nonlinear_active, N, m, F_a, ξ, l, fft_plan)

        elcoefficientmatrix, polarisationtioncoefficients = rk4_active(elcoefficientmatrix, Kmatrix, polarisationtioncoefficients, Dmatrix, dt, N, factors, K, μ, λ, y, active_contribution, nonlinear_active, BC0, derivspolar, derivselastic, double_derivsel)


        ux_temp, uy_temp, ϕ_temp = calculateuxuy_active(elcoefficientmatrix, polarisationtioncoefficients, N, Lx, x, m, ux_temp, uy_temp, ϕ_temp, l, i, dt)



        if i % Loggingstep == 0

            idx = i÷Loggingstep +1

            @views begin

                ux[:, :, idx] = ux_temp

                uy[:, :, idx] = uy_temp

                ϕ[:, :, idx] = ϕ_temp



            end

        end

    end

    return ux, uy, ϕ, x, y, m, l

end

