prefactor(L::Int, n::Int) = (2*π*n)/L


function prefactors(Lx::Int,N::Int)
    
    return [prefactor(Lx, n) for n in 0:N]

end



function orderparameter(ϕ::Array{Float64, 3}, T)

    orderparameter_t = zeros(Float64, T)

    points = size(ϕ, 1)*size(ϕ, 2)

    for t in 1:T 
       
        nx = cos.(ϕ[:, :, t])

        ny = sin.(ϕ[:, :, t])
    
        nx_avg = sum(nx) / points

        ny_avg = sum(ny) / points
       
        orderparameter_t[t] = sqrt(nx_avg^2 + ny_avg^2)
    end

    return orderparameter_t

end


function savedata(u_x::Array{Float64, 3}, u_y::Array{Float64, 3}, ϕ::Array{Float64, 3}, x::StepRangeLen, y::StepRangeLen, constants::Vector, flilename::String)

    save(flilename, Dict("u_x" => u_x, "u_y" => u_y, "ϕ" => ϕ, "x" => x, "y" => y, "constants" => constants))

end



function circularmean(ϕ::Array{Float64, 3}, t)

    ϕx = sum(sum(row) for row in cos.(ϕ[:, :, t]))

    ϕy = sum(sum(row) for row in sin.(ϕ[:, :, t]))

    return atan(ϕy, ϕx)

end

function circularmean_t(ϕ::Array{Float64, 3}, T)

    circularmean_t = zeros(Int(T))

    for t in 1:Int(T)

        circularmean_t[t] = circularmean(ϕ, t)

    end

    return circularmean_t

end

function getendangel(N::Int, dy::Float64, dx::Float64, Lx::Int, Ly::Int, F::Float64, T::Int, dt::Float64, K::Float64, μ::Float64, λ::Float64, F_a::Float64, ξ::Float64 ,BC0::Bool, filename::String)

 endangle = zeros(200)

 startangle = zeros(200)

    Threads.@threads for i in 1:200

        t = @elapsed ux, uy, ϕ, x, y, m, l, = run_activesystem(N, dy, dx, Lx, Ly, F, T, dt, K, μ, λ, F_a, ξ ,BC0)
        @info "Thread $(Threads.threadid()) iteration $i took $t seconds"

        θ = circularmean(ϕ, T)

        Θ = circularmean(ϕ, 2)

        startangle[i] = Θ

        # if θ < 0

        #     θ += 2*π

        # end

        endangle[i] = θ
        GC.gc()
    end
        

    save(filename, Dict("endangle" => endangle, "startangle" => startangle))z



end

function monitor_velocity(dt_u::Array{Float64, 3})
    # Maximum absolute value
    max_val = maximum(abs.(dt_u))
    
    # Average magnitude
    avg_val = mean(sqrt.(dt_u[:,:,1].^2 + dt_u[:,:,2].^2))
    
    return max_val, avg_val
end

function monitor_fourier_modes(elcoefficientmatrix::Array{Float64, 3}, polarisationtioncoefficients::Array{Float64, 3})
    N = size(elcoefficientmatrix, 3)  # N+1 modes
    max_elastic = zeros(N)
    max_polar = zeros(N)
    
    for k in 1:N
        
        # max_elastic[k] = sum(abs2.(view(elcoefficientmatrix, :, :, k)))

        # max_polar[k] = sum(abs2.(view(polarisationtioncoefficients, :, :, k)))
        
        max_elastic[k] = sqrt(mean(abs2.(view(elcoefficientmatrix, :, :, k))))
        max_polar[k] = sqrt(mean(abs2.(view(polarisationtioncoefficients, :, :, k))))
    end

    # max_elastic = max_elastic ./ sum(max_elastic[1])
    # max_polar = max_polar ./ sum(max_polar[1])
    
    return max_elastic, max_polar
end

function velocity_calc(ux::Array{Float64, 3}, uy::Array{Float64, 3}, x::StepRangeLen, y::StepRangeLen, μ::Float64, K::Float64, ϕ::Array{Float64, 3}, F_a::Float64)
    
    dt_ux =  deepcopy(ux)

    dt_uy = deepcopy(uy)

    for i in 1:size(ux, 3)
    
        @views begin
            ux_periodic = hcat(ux[:, end, i], ux[:, :, i], ux[:, 1, i])
            uy_periodic = hcat(uy[:, end, i], uy[:, :, i], uy[:, 1, i])
        end
        
        y_periodic = vcat(y[1] - 1, y, y[end] + 1)
        
    
        @views begin
            ux_spline = Spline2D(x, y_periodic, ux_periodic)
            uy_spline = Spline2D(x, y_periodic, uy_periodic)
            
    
            ux_xxderiv = similar(ux, size(ux, 1), size(ux, 2))
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
            dt_ux[:, :, i] = μ * (ux_xxderiv + ux_yyderiv) + K * (ux_xxderiv + uy_xyderiv) + F_a*cos.(ϕ[:, :, i])
            dt_uy[:, :, i] = μ * (uy_xxderiv + uy_yyderiv) + K * (ux_xyderiv + uy_yyderiv) + F_a*sin.(ϕ[:, :, i])
        end
    end
    return dt_ux, dt_uy
end

function avg_velocity(ux::Array{Float64, 3}, uy::Array{Float64, 3}, x::StepRangeLen, y::StepRangeLen, μ::Float64, K::Float64, ϕ::Array{Float64, 3}, F_a::Float64)

    avgvelocity = zeros(Float64, size(ux, 3))

    dtux, dtuy = velocity_calc(ux, uy, x, y, μ, K, ϕ, F_a)

    for i in 1:size(ux, 3)

        avgvelocity[i] = mean(sqrt.(dtux[:, :, i].^2 + dtuy[:, :, i].^2))

    end

    return avgvelocity

end


function ϕ_dt(ϕ::Array{Float64, 3}, T::Int, Loggingstep::Int, l::Int, m::Int, dt::Float64)

    t = T ÷ Loggingstep

    ϕ_dt = zeros(Float64, t)

    for i in 2:t

        ϕ_dt[i] = (ϕ[floor(Int, l/2), floor(Int, m/2), i] - ϕ[floor(Int, l/2), floor(Int, m/2), i-1])/dt

    end

    return ϕ_dt

end
