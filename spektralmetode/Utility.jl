prefactor(L::Int, n::Int) = (2*π*n)/L


function prefactors(Lx::Int,N::Int)
    
    return [prefactor(Lx, n) for n in 0:N]

end

function divergence(ux, uy, x, y, T, K , μ,  energy) #må skrives om for å sørge for at energien blir riktig

    divergence_t = deepcopy(ux)

    curl_t = deepcopy(ux)

    if energy == true

        defenergy_t = deepcopy(ux)

    end

    for i in 1:T


        ux_spline = Spline2D(x, y, ux[:, :, i])

        uy_spline = Spline2D(x, y, uy[:, :, i])

        ux_deriv = derivative(ux_spline, x, y, nux = 1, nuy = 0)

        uy_deriv = derivative(uy_spline, x, y, nux = 0, nuy = 1)

        ux_yderiv = derivative(ux_spline, x, y, nux = 0, nuy = 1)

        uy_xderiv = derivative(uy_spline, x, y, nux = 1, nuy = 0)

        divergence_t[:, :, i] = ux_deriv + uy_deriv

        curl_t[:, :, i] = (uy_xderiv - ux_yderiv)

        


        if energy == true

            ux_yderiv = derivative(ux_spline, x, y, nux = 0, nuy = 1)

            uy_xderiv = derivative(uy_spline, x, y, nux = 1, nuy = 0)

            defenergy_t[:, :, i] = (1/2)* K * (ux_deriv + uy_deriv).^2 + μ *((1/2)*(ux_yderiv+uy_xderiv) - (ux_deriv + uy_deriv)).^2

        end


    end


    if energy == true

        return divergence_t, curl_t, defenergy_t

    else

        return divergence_t, curl_t

    end

end

function orderparameter(ϕ::Array{Float64, 3}, N_p, T)

    orderparameter_t = zeros(Float64, T)

    for t in 1:T 
       
        nx = cos.(ϕ[:, :, t])

        ny = sin.(ϕ[:, :, t])
    
        nx_avg = sum(nx) / N_p

        ny_avg = sum(ny) / N_p
       
        orderparameter_t[t] = sqrt(nx_avg^2 + ny_avg^2)
    end

    return orderparameter_t

end


function savedata(u_x::Array{Float64, 3}, u_y::Array{Float64, 3}, ϕ::Array{Float64, 3}, x::StepRangeLen, y::StepRangeLen, constants::Vector, flilename::String)

    save(flilename, Dict("u_x" => u_x, "u_y" => u_y, "ϕ" => ϕ, "x" => x, "y" => y, "constants" => constants))

end

function xisweep(N::Int, dy::Float64, dx::Float64, Lx::Int, Ly::Int, F::Float64, T::Int, dt::Float64, K::Float64, μ::Float64, λ::Float64, F_a::Float64, BC0::Bool)
    
    file_lock = ReentrantLock()

    
    Threads.@threads for i in 10:5:20
        t = @elapsed begin
            ux, uy, ϕ, x, y, m, l = run_activesystem(N, dy, dx, Lx, Ly, F, T, dt, K, μ, λ, F_a, Float64(i), BC0)

            flilename = "simulation2_xi$(i).jld2"
            constants = [N, m, l, Lx, Ly, F, T, dt, K, μ, λ, F_a, i]


            lock(file_lock) do
                savedata(ux, uy, ϕ, x, y, constants, flilename)
            end
        end
        @info "Thread $(Threads.threadid()) completed ξ=$i in $t seconds"
    end
end

function openandplot()

    for i in 10:5:20

        flilename = "simulation2_xi$(i).jld2"
        
        data = load(flilename)

        m = data["constants"][2]

        l = data["constants"][3]

        orderparameter_t = orderparameter(data["ϕ"], m*l, data["constants"][7])

        plot_orderparameter(orderparameter_t, Int(data["constants"][7]), "orderparameterrun2_xi$(i).png")

        animation_deformation(data["u_x"], data["u_y"], data["x"], data["y"], Int(data["constants"][7]), 30, "deformationrun2_xi$(i).mp4")

        animation_polarisation(data["ϕ"], data["x"], data["y"], Int(data["constants"][7]), 30, "polarisationrun2_xi$(i).mp4")

        animation_heatmap(data["ϕ"], data["x"], data["y"], Int(data["constants"][7]), 30, "heatrun2_xi$(i).mp4")

    end

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
