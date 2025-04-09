
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

    orderparameter_t = zeros(Int(T))

    for t in 1:Int(T) 

        Px = sum(sum(row) for row in cos.(ϕ[:, :, t]))

        Py = sum(sum(row) for row in sin.(ϕ[:, :, t]))

        orderparameter_t[t] = 1/N_p * sqrt(Px^2 + Py^2)
    end

    
    return orderparameter_t
      
end


function savedata(u_x::Array{Float64, 3}, u_y::Array{Float64, 3}, ϕ::Array{Float64, 3}, x::StepRangeLen, y::StepRangeLen, constants::Vector, flilename::String)

    save(flilename, Dict("u_x" => u_x, "u_y" => u_y, "ϕ" => ϕ, "x" => x, "y" => y, "constants" => constants))

end

function xisweep(N::Int, dy::Float64, dx::Float64, Lx::Int, Ly::Int, F::Float64, T::Int, dt::Float64, K::Float64, μ::Float64, λ::Float64, F_a::Float64, BC0::Bool)
    # Create a lock for file operations
    file_lock = ReentrantLock()

    # Run simulations in parallel
    Threads.@threads for i in 0:5:30
        t = @elapsed begin
            ux, uy, ϕ, x, y, m, l = run_activesystem(N, dy, dx, Lx, Ly, F, T, dt, K, μ, λ, F_a, Float64(i), BC0)

            flilename = "simulation1_xi$(i).jld2"
            constants = [N, m, l, Lx, Ly, F, T, dt, K, μ, λ, F_a, i]

            # Lock file operations to prevent race conditions
            lock(file_lock) do
                savedata(ux, uy, ϕ, x, y, constants, flilename)
            end
        end
        @info "Thread $(Threads.threadid()) completed ξ=$i in $t seconds"
    end
end

function openandplot()

    for i in 0:5:30

        flilename = "simulation1_xi$(i).jld2"
        
        data = load(flilename)

        m = data["constants"][2]

        l = data["constants"][3]

        orderparameter_t = orderparameter(data["ϕ"], m*l, data["constants"][7])

        plot_orderparameter(orderparameter_t, Int(data["constants"][7]), "orderparameterrun1_xi$(i).png")

        animation_deformation(data["u_x"], data["u_y"], data["x"], data["y"], Int(data["constants"][7]), 30, "deformationrun1_xi$(i).mp4")

        animation_polarisation(data["ϕ"], data["x"], data["y"], Int(data["constants"][7]), 30, "polarisationrun1_xi$(i).mp4")

        animation_heatmap(data["ϕ"], data["x"], data["y"], Int(data["constants"][7]), 30, "heatrun1_xi$(i).mp4")

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

    Threads.@threads for i in 1:200

        t = @elapsed ux, uy, ϕ, x, y, m, l, = run_activesystem(N, dy, dx, Lx, Ly, F, T, dt, K, μ, λ, F_a, ξ ,BC0)
        @info "Thread $(Threads.threadid()) iteration $i took $t seconds"

        θ = circularmean(ϕ, T)

        # if θ < 0

        #     θ += 2*π

        # end

        endangle[i] = θ
    end
        

    save(filename, Dict("endangle" => endangle))

end
