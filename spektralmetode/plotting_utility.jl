
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






function animation_2d(u_x, u_y, ϕ, x, y, numiter, framerate)

    fig = Figure(size = (800, 800))

    ax1 = Axis(fig[1, 1])

    iplot = Observable(1)

    ux = @lift(u_x[:, :, $iplot])

    uy = @lift(u_y[:, :, $iplot])

    # px = @lift(cos.(ϕ[:, :, $iplot]))

    # py = @lift(sin.(ϕ[:, :, $iplot]))

    GLMakie.arrows!(x, y, ux, uy, arrowsize = 10)
    #GLMakie.arrows!(x, y, px, py, arrowsize = 10, color = :red)

    display(fig)



    # for i in 1:numiter

    #     iplot[] = i

    #     sleep(0.05)

    # end 



    record(fig, "animation4.mp4", 1:10:numiter; framerate = framerate) do i
      iplot[] = i

      sleep(0.05)

    end


end


function animation_1d(u, x, y, numiter, framerate)

    fig2 = Figure(size = (800, 800))

    ax2 = Axis(fig2[1, 1], aspect=  1)

    iplot1d = Observable(1)

    min = minimum(u)

    max = maximum(u)

    uplot = @lift(u[:, :, $iplot1d])

    GLMakie.heatmap!(x, y, uplot, colorrange = (min, max), colormap = :jet1)

    GLMakie.Colorbar(fig2[1, 2], limits = (min, max), colormap = :jet1)

    display(fig2)

    # for i in 1:numiter

    #     iplot1d[] = i

    #     sleep(0.05)
    
    # end

    record(fig2, "animation11.mp4", 1:10:numiter; framerate = framerate) do i
     iplot1d[] = i

     sleep(0.05)

    end

end


function animation_arrows2(u, x, y, numiter, framerate)

    fig3 = Figure(size = (800, 800))

    ax3 = Axis(fig3[1, 1], aspect=  1)

    iplot1d = Observable(1)


    uplot3 = @lift(u[:, :, $iplot1d])

    GLMakie.arrows!(x, y, uplot3, arrowsize = 10)


    display(fig3)

    # for i in 1:numiter

    #     iplot1d[] = i

    #     sleep(0.05)
    
    # end

    record(fig3, "animation4.mp4", 1:10:numiter; framerate = framerate) do i
     iplot1d[] = i

     sleep(0.05)

    end

end

function savedata(u_x::Array{Float64, 3}, u_y::Array{Float64, 3}, ϕ::Array{Float64, 3}, x::StepRangeLen, y::StepRangeLen, constants::Vector, flilename::String)

    jdlsave(flilename, u_x, u_y, ϕ, x, y, constants)

end