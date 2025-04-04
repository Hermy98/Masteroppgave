

function animation_deformation(u_x, u_y, x, y, numiter, framerate, filename) 

    fig = Figure(size = (800, 800))

    ax1 = Axis(fig[1, 1])

    iplot = Observable(1)

    ux = @lift(u_x[:, :, $iplot])

    uy = @lift(u_y[:, :, $iplot])

    GLMakie.arrows!(x, y, ux, uy, arrowsize = 10)

    # for i in 1:numiter

    #     iplot[] = i

    #     sleep(0.05)

    # end 



    record(fig, filename, 1:10:numiter; framerate = framerate) do i
      iplot[] = i

      sleep(0.05)

    end


end

function animation_polarisation(ϕ, x, y, numiter, framerate, filename)

    fig = Figure(size = (800, 800))

    ax = Axis(fig[1, 1])

    iplot = Observable(1)

    px = @lift(cos.(ϕ[:, :, $iplot]))

    py = @lift(sin.(ϕ[:, :, $iplot]))

    GLMakie.arrows!(x, y, px, py, arrowsize = 10, color = :red)


    record(fig, filename, 1:10:numiter; framerate = framerate) do i
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


    # for i in 1:numiter

    #     iplot1d[] = i

    #     sleep(0.05)
    
    # end

    record(fig2, "animation.mp4", 1:10:numiter; framerate = framerate) do i
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


    # for i in 1:numiter

    #     iplot1d[] = i

    #     sleep(0.05)
    
    # end

    record(fig3, "animation1.mp4", 1:10:numiter; framerate = framerate) do i
     iplot1d[] = i

     sleep(0.05)

    end

end

function plot_orderparameter(orderparameter_t::Vector{Float64}, T::Int, name::String)

    fig4 = Figure(size = (800, 800))

    ax4 = Axis(fig4[1, 1])

    GLMakie.lines!(1:T, orderparameter_t)

    save(name, fig4)

end


