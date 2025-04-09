

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



    record(fig, filename, 1:20:numiter; framerate = framerate) do i
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


    record(fig, filename, 1:20:numiter; framerate = framerate) do i
      iplot[] = i

      sleep(0.05)

    end

end


function animation_heatmap(u, x, y, numiter, framerate, filename)

    fig2 = Figure(size = (800, 800))

    ax2 = Axis(fig2[1, 1], aspect=  1)

    iplot1d = Observable(1)

    uplot = @lift(u[:, :, $iplot1d])

    GLMakie.heatmap!(x, y, uplot, colorrange = (-π, π), colormap = :hsv)

    GLMakie.Colorbar(fig2[1, 2], limits = (-π, π), colormap = :hsv)


    # for i in 1:numiter

    #     iplot1d[] = i

    #     sleep(0.05)
    
    # end

    record(fig2, filename, 1:10:numiter; framerate = framerate) do i
     iplot1d[] = i

     sleep(0.05)

    end

end




function plot_orderparameter(orderparameter_t::Vector{Float64}, T::Int, name::String)

    fig4 = Figure(size = (800, 800))

    ax4 = Axis(fig4[1, 1], xlabel = "timestep", ylabel = "Φ") 

    GLMakie.lines!(1:T, orderparameter_t)

    save(name, fig4)

end

function angulardistribution(figname::String, filename::String)

    data = load(filename)

    f = Figure(size = (800, 800))

    ax1 = Axis(f[1, 1], aspect = 1)

    ax2 = Axis(f[1, 2])

    GLMakie.scatter!(ax1,cos.(data["endangle"]), sin.(data["endangle"]), color = :blue)

    GLMakie.hist!(ax2, data["endangle"], bins = 10)

    save(figname, f)

end


function plot_heatmap(u, x, y, name)

    fig3 = Figure(size = (800, 800))

    ax3 = Axis(fig3[1, 1])

    GLMakie.heatmap!(x, y, u, colorrange = (-π, π), colormap = :hsv)

    GLMakie.Colorbar(fig3[1, 2], limits = (-π, π), colormap = :hsv)

    save(name, fig3)

end




