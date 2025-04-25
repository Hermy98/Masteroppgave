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

function plot_system_analysis(orderparameter::Vector{Float64}, max_velocity::Vector{Float64}, avg_velocity::Vector{Float64}, elastic_modes::Matrix{Float64}, polar_modes::Matrix{Float64}, T::Int, filename::String)
    fig = Figure(size=(1200, 1000))

    # Order parameter
    ax1 = Axis(fig[1,1], title="Order Parameter")
    lines!(ax1, 1:T, orderparameter)
    ax1.xlabel = "Time"
    ax1.ylabel = "Order Parameter"

    # Velocity statistics
    ax2 = Axis(fig[1,2], title="Velocity Statistics", yscale=log10)
    lines!(ax2, 1:T, max_velocity, label="Max |dt_u|")
    lines!(ax2, 1:T, avg_velocity, label="Avg |dt_u|")
    axislegend(ax2)
    ax2.xlabel = "Time"
    ax2.ylabel = "Velocity"

    # Elastic modes
    ax3 = Axis(fig[2,1], title="Elastic Modes", yscale=log10)
    for k in 1:size(elastic_modes, 1)
        lines!(ax3, 1:T, elastic_modes[k,:], label="k=$(k-1)")
    end
    axislegend(ax3)
    ax3.xlabel = "Time"
    ax3.ylabel = "Mode Amplitude"

    # Polarization modes
    ax4 = Axis(fig[2,2], title="Polarization Modes", yscale=log10)
    for k in 1:size(polar_modes, 1)
        lines!(ax4, 1:T, polar_modes[k,:], label="k=$(k-1)")
    end
    axislegend(ax4)
    ax4.xlabel = "Time"
    ax4.ylabel = "Mode Amplitude"

    save(filename, fig)
    return fig
end
