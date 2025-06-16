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

    record(fig2, filename, 1:5:numiter; framerate = framerate) do i
     iplot1d[] = i

     sleep(0.05)

    end

end




function plot_orderparameter(orderparameter_t::Vector{Float64}, T::Int, dt::Float64,name::String)

    n_steps = length(orderparameter_t)

    time = range(0, T*dt, length = n_steps)

    fig4 = Figure(size = (800, 800))

    ax4 = Axis(fig4[1, 1], xlabel = "time", ylabel = "Π") 

    GLMakie.lines!(time, orderparameter_t)

    GLMakie.save(name, fig4)

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
    ax1.xlabel = "t"
    ax1.ylabel = "Π"

    # Velocity statistics
    ax2 = Axis(fig[1,2], title="Velocity Statistics", yscale=log10)
    lines!(ax2, 1:T, max_velocity, label="Max |dt_u|")
    lines!(ax2, 1:T, avg_velocity, label="Avg |dt_u|")
    axislegend(ax2)
    ax2.xlabel = "t"
    ax2.ylabel = "Velocity"

    # Elastic modes
    ax3 = Axis(fig[2,1], title="Elastic Modes", yscale=log10)
    for k in 1:size(elastic_modes, 1)
        lines!(ax3, 1:T, elastic_modes[k,:], label="k=$(k-1)")
    end
    axislegend(ax3)
    ax3.xlabel = "t"
    ax3.ylabel = "Mode Amplitude"

    # Polarization modes
    ax4 = Axis(fig[2,2], title="Polarization Modes", yscale=log10)
    for k in 1:size(polar_modes, 1)
        lines!(ax4, 1:T, polar_modes[k,:], label="k=$(k-1)")
    end
    axislegend(ax4)
    ax4.xlabel = "t"
    ax4.ylabel = "Mode Amplitude"

    save(filename, fig)
    return fig
end

function plot_system_analysis2(filename::String, figname::String)

    data = load(filename)

    ux = data["u_x"]

    uy = data["u_y"]

    ϕ = data["ϕ"]

    x = data["x"]

    y = data["y"]

    Loggingstep = Int(data["constants"][1])

    T = Int(data["constants"][8])

    dt = data["constants"][9]

    K = data["constants"][10]

    μ = data["constants"][11]

    F_a = data["constants"][13]

    N_y = Int(data["constants"][6])

    N_x = Int(data["constants"][5])

    orderparameter_t = orderparameter(ϕ, Int(T/Loggingstep))

    avg_velocity_t = avg_velocity(ux, uy, x, y, μ, K, ϕ, F_a)

    t = (1:Loggingstep:T) * dt

    fig = CairoMakie.Figure(size = (1000, 1000))

    ax1 = CairoMakie.Axis(fig[1, 1], aspect = 1)
    Label(fig[1, 1, TopLeft()], "a)", padding = (0, 5, 5, 0))
    ax1.xlabel = "t"
    ax1.ylabel = "Π"  
    CairoMakie.lines!(ax1, t, orderparameter_t)

    ax2 = CairoMakie.Axis(fig[1, 2], aspect = 1)
    Label(fig[1, 2, TopLeft()], "b)", padding = (0, 5, 5, 0))
    ax2.xlabel = "t"
    ax2.ylabel = "Average velocity"
    CairoMakie.lines!(ax2, t, avg_velocity_t)
    

    CairoMakie.save(figname, fig)

end

function plot_phasesnapshots_hsv(filename::String, figname::String)

    data = load(filename)

    ϕ = data["ϕ"]

    x = data["x"]

    y = data["y"]

    ϕ = mod2pi.(ϕ)

    fig = CairoMakie.Figure(size= (1000, 750), aspect = 1)

    ax1 = CairoMakie.Axis(fig[1, 1], aspect = 1)
    Label(fig[1, 1, TopLeft()], "a)", padding = (0, 5, 5, 0))
    ax1.xlabel = "x"
    ax1.ylabel = "y"
    CairoMakie.arrows!(x, y, cos.(ϕ[:,:,1]), sin.(ϕ[:,:,1]), arrowsize = 5)

    ax2 = CairoMakie.Axis(fig[1, 2], aspect = 1)
    Label(fig[1, 2, TopLeft()], "b)", padding = (0, 5, 5, 0))
    ax2.xlabel = "x"
    ax2.ylabel = "y"
    CairoMakie.arrows!(x, y, cos.(ϕ[:,:,1500]), sin.(ϕ[:,:,1500]), arrowsize = 5)

    ax3 = CairoMakie.Axis(fig[1, 3], aspect = 1)
    Label(fig[1, 3, TopLeft()], "c)", padding = (0, 5, 5, 0))
    ax3.xlabel = "x"
    ax3.ylabel = "y"
    CairoMakie.arrows!(x, y, cos.(ϕ[:,:,end]), sin.(ϕ[:,:,end]), arrowsize = 5)
    

    ax4 = CairoMakie.Axis(fig[2, 1], aspect = 1)
    Label(fig[2, 1, TopLeft()], "e)", padding = (0, 5, 5, 0))
    ax4.xlabel = "x"
    ax4.ylabel = "y"
    CairoMakie.heatmap!(x, y, ϕ[:,:,1], colorrange = (π/2 *0.9 -0.3, π/2 *0.9 +0.3), colormap = :hsv)
    
    ax5 = CairoMakie.Axis(fig[2, 2], aspect = 1)
    Label(fig[2, 2, TopLeft()], "f)", padding = (0, 5, 5, 0))
    ax5.xlabel = "x"
    ax5.ylabel = "y"
    CairoMakie.heatmap!(x, y, ϕ[:,:,1500], colorrange = (π/2 *0.9 -0.3, π/2 *0.9 +0.3), colormap = :hsv)

    ax6 = CairoMakie.Axis(fig[2, 3], aspect = 1)
    Label(fig[2, 3, TopLeft()], "g)", padding = (0, 5, 5, 0))
    ax6.xlabel = "x"
    ax6.ylabel = "y"
    CairoMakie.heatmap!(x, y, ϕ[:,:,end], colorrange = (π/2 *0.9 -0.3, π/2 *0.9 +0.3), colormap = :hsv)

    CairoMakie.Colorbar(fig[2, 4], limits = (π/2 *0.9 -0.3, π/2 *0.9 +0.3), colormap = :hsv,
    ticks=([π/2 *0.9- 0.3, π/2 *0.9+0.3], ["π/2 *0.9 -0.3", "π/2 *0.9 +0.3"]))

    CairoMakie.save(figname, fig)

end



    
function ϕ_amp(figname::String)

    data = load("stick_K1.jld2")

    ϕ = data["ϕ"]

    constants = data["constants"]

    T = Int(constants[8])

    Loggingstep = Int(constants[1])

    dt = constants[9]

    data2 = load("stick_K2.jld2")

    ϕ2 = data2["ϕ"]

    data3 = load("stick_K3.jld2")

    ϕ3 = data3["ϕ"]

    data4 = load("stick_K4.jld2")

    ϕ4 = data4["ϕ"]

    data5 = load("stick_K5.jld2")

    ϕ5 = data5["ϕ"]

    fig = CairoMakie.Figure(size = (1000, 1000/2))

    ax1 = CairoMakie.Axis(fig[1, 1])
    Label(fig[1, 1, TopLeft()], "a)", padding = (0, 5, 5, 0))
    ax1.xlabel = "time"
    ax1.ylabel = "ϕ"

    ax2 = CairoMakie.Axis(fig[1, 2])
    Label(fig[1, 2, TopLeft()], "b)", padding = (0, 5, 5, 0))
    ax2.xlabel = "time"
    ax2.ylabel = "ϕ"

    ax3 = CairoMakie.Axis(fig[1, 3])
    Label(fig[1, 3, TopLeft()], "c)", padding = (0, 5, 5, 0))
    ax3.xlabel = "time"
    ax3.ylabel = "ϕ"

    ax4 = CairoMakie.Axis(fig[2, 1])
    Label(fig[2, 1, TopLeft()], "d)", padding = (0, 5, 5, 0))
    ax4.xlabel = "time"
    ax4.ylabel = "ϕ"

    ax5 = CairoMakie.Axis(fig[2, 2])
    Label(fig[2, 2, TopLeft()], "e)", padding = (0, 5, 5, 0))
    ax5.xlabel = "time"
    ax5.ylabel = "ϕ"

    CairoMakie.lines!(ax1, (1:T/Loggingstep)*dt , ϕ[Int(floor(21/2)), Int(floor(20/2)), :])
    CairoMakie.lines!(ax2, (1:T/Loggingstep)*dt , ϕ2[Int(floor(21/2)), Int(floor(20/2)), :])
    CairoMakie.lines!(ax3, (1:T/Loggingstep)*dt , ϕ3[Int(floor(21/2)), Int(floor(20/2)), :])
    CairoMakie.lines!(ax4, (1:T/Loggingstep)*dt , ϕ4[Int(floor(21/2)), Int(floor(20/2)), :])
    CairoMakie.lines!(ax5, (1:T/Loggingstep)*dt , ϕ5[Int(floor(21/2)), Int(floor(20/2)), :])
    
    CairoMakie.save(figname, fig)


end

function plot_orderparameters(figname::String)

    data = load("stick_K1.jld2")

    ϕ = data["ϕ"]

    constants = data["constants"]

    T = Int(constants[8])

    Loggingstep = Int(constants[1])

    dt = constants[9]

    data2 = load("stick_K2.jld2")

    ϕ2 = data2["ϕ"]

    data3 = load("stick_K3.jld2")

    ϕ3 = data3["ϕ"]

    data4 = load("stick_K4.jld2")

    ϕ4 = data4["ϕ"]

    data5 = load("stick_K5.jld2")

    ϕ5 = data5["ϕ"]

    orderparameter_t = orderparameter(ϕ, Int(T/Loggingstep))

    orderparameter_t2 = orderparameter(ϕ2, Int(T/Loggingstep))

    orderparameter_t3 = orderparameter(ϕ3, Int(T/Loggingstep))

    orderparameter_t4 = orderparameter(ϕ4, Int(T/Loggingstep))

    orderparameter_t5 = orderparameter(ϕ5, Int(T/Loggingstep))

    fig = CairoMakie.Figure(size = (1000, 1000/2))

    ax1 = CairoMakie.Axis(fig[1, 1])
    Label(fig[1, 1, TopLeft()], "a)", padding = (0, 5, 5, 0))
    ax1.xlabel = "time"
    ax1.ylabel = "Π"

    ax2 = CairoMakie.Axis(fig[1, 2])
    Label(fig[1, 2, TopLeft()], "b)", padding = (0, 5, 5, 0))
    ax2.xlabel = "time"
    ax2.ylabel = "Π"

    ax3 = CairoMakie.Axis(fig[1, 3])
    Label(fig[1, 3, TopLeft()], "c)", padding = (0, 5, 5, 0))
    ax3.xlabel = "time"
    ax3.ylabel = "Π"

    ax4 = CairoMakie.Axis(fig[2, 1])
    Label(fig[2, 1, TopLeft()], "d)", padding = (0, 5, 5, 0))
    ax4.xlabel = "time"
    ax4.ylabel = "Π"

    ax5 = CairoMakie.Axis(fig[2, 2])
    Label(fig[2, 2, TopLeft()], "e)", padding = (0, 5, 5, 0))
    ax5.xlabel = "time"
    ax5.ylabel = "Π"

    CairoMakie.lines!(ax1, (1:T/Loggingstep)*dt , orderparameter_t)
    CairoMakie.lines!(ax2, (1:T/Loggingstep)*dt , orderparameter_t2)
    CairoMakie.lines!(ax3, (1:T/Loggingstep)*dt , orderparameter_t3)
    CairoMakie.lines!(ax4, (1:T/Loggingstep)*dt , orderparameter_t4)
    CairoMakie.lines!(ax5, (1:T/Loggingstep)*dt , orderparameter_t5)
    
    CairoMakie.save(figname, fig)

end

    

    
function plot_ϕ_dt(figname::String)

    data = load("stick_K1.jld2")

    ϕ = data["ϕ"]

    constants = data["constants"]

    T = Int(constants[8])

    Loggingstep = Int(constants[1])

    dt = constants[9]

    data2 = load("stick_K2.jld2")

    ϕ2 = data2["ϕ"]

    data3 = load("stick_K3.jld2")

    ϕ3 = data3["ϕ"]

    data4 = load("stick_K4.jld2")

    ϕ4 = data4["ϕ"]

    data5 = load("stick_K5.jld2")

    ϕ5 = data5["ϕ"]

    ϕ_dt1 = ϕ_dt(ϕ, T, Loggingstep, 21, 20, dt)

    ϕ_dt2 = ϕ_dt(ϕ2, T, Loggingstep, 21, 20, dt)

    ϕ_dt3 = ϕ_dt(ϕ3, T, Loggingstep, 21, 20, dt)

    ϕ_dt4 = ϕ_dt(ϕ4, T, Loggingstep, 21, 20, dt)

    ϕ_dt5 = ϕ_dt(ϕ5, T, Loggingstep, 21, 20, dt)

    
    fig = CairoMakie.Figure(size = (1000, 1000/2))

    ax1 = CairoMakie.Axis(fig[1, 1])
    Label(fig[1, 1, TopLeft()], "a)", padding = (0, 5, 5, 0))
    ax1.xlabel = "time"
    ax1.ylabel = "ϕ time derivative"
    CairoMakie.lines!(ax1, (1:T/Loggingstep)*dt , ϕ_dt1)
    
    ax2 = CairoMakie.Axis(fig[1, 2])
    Label(fig[1, 2, TopLeft()], "b)", padding = (0, 5, 5, 0))
    ax2.xlabel = "time"
    ax2.ylabel = "ϕ time derivative"
    CairoMakie.lines!(ax2, (1:T/Loggingstep)*dt , ϕ_dt2)
    
    ax3 = CairoMakie.Axis(fig[1, 3])
    Label(fig[1, 3, TopLeft()], "c)", padding = (0, 5, 5, 0))
    ax3.xlabel = "time"
    ax3.ylabel = "ϕ time derivative"
    CairoMakie.lines!(ax3, (1:T/Loggingstep)*dt , ϕ_dt3)
    
    ax4 = CairoMakie.Axis(fig[2, 1])
    Label(fig[2, 1, TopLeft()], "d)", padding = (0, 5, 5, 0))
    ax4.xlabel = "time"
    ax4.ylabel = "ϕ time derivative"
    CairoMakie.lines!(ax4, (1:T/Loggingstep)*dt , ϕ_dt4)
    
    ax5 = CairoMakie.Axis(fig[2, 2])
    Label(fig[2, 2, TopLeft()], "e)", padding = (0, 5, 5, 0))
    ax5.xlabel = "time"
    ax5.ylabel = "ϕ time derivative"
    CairoMakie.lines!(ax5, (1:T/Loggingstep)*dt , ϕ_dt5)
    
    CairoMakie.save(figname, fig)
    
end
    

    

    