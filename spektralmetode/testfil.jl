using LinearAlgebra, CairoMakie, SparseArrays, Integrals, FFTW, Dierckx, GLMakie, FileIO, Statistics

include("plottingfunctions.jl")

include("Utility.jl")

#include("elasticsystem.jl")

#include("polarization.jl")

include("Activesystem.jl")


xi = [1., 2., 3., 4., 5.]

for i in eachindex(xi)

    @time ux1, uy1, ϕ1, x1, y1, m1, l1 = run_activesystem(10, 10,  1., 1., 21, 20, 1., 200000, 0.01, 1., 1., 1., 0.2,  1., true)

    savedata(ux1, uy1, ϕ1, x1, y1, [10, 10, 1., 1., 21, 20, 1., 200000, 0.01, 1., 1., 1., 0.2,  1.], "stick_xi$i.jld2")

    orderparameter_t = orderparameter(ϕ1, 20000)

    plot_orderparameter(orderparameter_t, 20000, 0.01, "orderparameter_stick_xi$i.png")

# animation_deformation(ux, uy, x, y, 10000, 30, "test_def2.mp4")

# animation_polarisation(ϕ, x, y, 10000, 30, "test_polar2.mp4")

# # data = load("RunKsweep1/stick_K1.jld2")

# # ϕ = data["ϕ"]

# # ux = data["u_x"]

# # uy = data["u_y"]

# # x = data["x"]

# # y = data["y"]

# # constants = data["constants"]

# # μ = constants[11]

# # K = constants[10]

# # F_a = constants[13]

# # logstep = constants[1]

# # T = constants[8]

# # dt = constants[9]

# # f_el = elastic_force(ux, uy, x, y, μ, K, ϕ, F_a)

# # circularmean1= circularmean_t(ϕ, Int(T/logstep))   

# # fig = CairoMakie.Figure(size= (1000, 1000), aspect = 1, fontsize = 20)

# # ax = CairoMakie.Axis(fig[1, 1])

# # CairoMakie.lines!(ax, (1:logstep:T)*dt, f_el)
# # CairoMakie.lines!(ax, (1:logstep:T)*dt, F_a*ones(Int(T/logstep)), color = :red)

# # CairoMakie.save("test.png", fig)


# # test = f_el[10000:end]

# # println(maximum(test))
# # println(minimum(test))
# # # println(mean(test))

CairoMakie.activate!()

CairoMakie.set_theme!(theme_latexfonts())
plot_limitcycle("limitcycle.pdf")

#plot_parabola("parabola.pdf")