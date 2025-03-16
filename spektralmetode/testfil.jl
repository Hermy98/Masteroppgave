using LinearAlgebra, CairoMakie, SparseArrays, Integrals, FFTW, Dierckx, GLMakie, JLD2

include("plotting_utility.jl")

#include("elasticsystem.jl")

#include("polarization.jl")

include("Activesystem.jl")


@time ux, uy, ϕ, x, y, m = run_activesystem(10, 1., 1.,  21, 40, 1., 60000, 0.001, 5., 4., 1., 0.75, 4., true)


animation_2d(ux, uy, ϕ, x, y, 60000, 30)


divu, curl_u = divergence(ux, uy, x, y, 60000, 1., 1., false)

animation_1d(curl_u, x, y, 60000, 30)

#animation_1d(ϕ, x, y, 10000, 30)

# fig  = Figure(size = (800, 800))

# ax1 = Axis(fig[1, 1])

# ax2 = Axis(fig[1, 2])

# arrows!(ax1, x, y, ux[:, :, 2], uy[:, :, 2])

# arrows!(ax2, x, y,ux[:, :, 3], uy[:, :, 3], color = :red)

# display(fig)



