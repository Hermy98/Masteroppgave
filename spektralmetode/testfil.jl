using LinearAlgebra, CairoMakie, SparseArrays, Integrals, FFTW, Dierckx, GLMakie

include("plotting_utility.jl")

#include("elasticsystem.jl")

#include("polarization.jl")

include("Activesystem.jl")


@time ux, uy, ϕ, x, y, m = run_activesystem(10, 1., 1.,  21, 40, 1., 50000, 0.001, 5., 2., 1., 0.75, 4.)


animation_2d(ux, uy, ϕ, x, y, 50000, 30)

#animation_1d(ϕ, x, y, 10000, 30)

# fig  = Figure(size = (800, 800))

# ax1 = Axis(fig[1, 1])

# ax2 = Axis(fig[1, 2])

# arrows!(ax1, x, y, ux[:, :, 2], uy[:, :, 2])

# arrows!(ax2, x, y,ux[:, :, 3], uy[:, :, 3], color = :red)

# display(fig)



