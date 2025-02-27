using LinearAlgebra, CairoMakie, SparseArrays, Integrals, FFTW, Dierckx, GLMakie

include("plotting_utility.jl")

include("elasticsystem.jl")

include("polarization.jl")

include("Activesystem.jl")


@time ux, uy, ϕ, x, y, m = run_activesystem(10, 1., 21, 20, 1., 100000, 0.001, 5., 2., 0.1, 10., 1.)

animation_2d(ux, uy, x, y, 100000, 30)

animation_1d(ϕ, x, y, 100000, 30)



