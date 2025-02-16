using LinearAlgebra, CairoMakie, SparseArrays, Integrals, FFTW, Dierckx, GLMakie

include("plotting_utility.jl")

include("elasticsystem.jl")

include("polarization.jl")

include("Activesystem.jl")



ux, uy, ϕ, x, y = run_activesystem(10, 1., 40, 20, 10., 10000, 0.001, 5., 2., 10., 1.)





animation_2d(ux, uy, x, y, 10000, 30)


