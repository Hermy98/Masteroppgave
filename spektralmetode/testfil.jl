using LinearAlgebra, CairoMakie, SparseArrays, Integrals, FFTW, Dierckx, GLMakie

include("plotting_utility.jl")

include("elasticsystem.jl")

include("polarization.jl")

include("Activesystem.jl")




@time ux, uy, ϕ, x, y = run_activesystem(10, 1., 40, 20, 10., 100000, 0.001, 0.1, 5., 2., 0.1)

animation_2d(ux, uy, x, y, 100000, 30)