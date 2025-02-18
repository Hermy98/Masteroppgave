using LinearAlgebra, CairoMakie, SparseArrays, Integrals, FFTW, Dierckx, GLMakie

include("plotting_utility.jl")

include("elasticsystem.jl")

include("polarization.jl")

include("Activesystem.jl")



@time x, y, ϕ, factors, m = runsystem(10, 1., 40, 20, 1., 20000, 0.001, 2.)

#animation_1d(ϕ ,x, y, 20000, 30)

animation_2d(cos.(ϕ), sin.(ϕ), x, y, 20000, 30)