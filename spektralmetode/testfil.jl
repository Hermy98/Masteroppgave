using LinearAlgebra, CairoMakie, SparseArrays, Integrals, FFTW, Dierckx, GLMakie, FileIO

include("plottingfunctions.jl")

include("Utility.jl")

#include("elasticsystem.jl")

#include("polarization.jl")

include("Activesystem.jl")


# @time ux, uy, ϕ, x, y, m , l = run_activesystem(10, 1., 1.,  21, 40, 1., 30000, 0.001, 5., 4., 2., 0.75, 3., false)

# circularmeant = circularmean_t(ϕ, 30000)

# plot_orderparameter(circularmeant, 30000, "circularmean.png")

# animation_polarisation(ϕ, x, y, 30000, 30, "polarisation6.mp4")


# print(rad2deg(circularmeant[30000]))



angle = getendangel(10, 1., 1.,  21, 40, 1., 20000, 0.001, 5., 4., 2., 0.75, 3., false)

print(rad2deg.(angle))

# animation_deformation(ux, uy, x, y, 30000, 30, "deformation6.mp4")

# orderparameter_t = orderparameter(ϕ, m*l, 50000)

# plot_orderparameter(orderparameter_t, 50000, "orderparameter6.png")



# #fig  = Figure(size = (800, 800))

#xisweep(10, 1., 1.,  21, 40, 1., 100000, 0.001, 5., 4., 2., 0.75, true)









