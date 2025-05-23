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

    animation_polarisation(ϕ1, x1, y1, 20000, 30, "animation_polarisation_stick_xi$i.mp4")

    animation_deformation(ux1, uy1, x1, y1, 20000, 30, "animation_deformation_stick_xi$i.mp4")

end

