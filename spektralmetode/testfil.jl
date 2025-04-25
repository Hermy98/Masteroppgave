using LinearAlgebra, CairoMakie, SparseArrays, Integrals, FFTW, Dierckx, GLMakie, FileIO, Statistics

include("plottingfunctions.jl")

include("Utility.jl")

#include("elasticsystem.jl")

#include("polarization.jl")

include("Activesystem.jl")





# Arrays to store monitoring data
max_velocity = zeros(50000)
avg_velocity = zeros(50000)

@time ux1, uy1, ϕ1, x1, y1, m1, l1, max_velocity1, avg_velocity1, elastic_modes, polar_modes = run_activesystem(10, 1., 1., 40, 80, 1., 50000, 0.001, 5., 1., 0., 0.1, 1., true)

animation_deformation(ux1, uy1, x1, y1, 50000, 30, "deformation5.mp4")

animation_polarisation(ϕ1, x1, y1, 50000, 30, "polarisation5.mp4")





# orderparameter_t1 = orderparameter(ϕ1, m1*l1, 50000)

# # Create a figure with 4 subplots
# fig = Figure(size=(1200, 1000))

# # Order parameter
# ax1 = Axis(fig[1,1], title="Order Parameter")
# lines!(ax1, 1:50000, orderparameter_t1)
# ax1.xlabel = "Time"
# ax1.ylabel = "Order Parameter"

# # Velocity statistics
# ax2 = Axis(fig[1,2], title="Velocity Statistics")
# lines!(ax2, 1:50000, max_velocity1, label="Max |dt_u|")
# lines!(ax2, 1:50000, avg_velocity1, label="Avg |dt_u|")
# axislegend(ax2)
# ax2.xlabel = "Time"
# ax2.ylabel = "Velocity"

# # Elastic modes
# ax3 = Axis(fig[2,1], title="Elastic Modes")
# for k in 2:size(elastic_modes, 1)
#     lines!(ax3, 1:50000, elastic_modes[k,:], label="k=$(k-1)")
# end
# axislegend(ax3)
# ax3.xlabel = "Time"
# ax3.ylabel = "Mode Amplitude"

# # Polarization modes
# ax4 = Axis(fig[2,2], title="Polarization Modes")
# for k in 2:size(polar_modes, 1)
#     lines!(ax4, 1:50000, polar_modes[k,:], label="k=$(k-1)")
# end
# axislegend(ax4)
# ax4.xlabel = "Time"
# ax4.ylabel = "Mode Amplitude"

# save("system_analysis_long5.png", fig)

# # # Also create snapshots of the system state at key times
# # # Just before breakdown (t=19000), during breakdown (t=20000), and after (t=21000)
# # for t in [19000, 20000, 21000]
# #     fig_snap = Figure(size=(800, 400))
# #     ax_pol = Axis(fig_snap[1,1], title="Polarization t=$t")
# #     ax_vel = Axis(fig_snap[1,2], title="Velocity t=$t")
    
# #     # Polarization field
# #     arrows!(ax_pol, x1, y1, cos.(ϕ1[:,:,t]), sin.(ϕ1[:,:,t]))
    
# #     # Velocity field
# #     arrows!(ax_vel, x1, y1, ux1[:,:,t], uy1[:,:,t])
    
# #     save("snapshot_t$(t)4.png", fig_snap)
# # end


# animation_deformation(ux1, uy1, x1, y1, 50000, 30, "deformation5.mp4")

# animation_polarisation(ϕ1, x1, y1, 50000, 30, "polarisation5.mp4")


