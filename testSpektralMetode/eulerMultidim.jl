using LinearAlgebra, Plots, SparseArrays

function multidimmatrix(N::Int, M::Int, α::Float64, β::Float64, γ::Float64, r)


    S = N*M

    diag = ones(S)*(1 -2* α)

    off_diagup = ones(S-1)*(γ)
    off_diaglo = ones(S-1)*(γ)

    m_diaglo = ones(S-M)*(α)

    m_diagup = ones(S-M)*(α)

    periodic_intup = zeros(S-(M-1)) # migth not work for all

    periodic_intlo = zeros(S-(M-1)) # migth not work for all
    diag[(N*(M-1)):end] .= 1

    for i in M:M:S

        periodic_intlo[i-(M-1)] = 1

        periodic_intup[i-(M-1)] = 1


        if i < S

            off_diagup[i] = 0

            off_diaglo[i] = 0

        end
    
    end

    diag[(N*(M-1)):end] .= 1

    return A


end



A = multidimmatrix(3, 4, 1., 1., 1., 0:10:1)


heatmap(A)