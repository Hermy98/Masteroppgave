using LinearAlgebra, Plots, SparseArrays

function multidimmatrix(N::Int, M::Int, α::Float64, β::Float64, γ::Float64, r)


    S = N*M


    beta_array = β*ones(S)./r

    gamma_array = γ*ones(S)./(r)^2

    diag = ones(S)*(1 -2* α) .+ 2 .*gamma_array

    off_diagup = ones(S-1).*gamma_array[1:end-1]
    off_diaglo = ones(S-1).*gamma_array[2:end]

    m_diaglo = ones(S-M)*(α)

    m_diagup = ones(S-M)*(α)

    periodic_intup = zeros(S-(M-1))

    periodic_intlo = zeros(S-(M-1)) 


    for i in M:M:S

        periodic_intlo[i-(M-1)] = 1

        periodic_intup[i-(M-1)] = 1


        if i < S

            off_diagup[i] = 0

            off_diaglo[i] = 0

        end
    
    end

    diag[S-M+1:end] .= 1

    off_diagup[S-(M-1):end] .= 0

    off_diaglo[S-(M-1):end] .= 0

    m_diagup[M+1:end] .= 0

    m_diaglo[M+1:end] .= 0

    periodic_intlo[end] = 0

    periodic_intup[end] = 0

    print(m_diaglo)

    
    #print lenth of all arrays

    print(length(diag), length(off_diagup), length(off_diaglo), length(m_diagup), length(m_diaglo), length(periodic_intup), length(periodic_intlo))

    


    return spdiagm(0 => diag, 1 => off_diagup, -1 => off_diaglo, M => m_diagup, -(M) => m_diaglo, -(M-1) => periodic_intup, M-1 => periodic_intlo)


end



A = multidimmatrix(3, 4, 1., 1., 1., 0:10:1)



heatmap(A)