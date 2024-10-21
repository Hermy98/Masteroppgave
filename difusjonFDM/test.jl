using Plots, LinearAlgebra, Statistics


function get_matrix(N, α, β, r, d_t)


    beta_array = β*ones(N)./(2 .*r)

    
    off_diagup = -1*ones(N-1)*(α/2)
    diag = ones(N)*(1 + α)  # N for interpolasjon
    off_diaglo = -1*ones(N-1)*(α/2)#N-1 for interpolasjon

    off_diaglo .+= beta_array[2:end]

    off_diagup .-= beta_array[1:end-1]

    diag[2:end] = diag[2:end] .-(4*d_t ./(r[2:end].^2))


    M = Tridiagonal(off_diaglo, diag, off_diagup)


    return M

end

function initialcondition(N)

    u_1 = zeros(N)

    for i in 1:N
        u_1[i] = 5+((i-1)*(i-100))/2000
    end

    return  u_1

end


function initialcondition2(r, θ)

    u1 = zeros(length(2*r), length(θ/2))

    r1 = range(0, 2*pi, length=length(2*r))

    for j in 1:length(θ/2)

        u1[:, j] = cos.(r1)

    end

    return  u1

end


function initialcondition3(r, θ)

    u2 = zeros(length(r), length(θ))

    r1 = range(0, 6*pi, length=length(r))

    #u2[end, : ] = cos.(r1)

    for i in eachindex(θ)

        u2[:, i] = (r.^2 .+4).*cos.(θ[i])

    end

    return  u2

end

origoInterpolasjon(u) = (1/3)*(4*mean(u[2, :]) - mean(u[3, :]))

#u_originLaplace(u, d_r) = u[2] - ((d_r^2)/4)*u[1] 



function initialilize(N, J, d_t, D)

    r = range(0, 1, length=N)

    Θ = range(0, 2*π, length=J)

    d_r = r[2] - r[1]

    α = (D*d_t)/(d_r^2)

    β = (D*d_t)/(2*d_r)


    A = get_matrix(N, α, β, r, d_t)

    B = get_matrix(N, -α, -β, r, d_t)

    #setting the boundary conditions

  #  A[1,1] *= 2

   # B[1,1] *= 2

    A[1,2] = -α/2

    B[1,2] = α/2

    A[end, end] = 1
    A[end, end-1] = 0

    B[end, end] = 1
    B[end, end-1] = 0

    return r, Θ, A, B, d_r


end


function main(N, J, T, d_t, D)

    r, Θ, A, B , d_r = initialilize(N, J, d_t, D)

    u = zeros(N, J, T)

    u[:, :, 1] = initialcondition3(r, Θ)


    for i in 1:T-1
        for j in 1:J
            u[:, j, i+1] = (A\B)*u[:, j ,i]

            #u[1, :,  i+1] .= mean(u[1, :, i+1])

            #u[1, :,  i+1] .= origoInterpolasjon(u[:, :, i+1])
        end

    end

    return u, r, Θ
 

end


function coordinatetransform(Θ, r, N, J)
    

    X = r'.*cos.(Θ)

    Y = r'.*sin.(Θ)


    return X, Y

end
@time u , r, Θ = main(100, 100, 1000, 1E-5, 1)


x, y = coordinatetransform(Θ, r, 100, 100)
display(surface(x, y, transpose(u[:, :, 2])))
display(surface(x, y, transpose(u[:, :, end])))

